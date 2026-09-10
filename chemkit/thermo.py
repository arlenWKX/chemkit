"""chemkit.thermo —— 独立温度模块（非恒温热效应计算与绝热耦合求解）。

设计纪律（与平衡求解完全解耦）：
    本模块只消费 judge() 的结果。事后分析（analyze）从净变换（初态→终态
    +逸出）出发，用焓的状态函数性质计算反应热与溶液温升；耦合求解
    （coupled）把 judge() 当黑盒在外层迭代温度不动点——平衡求解路径上
    不加任何判断语句，调用点只有 system.Reaction 包装层。

方法（焓是状态函数，无需逐步反应热）：
    1. 净变换 net[物种] = production − consumption（production 已含逸出
       气体——它们是携带 ΔHf(g) 离开体系的真实产物，不减去 escaped；
       否则产气体系丢焓且 H/O 原子平衡必不闭合）+ H⁺/OH⁻ 账户差
       （H_excess 初/终值展开 + 规范化中和步的消耗——中和发生在初态
       快照之前，其 H⁺/OH⁻ 必须显式计入 net，否则纯中和体系报零热）
    2. 水量重构：整个变换（含逸出气体）原子守恒 ⇒ 由全部物种的
       H、O 原子量净差解出 Δn(H₂O)；两条独立路径（H 平衡/2 与 O 平衡）
       互为校验，不一致（>2% 或 1e-3 mol）即拒绝并报告。
    3. ΔH = Σ net·ΔHf°(物种) + Δn(H₂O)·ΔHf°(H₂O,l)；
       放热 q = −ΔH（kJ）。
    4. ΔT = q×1000 / (m·c_p)，m = V×1000 g（以水的比热容与密度为准，
       忽略溶质质量——教材级近似），c_p = 4.184 J/(g·K)。
    5. ΔHf° 单一数据源 data/thermo.json（裸名 = 账本态，"X(g)" 键 =
       逸出气相态）；净变换中含量 ≥1e-3 mol 的物种缺数据时返回
       heat=None（宁缺毋假），<1e-3 mol 的痕量物种忽略。

绝热耦合求解 coupled()（isothermal=False 的真实语义）：
    物理问题：溶液初温 T₀，反应绝热进行，终态 = 平衡(T_f) 且能量守恒：
        m·c_p·(T_f − T₀) = q(T_f)
    外层不动点（锚定初温——q 随 T 单调减弱（放热反应 van't Hoff），
    映射 T ↦ T₀ + q(T)/(m·c_p) 是收缩映射，实测 2-4 轮收敛）：
        1. judge() 在 T_cur 求化学平衡；
        2. analyze() 算 q → T_next = T₀ + q/(m·c_p)；
        3. |T_next − T_cur| < 0.25 K 收敛；否则 T_cur ← T_next 重解。
    相变平台：T_next 越出 273.15–373.15 K 时钳在边界（沸腾/凝固体系
    恒温于相变点，多余热量走潜热——线性外推不建模，如实给标志）。

约定与已知边界（诚实记录）：
    · 投料按"已溶解的溶质"计（引擎语义：可溶性强电解质投料即电离），
      固体投料的溶解热不计（如 NaOH(s) 溶解 −44.5 kJ/mol）。
    · SO₃ 等在规范化阶段与水反应的投料：H/O 平衡不闭合（投料端形态
      未入账），此时返回 None 并注明原因。
    · 逸出气体按同温离开体系计（不追踪气相显热）。
    · ΔCp≈0（Kirchhoff）近似：反应焓不随 T 变；热容只计水的 m·c_p。
    · T_final 超出 273.15–373.15 K 时给出 boils/freezes 标志——线性外推
      未计汽化/凝固潜热，仅提示物理上已达相变。
    · 逸出气相态拆分：escaped 部分按 "X(g)" 键计价（CO₂/SO₂/H₂S/
       NH₃ 有数据；无键的气体如 H₂/NO 本就同值或元素零），残留溶解
      态按水溶值——两者焓差即溶解热如实入账。
"""

from .core import elements_of
from .data import dhf_of

# ---- 水物性（教材级常数） ----
C_P_WATER = 4.184          # J/(g·K)
RHO_WATER = 1000.0         # g/L（25 °C，忽略溶质贡献）
DH_WATER = -285.83         # kJ/mol  H2O(l)
T_BOIL = 373.15
T_FREEZE = 273.15

# 绝热耦合参数
MAX_OUTER = 6              # 外层不动点最大轮数（收缩映射 2-4 轮收敛）
TOL_K = 0.25               # 温度收敛容差（K）

# ΔHf° 数据面：data.thermo.json（v0.3.5 起单一事实源；v0.3.4 及之前本
# 模块自带硬编码 DHF 副本，与数据表双源冲突——已统一，见 data.dhf_of）。
# 裸名 = 账本态（水溶/凝聚）；"X(g)" 键 = 逸出气相态。

# 净变换中缺 ΔHf 数据仍可忽略的摩尔上限（热贡献 ~0.5 J 量级）
TRACE_MOL = 1e-3


def net_transformation(res: dict) -> dict[str, float]:
    """净变换 {物种: mol}（正值=生成，负值=消耗；含逸出气体）。

    production − consumption（引擎的 production 已把逸出气体并入——
    它们是携带 ΔHf(g) 离开体系的产物，参与焓加和与原子平衡；减掉
    会使产气体系丢焓且 H/O 重构必不闭合）+ 游离 H⁺/OH⁻ 账户差
    （H_excess 初/终各按正负号展开）+ 规范化中和步的 H⁺/OH⁻ 消耗
    （中和发生在初态快照之前，不补入则纯中和体系净变换恒为空）。
    """
    prod = {e["name"]: e["mol"] for e in res.get("production", [])}
    cons = {e["name"]: e["mol"] for e in res.get("consumption", [])}
    net = {s: prod.get(s, 0.0) - cons.get(s, 0.0)
           for s in set(prod) | set(cons)}
    he0 = res.get("H_excess_initial", 0.0)
    hef = res.get("H_excess", 0.0)
    d_h = max(hef, 0.0) - max(he0, 0.0)          # 游离 H⁺ 净变化
    d_oh = max(-hef, 0.0) - max(-he0, 0.0)       # 游离 OH⁻ 净变化
    if d_h:
        net["H^+"] = net.get("H^+", 0.0) + d_h
    if d_oh:
        net["OH^-"] = net.get("OH^-", 0.0) + d_oh
    # 规范化阶段的强酸强碱中和（初态快照之前）：H⁺ + OH⁻ → H₂O
    # 的反应物端显式计入净变换（水由原子平衡重构，不在此加）
    for st in res.get("steps", []):
        if st.get("kind") == "neutralize":
            e = st.get("extent", 0.0)
            if e > 0:
                net["H^+"] = net.get("H^+", 0.0) - e
                net["OH^-"] = net.get("OH^-", 0.0) - e
    return {s: v for s, v in net.items() if abs(v) > 1e-9}


def water_formed(net: dict[str, float], res: dict) -> float | None:
    """净变换生成的水量（mol）。

    H 平衡：Σ net·n_H + 2·Δn(H₂O) = 0
    O 平衡：Σ net·n_O + Δn(H₂O) = 0
    （net 含逸出气体与中和步反应物端，原子平衡自动闭合；
    两路一致（容差 max(2%, 1e-3)）才返回；否则 None——变换未闭合，
    典型为 SO₃ 类投料在规范化阶段与水反应。）
    """
    dH_atoms = 0.0
    dO_atoms = 0.0
    for s, v in net.items():
        n = elements_of(s)
        dH_atoms += v * n.get("H", 0)
        dO_atoms += v * n.get("O", 0)
    w_h = -dH_atoms / 2.0
    w_o = -dO_atoms
    if abs(w_h - w_o) > max(0.02 * max(abs(w_h), abs(w_o)), 1e-3):
        return None
    return (w_h + w_o) / 2.0


def analyze(res: dict, T) -> dict | None:
    """非恒温分析：从 judge() 结果与数据表计算 (heat_kJ, dT_K, T_final_K)。

    ΔHf 取自数据表（v0.3.5 起单一事实源；T 为 load_tables() 的 Tables）。
    **逸出气相拆分**：净变换中的气体全部按账本态（水溶）计价，逸出部分
    再按 "X(g)" 气相键校正（溶解→逸出的焓差：CO₂ 20.3 / NH₃ 34.4
    kJ/mol）——比 v0.3.4 的「全部按气相值」更准确（残留溶解态不再
    多计）。

    返回 dict（可能含 None 字段）：
        heat_kJ     放热为正（= −ΔH；None = 数据不足，宁缺毋假）
        dT_K        温升（负=降温；None 同上）
        T_final_K   终温（数据不足时 = 初温）
        water_mol   净生成水量
        mass_g      供热容用的水质量（V×1000 g）
        missing     缺 ΔHf 数据的显著物种（|net| ≥ 1e-3）
        flags       ["boils"/"freezes", "h_o_open"...]
        reason      失败原因（数据不足/H-O 不闭合/无变换）
    """
    cond = res.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    T_K = float(cond.get("T_K", 298.15))
    out = {"heat_kJ": None, "dT_K": None, "T_final_K": T_K,
           "water_mol": None, "mass_g": V * RHO_WATER,
           "missing": [], "flags": [], "reason": None}
    net = net_transformation(res)
    if not net:
        out["reason"] = "无净变换（无反应或纯形态变化）"
        out["water_mol"] = 0.0
        out["heat_kJ"] = 0.0
        out["dT_K"] = 0.0
        return out
    w = water_formed(net, res)
    if w is None:
        out["flags"].append("h_o_open")
        out["reason"] = "H/O 原子不闭合（投料在规范化阶段转变，如 SO3+H2O）"
        return out
    out["water_mol"] = w
    missing = sorted(s for s, v in net.items()
                     if dhf_of(T, s) is None and abs(v) >= TRACE_MOL)
    if missing:
        out["missing"] = missing
        out["reason"] = f"缺 ΔHf° 数据：{', '.join(missing[:6])}"
        return out
    # 逸出气相拆分校正：escaped ⊆ production；gas 键存在时逸出部分按
    # ΔHf(g) 计价（账本态基线 + (g−aq)·esc）
    esc = {e["name"]: e["mol"] for e in res.get("escaped", [])}
    dH = 0.0
    for s, v in net.items():
        dH += v * dhf_of(T, s)
        e = esc.get(s)
        if e:
            g = T.thermo.get(s + "(g)")
            if g is not None:
                dH += e * (g - dhf_of(T, s))
    dH += w * DH_WATER
    q = -dH                       # kJ，放热为正
    out["heat_kJ"] = round(q, 3)
    dT = q * 1000.0 / (out["mass_g"] * C_P_WATER)
    out["dT_K"] = round(dT, 3)
    out["T_final_K"] = round(T_K + dT, 2)
    if out["T_final_K"] > T_BOIL:
        out["flags"].append("boils")
        out["reason"] = ("终温超过沸点（线性外推未计汽化潜热；"
                         "实际溶液沸腾恒温于 373.15 K）")
    elif out["T_final_K"] < T_FREEZE:
        out["flags"].append("freezes")
        out["reason"] = "终温低于冰点（线性外推未计凝固潜热）"
    return out


# ========================================================== 绝热耦合求解

def _clamp_phase(T_K: float) -> tuple[float, list[str]]:
    """液态水域钳制（相变平台）：越界停在边界并给出标志。"""
    flags = []
    if T_K > T_BOIL:
        T_K = T_BOIL
        flags.append("boils")
    elif T_K < T_FREEZE:
        T_K = T_FREEZE
        flags.append("freezes")
    return T_K, flags


def coupled(substances: list[dict], conditions: dict, T) -> dict:
    """非恒温（绝热）耦合求解：化学平衡 ⇌ 能量平衡的外层温度不动点。

    参数同 judge()（substances/conditions/tables）。返回 judge() 结果 dict，
    附加 "thermal" 键：
        heat_kJ     放热为正（终态平衡下的净反应热；None=数据不足）
        dT_K        溶液温升 = T_final − T₀（锚定初温，能量守恒口径）
        T_final_K   终温（末轮能量平衡估计 T₀ + q/(m·c_p) 的钳制值；
                    化学解在 T_cur 求得，与 T_final 相差 < 收敛容差 0.25 K）
        trace       逐轮 [{T_K, heat_kJ, dT_K, clamped}]（含第一轮初温解）
        converged   外层不动点是否收敛（|T_next − T_cur| < 0.25 K）
        flags       boils/freezes（末轮相变平台钳制）
        其余字段同 analyze()（water_mol/mass_g/missing/reason）

    语义：
        · 热效应 |ΔT| < 0.25 K 或数据不足 → 单遍直出（与事后分析等价）；
        · 化学结果取**收敛温度下的平衡**（K(T) 经 van't Hoff 随 T 移动，
          温度反馈真实进入求解——非"初温解完再补一个 ΔT 数字"的事后
          假变温）；
        · 收敛判据锚定 T₀：T_next = T₀ + q(T_cur)/(m·c_p)，放热反应
          q 随 T 单调减弱 ⇒ 收缩映射（相对更新 T += dT 会收敛到 q=0 的
          错误不动点——那是"把反应热耗干"而不是"绝热升温"）。
    """
    from .engine import judge        # 延迟导入（engine 不依赖本模块，无环）
    cond = dict(conditions or {})
    T0 = float(cond.get("T_K", 298.15))
    res = judge(substances, cond, T)
    th = _safe_analyze(res, T)
    dT = th.get("dT_K")
    trace = [{"T_K": round(T0, 2), "heat_kJ": th.get("heat_kJ"),
              "dT_K": dT, "clamped": False}]
    if dT is None or abs(dT) < TOL_K:
        # 近零热效应或数据不足：初温解即自洽绝热解（终温按能量平衡口径，
        # 与化学解温度相差 < 容差）
        T_fin, ph = _clamp_phase(T0 + (dT or 0.0))
        res["thermal"] = _pack(th, trace, True, ph, T0, T_fin)
        return res
    T_cur = T0
    flags: list[str] = []
    converged = False
    for _ in range(MAX_OUTER):
        T_next, ph = _clamp_phase(T0 + dT)   # 锚定初温（能量守恒对 T₀ 结算）
        if abs(T_next - T_cur) < TOL_K:
            converged = True
            flags = ph
            break
        T_cur = T_next
        flags = ph
        res = judge(substances, dict(cond, T_K=T_cur), T)
        th = _safe_analyze(res, T)
        trace.append({"T_K": round(T_cur, 2), "heat_kJ": th.get("heat_kJ"),
                      "dT_K": th.get("dT_K"), "clamped": bool(ph)})
        dT = th.get("dT_K")
        if dT is None:
            # 新温度下的平衡出现缺 ΔHf 物种（形态随温度改变）：保留当前
            # 化学解与已有热数据，如实标注未收敛
            res["thermal"] = _pack(th, trace, False, flags, T0, T_cur)
            return res
    # 末轮能量平衡估计（收敛时与 T_cur 相差 < 容差；轮数耗尽时为
    # 最好估计并如实报未收敛）
    T_fin, ph = _clamp_phase(T0 + dT)
    flags = ph if converged else sorted(set(flags) | set(ph))
    res["thermal"] = _pack(th, trace, converged, flags, T0, T_fin)
    return res


def _pack(th: dict, trace: list, converged: bool, flags: list,
          T0: float, T_fin: float) -> dict:
    """组装 Reaction.thermal 视图：dT/T_final 以初温 T₀ 为基准。"""
    out = dict(th)
    out["trace"] = trace
    out["converged"] = converged
    phase = {"boils", "freezes"} & set(flags)
    out["flags"] = sorted(set(out.get("flags", [])) | set(flags))
    if phase:
        out["reason"] = ("终温达相变（耦合钳制在相变平台 373.15/273.15 K；"
                         "线性外推未计汽化/凝固潜热）")
    out["T_final_K"] = round(T_fin, 2)
    out["dT_K"] = round(T_fin - T0, 3)
    return out


def _safe_analyze(res: dict, T) -> dict:
    """analyze 的防御包装：异常不拖垮主结果（热层纪律）。"""
    try:
        th = analyze(res, T)
    except Exception:                     # noqa: BLE001
        th = {"heat_kJ": None, "reason": "温度模块异常"}
    return th or {"heat_kJ": None, "reason": "温度模块异常"}
