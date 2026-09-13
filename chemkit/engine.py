"""chemkit.engine：水溶液反应判定主引擎（预测 walk 与 pH 估计）。

六条核心：
  ① 候选即平衡（logK 来自电对/pKa/pKsp/logβ 或其 Hess 精确加和）；
  ② pH 连续，由账本现场推算；
  ③ S = logK − logQ（实际活度），永远生效；
  ④ 执行 = 二分求解平衡程度 x*（S(x*)=0 或计量上限）；
  ⑤ H_excess 单一带符号账本，方程一律 H+ 正则形（无 OH-）；
  ⑥ 程度判定按极限试剂转化率。

动力学层（与热力学显式分离）：gate 闸门 / slow 标注 / 膜 blocked /
溶剂优先（仅此一条界面动力学规则）/ OVERRIDE 逃生舱。

v0.3.4 起模块族分工（单向依赖链 core←data←candidates←normalize←
speciation/templates←engine←system←interfaces）：
    candidates  Cand/logK_T/配平/兜底电离/派生候选/全局常数
    normalize   投料规范化（投料→账本）
    speciation  pH 机器（滴定/形态分布/强酸再电离）
    templates   候选宇宙（红氧模板四级缓存 + 枚举）
    engine      本模块：S_of/solve_extent/blocked_extent/judge walk
下方再导出保持 v0.3.3 前经由 engine 访问全部符号的兼容路径。
"""
from __future__ import annotations
import os as _os
from math import log10

from .core import (elements_of, pKw_of, _vant, k_nernst,
                   K_NERNST_298, PKW_298, charge_of, balance)
from .acidbase import ledger_charge as _ledger_charge
from .data import Tables, henry_of
from .candidates import (Cand, logK_T, build_derived, _bal, _bal_fast,
                         _half_pairs, _redox_mix_ok, _half_scale,
                         _split_salt_general, _ion_name,
                         WATER, H_ION, X_MIN, ACT_FLOOR, A_GAS,
                         P_EXT_KPA, P_STD_KPA, P_RES_KPA, SAT_SKIP,
                         SOLVENT_FIRST, BLOCKED_EXTENT, ANN_MIN_EXTENT,
                         DEGREE_COMPLETE, DEGREE_PARTIAL, MAX_ITER,
                         WATER_MOL_PER_L, STRONG_MOLECULAR_ACIDS,
                         WATER_FIRST_METALS, NONMETAL_SOLIDS,
                         HALATE_DISP_T, HALATE_BASE_PH, _NONMETAL_ELEMS,
                         _STRONG_ACID, _HALF_CACHE, _PROTON_ROOTS,
                         _set_proton_roots)
from .normalize import (normalize, _mol_fraction, _split_acid,
                        _salt_ksp_cell, _ionize_map, _split_salt)
from . import speciation as _speciation
from .speciation import (_buffer_titration, estimate_pH, estimate_state,
                         _respeciate_strong_acids, _full_speciation,
                         _RESPECIATE_ACIDS, _pksp, closed_pH,
                         weak_species_set, exact_proton_pH)
from .joint import joint_solve, _solve_ph, JOINT_MAX_M, JOINT_MIN_M
from .templates import (_redox_pair_static, _redox_templates,
                        _oxide_dissolve_info, _build_static_cands,
                        _gate_species, enumerate_candidates, _gate_check,
                        _acid_conc, _solid_acid_pka, _ksp_xy, _ksp_asat,
                        _product_form_ok)

_TRACE = bool(_os.environ.get("CHEM_TRACE"))

# ---- 确定性性能计数器（v0.5.0 性能审计，architecture §7 W-4）----------
# 墙钟在共享沙箱里逐轮抖动 ±10–15%（同一份代码三次全量：均值
# 29.1 / 32.5 / 33.9 ms），**不足以裁决"求根改动是否真的省了"**。
# `S_of` 调用次数是确定性的：同一份代码 ⇒ 同一计数，且它正是走步每步
# 的钱主要花在哪（§7 U-1 剖面：S_of 33% + pH 机器 20%，每次 f 各一遍）。
# 默认关闭（每次调用一次全局判空，开销可忽略）；`converg.dump()` 临时置位。
SOF_CALLS: list[int] | None = None

# ---- 求根审计钩子（仅 tools/roots.py 启用；生产路径恒为 None）---------
# 待裁决的问题（§7 W-5 之后剩下的最大性能杠杆）：二分在主求解分支固定
# ~52 次迭代（容差落地后 ~35 次求值/步），而**保括号的假位法**（Illinois）
# 在光滑单调括号上 ~10 次即可到同一精度 ⟹ 每步成本还能再降一大截。
# 唯一的语义风险是**根选定**：括号内 f 非单调（多根）时，快方法可能收敛到
# 另一个根，而"选中哪个根"就是"实际执行多少"就是化学——不能靠断言兜底
# （断言锁化学，但换根 = 换化学，是**必须**避免的位移）。故先审计：
# 在同一括号上同时跑两种方法，记录根差、求值次数与网格符号翻转数。
ROOT_AUDIT: dict | None = None

# 二分收敛容差（相对 x_max 的绝对值下限见 solve_extent）：见 §7 W-4 的
# 论证与实测。置 0.0 可复现"跑满迭代到浮点饱和"的旧路径（实验用）。
_EXTENT_TOL_REL = 1e-11

# ========================================================== ③ S = logK − logQ

def _logc_of(s: str, ledger: dict, V: float, logc: dict | None) -> float:
    """log10(max(c/V, ACT_FLOOR))——同一表达式集中一处；logc 非空时按物种
    记忆（账本不变期间跨候选/跨二分点复用，值与逐次计算 bit 级一致，
    由调用方在账本物种变动时失效对应条目）。"""
    if logc is None:
        return log10(max(ledger.get(s, 0.0) / V, ACT_FLOOR))
    v = logc.get(s)
    if v is None:
        v = log10(max(ledger.get(s, 0.0) / V, ACT_FLOOR))
        logc[s] = v
    return v


def S_of(c: Cand, ledger: dict, V: float, pH: float, T_K: float, T,
         gsup: frozenset = frozenset(), p_ext_kpa: float = P_EXT_KPA,
         gas_escape: bool = True, logc: dict | None = None) -> float:
    """计算候选反应的亲和势 S = logK − logQ。

    p_ext_kpa: 外界气相总压（kPa）。低于常压时气体更易逸出（泡点降低），
    高于常压时更多气体留在溶液。默认 P_EXT_KPA（101.3 kPa 常压）。
    gas_escape: False = 闭口体系——自产气体不逸出，活度按溶解态浓度
    c/c° 计（产物积累抑制反应，勒沙特列效应），与投料气体同待遇。
    logc: 浓度对数缓存（可选；绑定当前 ledger，由调用方维护失效）。

    实现：项执行计划按 Cand 缓存（物种→项类型分类只依赖数据表），
    求值只做查表与算术；项序与原逐项分支一致（r 全项后 pr 全项），
    logQ 累加顺序不变 ⇒ 数值 bit 级等价。
    """
    if SOF_CALLS is not None:
        SOF_CALLS[0] += 1
    tid = id(T)
    plan = c._plan
    if plan is None or c._plan_tid != tid:
        terms = []
        for s, nu in c.r.items():
            if s == H_ION:
                terms.append((0, s, nu))
            elif s == WATER or s in T.solids:
                pass
            else:
                terms.append((2, s, nu))
        for s, nu in c.pr.items():
            if s == H_ION:
                terms.append((1, s, nu))
            elif s == WATER or s in T.solids:
                pass
            elif s in T.gases:
                terms.append((4, s, nu))
            else:
                terms.append((3, s, nu))
        plan = tuple(terms)
        c._plan = plan
        c._plan_tid = tid
    logQ = 0.0
    for typ, s, nu in plan:
        if typ == 2:
            logQ -= nu * _logc_of(s, ledger, V, logc)
        elif typ == 3:
            logQ += nu * _logc_of(s, ledger, V, logc)
        elif typ == 0:
            logQ += nu * pH
        elif typ == 1:
            logQ -= nu * pH
        elif gas_escape and s not in gsup:
            # 气体产物活度：外加供给的气体按账本浓度（持续供给维持）；
            # 自产气体超过泡点（c > H(T)·p_ext）的部分已被扫气移入账本
            # （见 _sweep_gases，逸出即消失），账上浓度恒 ≤ 饱和浓度，
            # 活度 a = c/c°（饱和时 c = H(T)·p_ext，等价 p/p° = p_ext/p°），
            # 下限 A_GAS 为惰性环境残余分压约定（避免浓度地板制造虚假驱动）；
            # 无 Henry 数据的物种回退固定 A_GAS（视为全逸出）
            H = henry_of(T, s, T_K)
            c_g = ledger.get(s, 0.0) / V
            if H is not None:
                a = max(min(c_g, H * p_ext_kpa), A_GAS)
            else:
                a = A_GAS
            logQ += nu * log10(a)
        else:
            # H+ 项之外的兜底（gas_escape=False 的气体 / 持续供给气体）：
            # 活度一律按溶解态浓度计
            logQ += nu * _logc_of(s, ledger, V, logc)
    return logK_T(c, T_K) - logQ


# ========================================================== ④ 平衡程度求解

def _sweep_gases(ledger: dict, escaped: dict, gsup: frozenset,
                 V: float, T_K: float, p_ext_kpa: float, T,
                 gas_escape: bool = True) -> None:
    """泡点扫气：自产气体（非持续供给）溶解浓度超过 c_sat = H(T)·p_ext 时，
    超额部分鼓泡逸出，移入 escaped 账户——逸出即消失，不再参与任何后续
    反应（取代旧设计"只压活度不离账"导致的隔空反应与 rev_gate 补丁）。
    低于泡点的气体保持溶解（惰性环境下不强制脱气）。无 Henry 数据的物种
    不扫（由 S_of 的 A_GAS 回退承载其逸出驱动）。"""
    if not gas_escape:
        return    # 闭口体系：不扫气，自产气体保留在溶液账本
    # sorted：T.gases 是 set，迭代序随 PYTHONHASHSEED 变化——各气体扫气
    # 相互独立（量值与序无关），但 escaped 字典的插入序（及下游 esc_list/
    # production 呈现序）会跨进程翻转（N08/H25 的 SO2/CO2 顺序掷骰子）。
    # 固定字典序后结果呈现跨进程确定
    for g in sorted(T.gases):
        if g in gsup:
            continue
        m = ledger.get(g, 0.0)
        if m <= X_MIN:
            continue
        H = henry_of(T, g, T_K)
        if H is None:
            continue
        cap = H * p_ext_kpa * V
        if m > cap:
            escaped[g] = escaped.get(g, 0.0) + (m - cap)
            ledger[g] = cap


def _audit_bracket(f, c, direction: int, x_max: float, x_bis: float,
                   n_bis: int, f_hi: float, micro: bool) -> float | None:
    """求根审计（仅 `tools/roots.py` 启用；生产路径不调用）。

    在**同一个括号** [0, x_max] 上记录三件事：

      ① 生产二分的根与求值次数（`x_bis`/`n_bis`，由调用方传入）；
      ② **Illinois 保括号假位法**的根与求值次数（`x_ill`/`n_ill`）——
         它收敛到同一精度的求值次数就是这块性能的量级；
      ③ f 在 33 点网格上的**符号翻转次数**：>1 ⟹ 括号内多根，快方法
         可能选中另一个根。这不是"测试被锁"的问题，而是"换根 = 换
         实际执行量 = 换化学"，必须避免。

    返回 `x_ill`（供调用方在"换根"时补一份 pH 轨迹细扫）。
    """
    n = 33
    signs = [f(x_max * k / (n - 1)) > 0 for k in range(n)]
    flips = sum(1 for i in range(1, n) if signs[i] != signs[i - 1])
    tol = _EXTENT_TOL_REL * max(1.0, x_max)
    f0 = f(0.0)
    n_ill = 1
    x_ill = 0.0
    if f0 > 0:
        lo, hi = 0.0, x_max
        flo, fhi = f0, f_hi
        last = ""
        while hi - lo > tol and n_ill < 200:
            den = fhi - flo
            x = (lo * fhi - hi * flo) / den if den != 0 else 0.5 * (lo + hi)
            if not (lo < x < hi):
                x = 0.5 * (lo + hi)
            vx = f(x)
            n_ill += 1
            if vx > 0:
                lo, flo = x, vx
                if last == "lo":
                    fhi *= 0.5        # Illinois 端点衰弱
                last = "lo"
            else:
                hi, fhi = x, vx
                if last == "hi":
                    flo *= 0.5
                last = "hi"
        x_ill = lo
    ROOT_AUDIT["rec"].append({
        "case": ROOT_AUDIT.get("case"),
        "eq": (" + ".join(f"{_fmt(nu)}{s}" for s, nu in c.r.items() if s != WATER)
               + " -> "
               + " + ".join(f"{_fmt(nu)}{s}" for s, nu in c.pr.items()
                            if s != WATER)),
        "kind": c.kind, "dir": direction, "x_max": x_max,
        "x_bis": x_bis, "x_ill": x_ill, "n_bis": n_bis, "n_ill": n_ill,
        "flips": flips, "f0": f0, "micro": micro, "f": f,
    })
    return x_ill if f0 > 0 else None


def solve_extent(c: Cand, direction: int, ledger: dict, H_excess: float,
                 V: float, T_K: float, T,
                 gsup: frozenset = frozenset(), iters: int = 60,
                 p_ext_kpa: float = P_EXT_KPA,
                 gas_escape: bool = True,
                 micro_rel: float | None = None) -> tuple[float, float]:
    """返回 (x*, x_max)。x* 为使 S 归零的程度；若全程为正则取计量上限。

    p_ext_kpa/gas_escape 透传给 S_of（泡点判据/闭口体系）。
    micro_rel：微步快通道（仅主循环 pick 位调用）——当二分区间上界
    hi ≤ max(X_MIN, micro_rel·x_max) 时提前止损：完整二分的返回值
    lo 恒 < hi ≤ 阈值，只进主循环的微步废弃分支（disable/idle），
    返回值本身不被消费，提前退出与跑满送代在 walk 语义上完全等价；
    探针数从 ~47 降到 ~14（log2(1/1e-4)），微步占六成以上的多平衡
    耦合体系（FeCl3+SCN 等）提速显著。调用方需保证 osc 外推块
    不会引用该候选的 ext（nk 不在近期 hist 键中）且 kind 非
    dissolve/precip（微溶盐微步会真实执行，需精确值）。"""
    rr = c.r if direction > 0 else c.pr
    pp = c.pr if direction > 0 else c.r
    limited = [(s, nu) for s, nu in rr.items() if s not in (WATER, H_ION)]
    if not limited:
        return 0.0, 0.0
    x_max = min(ledger.get(s, 0.0) / nu for s, nu in limited)
    if x_max <= X_MIN:
        return 0.0, x_max
    nu_H = pp.get(H_ION, 0) - rr.get(H_ION, 0)
    if c.kind == "redox":
        # 游离强酸/强碱为定量资源：redox 净耗 H+ 受 H_excess 限制（耗尽后由
        # 固相/水氧化路径另行处理）；净产 H+（耗 OH-）受碱储备限制。
        # 例外：无固相电对的还原剂（Na/Ca/K/Fe2+/I- 等）——水供质子，
        # 耗 H+ 升 pH 由 f(x) 的 S 归零自限，整个 nu_H<0 分支都不受游离酸
        # 硬约束（原仅在 He<=0 时豁免，He 为微小正残差时仍被误顶死）
        if nu_H < 0:
            if c.meta.get("acid_limited", True):
                if H_excess <= 0.0:
                    return 0.0, 0.0
                x_max = min(x_max, H_excess / (-nu_H))
        elif nu_H > 0 and H_excess < 0.0:
            # 净产 H+ ≡ 耗 OH-：仅当体系确实呈碱性（有真实 OH- 储备）才限量。
            # 酸性/近中性时 He 的微小负残差是记账幻影（speciation 把质子编入
            # HNO2 等弱酸形态，estimate_pH 仍报酸），真实 OH- 储备 ~1e-10、
            # 无实际限量对象——S 随 pH 下降自限。否则 NO2 歧化类产酸通道被
            # 幻影残差顶成每轮 -He/nu_H 的微步爬行（T78 曾 1500 轮 6.6s）
            _cls0 = closed_pH(ledger, H_excess, V, T, T_K)
            _ph0 = _cls0[0] if _cls0 is not None else estimate_pH(
                ledger, H_excess, V, T, T_K)
            if _ph0 > 9.0:
                x_max = min(x_max, -H_excess / nu_H)
        if x_max <= X_MIN:
            return 0.0, x_max

    # 预计算变化物种及其每单位 x 的净增量（rr 消耗为负、pp 生成为正）。
    # 原实现每次 f(x) 调用都 dict(ledger) 全拷贝并逐项 led2.get(s, 0.0)，
    # 是 solve_extent 的最大开销（318k 次调用 × O(N) 拷贝）。改为复用单一
    # 工作账本：只更新变化物种，其余保持 ledger 原值——estimate_state /
    # _buffer_titration 内部都自行 dict(ledger) 拷贝，不会污染工作账本。
    changing: list[tuple[str, float, float]] = []  # (species, net_delta_per_x, orig_value)
    _seen = set()
    for s, nu in rr.items():
        if s == WATER or s == H_ION:
            continue
        changing.append((s, -float(nu), ledger.get(s, 0.0)))
        _seen.add(s)
    for s, nu in pp.items():
        if s == WATER or s == H_ION:
            continue
        if s in _seen:
            # 同物种出现在 rr/pp 两侧（理论上 _try_vec 后不会，防御性合并）
            for i, (sp, d, o) in enumerate(changing):
                if sp == s:
                    changing[i] = (sp, d + float(nu), o)
                    break
        else:
            changing.append((s, float(nu), ledger.get(s, 0.0)))
    led_work = dict(ledger)
    # _buffer_titration 堆条目缓存（solve_extent 级生命周期）：二分内
    # led_work 仅 changing 物种的量变化，键序恒定；touch = changing 物种
    # 名集（堆条目对这些物种按当前量重建，其余整条目复用——平局 cnt
    # 恒为首次构建序，与无缓存路径 bit 级等价，见 speciation.py 文档）
    _bt_cache: dict = {}
    _touch = frozenset(s for s, _, _ in changing)
    # 非 redox 候选且不含 H+ 时 S 与 pH 无关（S_of 中 pH 仅用于 H+ 项；
    # estimate_pH 不改账本），二分全程跳过缓冲滴定（每次 f 调用省一次
    # 全套件最大热点）。redox 即使无 H+ 也必须走 estimate_state——
    # 其副产物 led_v 做强酸形态重排（NO3-/HNO3 等），影响 S。
    _need_ph = (c.kind == "redox") or (H_ION in c.r) or (H_ION in c.pr)
    # 闭式 pH 快路径（v0.4.2 迭代 A，bit 级等价）：无弱组分体系的 pH
    # 有闭式解（He 直读/纯水 pKw/2），滴定虚拟账本 = 原账本身份（两侧
    # 堆空）。判定在二分外一次完成：非 changing 物种的量在二分中不变
    # （floor 翻转不可能），changing 物种可能从无到有（生成型）——一律
    # 视为弱组分在场（保守方向：宁可走完整路径）。闭式时 f(x) 跳过全套
    # estimate_state（滴定+分支扫描——二分的最大热点），led_v 直接取
    # led_work 身份（与完整路径无滴定时的返回身份一致，_logc 缓存衔接）。
    _ph_closed = None
    if _need_ph:
        _cls = closed_pH(ledger, H_excess, V, T, T_K)
        if _cls is not None:
            _weak = weak_species_set(T)
            if any(s in _weak for s, _, _ in changing):
                _cls = None
        if _cls is not None:
            _pKw_c = pKw_of(T_K)
            _ph_water = -log10(10.0 ** (-_pKw_c / 2))
            _V_c = V

            def _ph_closed(he: float) -> float:
                h = he / _V_c
                if h >= 1e-3:
                    return max(-1.0, -log10(h))
                if h <= -1e-3:
                    return min(_pKw_c + 1.0, _pKw_c + log10(-h))
                return _ph_water
    # f 求值间的浓度对数缓存：led_work 仅 changing 物种随 x 变化，
    # 其余物种的 log10(c/V) 在整个二分期间不变——按物种记忆、逐次失效
    # changing 条目（值与逐次计算 bit 级一致）
    _logc: dict = {}
    # 求根审计专用（生产路径 `_atrace is None` ⟹ 三处 append 全跳过）：
    # `_atrace` = f 的 (x, pH, 固相在场) 轨迹；`_solid0/_solid_chg` 供标注
    # "固相在场"——用于判定 pH 跳变是否恰好伴随固相出现/消失。
    _atrace: list | None = None
    _solid0 = False
    _solid_chg: list = []
    if ROOT_AUDIT is not None:
        _solid0 = any(s in T.solids and m > X_MIN for s, m in ledger.items())
        _solid_chg = [(s, d, o) for s, d, o in changing if s in T.solids]

    def _sol(x: float) -> bool:
        return _solid0 or any(o + d * x > X_MIN for s, d, o in _solid_chg)

    def f(x: float) -> float:
        for s, d, orig in changing:
            led_work[s] = orig + d * x
            if s in _logc:
                del _logc[s]
        if not _need_ph:
            if _atrace is not None:
                _atrace.append((x, None, False, "与 pH 无关"))
            return direction * S_of(c, led_work, V, 7.0, T_K, T, gsup,
                                    p_ext_kpa, gas_escape, _logc)
        if c.kind == "redox":
            if _atrace is not None:
                _speciation.PH_TAGS = []
            if _ph_closed is not None:
                pH_x, led_v = _ph_closed(H_excess + nu_H * x), led_work
            else:
                pH_x, led_v, _ = estimate_state(led_work, H_excess + nu_H * x,
                                                V, T, T_K, _bt_cache, _touch)
            if _atrace is not None:
                _tags = _speciation.PH_TAGS
                _atrace.append((x, pH_x, _sol(x),
                                _tags[-1] if _tags else "闭式"))
                _speciation.PH_TAGS = None
            # led_v 是滴定后的虚拟账本：与 led_work 同一对象时（无滴定）
            # 浓度缓存仍有效；新生成的 dict 必须回退逐项计算
            return direction * S_of(c, led_v, V, pH_x, T_K, T, gsup,
                                    p_ext_kpa, gas_escape,
                                    _logc if led_v is led_work else None)
        if _atrace is not None:
            _speciation.PH_TAGS = []
        pH_x = (_ph_closed(H_excess + nu_H * x) if _ph_closed is not None
                else estimate_pH(led_work, H_excess + nu_H * x, V, T, T_K,
                                 _bt_cache, _touch))
        if _atrace is not None:
            _tags = _speciation.PH_TAGS
            _he = H_excess + nu_H * x
            _atrace.append((x, pH_x, _sol(x), _tags[-1] if _tags else "闭式",
                            _ledger_charge(led_work) + _he))
            _speciation.PH_TAGS = None
        return direction * S_of(c, led_work, V, pH_x, T_K, T, gsup,
                                p_ext_kpa, gas_escape, _logc)

    f_hi = f(x_max)
    if f_hi > 0:
        return x_max, x_max
    if iters > 20:
        # 主求解（参与路径选择）：纯二分，迭代次数固定。
        # 鞍点/多根体系（NaClO+CO2、Cu+HNO3 的 NO2→NO 脱气伪解通道）中
        # f 非单调、存在多个过零点，**根选定**直接决定 walk 路径——即
        # "执行多少"就是化学。旧的约束写法是"被测试套件锁定"（§7 V），
        # 这条已按 §7 W 纠正：断言锁化学不锁轨迹，**真约束是根选定**
        # ——换更快的求根方法时，必须保证它选出同一个根（`tools/roots.py`
        # 的括号审计：网格符号翻转 + 两法根差）。
        # 微步快通道例外（见 micro_rel 文档）：lo 恒 < hi ≤ 阈值，与跑满
        # 送代在主循环的分支决策上 bit 级等价（非单调 f 同样成立：
        # 后续送代 lo=mid<hi ≤ 阈值的单调推理不依赖 f 的形态）
        _micro_thr = (max(X_MIN, micro_rel * x_max)
                      if micro_rel is not None else None)
        # 收敛容差（v0.5.0 性能，§7 U-3/W）：二分原本跑到**浮点饱和**——
        # 源码注释已记"60 次二分超出 double 分辨率 ~52 bit，剩余迭代 mid 与
        # 端点重合、f 重复求值同值，纯空转"。而走步对 ext 的实际需求远低于
        # double 分辨率：执行量最终以 1e-6 mol 呈现、判据阈值是
        # ANN_MIN_EXTENT(1e-3) 与微步线(0.02·x_max)。取 `1e-11·max(1,x_max)`
        # 绝对容差——比任何消费 ext 的判据小 5 个数量级以上。
        # 实测（SOF_CALLS 确定性计数，§7 W-4）：S_of 调用 −19%；
        # 代价是部分用例的走步轨迹改变（多为收敛变好：Z31 resid 0.02→0.0、
        # N23 8.28→0、Y03 10.2→1.8），断言侧锁轨迹的 4 条已按 §7 W 改为
        # 锁化学（张成判据 + 区间断言）⟹ 全量 1294 断言不降。
        _tol = _EXTENT_TOL_REL * max(1.0, x_max)
        _audit = None
        if ROOT_AUDIT is not None:
            _i = ROOT_AUDIT["n"]
            ROOT_AUDIT["n"] += 1
            if _i % ROOT_AUDIT["every"] == 0:
                _audit = [0]                      # 本括号的 f 求值计数
        lo, hi = 0.0, x_max
        for _ in range(iters):
            if hi - lo <= _tol:
                break   # 区间已到需求精度：返回的 lo 与真根差 ≤ _tol
            mid = (lo + hi) / 2
            if mid <= lo or mid >= hi:
                break   # 浮点饱和：区间已不可再分（数学上与跑完全部迭代等价）
            if _audit is not None:
                _audit[0] += 1
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
            if _micro_thr is not None and hi <= _micro_thr:
                break   # 微步快通道：根的上界已坍缩到主循环微步阈值之下
        if _audit is not None:
            _x_ill = _audit_bracket(f, c, direction, x_max, lo, _audit[0], f_hi,
                                    _micro_thr is not None)
            if _x_ill is not None:
                # pH(x) 细扫轨迹：判定"口袋是不是 pH 造成的"（§7 X-3 的假设）
                # 并普查 pH 的连续性（近似机器分支切换会让 pH 跳变）
                _atrace = []
                _hi_s = 1.3 * max(lo, _x_ill)
                for _k in range(61):
                    f(_hi_s * _k / 60)
                ROOT_AUDIT["rec"][-1]["trace"] = _atrace
                _atrace = None
        return lo, x_max
    # 粗精度求解（iters<=20，仅用于慢标注等布尔阈值判定，不影响路径）：
    # 5 次二分定盆 + Brent 抛光，~12 次求值达到足够精度
    lo, hi = 0.0, x_max
    for _ in range(5):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    xtol = 1e-6 * max(1.0, x_max)
    a, b = lo, hi
    fa, fb = f(a), f(b)
    if fa <= 0:
        return a, x_max
    cc, fcc = b, fb
    d = e = b - a
    for _ in range(iters):
        if (fb > 0) == (fcc > 0):
            cc, fcc = a, fa
            d = e = b - a
        if abs(fcc) < abs(fb):
            a, b, cc = b, cc, b
            fa, fb, fcc = fb, fcc, fb
        tol1 = 2.0 * xtol
        xm = 0.5 * (cc - b)
        if abs(xm) <= tol1 or fb == 0.0:
            break
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            bs = fb / fa
            if a == cc:
                bp = 2.0 * xm * bs
                bq = 1.0 - bs
            else:
                bq = fa / fcc
                br = fb / fcc
                bp = bs * (2.0 * xm * bq * (bq - br) - (b - a) * (br - 1.0))
                bq = (bq - 1.0) * (br - 1.0) * (bs - 1.0)
            if bp > 0:
                bq = -bq
            else:
                bp = -bp
            if 2.0 * bp < min(3.0 * xm * bq - abs(tol1 * bq), abs(e * bq)):
                e, d = d, bp / bq
            else:
                d = e = xm
        else:
            d = e = xm
        a, fa = b, fb
        b = b + d if abs(d) > tol1 else b + (tol1 if xm > 0 else -tol1)
        fb = f(b)
    return b, x_max


def _film_acid_soluble(fp: str, aOH: float, T, T_K: float = 298.15) -> bool:
    """氢氧化物膜的酸可溶性：其阳离子在当前 pH 的饱和活度 >1（可溶 >1 mol/L）
    → 膜无法积累。非氢氧化物膜（PbSO4 等）返回 False（仍按膜封锁判定）。"""
    cell = T.ksp_by_solid.get(fp)
    if cell is None or cell["pair"][1] != "OH^-":
        return False
    return _ksp_asat(cell, aOH, T_K) > 1.0


# ========================================================== 动力学层：blocked

def blocked_extent(pick: Cand, direction: int, evals: list, T_K: float, T,
                   ledger: dict, H_excess: float, V: float, x_max: float,
                   pH: float = 7.0, p_ext_kpa: float = P_EXT_KPA,
                   gas_escape: bool = True) -> tuple:
    # 返回 (extent|None, 膜电位 species 字典)；未受阻为 (None, None)
    rr = pick.r if direction > 0 else pick.pr
    pr = pick.pr if direction > 0 else pick.r
    if pick.kind != "redox":
        return None, None
    film_ps = [p for p in pr if T.ksp_by_solid.get(p, {}).get("film")]
    if not film_ps or not any(s in T.solids for s in rr):
        return None, None
    # 酸可溶的氢氧化物膜豁免封锁：膜阳离子在当前 pH 的饱和活度 >1 时
    # 膜生成即溶、无法积累（Zn(OH)2 在 8M HNO3 中——溶解候选只在膜在账时
    # 才枚举（§4.5 需求集），而决策点膜恒被上步清 0，通道检查结构性失效，
    # 故直接以饱和活度判定）。PbSO4/CaSO4 等非氢氧化物膜、中性水的
    # Mg(OH)2/Al(OH)3 膜不受影响
    # 豁免仅在酸性介质成立（pH<6 有酸库持续溶膜）；中性/碱中氢氧化物膜
    # 溶解会使界面 pH 自缓冲升高（Mg+冷水：Mg(OH)2 饱和界面 pH≈10.5，
    # a_sat 跌回 <1），膜依然封锁
    if pH < 6.0:
        aOH = 10.0 ** (pH - pKw_of(T_K))
        film_ps = [p for p in film_ps if not _film_acid_soluble(p, aOH, T, T_K)]
        if not film_ps:
            return None, None
    film_ps = [p for p in film_ps if T_K < T.ksp_by_solid[p].get("unlock_T", 1e9)]
    if not film_ps:
        return None, None
    # 溶解通道能溶走的膜量 ≥ 封锁允许程度 → 膜可穿，不封锁；
    # 否则通道只溶解微量（ACT_FLOOR 播种的假阳性），仍封锁
    allow = BLOCKED_EXTENT * max(x_max, 1e-12)
    for c, d, S in evals:
        reactants = c.r if d > 0 else c.pr
        if S > 0 and any(fp in reactants for fp in film_ps):
            ext_c, _ = solve_extent(c, d, ledger, H_excess, V, T_K, T,
                                    p_ext_kpa=p_ext_kpa,
                                    gas_escape=gas_escape)
            if ext_c >= allow:
                return None, None
    return BLOCKED_EXTENT, film_ps


# ========================================================== OVERRIDE（逃生舱）

def _swap_conserves(old: dict, new: dict) -> bool:
    """账本整体替换的守恒闸门：非 H/O 元素的量不得改变。

    H/O 走水与质子账本（`H_excess`），不纳入；其余元素必须守恒。

    **尺度只取溶质**：首版把溶剂水（55.6 mol）算进 scale，容差被抬到
    0.056，于是 0.0242 mol 的碳缺口照样放行——"把溶剂混进溶质尺度的
    判据"是这里最容易犯的错（净方程契约里也踩过同类）。判据尺度必须是
    **被约束对象的量级**。容差取 max(1e-6, 1e-3 × 溶质尺度)，1e-3 是
    "呈现上可忽略"的分辨率，与净方程契约同口径。
    """
    els: set = set()
    for sp in list(old) + list(new):
        if sp != WATER and not sp.startswith("__"):
            els |= set(elements_of(sp))
    els -= {"H", "O"}
    if not els:
        return True
    scale = max([abs(v) for s, v in list(old.items()) + list(new.items())
                 if s != WATER and not s.startswith("__")] or [1.0])
    tol = max(1e-6, 1e-3 * scale)
    for el in els:
        a = sum(elements_of(s).get(el, 0) * v for s, v in old.items()
                if not s.startswith("__"))
        b = sum(elements_of(s).get(el, 0) * v for s, v in new.items()
                if not s.startswith("__"))
        if abs(a - b) > tol:
            return False
    return True


def _virt_redox_gain(ledger: dict, vled: dict, T) -> bool:
    """滴定后的虚拟账本是否出现了账本上不存在（或已耗尽）的、且能真正配成
    新氧化还原候选的物种。仅"出现且可配对"方向触发惰性实现：
    - "出现"但无配对对象不触发（如 NH3+HCl 中和产生的 NH4+ 虽是 NO3-/NH4+
      电对的还原侧，但账上无 E0>0.88V 的氧化剂，实现滴定只会白白吞掉
      中和步骤的诚实记账）；
    - "消失"方向不触发（消失的形态是酸碱候选自己的记账职责）。
    配对用物种级最高/最低 E0 近似（存在性判定，假阳性无害——枚举层会再过滤）。"""
    for sp, m in vled.items():
        if m <= X_MIN or ledger.get(sp, 0.0) > X_MIN:
            continue
        e_ox = T.redox_ox_E.get(sp)    # 作为氧化剂（被还原）需账上有更弱氧化剂电对的还原剂
        if e_ox is not None:
            for sp2, m2 in vled.items():
                if m2 > X_MIN and sp2 != sp:
                    e_red = T.redox_red_E.get(sp2)
                    if e_red is not None and e_red < e_ox:
                        return True
        e_red = T.redox_red_E.get(sp)  # 作为还原剂（被氧化）需账上有更强氧化剂
        if e_red is not None:
            for sp2, m2 in vled.items():
                if m2 > X_MIN and sp2 != sp:
                    e_o2 = T.redox_ox_E.get(sp2)
                    if e_o2 is not None and e_o2 > e_red:
                        return True
    return False


def _match_override(substances: list[dict], V: float, T_K: float, T):
    names = {i["name"] for i in substances}
    # 按匹配物种数降序：最具体的规则优先（如 KO2+CO2 先于 KO2+水）
    for o in sorted(T.overrides, key=lambda x: -len(x["match"].get("species", []))):
        m = o["match"]
        if not set(m.get("species", [])) <= names:
            continue
        if m.get("conc") == "concentrated":
            if not any(i["name"] in T.conc and i["mol"] / V >= T.conc[i["name"]] for i in substances):
                continue
        if "T_min" in m and T_K < m["T_min"]:
            continue
        if "T_max" in m and T_K > m["T_max"]:
            continue
        return o
    return None


# ========================================================== 主循环

_COND_KEYS = {"V_L", "T_K", "T_C", "c_H", "c_OH", "pH", "p_kpa",
               "isothermal", "kinetics", "gas_escape"}
# isothermal/kinetics/gas_escape 为模式开关（bool）：isothermal 仅透传给
# 结果的 cond 元数据（由独立温度模块 thermo.py 事后消费，不影响求解）；
# kinetics 控制动力学层（slow/gate/膜封锁等标记是否生效）；
# gas_escape 控制自产气体逸出（False=闭口体系，气体保留在溶液账本）。


def judge(substances: list[dict], conditions: dict | None, T: Tables,
         _probe: dict | None = None) -> dict:
    conditions = dict(conditions or {})
    # 质子化族根表惰性初始化（v0.4.2 触发点③：社区质子交换族数计数）
    if not _PROTON_ROOTS:
        _set_proton_roots(T)
    bad = set(conditions) - _COND_KEYS
    if bad:
        # 条件键静默忽略是定义错位的温床（如 T_C 被当 298K 跑），宁可报错
        raise ValueError(f"未知条件键 {sorted(bad)}；支持 {sorted(_COND_KEYS)}")
    if "T_C" in conditions and "T_K" not in conditions:
        conditions["T_K"] = conditions["T_C"] + 273.15
    cond = {"V_L": 1.0, "T_K": 298.15, "c_H": None, "c_OH": None, "pH": None,
            "p_kpa": P_EXT_KPA}
    cond.update(conditions)
    V, T_K = cond["V_L"], cond["T_K"]
    p_ext_kpa = float(cond["p_kpa"])
    # 模式开关（经 conditions 传入；默认值与模块级语义一致）
    kinetics = bool(cond.get("kinetics", True))
    gas_escape = bool(cond.get("gas_escape", True))
    if not 273.15 <= T_K <= 373.15:
        # 常压液态水温度域：域外水的存在形式/活度约定全部失效（van't Hoff
        # ΔCp≈0 近似与 Hill 形态曲线也仅在此域内标定），宁可报错不静默
        raise ValueError(f"T_K={T_K} 超出液态水温度域 273.15–373.15 K")
    if p_ext_kpa <= 0.0:
        raise ValueError(f"p_kpa={p_ext_kpa} 必须为正（kPa）")

    ov = _match_override(substances, V, T_K, T)
    if ov is not None:
        return _apply_override(substances, cond, V, T_K, T, ov)

    origins: dict[str, set] = {}
    ledger, H_excess, steps, unknown = normalize(substances, cond, T, origins)
    H_excess0 = H_excess   # 初始质子账本（净离子方程式的 H+/OH- 净差来源）
    initial = dict(ledger)
    # 初始投料中的气体 = 持续供给（按账本活度）；反应自产气体按 A_GAS 逸出
    gsup = frozenset(sp for sp in initial if sp in T.gases and initial[sp] > X_MIN)
    escaped: dict[str, float] = {}   # 逸出气相账户（泡点扫气，逸出即消失）
    chem_net: dict = {}   # key -> [净程度, kind, 投料来源集]（reacted 判据）
    annotations: list[str] = []
    blocked_solids: dict[str, list] = {}   # 被膜封锁的金属 -> 膜固相列表
    disabled: dict = {}   # (key, d) -> 禁用时的账本签名；状态实质漂移后自动解禁
    slow_seen = False
    _exit_reason = "max-iter"   # 探针专用：主循环退出原因（只读诊断）
    _it_total = 0               # 探针专用：累计主循环迭代数（含 sweep 轮）

    def _sig():
        # 账本签名：用于 disabled 过期判定。纯震荡会回到原签名（禁用保持），
        # 状态漂移（其他通道推进）则签名改变 → 解禁，避免误杀平衡已移动的通道
        return (tuple(sorted((s, round(m, 5)) for s, m in ledger.items()
                             if s != WATER and m > X_MIN)), round(H_excess, 5))

    def _refresh_disabled():
        sig = _sig()
        for k, sig0 in list(disabled.items()):
            if sig != sig0:
                del disabled[k]
    # 净反应震荡冻结：同一净反应（忽略 H2O/H+）正反向净零空转经三种
    # 检测器（极限环签名复现 / 近窗程度对消 / 短周期重复）判定后永久冻结——
    # 阻断浓酸体系中同一反应因 H+ 挂侧不同生成多个配平形式互相逆转净零
    hist: list = []   # 已执行净反应键序列（用于循环震荡检测）
    frozen_perm: set = set()   # 实测/极限环/短周期震荡判定后永久冻结的净反应键
    # 慢标注采样节奏记忆：(key, d) -> 上次阈值下评的 hist 位点。
    # 原实现 slow_seen 恒 False 的体系（慢通道永不可达显著量）每迭代对
    # 全部慢候选重解（Fe31：104 迭代×12 候选=1209 次 solve_extent，占
    # 总求值 2/3）——低于阈值的负结论按“≥32 迭代或步量 ≥0.01 的实质步
    # 发生后”重评，收敛终局前统一复核一次，保证终态标注不漏报
    slow_ann: dict = {}
    _last_mat_step = -1
    # 爬行收敛加速（窗口几何外推）状态：每执行步的 (账本快照, He)；
    # 40 步为一窗，窗间 L1 范数等比衰减（非精确周期的不规则爬行——
    # 精确短周期已由上方 [osc] 外推覆盖）时，几何尾部 Σρ^k·D ≈ D·ρ/(1−ρ)
    # 一次补齐；外推置零负值钳制，失准由不动点两侧 S 变号自纠回。
    # D 是平衡步计量的线性组合，任意缩放仍原子/电荷守恒（钳制物种除外）
    snaps: list = []
    last_jump = -1
    _CRAWL_W = 40
    # pH 悬崖乒乓检测状态（v0.4.0 第五检测器）：每执行步的 H⁺ 参与标记
    # （走步级兜底 pH 机器的端点不连续，He 跨 1e-3 时 pH 3.0↔6 跳变——
    # 跨键喂食循环里各步的 S 评估用了不同相位快照，各自为正而联合为负。
    # O02：Zn(OH)₂ 沉淀放 H⁺ ↔ Fe 溶解耗 H⁺，净循环 logK=−11.3 不自发
    # 却乒乓 300+ 步把 Fe 溶掉 0.34）。pH 序列用 snaps 懒评估（drain
    # 循环内迭代 pH 不刷新——逐快照重估才可见摆动；每 8 步去抖）。
    hexec: list[bool] = []
    # 走步画像诊断计数（v0.4.2 探针扩展，纯诊断零行为影响——不进
    # digest，仅供 converg/慢例分析）：联立尝试/冻结事件/实质微步数；
    # CHEM_TRACE_WINDOWS=1 时每 32 步窗口导出 (drift, turnover) 标定数据
    _diag = {"joint_tries": 0, "joint_ok": 0, "freeze_events": 0,
             "micro_steps": 0, "windows": []}
    _WINDOWS = bool(_os.environ.get("CHEM_TRACE_WINDOWS"))

    def _netkey(c, direction):
        # 净反应键（忽略 H2O/H+ 的 sorted 物种对）——只依赖候选本身，
        # 按 (Cand, 方向) 缓存（原每次调用重排 64k 次/慢例）
        return c.netkey_fwd if direction > 0 else c.netkey_rev

    # ---- 步执行记账（v0.3.8 从主循环体内提出，walk 步与联立跳步共用） ----
    # 账本落实 + steps/chem_net/origins/hist/snaps 全套记账。纯剪切
    # （语句序与数值路径不变，bit 级等价）；H_excess/_last_mat_step 为
    # nonlocal 绑定。
    def _exec(pick, d, ext, x_max, S, nk):
        nonlocal H_excess, _last_mat_step
        _diag["micro_steps"] += 1
        rr = pick.r if d > 0 else pick.pr
        pp = pick.pr if d > 0 else pick.r
        for s, nu in rr.items():
            if s == H_ION:
                H_excess -= nu * ext
            elif s != WATER:
                ledger[s] = ledger.get(s, 0.0) - nu * ext
        for s, nu in pp.items():
            if s == H_ION:
                H_excess += nu * ext
            elif s != WATER:
                ledger[s] = ledger.get(s, 0.0) + nu * ext
        eq = " + ".join(f"{_fmt(nu)}{s}" for s, nu in rr.items() if s != WATER)
        eq += " -> " + " + ".join(f"{_fmt(nu)}{s}" for s, nu in pp.items() if s != WATER)
        lims = [initial.get(s, 0.0) / nu for s, nu in rr.items()
                if s not in (WATER, H_ION) and initial.get(s, 0.0) > 0]
        lim0 = min(lims) if lims else 0.0
        abs_conv = min(1.0, ext / lim0) if lim0 > 0 else (1.0 if ext > 0 else 0.0)
        steps.append({"kind": pick.kind, "equation": eq,
                      "logK": round(logK_T(pick, T_K), 2),
                      "S": round(S, 2), "extent": round(ext, 6),
                      "conversion": round(min(1.0, ext / x_max), 4) if x_max > 0 else 1.0,
                      "abs_conv": round(abs_conv, 4)})
        # 来源追踪：产物继承反应物的投料来源集（中间体由此携带血统，
        # 供 reacted（狭义化学反应）的"跨投料相互作用"判据使用）；并按候选
        # 累计净程度（正反向对消），判据在收敛后统一应用
        _fs: set = set()
        for s in rr:
            if s != WATER:
                _fs |= origins.get(s, set())
        # H+ 参与（任一侧）时并入质子账本来源：沉淀/溶解表观只写 Fe3+ -> Fe(OH)3 + 3H+，
        # 驱动它的 OH- 来自另一投料（FeCl3+NaOH），跨投料相互作用借 H+ 通道显现。
        # 动态传播：凡触碰质子账本的投料都登记进 origins[H_ION]——固相碱
        # （石灰乳 Ca(OH)2）在规范化时保持固相不登记，其溶解放酸驱动他种沉淀
        # （海水提镁 Mg2+ + Ca(OH)2）也属跨投料相互作用
        if H_ION in rr or H_ION in pp:
            origins.setdefault(H_ION, set()).update(_fs)
            _fs |= origins[H_ION]
        for s in pp:
            if s not in (WATER, H_ION):
                origins.setdefault(s, set()).update(_fs)
        _cn = chem_net.get(pick.key)
        if _cn is None:
            chem_net[pick.key] = [0.0, pick.kind, _fs]
        else:
            _cn[2] |= _fs
        chem_net[pick.key][0] += d * ext
        hist.append((nk, ext))
        hexec.append(H_ION in rr or H_ION in pp)
        if ext >= 0.01:
            _last_mat_step = len(hist)   # 实质步标记（慢标注负结论的失效钩子）
        snaps.append((dict(ledger), H_excess))

    # ---- 联立求解加速（v0.3.8 深水区；详见 joint.py 模块头） ----
    # 触发点①（爬行检测）：it≥64 且近 24 步全微步且 ≥2 个不同规范净键
    # → 只联立解实际在循环的平衡（Gauss-Seidel 循环坐标的 Newton 加速；
    # 不全收——候选集存在跨数据源 Hess 互斥副本，全收无解）。成功 →
    # 走步同款记账落实 + disabled.clear()，walk 全量重评估自校验；
    # "boundary"（联立不动点在物理域外——数据张力型）→ 提前冻结循环键
    # （Ag32 型从 ~1786 步提前到 ~101 步）；失败 = 严格无操作 + 黑名单。
    # （触发点② idle 退出前精修已试已回退——全量差分否决：idle 点上走步
    # 仲裁已完成，活动子集的联立不动点 ≠ 仲裁点，跳步=语义翻案，
    # Co32/Co33/T34/Fe33 四例翻车；J06 型欠收敛真因是 pH 滴定端点悬崖。）
    _joint_crawl_last = -32
    _joint_idle_sig = None
    _joint_fires = 0
    _joint_blacklist: set = set()   # 失败循环签名（canonical keys 排序元组）

    def _joint_collect(cycle_keys=None):
        # 活性平衡集：两侧在场、未冻结、非 slow/deferred、未被膜封锁。
        # 全收（含近平衡者——它们是跳步的耦合约束方程：排除会让跳步把它们
        # 扰离平衡、走步反向执行拆掉跳步，J06 曾三跳三拆零净效果）；驱动
        # 方向由 S 符号定，|S| 排序截前 JOINT_MAX_M 个。返回 (Cand, 方向, S)
        # 三元组。cycle_keys 非空时仅收循环键匹配者（触发点①口径）
        seen: dict = {}
        for c in cands:
            if kinetics and (c.meta.get("slow") or c.meta.get("deferred")):
                continue
            if (c.kind == "redox" and c.meta.get("ox_couple") == H_ION
                    and blocked_solids
                    and any(s in blocked_solids
                            for s in list(c.r) + list(c.pr)
                            if s not in (WATER, H_ION))):
                continue
            nk_c = (c.netkey_fwd if c.netkey_fwd <= c.netkey_rev
                    else c.netkey_rev)   # 规范净键：fwd/rev 归一
            if nk_c in frozen_perm:   # 冻结集含双方向，规范键必在其一
                continue
            if cycle_keys is not None and nk_c not in cycle_keys:
                continue
            ps = c.pres_specs
            if not all(ledger.get(s, 0.0) > X_MIN for s in ps[0] + ps[1]):
                continue
            S_f = S_of(c, ledger, V, pH, T_K, T, gsup, p_ext_kpa, gas_escape)
            prev = seen.get(nk_c)
            if prev is None or abs(S_f) > abs(prev[1]):
                seen[nk_c] = (c, S_f)
        out = [(c, 1 if S > 0 else -1, S) for c, S in seen.values()]
        out.sort(key=lambda e: -abs(e[2]))
        return out[:JOINT_MAX_M]

    def _joint_fire(cycle_keys=None, freeze_on_boundary=False) -> bool:
        """返回 True = 状态实质移动（"ok" 跳步落实）。boundary 时可选冻结
        循环键（提前宣告平衡——与走步晚期 freeze 同语义）；fail 严格无操作。"""
        nonlocal _joint_fires
        actives = _joint_collect(cycle_keys)
        if len(actives) < JOINT_MIN_M:
            return False
        if max(abs(S) for _c, _d, S in actives) < 0.02:
            return False   # 无驱动者（全近平衡）：无事可做
        status, x, _res = joint_solve(ledger, H_excess,
                                      [(c, d) for c, d, _S in actives],
                                      V, T_K, T, gsup, p_ext_kpa, gas_escape,
                                      S_of, H_ION, WATER)
        _diag["joint_tries"] += 1
        if status == "boundary":
            if not freeze_on_boundary or cycle_keys is None:
                return False
            # 联立不动点在物理域外（数据张力型爬行）：提前冻结循环键
            # （Ag32 型的正确仲裁从 ~1786 步提前到 ~80 步）
            for k in cycle_keys:
                frozen_perm.add(k)
                frozen_perm.add((k[1], k[0]))
            _diag["freeze_events"] += 1
            if _TRACE:
                print(f'  [joint-boundary-freeze] {len(cycle_keys)} keys '
                      f'resid={_res:.3f}')
            return False
        if status != "ok":
            return False
        _joint_fires += 1
        _diag["joint_ok"] += 1
        if _TRACE:
            print(f'  [joint] m={len(actives)} resid={_res:.3f} '
                  f'x=[{" ".join(f"{xj:+.4g}" for xj in x)}]')
        # 逐反应落实（走步同款记账）；x_j<0 = 净反向，翻方向记正向步；
        # x_max/步 S 均按跳步前账本计（同一批的联合语义）
        for (c, dj, S_f), xj in zip(actives, x):
            if abs(xj) < X_MIN:
                continue
            dd = dj if xj > 0 else -dj
            rr_j = c.r if dd > 0 else c.pr
            x_max_j = min((ledger.get(s, 0.0) / nu for s, nu in rr_j.items()
                           if s not in (WATER, H_ION)), default=0.0)
            _exec(c, dd, abs(xj), x_max_j, S_f if xj > 0 else -S_f,
                  _netkey(c, dd))
        disabled.clear()   # 状态实质移动：全量解禁让 S 重验（既有自校正）
        return True

    def _joint_cycle_keys():
        # hist 尾部净键的规范集（fwd/rev 归一）：实际在循环的平衡
        out = set()
        for k, _e in hist[-24:]:
            out.add(k if k <= (k[1], k[0]) else (k[1], k[0]))
        return out

    # ---- 触发点③（v0.4.2 迭代 D）：社区级 pH 一致化联立 ----
    # 微步窗同①，但社区口径（_joint_collect(None) 全收活性平衡，含近
    # 平衡者——跳步的耦合约束方程）+ pH 提升为联立变量（_solve_ph 的
    # 自洽闭合行——快照 pH 机器锚定联立 pH）。张力族（H46 Ag₂O 极限
    # 环 / N34 草酸阶梯）的"一步数学"。族闸：社区质子交换族数 ≥2 时
    # pH 闭包只是启发式（E35 型假不动点一枪跳 0.589 mol 凌驾仲裁的
    # 翻案教训）→ 硬失败，绝不 fall-through 到 legacy joint_solve
    # （对同样行做同样的坏跳）。boundary 冻结循环键（Ag32 仲裁提前）。
    # 独立冷却 32 迭代 + 失败签名黑名单（与①互不干扰）。
    _joint_ph_last = -32
    _joint_ph_blacklist: set = set()

    def _joint_families(actives) -> set:
        """社区质子交换族数（族闸）：H⁺ 参与的平衡经其 pKa 共轭酸碱族
        （真弱酸碱形态——H₂S/HS⁻/S²⁻、NH₃/NH₄⁺ 等）计数。金属/沉淀/
        配合物物种不计（非质子交换族——H46 的 Ag⁺/Ag₂O 挂侧不是
        pKa 共轭语义；E35 的判别子是 H₂S+NH₃ 两个真弱酸碱族）。"""
        fam: set = set()
        for c, _d, _S in actives:
            if H_ION in c.r or H_ION in c.pr:
                for s in list(c.r) + list(c.pr):
                    if s in (H_ION, WATER):
                        continue
                    root = _PROTON_ROOTS.get(s)
                    if root is not None:
                        fam.add(root)
        return fam

    def _joint_fire_ph() -> bool:
        """触发点③：社区级 pH 一致化联立。返回 True = 状态实质移动。"""
        nonlocal _joint_fires
        actives = _joint_collect(None)   # 社区口径：全活性（含近平衡约束）
        if len(actives) < JOINT_MIN_M:
            return False
        if max(abs(S) for _c, _d, S in actives) < 0.02:
            return False   # 无驱动者（全近平衡）：无事可做
        fam = _joint_families(actives)
        if len(fam) >= 2:
            # 多族 pH 闭包 = 启发式（E35 翻案教训）：硬失败——不跳、
            # 不 fall-through 到 legacy（同样的坏跳）。黑名单 + 冷却
            # 由调用方处理（这里只报闸门拦截）
            _diag["joint_tries"] += 1
            if _TRACE:
                print(f'  [joint-ph-gate] {len(fam)} proton families '
                      f'blocked (community m={len(actives)})')
            return False
        status, x, _res = _solve_ph(ledger, H_excess,
                                    [(c, d) for c, d, _S in actives],
                                    V, T_K, T, gsup, p_ext_kpa, gas_escape,
                                    S_of, H_ION, WATER)
        _diag["joint_tries"] += 1
        if _TRACE and status != "ok":
            print(f'  [joint-ph-{status}] m={len(actives)} resid={_res:.3f} '
                  f'families={sorted(_joint_families(actives))}')
        if status == "boundary":
            _ck = _joint_cycle_keys()
            for k in _ck:
                frozen_perm.add(k)
                frozen_perm.add((k[1], k[0]))
            _diag["freeze_events"] += 1
            if _TRACE:
                print(f'  [joint-ph-boundary] {len(_ck)} keys '
                      f'resid={_res:.3f}')
            return False
        if status != "ok":
            return False
        _joint_fires += 1
        _diag["joint_ok"] += 1
        if _TRACE:
            print(f'  [joint-ph] m={len(actives)} resid={_res:.3f} '
                  f'x=[{" ".join(f"{xj:+.4g}" for xj in x)}]')
        for (c, dj, S_f), xj in zip(actives, x):
            if abs(xj) < X_MIN:
                continue
            dd = dj if xj > 0 else -dj
            rr_j = c.r if dd > 0 else c.pr
            x_max_j = min((ledger.get(s, 0.0) / nu for s, nu in rr_j.items()
                           if s not in (WATER, H_ION)), default=0.0)
            _exec(c, dd, abs(xj), x_max_j, S_f if xj > 0 else -S_f,
                  _netkey(c, dd))
        disabled.clear()   # 状态实质移动：全量解禁让 S 重验（既有自校正）
        return True

    for _sweep_round in range(20):
      seen_sig: dict = {}
      idle = 0   # 连续零执行迭代计数：签名不变 ⇒ disabled/frozen 永不刷新，
               # 空转只会无限重复同一评估，达阈值即宣告收敛
      enum_memo: dict = {}   # 枚举结果缓存（present/pH 桶不变时跳过 30k 模板重扫）
      for it in range(MAX_ITER):
        _it_total += 1
        _refresh_disabled()
        H_excess = _respeciate_strong_acids(ledger, H_excess, V, T)
        # 闭式 pH 快路径（v0.4.2 迭代 A）：无弱组分时 pH/vled/He_v 三合一
        # 闭式直出（vled = ledger 原身份 → 惰性实现检查天然跳过，与完整
        # 路径无滴定时的身份语义一致）
        _cls = closed_pH(ledger, H_excess, V, T, T_K)
        if _cls is not None:
            pH, vled, He_v = _cls[0], ledger, _cls[2]
        else:
            pH = estimate_pH(ledger, H_excess, V, T, T_K)
            vled, He_v = _full_speciation(ledger, H_excess, V, T, T_K)
        if vled is not ledger and _virt_redox_gain(ledger, vled, T):
            # 惰性实现酸碱平衡：质子转移远快于氧化还原，强酸/强碱下的自由形态
            # （如 VO3-→VO2+）应直接参与氧化还原竞争——否则 H+/Zn 会在 V(V)
            # 未现身时抢跑耗尽 Zn。仅在"多级形态分布出现了账本上不存在的、
            # 且能配成新电对的氧化还原物种"时才落实为真实账本：全局落实会让
            # 酸碱候选失去诚实的 He 累积（Al3+ 水解酸性被吸回、NH4Ac 解离
            # 幻影循环、半中和 Henderson 失效），故保持惰性。pH 重算保持一致。
            #
            # **守恒闸门（v0.5.0）**：换账本 = 用虚拟形态分布整体替换真实
            # 账本，必须**元素守恒**才许换。`_buffer_titration` 的多级堆在
            # 某些体系会漏掉未滴定的族尾（实测 FeCl3+Na2CO3：虚拟账本把
            # 0.024234 mol CO3^2- 整个丢掉，碳 3.0 → 2.975767），于是净差
            # 向量本身不平、美化器永远找不到干净整数式（E42/K02 的
            # `375.187CO3^2- + …` 即此）。非 H/O 元素缺失即拒绝替换——
            # 留在真实账本上是诚实的质量态，虚拟账本只是形态分布优化。
            if _swap_conserves(ledger, vled):
                if _TRACE: print('  [realize]', [sp for sp, m in vled.items()
                                  if m > X_MIN and ledger.get(sp, 0.0) <= X_MIN])
                ledger, H_excess = vled, He_v
                pH = estimate_pH(ledger, H_excess, V, T, T_K)
            elif _TRACE:
                print('  [realize-rejected] 虚拟账本不守恒，保留真实账本')
        cands = enumerate_candidates(ledger, H_excess, pH, V, T_K, T, kinetics,
                                     memo=enum_memo)
        evals = []
        deferred_evals = []   # 让位档（晶格氧化等慢氧化通道）：无快候选时才出手
        slow_now = False
        # 候选评估循环中账本/pH 固定不变：浓度对数按物种记忆，跨全部
        # 候选复用（与逐次计算 bit 级一致；步后由下轮重建）
        logc: dict = {}
        for c in cands:
            if (c.key, 1) in disabled and (c.key, -1) in disabled:
                continue
            # 膜封锁只抑制溶剂（H+/H2O）氧化通道：致密膜（Mg(OH)2、Zn(OH)2、
            # 钝化 Al2O3…）阻断的是金属-水界面腐蚀；溶液中更强的氧化剂
            # （Fe2+、Cu2+…）经膜缺陷/置换路径仍可反应（镀锌层牺牲、置换沉积
            # 均为实验事实），不因自腐蚀成膜而冻结
            if c.kind == "redox" and c.meta.get("ox_couple") == H_ION and any(
                    s in blocked_solids for s in list(c.r) + list(c.pr)
                    if s not in (WATER, H_ION)):
                continue
            _ps = c.pres_specs
            pres_r = all(ledger.get(s, 0.0) > X_MIN for s in _ps[0])
            pres_p = all(ledger.get(s, 0.0) > X_MIN for s in _ps[1])
            if not pres_r and not pres_p:
                continue
            S_fwd = S_of(c, ledger, V, pH, T_K, T, gsup, p_ext_kpa,
                         gas_escape, logc)
            if pres_r and S_fwd > 0 and (c.key, 1) not in disabled:
                d, S = 1, S_fwd
            elif pres_p and S_fwd < 0 and (c.key, -1) not in disabled:
                d, S = -1, -S_fwd
            else:
                continue
            if kinetics and (c.meta.get("slow") or d in c.meta.get("slow_dirs", ())):
                # 仅当慢反应可达显著程度才标注（忽略 ACT_FLOOR 量级通道）
                # 慢标注只需判定"能否达 1e-3 量级"：粗精度 12 次二分足够（分辨率 ~2.4e-4，对 1e-3 阈值充分）
                # （60 次高精度为鞍点路径敏感场景保留，此处无需）；
                # 标注是布尔终态——一旦确认过，后续迭代不再重复求解。
                # 负结论（低于阈值）按节奏重评：距上次评估 ≥32 迭代、或
                # 其间发生步量 ≥0.01 的实质步（慢通道平衡程度随状态演化）
                if not slow_seen:
                    pos = slow_ann.get((c.key, d))
                    if (pos is None or len(hist) - pos >= 32
                            or _last_mat_step >= pos):
                        ext_s, _ = solve_extent(c, d, ledger, H_excess, V, T_K, T, gsup,
                                                iters=12, p_ext_kpa=p_ext_kpa)
                        if ext_s >= ANN_MIN_EXTENT:
                            slow_now = True
                        else:
                            slow_ann[(c.key, d)] = len(hist)
                continue
            # 平衡冻结在候选评估阶段排除（而非 pick 后跳过）：
            # 已宣告平衡的净反应不参与竞争，让次优通道（如 ksp_beta）接手
            nk_c = _netkey(c, d)
            if nk_c in frozen_perm:
                if _TRACE: print('  [frozen]', c.key, d)
                continue
            if c.meta.get("deferred"):
                # 让位档：先收集；本轮若有任何快候选（evals 非空）则让位，
                # 快通道耗尽后才执行（如 PbS+H2O2 无快通道竞争，照常氧化）
                deferred_evals.append((c, d, S))
                continue
            evals.append((c, d, S))
        slow_seen = slow_seen or slow_now
        if not evals:
            if deferred_evals:
                # 无快候选：让位档出手，并按可观察慢反应标注
                evals = deferred_evals
                slow_seen = True
            else:
                # 终局慢通道复核：收敛点上的慢标注不受采样节奏影响（节奏
                # 缓存的负结论可能是收敛前的陈旧状态，此处统一重评一次）
                if kinetics and not slow_seen:
                    for c2 in cands:
                        _sl = (c2.meta.get("slow") or
                               1 in c2.meta.get("slow_dirs", ()) or
                               -1 in c2.meta.get("slow_dirs", ()))
                        if not _sl or (c2.key, 1) in disabled and (c2.key, -1) in disabled:
                            continue
                        _ps2 = c2.pres_specs
                        for d2 in (1, -1):
                            if d2 not in c2.meta.get("slow_dirs", ()) and not c2.meta.get("slow"):
                                continue
                            side = _ps2[0] if d2 > 0 else _ps2[1]
                            if not side or not all(ledger.get(s, 0.0) > X_MIN for s in side):
                                continue
                            S2f = S_of(c2, ledger, V, pH, T_K, T, gsup,
                                       p_ext_kpa, gas_escape, logc)
                            if d2 > 0 and S2f <= 0:
                                continue
                            if d2 < 0 and S2f >= 0:
                                continue
                            ext_s, _ = solve_extent(c2, d2, ledger, H_excess, V,
                                                    T_K, T, gsup, iters=12,
                                                    p_ext_kpa=p_ext_kpa)
                            if ext_s >= ANN_MIN_EXTENT:
                                slow_seen = True
                                break
                _exit_reason = "no-cands"   # 快慢候选全部耗尽（自然收敛点）
                break

        # ---- pick 与求解（主路径 / 微步排空共用）----
        # 微步排空（drain）：pick 被判微步禁用后，若账本/He 冻结（无强酸
        # 分子影子 → 顶部 respeciate 结构性空转；realize/枚举/pH 全为纯
        # 函数），本轮评估结果 evals 仍然完全有效——直接剔除已禁用者重挑
        # 次优，免去每个空转迭代一整轮 pH/形态/枚举/评估重算（多平衡体系
        # 空转可达数百迭代）。任何真实步/集合耗尽/idle 达上限即退回外层。
        _drain = False
        _drain_top = False
        _outer_break = False
        _joint_refined = False
        while True:
          # blocked 再验证（每轮）：膜已消失，或溶解通道能力 ≥ 现存膜量 → 解除封锁。
          # 封锁只在"通道溶不动膜"时维持（如 Mg/冷水、Zn/纯水）
          for m, films in list(blocked_solids.items()):
            alive = [fp for fp in films if ledger.get(fp, 0.0) > X_MIN]
            if not alive:
                del blocked_solids[m]
                continue
            for c2, d2, S2 in evals:
                reactants2 = c2.r if d2 > 0 else c2.pr
                targets = [fp for fp in alive if fp in reactants2]
                if S2 > 0 and targets:
                    ext_c, _ = solve_extent(c2, d2, ledger, H_excess, V, T_K, T, gsup,
                                            p_ext_kpa=p_ext_kpa)
                    if ext_c >= min(ledger.get(fp, 0.0) for fp in targets):
                        del blocked_solids[m]
                        break

          # 溶剂优先（界面动力学规则，仅此一条）：无游离强酸且**无更强氧化剂在账**
          # 时水还原优先（防止裸离子路径虚报；强氧化剂在场时按 S 正常竞争）
          strong_ox = [e for e in evals if e[0].kind == "redox"
                       and e[0].meta.get("ox_couple") != H_ION and e[2] >= SOLVENT_FIRST]
          water_c = [] if strong_ox else [
              e for e in evals if e[0].meta.get("ox_couple") == H_ION
              and e[2] >= SOLVENT_FIRST and H_excess / V < 1e-3]
          pick, d, S = max(water_c, key=lambda e: e[2]) if water_c else max(evals, key=lambda e: e[2])
          if _TRACE: print('  [pick]', pick.key, d, round(S,2))

          nk = _netkey(pick, d)

          # 微步快通道启用条件（与 solve_extent 的 micro_rel 契约配套）：
          # ① kind 非 dissolve/precip——微溶盐微步（X_MIN < ext < 1e-4·x_max）
          #   会真实执行并改账本，需精确 ext；
          # ② nk 不在最近 12 个已执行步的键中——振荡外推块以 hist 尾部
          #   周期模式 + keys[0]==nk 为前置，ext 参与跳幅计算；此守卫保证
          #   外推块对本 pick 恒为空转；
          # ③ x_max 计量上界 * 1e-4 < BLOCKED_EXTENT(0.02)——膜封锁比较
          #   ext > bext 在微步量级恒不触发（双保险，防极端大投料）
          _micro_rel = None
          if pick.kind not in ("dissolve", "precip"):
            _rr_pk = pick.r if d > 0 else pick.pr
            _xm = min((ledger.get(s, 0.0) / nu for s, nu in _rr_pk.items()
                       if s not in (WATER, H_ION)), default=0.0)
            if (1e-4 * _xm < BLOCKED_EXTENT
                    and nk not in {k for k, _e in hist[-12:]}):
                _micro_rel = 1e-4
          ext, x_max = solve_extent(pick, d, ledger, H_excess, V, T_K, T, gsup,
                                    p_ext_kpa=p_ext_kpa, micro_rel=_micro_rel)
          if _TRACE: print('    [ext]', round(ext,5), 'pH', round(pH,2), 'He', round(H_excess,4))
          bext, film_ps = (blocked_extent(pick, d, evals, T_K, T, ledger, H_excess, V,
                                          x_max, pH, p_ext_kpa=p_ext_kpa)
                           if kinetics else (None, []))
          if bext is not None and ext > bext:
            ext = bext
            if "blocked" not in annotations:
                annotations.append("blocked")
            rr0 = pick.r if d > 0 else pick.pr
            for s in rr0:
                if s in T.solids:
                    blocked_solids[s] = film_ps
          # 振荡外推加速：短周期（2/3 通道）等比衰减爬行（沉淀↔逆转化乒乓、
          # 溶解↔双沉淀三循环等，步长比 ρ→1）时，几何级数剩余工作量
          # ≈ ext·ρ/(1−ρ)，一次补齐直抵不动点附近，避免数百步微步爬行；
          # 补齐受 x_max 钳制，外推失准由 S<0 反向步自动纠回
          # （不动点两侧 S 变号是引擎的既有自校正）
          for _p in range(2, 7):
            if len(hist) < 2 * _p or ext <= 0:
                continue
            tail = hist[-2 * _p:]
            keys = [k for k, _ in tail]
            if keys[:_p] != keys[_p:] or keys[0] != nk:
                continue
            if len({k for k in keys}) < 2:
                continue
            e0, e1_ = tail[0][1], tail[_p][1]
            if e0 <= 0 or e1_ <= 0:
                continue
            rho = e1_ / e0
            if 0.5 < rho < 0.999:
                jump = min(x_max, ext / (1 - rho)) - ext
                if jump > 10 * ext:
                    if _TRACE: print('    [osc]', _p, round(rho, 3), '+', round(jump, 4))
                    ext += jump
            break
          if ext <= max(X_MIN, 1e-4 * x_max):
            # 微步（已达平衡或程度可忽略）：双向禁用该平衡，等待状态实质改变。
            # 程度绝对可辨（>X_MIN）时先执行这一次再禁用——微溶盐（AgCl、
            # BaSO4 等）的溶解度就是这一步到位的热力学终态，直接跳过会把
            # 微溶事实整体漏报；真正的零推进由签名机制在下一轮挡下
            sig0 = _sig()
            disabled[(pick.key, d)] = sig0
            disabled[(pick.key, -d)] = sig0
            # 仅溶解/沉淀类微步执行（微溶盐终态）；质子/氧化还原微步仍跳过——
            # 后者执行会经签名变化逐对渗漏（NH4Ac 双水解曾被渗到 pH 9.4）
            if ext <= X_MIN or pick.kind not in ("dissolve", "precip"):
                idle += 1
                if idle >= 8:
                    # （触发点② idle 退出前精修已回退——v0.3.8 全量差分否决：
                    # Co32/Co33/T34/Fe33 四例翻车。根因：idle 点上走步仲裁
                    # （freeze/max-S 竞争）已完成，活动子集的联立不动点 ≠
                    # 走步仲裁点（候选网络跨数据源 Hess 不一致，走步的妥协
                    # 点才是被 1247 例锁定的语义）；跳步凌驾仲裁=语义翻案。
                    # J06 型欠收敛的真根因另查：pH 滴定端点悬崖（He=1e-3
                    # 恰为缓冲容量时 pH 3.0→6.2 突变，见 tests.json 注记），
                    # 非 walk 收敛缺陷，联立解无法跨越不连续。）
                    _exit_reason = "idle"
                    _outer_break = True
                    break
                if not _drain:
                    # 入排空前提：respeciate 结构性空转（无强酸分子影子——
                    # 否则顶部再电离可能改账本，evals 会陈旧）；
                    # realize/枚举/pH 均为账本纯函数，冻结即不变
                    if (ledger.get("__tot_HNO_3", 0.0) <= 0.0
                            and ledger.get("__tot_H_2SO_4", 0.0) <= 0.0):
                        _drain = True
                    else:
                        _drain_top = True
                        break
                # 剔除已禁用 pick 后重挑（候选集只缩不涨；双向禁用与
                # 顶部评估循环的跳过条件等价）
                evals = [e for e in evals if e[0] is not pick]
                if not evals:
                    # 候选集耗尽 → 外层顶（evals 重建为空 → deferred/
                    # 终局慢复核/break 路径，与原实现逐迭代一致）
                    _drain_top = True
                    break
                continue
            # 微溶盐溶解/沉淀微步：执行这一次（热力学终态），双向禁用已在上方
            # 登记——签名机制保证下一轮挡下，无需额外停滞语义
          else:
            disabled.clear()   # 状态将发生实质改变，解禁全部（签名机制双保险）
          break
        if _outer_break:
            break
        if _drain_top or _joint_refined:
            continue   # 联立精修刚移动状态：重评估（本微步不执行）

        idle = 0
        _exec(pick, d, ext, x_max, S, nk)
        # 极限环检测（仅在实质步后判定）：账本签名精确复现 ⇒ 确定性求解器
        # 进入零净推进循环，永久冻结窗口内全部净反应让其他通道接手；
        # 窗口为空（微步原地）则不动作，交由 disabled/微步机制处理
        sig_now = _sig()
        if sig_now in seen_sig:
            window = {k for k, _ in hist[seen_sig[sig_now]:]}
            window -= {k for k in window if k in frozen_perm}
            if window:
                for k in window:
                    frozen_perm.add(k)
                    frozen_perm.add((k[1], k[0]))
                _diag["freeze_events"] += 1
                if _TRACE: print('  [freeze-limit-cycle]', window)
        else:
            seen_sig[sig_now] = len(hist)
        # 实测震荡判定（近窗）：最近 10 步内正反向执行程度接近抵消
        # （|净| < 10% 总量）且总量显著 → 原地空转，永久冻结。
        # 限近窗是因为跨长程的"抵消"往往是平衡被其他通道移动后的正常演化
        rev = (nk[1], nk[0])
        recent = hist[-10:]
        ext_fwd = sum(e for k, e in recent if k == nk)
        ext_rev = sum(e for k, e in recent if k == rev)
        gross = ext_fwd + ext_rev
        if gross >= 0.05 and abs(ext_fwd - ext_rev) <= 0.1 * gross:
            frozen_perm.add(nk)
            frozen_perm.add(rev)
            _diag["freeze_events"] += 1
            if _TRACE: print('  [freeze-perm]', nk, 'gross', round(gross, 3))
        # 循环震荡检测：最近若干步为同一短周期（长度 2 或 3）反复且总推进量
        # 低于显著阈值 → 冻结该周期涉及的全部净反应（E35 类阶梯每周期有实质
        # 推进，不受影响；34 类 A->B->A->B 原地空转被捕获）
        for period in (2, 3):
            w = 3 * period
            if len(hist) >= w:
                tail = [k for k, _ in hist[-w:]]
                unit = tail[:period]
                recent_ext = [e for _, e in hist[-w:]]
                if tail == unit * 3 and sum(recent_ext) < 1e-3:
                    for k in set(unit):
                        frozen_perm.add(k)
                        frozen_perm.add((k[1], k[0]))
                    _diag["freeze_events"] += 1
                    if _TRACE: print('  [freeze-cycle]', unit)

        # 物种级周转冻结（v0.3.8）：多平衡乒乓（不同净键互为往返——extent
        # 层净/毛恒 1 但物种账本原地踏步；Ag32 的配位-氧化银 876 次乒乓：
        # R1 产 Ag2O、R2 溶回，数据张力型无联立不动点，走步 freeze 是唯一
        # 正解但触发太晚）。判据：近窗（40 步）物种总周转（逐步 L1 累计）
        # 显著 ≥0.05 而账本净移 ≤2% → 空转实锤，冻结窗内全部净键。
        # 与既有检测器同族；genuine 爬行（净/毛 > 2%，ρ<0.995）不触发。
        # 每 20 步评估一次（O(窗×物种)，热路径成本可忽略）
        if (len(snaps) >= 81 and len(snaps) % 20 == 0):
            _FW = 40
            s_a, _ = snaps[-1 - _FW]
            turnover = 0.0
            for _i in range(len(snaps) - _FW, len(snaps) - 1):
                s_0, _ = snaps[_i]
                s_1, _ = snaps[_i + 1]
                turnover += sum(abs(m - s_0.get(s, 0.0))
                                for s, m in s_1.items()
                                if not s.startswith("__") and s != WATER)
            drift = sum(abs(m - s_a.get(s, 0.0)) for s, m in snaps[-1][0].items()
                        if not s.startswith("__") and s != WATER)
            if turnover >= 0.05 and drift <= 0.02 * turnover:
                window = {k for k, _ in hist[-_FW:]}
                window -= {k for k in window if k in frozen_perm}
                # v0.4.2 迭代 C：冻结范围升级为整个活性社区——多拼写
                # 乒乓（同一化学转化的 H⁺/NH₄⁺ 挂侧变体）逐键冻结太慢
                # （H46 死因②：冻结逐拼写进行），社区级一次止震。
                # 真爬行（净/毛 > 2%）不触发本检测器，天然不受影响
                for c, _d, _S in _joint_collect(None):
                    nk_c = (c.netkey_fwd if c.netkey_fwd <= c.netkey_rev
                            else c.netkey_rev)
                    if nk_c not in frozen_perm:
                        window.add(nk_c)
                if window:
                    for k in window:
                        frozen_perm.add(k)
                        frozen_perm.add((k[1], k[0]))
                    _diag["freeze_events"] += 1
                    if _TRACE:
                        print(f'  [freeze-turnover] {len(window)} keys '
                              f'drift/turnover={drift / max(turnover, 1e-9):.4f}')

        # ---- 窗口净移/毛周转标定导出（CHEM_TRACE_WINDOWS=1，纯诊断）----
        # 每 32 步窗口的 (drift, turnover)：真爬行（净/毛 > 0.5）、螺旋
        # （~0.2）、纯乒乓（<0.05）的判别数据源——不进 digest、不导出到
        # converg dump（env 专属，走 probe["windows"] 由单独脚本消费）
        if _WINDOWS and len(snaps) >= 33 and len(snaps) % 32 == 0:
            _s0, _ = snaps[-33]
            _to = 0.0
            for _i in range(len(snaps) - 32, len(snaps) - 1):
                _sa, _ = snaps[_i]
                _sb, _ = snaps[_i + 1]
                _to += sum(abs(m - _sa.get(s, 0.0)) for s, m in _sb.items()
                           if not s.startswith("__") and s != WATER)
            _dr = sum(abs(m - _s0.get(s, 0.0))
                      for s, m in snaps[-1][0].items()
                      if not s.startswith("__") and s != WATER)
            _diag["windows"].append((round(_dr, 9), round(_to, 9)))

        # ---- pH 悬崖乒乓冻结（v0.4.0 第五检测器）----
        # J06 型失稳的走步级兜底：pH 机器在 He 跨 1e-3 时跳分支，跨键 H⁺
        # 喂食循环里各步的 S 评估在 bisection 探针内跨悬崖相位——各自为
        # 正而联合为负，走步把伪过程当真爬山（净推进 2-20% 的净键循环，
        # 四种既有检测器均不覆盖；主迭代 pH 被固相滴定掩盖恒 5.8，摆动
        # 只在探针内——改用步后 He 序列的符号交替率做喂食零和签名）。
        # 判据：近窗全微步 + H⁺ 参与步过半 + He 符号交替率 ≥ 2/3——真实
        # 酸碱过程 He 单向变化不满足。冻结窗内净键（平衡止震语义：
        # 真实化学的 H⁺ 是连续分布，乒乓是 pH 估计不连续的数值实现）。
        # O02 锚点：Zn+FeCl₂ 置换 ext=1 完全后悬崖乒乓首检出 ≤32 步
        # （Fe 损失 ≤0.02）；全量差分审误杀面。
        _PHW = 24
        if (len(hist) >= _PHW and len(hist) % 8 == 0
                and max(e for _, e in hist[-_PHW:]) < 0.02
                and sum(hexec[-_PHW:]) >= _PHW // 2
                and len(snaps) >= _PHW):
            _he_seq = [_h for _s, _h in snaps[-_PHW:]]
            _alt = sum(1 for _i in range(1, _PHW)
                       if _he_seq[_i - 1] * _he_seq[_i] < 0.0)
            if (_alt * 3 >= 2 * (_PHW - 1)
                    and max(abs(min(_he_seq)), abs(max(_he_seq))) <= 3e-3):
                # 悬崖带判据：He 贴 ±1e-3 阈值摆（pH 分支界上的乒乓），
                # 幅度 ≤3e-3——E23 型真实缓冲震荡（He ±1.3e-2，Ag⁺ 水解
                # Henderson 吸收释放）不触发，只有悬崖乒乓贴地摆动
                # 只冻结 H⁺ 喂食步的键（E23 教训：真实置换键与喂食步
                # 同窗交织时，全窗冻结会误杀置换——旁观者不冻）
                window = {k for _i, (k, _e) in enumerate(hist[-_PHW:])
                          if hexec[-_PHW:][_i]}
                window -= {k for k in window if k in frozen_perm}
                if window:
                    for k in window:
                        frozen_perm.add(k)
                        frozen_perm.add((k[1], k[0]))
                    _diag["freeze_events"] += 1
                    if _TRACE:
                        print(f'  [freeze-ph-cliff] {len(window)} keys '
                              f'He-alt {_alt}/{_PHW - 1} '
                              f'[{min(_he_seq):.2e},{max(_he_seq):.2e}]')

        # ---- 爬行收敛加速（窗口几何外推，非精确周期的慢收敛体系） ----
        # 触发条件：深度入局（≥2 窗历史）+ 冷却期满（上次外推后 ≥1 窗新步）
        # + 近窗全部微步（<0.02 mol，排除中途大步使外推失义）。
        # 守卫三重：窗差 D 的 L1 范数相对前窗等比衰减 ρ∈(0.25, 0.85)；
        # 净/毛比 n1 ≥ 0.4×窗内步量绝对值和（纯爬行逐步向前、净差与步量
        # 同量级；正逆步对消的震荡窗净差≪步量和——外推方向≠不动点方向，
        # Ag32 曾因 ρ=0.9 的 9× 外推沿震荡方向冲出物料守恒）；外推倍率
        # 封顶 4×（ρ→1 时几何尾部爆炸，宁可分多窗跳）。
        if (len(snaps) >= 2 * _CRAWL_W + 1
                and len(hist) - last_jump >= _CRAWL_W
                and max(e for _, e in hist[-_CRAWL_W:]) < 0.02):
            led1, he1 = snaps[-1]
            led2, he2 = snaps[-1 - _CRAWL_W]
            led3, he3 = snaps[-1 - 2 * _CRAWL_W]
            D = {}
            n1 = 0.0
            for s, m in led1.items():
                d1 = m - led2.get(s, 0.0)
                if abs(d1) > 1e-12:
                    D[s] = d1
                    n1 += abs(d1)
            d_he1 = he1 - he2
            n1 += abs(d_he1)
            n2 = sum(abs(m - led3.get(s, 0.0)) for s, m in led2.items())
            n2 += abs(he2 - he3)
            gross = sum(abs(e) for _, e in hist[-_CRAWL_W:])
            if (n1 > 1e-9 and n2 > 1e-9 and gross > 1e-9
                    and 0.25 < n1 / n2 < 0.85 and n1 >= 0.4 * gross):
                rho = n1 / n2
                scale = min(rho / (1.0 - rho), 4.0)
                # v0.4.0 守恒守卫：D 的线性缩放保持原子守恒的前提是不把
                # 任何物种推负——v0.3.8 的"钳零"路径在弱配+沉淀强耦合体系
                # （N34 的 decomplex↔derived 配对交换循环）上把 Zn-Cl 池
                # 成员推负后钳零，Zn 总量凭空 +0.26、净方程系数 16091。
                # 净减物种的安全倍率上界（0.95 折扣留收敛余量），倍率降到
                # 线性域内——外推加速保留，守恒硬保证
                _safe = 4.0
                for _s, _d in D.items():
                    if _d < -1e-12 and not _s.startswith('__'):
                        _cap = 0.95 * ledger.get(_s, 0.0) / (-_d)
                        if _cap < _safe:
                            _safe = _cap
                scale = min(scale, _safe)
                if scale >= 0.05:
                    for s, d in D.items():
                        if s.startswith("__"):
                            continue   # 影子库存由 _respeciate 每轮重建
                        tgt = ledger.get(s, 0.0) + d * scale
                        ledger[s] = tgt if tgt > 0.0 else 0.0
                    H_excess += d_he1 * scale
                    disabled.clear()   # 状态实质改变，解禁全部让 S 重验
                    if _TRACE:
                        print(f'  [crawl-jump] rho={rho:.3f} +{scale:.1f}x '
                              f'D[{len(D)}] n1/gross={n1/max(gross,1e-9):.2f} '
                              f'HeΔ={d_he1:.2e}')
                last_jump = len(hist)

        # ---- 联立求解加速·触发点①（爬行检测，迭代型）----
        # it≥64 且近 24 步全微步（<0.02）且 ≥2 个不同规范净键（耦合循环，
        # 非单平衡停滞）；冷却 32 迭代防反复失败；失败签名黑名单（同循环
        # 结构重试必然同果——RX13 曾 6 次尝试 4 次白费）。E23 型（249 迭代
        # 仅 ~20 执行步）与 Ag32 型（乒乓）在此命中；只联立解实际循环的
        # 平衡——候选集存在 Hess 互斥副本，全收会因数据不一致无解
        if (it >= 64 and len(hist) >= 8
                and it - _joint_crawl_last >= 32
                and len({k for k, _ in hist[-24:]}) >= 2
                and max((e for _, e in hist[-24:]), default=1.0) < 0.02):
            _ck = _joint_cycle_keys()
            if tuple(sorted(_ck)) not in _joint_blacklist:
                _joint_crawl_last = it
                if not _joint_fire(_ck, freeze_on_boundary=True):
                    _joint_blacklist.add(tuple(sorted(_ck)))

        # ---- 联立求解加速·触发点③（社区级 pH 一致化，v0.4.2 迭代 D）----
        # 与①同微步窗（独立冷却），社区口径 pH 一致化联立（族闸内建：
        # ≥2 质子族硬失败不跳）。失败签名黑名单防同结构重试
        if (it >= 64 and len(hist) >= 8
                and it - _joint_ph_last >= 32
                and len({k for k, _ in hist[-24:]}) >= 2
                and max((e for _, e in hist[-24:]), default=1.0) < 0.02):
            _sig_ph = tuple(sorted(_joint_cycle_keys()))
            if _sig_ph not in _joint_ph_blacklist:
                _joint_ph_last = it
                if not _joint_fire_ph():
                    _joint_ph_blacklist.add(_sig_ph)

      # 不动点收敛后扫气一轮：逸出离账会移动平衡（Le Chatelier），
      # 有新增逸出则再跑一轮不动点；无新增即全局收敛
      _n0 = sum(escaped.values())
      _sweep_gases(ledger, escaped, gsup, V, T_K, p_ext_kpa, T)
      if sum(escaped.values()) - _n0 <= 1e-9:
        break

    if _probe is not None:
        _probe_exit(_probe, ledger, H_excess, escaped, gsup, V, T_K, T,
                    kinetics, gas_escape, p_ext_kpa, disabled, frozen_perm,
                    _exit_reason, _it_total, len(hist), len(steps),
                    blocked_solids)
        # 走步画像诊断（v0.4.2 探针扩展）：纯诊断字段，不进
        # _result_digest（digest 键集不含它们），零行为影响
        _probe["joint_tries"] = _diag["joint_tries"]
        _probe["joint_ok"] = _diag["joint_ok"]
        _probe["freeze_events"] = _diag["freeze_events"]
        _probe["micro_steps"] = _diag["micro_steps"]
        if _diag["windows"]:
            _probe["windows"] = _diag["windows"]
    if slow_seen:
        annotations.append("slow")
    return _finalize_result(ledger, initial, H_excess, H_excess0, escaped, steps,
                            chem_net, annotations, unknown, cond, V, T_K, T)


def _presentation_He(ledger: dict, H_excess: float, V: float, T,
                     T_K: float) -> float:
    """呈现层质子账本自洽闸（v0.4.1）：强碱分支的阳离子-氢氧化物
    过饱和检查。

    欠收敛残余的幻影碱（RX13 型：Ag⁺/NO₃⁻/HNO₂ 硝酸循环 idle 退出的
    负 He 恰落在 −1e-3 强碱分支界，字面读数 pH 11+）与在账阳离子在
    [OH⁻] = −He/V 下的离子积超 Ksp ≥ 3 个数量级时——热力学上不可能
    共存（真实碱体系的走步已把阳离子沉淀至 Q≈Ksp，本检查天然通过），
    判定“自由强碱”为欠收敛幻影：He 呈现归零（被沉淀概念吸收），
    pH 回落分支 4（缓冲/阳离子水解主导——RX13 的 Cu²⁺ 水解 ~4.7）。
    只改终态呈现（final_pH / H_excess / 收敛探针画像），走步机器
    零接触（在环 pH 估计/酸闸门/探针 f(x) 全部维持原语义）。

    分支镜像：先过缓冲滴定（与 estimate_state 同口径）——可被在账
    弱酸吸收的碱是缓冲化学（SE03 的 NH₄⁺ Henderson 7.5），非分支 3
    场景原样放行；只有滴定无解且残余 He ≤ −1e-3（真会走强碱分支）
    才做阳离子过饱和检查。"""
    pKw = pKw_of(T_K)
    tit, He_res, _led_t = _buffer_titration(ledger, H_excess, V, T, pKw,
                                            T_K=T_K)
    if tit is not None or He_res > -1e-3:
        return H_excess   # 缓冲有解 / 残余不到分支 3：原样放行
    oh = (-He_res) / V
    for e in T.ksp:
        cat, an = e["pair"]
        if an != "OH^-":
            continue
        m = ledger.get(cat, 0.0)
        if m <= X_MIN:
            continue
        x, y = _ksp_xy(e)
        q = (m / V) ** x * oh ** y
        if q > 1e3 * 10.0 ** (-_pksp(e, T_K)):
            return 0.0    # 幻影碱：被阳离子氢氧化物吸收（呈现层）
    return H_excess


def presentation_pH(ledger: dict, H_excess: float, V: float, T,
                    T_K: float) -> float:
    """**呈现层 pH**（冷路径，每次判定只算一次）：B3 档走精确质子条件，
    其余档走 pH 机器（architecture §7 O-4 的接线决定）。

    为什么只在这里接：分支 4 的两性中点式/缓冲对加权在 B3 档有实测误差
    （AB03 10.50 vs 精确 9.84；Q05 的 1:1 H₂S/HS⁻ 缓冲对被当成纯两性盐
    给 10.50，真值 = pKa₁ = 7.00；H88 7.246 vs 7.005），而精确解进二分
    探针要 +16% 墙钟（§7 O-2/O-3）。呈现层不在探针里 ⟹ 零成本，
    且改的正是用户直接看到的那个数（`final_pH` / 探针画像 / 质量口径）。
    储库档（固相 / 水解阳离子）与账本不自洽的档由 `exact_proton_pH`
    的判据自动排除，机器的既有行为在那两档原样保留。"""
    p = exact_proton_pH(ledger, H_excess, V, T, T_K)
    return estimate_pH(ledger, H_excess, V, T, T_K) if p is None else p


def _probe_exit(probe: dict, ledger: dict, H_excess: float, escaped: dict,
                gsup: frozenset, V: float, T_K: float, T, kinetics: bool,
                gas_escape: bool, p_ext_kpa: float, disabled: dict,
                frozen_perm: set, exit_reason: str, it_total: int,
                hist_len: int, steps_len: int,
                blocked_solids: set = frozenset()) -> None:
    """收敛质量探针（只读诊断，v0.3.8）：walk 退出点上的平衡残差画像。

    对最终账本重新枚举候选，记录"两侧均在场"的平衡的 S = logK − logQ
    （热力学上应全部为零；|S| 大 = 欠收敛）。联立求解器的验收基准
    （converg.py dump/diff）与触发判据的数据源。不改任何求解行为。

    v0.5.0：每条 active 记录带**逐方向 disabled**、`slow`、`blocked`
    （膜封锁）三个"引擎既定语义"标记。质量口径必须只统计 walk **真正
    会考虑**的平衡——否则指标会去追引擎明确拒绝的方向（实测：TS04 的
    `resid_live` 长期被幻影硫酸盐慢通道顶到 17.975，E24 被膜封锁的
    金属-水通道顶到 156.12）。"""
    from .templates import enumerate_candidates as _enum
    led = dict(ledger)   # 副本隔离：respeciate 会原地改账本（探针零副作用）
    H_excess = _respeciate_strong_acids(led, H_excess, V, T)
    H_excess = _presentation_He(led, H_excess, V, T, T_K)
    pH_f = presentation_pH(led, H_excess, V, T, T_K)
    cands_f = _enum(led, H_excess, pH_f, V, T_K, T, kinetics)
    active = []
    logc: dict = {}

    def _peq(c) -> str:
        rr = " + ".join(f"{_fmt(nu)}{s}" for s, nu in c.r.items() if s != WATER)
        pp = " + ".join(f"{_fmt(nu)}{s}" for s, nu in c.pr.items() if s != WATER)
        return f"{rr} -> {pp}"

    for c in cands_f:
        ps = c.pres_specs
        pres_r = all(led.get(s, 0.0) > X_MIN for s in ps[0])
        pres_p = all(led.get(s, 0.0) > X_MIN for s in ps[1])
        if not (pres_r or pres_p):
            continue
        S_f = S_of(c, led, V, pH_f, T_K, T, gsup, p_ext_kpa,
                   gas_escape, logc)
        # 逐方向 disabled：walk 的评估循环按"该方向未禁用"才入选
        # （engine L903/905），只看双向同时禁用会漏报——TS04 的三条
        # S≈+17.6 候选正是**单向**禁用，旧口径报 dis=0 造成"walk 看不见
        # 仍有驱动的候选"的假象。
        df, dr = (c.key, 1) in disabled, (c.key, -1) in disabled
        sl = bool(c.meta.get("slow") or 1 in c.meta.get("slow_dirs", ())
                  or -1 in c.meta.get("slow_dirs", ()))
        # 膜封锁（与 walk L892 同口径）：致密膜只抑制溶剂氧化通道
        bl = bool(c.kind == "redox" and c.meta.get("ox_couple") == H_ION
                  and any(s in blocked_solids
                          for s in list(c.r) + list(c.pr)
                          if s not in (WATER, H_ION)))
        # 可达程度上界（化学计量上限，无需求解）：|S| 大 ≠ 会反应——
        # 痕量物种对可以有巨大的 log 驱动力而只能走 ~1e-6 mol。
        # D38 即此：PbO₂ 2.9e-06 / Cr³⁺ 1.9e-06 mol 给出 S=+33.4。
        # 引擎自身的"显著程度"判据是 ANN_MIN_EXTENT（slow 标注同用），
        # 质量口径据此区分"欠收敛"与"无关的痕量方向"。
        _dr = c.r if S_f > 0 else c.pr
        _lim = [led.get(s, 0.0) / nu for s, nu in _dr.items()
                if s not in (WATER, H_ION) and nu > 0]
        ex = min(_lim) if _lim else float("inf")
        if not (pres_r and pres_p):
            # 单侧在场：驱动属正常（反应物耗尽/产物未生）；仅记录不判残差
            active.append({"kind": c.kind, "eq": _peq(c),
                           "S": round(S_f, 3), "two_sided": False,
                           "frozen": c.netkey_fwd in frozen_perm
                                     or c.netkey_rev in frozen_perm,
                           "slow": sl, "blocked": bl, "ext_max": round(ex, 9),
                           "dis_fwd": df, "dis_rev": dr,
                           "disabled": df and dr})
            continue
        active.append({"kind": c.kind, "eq": _peq(c), "S": round(S_f, 3),
                       "two_sided": True,
                       "frozen": c.netkey_fwd in frozen_perm
                                 or c.netkey_rev in frozen_perm,
                       "slow": sl, "blocked": bl, "ext_max": round(ex, 9),
                       "dis_fwd": df, "dis_rev": dr,
                       "disabled": df and dr})
    probe.clear()
    probe.update({
        "exit": exit_reason, "iters": it_total, "hist": hist_len,
        "steps_n": steps_len, "pH": round(pH_f, 3),
        "H_excess": round(H_excess, 9),
        "ledger": {s: round(m, 9) for s, m in led.items() if m > 0.0},
        "active": active,
        "max_abs_S": max((abs(a["S"]) for a in active if a["two_sided"]),
                         default=0.0),
        "frozen_n": sum(1 for a in active if a["frozen"]),
        "disabled_n": sum(1 for a in active if a["disabled"]),
    })



def _apply_override(substances, cond, V, T_K, T, ov) -> dict:
    """OVERRIDE 路径（v0.3.7 从 judge() 提出）：注册反应事实的直出结果，
    化学计量按限量试剂缩放（H+/OH- 按游离酸碱储备限量——酸碱不足时
    override 不许过量产出）。"""
    res = dict(ov["result"])
    # 按限量试剂缩放化学计量（result 中的 mol 为每单位反应式的量）
    led0, _He0, _st0, _un0 = normalize(substances, cond, T)
    scale = float("inf")
    for c0 in res.get("consumption", []):
        # H2O 为溶剂不记账；H+/OH- 以 He 记账——按游离强酸/强碱储备限量
        # （原一律跳过限量，酸/碱不足时 override 会过量产出）
        if c0["name"] == WATER:
            continue
        if c0["name"] == H_ION:
            scale = min(scale, max(_He0, 0.0) / c0["mol"])
            continue
        if c0["name"] == "OH^-":
            scale = min(scale, max(-_He0, 0.0) / c0["mol"])
            continue
        scale = min(scale, led0.get(c0["name"], 0.0) / c0["mol"])
    if scale == float("inf"):
        scale = 1.0
    scale = max(scale, 0.0)
    res["consumption"] = [dict(c0, mol=round(c0["mol"] * scale, 6))
                         for c0 in res.get("consumption", [])]
    res["production"] = [dict(c0, mol=round(c0["mol"] * scale, 6))
                        for c0 in res.get("production", [])]
    res["steps"] = []
    res["final"] = []
    res["unknown"] = []
    res["escaped"] = []
    res["override"] = ov["id"]
    # OVERRIDE 注册的是反应事实：真实转化的（如白磷歧化）reacted=True，
    # 钝化类（changed=False，仅成膜保护）无化学反应
    res["reacted"] = bool(res.get("changed"))
    res["final_pH"] = None
    res["H_excess"] = 0.0
    res["H_excess_initial"] = round(_He0, 6)
    # 初态：post-normalize ledger（仅 strong electrolyte 已电离 + 中和已记账）
    init_dict: dict[str, float] = {s: round(m, 6) for s, m in led0.items()
                                   if s != WATER and not s.startswith("__")
                                   and m > 1e-6}
    if _He0 > 1e-6:
        init_dict[H_ION] = round(_He0, 6)
    elif _He0 < -1e-6:
        init_dict["OH^-"] = round(-_He0, 6)
    res["initial"] = [{"name": s, "mol": m} for s, m in init_dict.items()]
    res["ionize"] = _ionize_map(
        [c0["name"] for c0 in res.get("consumption", [])
         + res.get("production", [])], T)
    # 条件元数据透传（独立温度模块事后消费；不影响求解）
    res["cond"] = {"V_L": round(V, 4), "T_K": T_K,
                   "isothermal": bool(cond.get("isothermal", False)),
                   "kinetics": bool(cond.get("kinetics", True)),
                   "gas_escape": bool(cond.get("gas_escape", True))}
    return res


def _finalize_result(ledger, initial, H_excess, H_excess0, escaped, steps,
                     chem_net, annotations, unknown, cond, V, T_K, T) -> dict:
    """judge() 收尾（v0.3.7 从 654 行单体中提出的纯后处理段）：

    账本净差 → consumption/production/final/escaped 列表（1e-6 报告阈）、
    A1 changed 三通道事件判定、B4 reacted 狭义化学判定、degree 程度、
    初态重建与结果字典装配。只读输入状态，无走步副作用——判定语义的
    唯一居所（评审入口）。
    """
    # 呈现层自洽闸（v0.4.1）：幻影碱不进 final_pH / H_excess 呈现
    # （阳离子-氢氧化物过饱和检查；方程构建器同吃这两个字段，
    # 净方程的 H⁺/OH⁻ 净差线随之自洽）
    H_excess = _presentation_He(ledger, H_excess, V, T, T_K)
    consumption = [{"name": s, "mol": round(initial[s] - ledger.get(s, 0.0), 6)}
                  for s in initial if s != WATER and not s.startswith("__")
                  and initial[s] - ledger.get(s, 0.0) > 1e-6]
    production = [{"name": s, "mol": round(ledger.get(s, 0.0) - initial.get(s, 0.0), 6)}
                 for s in ledger if s != WATER and not s.startswith("__")
                 and ledger.get(s, 0.0) - initial.get(s, 0.0) > 1e-6]
    # 逸出气体计入 production（反应事实），但在 final 中只剩溶解态
    esc_list = [{"name": s, "mol": round(m, 6)} for s, m in escaped.items() if m > 1e-6]
    if esc_list:
        have = {e["name"] for e in production}
        for e in esc_list:
            if e["name"] in have:
                for pe in production:
                    if pe["name"] == e["name"]:
                        pe["mol"] = round(pe["mol"] + e["mol"], 6)
            else:
                production.append(dict(e))
    final = [{"name": s, "mol": round(m, 6)} for s, m in ledger.items()
             if s != WATER and not s.startswith("__") and m > 1e-6]
    main_steps = [st for st in steps if st.get("extent", 0) >= ANN_MIN_EXTENT
                  and st["kind"] != "neutralize"]
    # A1 净反应事件判定（v0.3.7）：changed = 存在显著净变化，按三条正交通道：
    #   (a) 溶液相重分布：物种账本差 ÷ 该物种在净显著步中的最大计量系数
    #       ν̄（事件量口径）。直除账本差（v0.3.6 口径）被配位/氧化还原计量
    #       系数放大——ν=4 弱配位金属仅 0.86% 转化，配体账本差 4×=3.4 mmol
    #       即触发阈值，在 NR 盐混合中制造 changed=True 幻影，迫使 v0.3.6
    #       撤回真实存在的弱卤配合物数据（违反"宁缺毋假≠删真数据"原则）。
    #       事件口径双向守恒：HAc 电离 1.3 mmol 事件（pH 可测）仍显著；
    #       NiCl_4 痕量配位 0.86 mmol 事件不显著。
    #   (b) 新相生成：固相产物量 ≥ 阈值（不折扣——可见浑浊量与事件量纲
    #       不同：KIN06 缓氧化 0.5 mmol 事件生成 1.6 mmol Fe(OH)_3 浑浊；
    #       RX04 逆向"不反应"体系中 Fe^{3+} 水解沉淀 6 mmol——教材口径
    #       FeCl_3 溶液必须加盐酸配制抑制的水解事实）。
    #   (c) 质子账本：游离 H+/OH- 不入 consumption 账本，其中和以
    #       neutralize 步显著程度判定；气体逸出（非步进通道）单独计显著。
    #   pH 实现通道（W21 HClO+NaOH 整域酸碱转化不经步进记录）物种未被
    #   步进触碰，ν̄=1 直用账本差。AlCl_3 水解三步循环（derived→derived→
    #   precip 溶解）各键净量大但账本净差 ~0.3 mmol：循环自抵不过 ν̄ 通道
    _nu_max: dict[str, int] = {}
    for _r, _pr in chem_net:
        for _sp, _n in _r + _pr:
            if _n > _nu_max.get(_sp, 1):
                _nu_max[_sp] = _n
    # v0.4.0 弱池中心守恒口径：池成员的账本差不直接参与通道 (a)（池内
    # 再分布是形态变化非事件——B4 判据的 A1 延伸），改用池总量
    # （Σ成员）初→终变化判定事件（ν=1 中心转移计量）。NR19 盐混合
    # 78/22 池内分布不再制造 changed=True 幻影；O02 金属→池 0.95 仍事件。
    # 配体守恒配套：池内配合物净增吸入的配体（Cl⁻→[ZnCl]⁺ 摄入）从
    # 配体的通道 (a) 消耗中抵扣——NR19 的 Cl⁻ 0.42 全部进池内分布，
    # 事件消耗归零。
    _pool_evt = False
    _lig_intake: dict[str, float] = {}
    for _anchor, _members in (T.pool_members or {}).items():
        _p0 = sum(initial.get(_m, 0.0) for _m in _members)
        _p1 = sum(ledger.get(_m, 0.0) for _m in _members)
        if abs(_p1 - _p0) >= ANN_MIN_EXTENT:
            _pool_evt = True
        for _m in _members:
            if _m == _anchor:
                continue
            _b = T.beta_by_complex.get(_m)
            if _b:
                _d = ledger.get(_m, 0.0) - initial.get(_m, 0.0)
                if abs(_d) > 1e-12:
                    _lig = _b["ligand"]
                    _lig_intake[_lig] = _lig_intake.get(_lig, 0.0) + _b["nu"] * _d
    changed = (any(
                   (e["mol"] - max(_lig_intake.get(e["name"], 0.0), 0.0))
                   / _nu_max.get(e["name"], 1) >= ANN_MIN_EXTENT
                   for e in consumption
                   if e["name"] not in (T.pools or {}))
               or _pool_evt
               or any(e["mol"] >= ANN_MIN_EXTENT
                      for e in production if e["name"] in T.solids)
               or any(st["kind"] == "neutralize" and st.get("extent", 0) >= ANN_MIN_EXTENT
                      for st in steps)
               or any(e["mol"] >= ANN_MIN_EXTENT for e in esc_list))
    # B4 狭义化学反应：净显著步骤中存在 redox/中和/非常规候选，或某净步骤的
    # 反应物横跨 ≥2 个投料来源（复分解、沉淀、配位溶解等）；纯溶解/电离/
    # 自互变（proton/dissolve/complex/decomplex 且单一来源）只是形态变化。
    # 例外：单一来源的 precip 步骤 = 盐类水解成淀（Fe^{3+}+3H2O ⇌ Fe(OH)3+3H+
    # 生成新相新物质，教材标准可逆反应）——计为化学反应；单一来源的 proton
    # 纯形态分布（Na2CO3 溶液中 HCO3- 分率）无新相生成，仍视为形态变化
    _SPEC_KINDS = {"proton", "dissolve", "complex", "decomplex"}
    chemical = False
    for st in steps:   # 规范化阶段的酸碱中和（steps0，无 key 入账）
        if st["kind"] == "neutralize" and st.get("extent", 0) >= ANN_MIN_EXTENT:
            chemical = True
    for net, kind, feeds in chem_net.values():
        if abs(net) < ANN_MIN_EXTENT:
            continue
        if kind not in _SPEC_KINDS or len(feeds) >= 2:
            chemical = True
            break
    # degree 程度整数：2=完全反应 / 1=可逆（部分）反应 / 0=难反应或未反应
    if not main_steps:
        conv0 = max((st.get("abs_conv", st.get("conversion", 0.0)) for st in steps),
                    default=0.0)
        if conv0 >= DEGREE_COMPLETE:
            degree = 2
        elif conv0 >= DEGREE_PARTIAL:
            degree = 1
        else:
            degree = 0
    else:
        lead = max(main_steps, key=lambda st: st["extent"])   # 主步：程度最大者
        conv = lead.get("abs_conv", lead["conversion"])
        # 多通道补全：主步转化率略低于阈值，但某初始反应物跨通道总体耗尽
        # （≥DEGREE_COMPLETE）时仍判 2（反应已被驱动到底）
        if conv < DEGREE_COMPLETE and conv >= DEGREE_PARTIAL:
            exhausted = any(
                initial.get(s, 0.0) > ANN_MIN_EXTENT
                and (initial[s] - ledger.get(s, 0.0)) / initial[s] >= DEGREE_COMPLETE
                for s in initial if s != WATER)
            if exhausted:
                conv = 1.0
        degree = 2 if conv >= DEGREE_COMPLETE else (
            1 if conv >= DEGREE_PARTIAL else 0)
    if not changed:
        degree = 0   # 一致性：无显著净变化即未反应，程度与 changed 同口径
    # 初态：post-normalize 的 ledger 副本（强电解质已完成电离、气体如 SO3
    # 已与水反应为酸、酸碱中和已记账），加入 H+/OH-（来自 H_excess_initial）
    # 与 H2O（溶剂，不入）。这是反应前的化学初态。
    initial_dict: dict[str, float] = {s: round(m, 6) for s, m in initial.items()
                                      if s != WATER and not s.startswith("__")
                                      and m > 1e-6}
    if H_excess0 > 1e-6:
        initial_dict[H_ION] = round(H_excess0, 6)
    elif H_excess0 < -1e-6:
        initial_dict["OH^-"] = round(-H_excess0, 6)
    # 精确净差（只给方程式装配吃；公开 consumption/production 仍是 1e-6
    # 报告口径）。**为什么必须分开**：迹量反应的整条净差可以小于报告阈
    # （P10：Fe(OH)₃ 1.172e-6 mol），报告口径的 `> 1e-6` 阈会把真实项
    # （NH₃ 4e-7）整块删掉 ⟹ 净差本身电荷不平 ⟹ **任何**呈现都无法配平，
    # 最后印出 `4H^+ + Fe(OH)_3 -> Fe^{3+}`（电荷 +4≠+3）这种伪方程。
    # 方程式是化学事实的呈现，必须吃精确量；目录式摘要才吃报告口径。
    _nx_c = {s: initial[s] - ledger.get(s, 0.0) for s in initial
             if s != WATER and not s.startswith("__")
             and initial[s] - ledger.get(s, 0.0) > 0.0}
    _nx_p = {s: ledger.get(s, 0.0) - initial.get(s, 0.0) for s in ledger
             if s != WATER and not s.startswith("__")
             and ledger.get(s, 0.0) - initial.get(s, 0.0) > 0.0}
    for s, m in escaped.items():          # 逸出气体计入净生成（同 esc_list 口径）
        if m > 0.0:
            _nx_p[s] = _nx_p.get(s, 0.0) + m
    net_exact = {"c": _nx_c, "p": _nx_p, "He_i": H_excess0, "He_f": H_excess}
    return {
        "changed": changed, "reacted": chemical and changed,
        "degree": degree, "annotations": annotations,
        "consumption": consumption, "production": production, "final": final,
        "net_exact": net_exact,
        "initial": [{"name": s, "mol": m} for s, m in initial_dict.items()],
        "escaped": esc_list,
        "ionize": _ionize_map(
            [e["name"] for e in consumption + production], T),
        "pool_map": dict(T.pools) if T.pools else {},
        "pool_members": {a: list(m) for a, m in T.pool_members.items()},
        "pool_ligands": {a: sorted(l) for a, l in T.pool_ligands.items()},
        # 池 complex 的 (ligand, nu)——配体摄入抵扣用（净方程折叠的元素
        # 守恒配套：[ZnCl]⁺ 折到 Zn²⁺ 时其 Cl 摄入量从 Cl⁻ 消耗中扣除）
        "pool_nu": {e["complex"]: (e["ligand"], e["nu"])
                    for e in T.beta if e["complex"] in T.pools},
        "steps": steps, "unknown": unknown,
        "final_pH": round(presentation_pH(ledger, H_excess, V, T, T_K), 2),
        "H_excess": round(H_excess, 6),
        "H_excess_initial": round(H_excess0, 6),
        "override": None,
        # 条件元数据透传（独立温度模块事后消费；不影响求解）
        "cond": {"V_L": round(V, 4), "T_K": T_K,
                 "isothermal": bool(cond.get("isothermal", False)),
                 "kinetics": bool(cond.get("kinetics", True)),
                 "gas_escape": bool(cond.get("gas_escape", True))},
    }


def _fmt(nu: float) -> str:
    """方程式系数格式化：1 省略；整数显示整数；非整数定点小数（禁用科学
    计数法——'1e+05' 这类输出会被方程式解析器当成物种名的一部分）。"""
    if nu == 1:
        return ""
    if nu == int(nu):
        return str(int(nu))
    return f"{nu:.6f}".rstrip("0").rstrip(".")
