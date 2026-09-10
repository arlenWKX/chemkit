"""chemkit.data：数据表加载与 Hess 运行时派生。

数据库最小化原则：JSON 只存规范可查的标准常数（E0/pKa/pKsp/logβ/ΔHf°），
派生量（反应焓 dH、半反应配平系数等）一律在加载时按 Hess 定律计算，不落地。
默认数据目录为包内 ./data（含 tests.json 测试用例库）。
"""
from __future__ import annotations
import json
import math
from math import gcd
import os
from dataclasses import dataclass, field

from .core import elements_of, charge_of, K_NERNST_298

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


# 单质标准态（ΔHf° ≡ 0，无需入库）
_ELEMENTS = {"Ag", "Al", "Au", "Ba", "Be", "C", "Ca", "Cd", "Co", "Cr", "Cs",
             "Cu", "Fe", "Ga", "Hg", "K", "Li", "Mg", "Mn", "Na", "Ni", "Pb",
             "Pt", "Rb", "S", "Sc", "Se", "Sn", "Sr", "Ti", "V", "Zn"}


@dataclass
class Tables:
    couples: list[dict]
    pka: list[dict]
    ksp: list[dict]
    beta: list[dict]
    conc: dict[str, float]
    ex: dict[str, dict]
    overrides: list[dict]
    thermo: dict[str, float] = field(default_factory=dict)   # 物种 ΔHf° (kJ/mol)
    # 派生
    cations: set[str] = field(default_factory=set)
    anions: set[str] = field(default_factory=set)
    ksp_by_pair: dict[tuple[str, str], dict] = field(default_factory=dict)
    ksp_by_solid: dict[str, dict] = field(default_factory=dict)
    beta_by_pair: dict[tuple[str, str], dict] = field(default_factory=dict)
    beta_by_complex: dict[str, dict] = field(default_factory=dict)
    pka_acid: dict[str, list[dict]] = field(default_factory=dict)   # acid -> entries
    pka_base: dict[str, list[dict]] = field(default_factory=dict)   # base -> entries
    solids: set[str] = field(default_factory=set)
    gases: set[str] = field(default_factory=set)
    # Henry 定律常数 H（mol/(L·kPa)，298K）：p = c/H。来源：NIST WebBook /
    # Sander, Atmos. Chem. Phys. 15, 4399 (2015) 汇编值（25℃，换算自 M/atm）
    henry: dict[str, float] = field(default_factory=lambda: {
        "CO_2": 3.4e-4,   # 0.034 M/atm（纯物理溶解，不含水合/电离）
        "O_2": 1.3e-5,    # 0.0013
        "H_2": 7.8e-6,    # 0.00078
        "N_2": 6.5e-6,    # 0.00065
        "NH_3": 0.57,     # ~57 M/atm（极易溶）
        "H_2S": 1.0e-3,   # 0.10
        "SO_2": 1.2e-2,   # 1.2
        "Cl_2": 6.2e-4,   # 0.062（物理溶解部分）
        "NO": 1.9e-5,     # 0.0019
        "N_2O": 2.5e-4,   # 0.025
        "CH_4": 1.4e-5,   # 0.0014
        "C_2H_2": 4.1e-4, # 0.041
        "PH_3": 8.0e-5,   # 0.008
        "H_2Se": 8.4e-4,  # 0.084
    })
    # 气体溶解焓 ΔH_sol（kJ/mol，放热为负；来源同 henry 汇编）。
    # H(T) = H_298 · exp(−ΔH_sol/R·(1/T − 1/298.15))：升温溶解度下降。
    henry_dH: dict[str, float] = field(default_factory=lambda: {
        "CO_2": -19.4, "O_2": -12.1, "H_2": -4.2, "N_2": -10.6,
        "NH_3": -30.6, "H_2S": -18.7, "SO_2": -25.5, "Cl_2": -19.2,
        "NO": -11.2, "N_2O": -16.5, "CH_4": -13.0, "C_2H_2": -17.5,
        "PH_3": -14.6, "H_2Se": -18.0,
    })
    redox_species: set[str] = field(default_factory=set)   # 所有电对的 ox/red 物种
    redox_ox_E: dict[str, float] = field(default_factory=dict)    # 作为氧化剂的最高 E0
    redox_red_E: dict[str, float] = field(default_factory=dict)   # 作为还原剂的最低 E0
    # 弱配形态池（v0.4.0）：累积 logβ 全部 < _POOL_LOG_BETA 的 (center, ligand)
    # 族是"部分混合"形态（NR19 盐混合 78/22 再分布——B4 判据既有认定：单源
    # complex/decomplex = 形态变化非化学反应），区别于"完全转化"事件（银氨
    # β2=7.2、PbCl3- 2.0）。事件口径三层（A1 changed / 净方程 / 断言量值）
    # 池折叠到 anchor（center），物种账本层如实逐物种（DB13 原则）。
    pools: dict[str, str] = field(default_factory=dict)           # 成员 → anchor
    pool_members: dict[str, list[str]] = field(default_factory=dict)  # anchor → 排序成员
    pool_ligands: dict[str, set[str]] = field(default_factory=dict)   # anchor → 弱族游离配体集合


def _load(name: str):
    with open(os.path.join(DATA_DIR, name), encoding="utf-8") as f:
        return json.load(f)


_R_GAS = 8.314e-3   # kJ/(mol·K)
_POOL_LOG_BETA = 1.0   # 弱配形态池阈值：族内累积 logβ 全部低于此值 ⟹ 池折叠


def henry_of(T: "Tables", gas: str, T_K: float) -> float | None:
    """Henry 常数 H（mol/(L·kPa)）在 T_K 的值；van't Hoff 温度修正。
    无数据的物种返回 None（调用方回退为完全逸出处理）。
    非 298.15 K 的修正含 exp 求值——按 (gas, T_K) 缓存于数据表
    （值纯函数确定，缓存与逐次计算一致；热路径 S_of 每次求值都经过这里）。"""
    H298 = T.henry.get(gas)
    if H298 is None:
        return None
    dH = T.henry_dH.get(gas)
    if dH is None or T_K == 298.15:
        return H298
    cache = T.__dict__.get("_henry_c")
    if cache is None:
        cache = T._henry_c = {}
    k = (gas, T_K)
    v = cache.get(k)
    if v is None:
        v = H298 * math.exp(-dH / _R_GAS * (1.0 / T_K - 1.0 / 298.15))
        cache[k] = v
    return v


def _expand_kinetics(e: dict) -> dict:
    """统一动力学字段 -> 引擎内部键（数据层单一事实源，引擎语义不变）。

    couples.json 中每个电对至多一个 "kinetics" 对象：
      ox_closed        作氧化剂方向动力学封闭（原 ox_inert；如 SO4^2- 稀溶液无氧化性）
      red_closed       作还原剂方向封闭（原 slow_as_reductant；如 Mn2+->MnO4-）
      below_T          低于该温度双向封闭（原 slow_below / gate.T_min）
      below_T_only_red 温度闸门只对指定还原剂生效（原 gate.only_vs_red）
      closed_with_red / closed_with_ox / closed_except_red
                       与特定配对封闭（原 slow_with_red/ox、slow_except_red）
      red_pH_min       还原通道 pH 下限
      rev_gate         逆向闸门（产物选择性）
    均为"无限时间也不发生"的化学硬事实，与热力学无关。
    """
    k = e.pop("kinetics", None)
    if not k:
        return e
    if k.get("ox_closed"):
        e["ox_inert"] = True
    if k.get("red_closed"):
        e["slow_as_reductant"] = True
    if "below_T" in k:
        if k.get("below_T_only_red"):
            e.setdefault("gate", {})["T_min"] = k["below_T"]
            e["gate"]["only_vs_red"] = k["below_T_only_red"]
        else:
            e["slow_below"] = k["below_T"]
    if k.get("closed_with_red"):
        e["slow_with_red"] = k["closed_with_red"]
    if k.get("closed_with_ox"):
        e["slow_with_ox"] = k["closed_with_ox"]
    if k.get("closed_except_red"):
        e["slow_except_red"] = k["closed_except_red"]
    if "red_pH_min" in k:
        e["red_pH_min"] = k["red_pH_min"]
    if "h2o_red_oh_min" in k:
        e["h2o_red_oh_min"] = k["h2o_red_oh_min"]
    if "ox_pH_max" in k:
        e["ox_pH_max"] = k["ox_pH_max"]
    if k.get("rev_gate"):
        e["rev_gate"] = k["rev_gate"]
    return e


# 数据表缓存：按 data_dir 键缓存，确保同一数据目录在程序生命周期内只加载一次。
# import chemkit 时 system.TABLES = load_tables() 首次加载；后续任何
# load_tables() 调用（如 testsuit.main）直接返回缓存实例，零磁盘 I/O。
_TABLES_CACHE: dict[str, Tables] = {}


def load_tables(data_dir: str | None = None) -> Tables:
    global DATA_DIR
    if data_dir:
        DATA_DIR = data_dir
    cache_key = DATA_DIR
    cached = _TABLES_CACHE.get(cache_key)
    if cached is not None:
        return cached
    ex = _load("substance_ex.json")
    # 浓溶液阈值（mol/L）作为物质属性并入 substance_ex.json 的 "conc_M" 字段，
    # 此处派生（原独立 conc.json 已废除——4 个阈值单建文件不值，且与
    # conc_forms 分居两处易改漏）
    conc = {name: e["conc_M"] for name, e in ex.items() if "conc_M" in e}
    t = Tables(
        couples=[_expand_kinetics(e) for e in _load("couples.json")],
        pka=_load("pka.json"),
        ksp=_load("ksp.json"),
        beta=_load("beta.json"),
        conc=conc,
        ex=ex,
        overrides=_load("overrides.json"),
        thermo=_load("thermo.json"),
    )
    for e in t.ksp:
        cat, an = e["pair"]
        t.ksp_by_pair[(cat, an)] = e
        t.ksp_by_solid[e["solid"]] = e
        t.cations.add(cat)
        t.anions.add(an)
        t.solids.add(e["solid"])
    for e in t.beta:
        t.beta_by_pair[(e["center"], e["ligand"])] = e
        t.beta_by_complex[e["complex"]] = e
        (t.cations if charge_of(e["center"]) > 0 else t.anions).add(e["center"])
        (t.cations if charge_of(e["ligand"]) > 0 else t.anions).add(e["ligand"])
        (t.cations if charge_of(e["complex"]) > 0 else t.anions).add(e["complex"])
    # ---- 弱配形态池派生（v0.4.0）----
    # 折叠条件（两要件同时满足）：
    #   ① 族内累积 logβ 全部 < 1.0（弱配 = 部分混合形态）；
    #   ② 族条目数 ≥ 2（逐级 β 谱 = 形态分布族，"同一份溶解盐"语义成立）。
    #   单条目族（Cu-Cl/Cu-Br 的 β4-only、PbCl₃⁻、[AgI₂]⁻）是条件性事件
    #   配合物（浓介质生成/稀介质分解），教学口径锁定为事件呈现
    #   （FeCl₃+Cu 刻蚀 → [CuCl₄]²⁻ 族六例），不折叠。
    # 以 center 为锚（中心金属守恒），池成员 = center + 该 center 全部
    # 弱族的 complexes。排序固定序（跨进程确定性纪律，v0.3.9 教训）。
    # 边界例：Ni-SCN/Zn-Br/Zn-SCN 的 logβ4=1.0 恰在阈值上不折叠。
    fam: dict[tuple[str, str], list[dict]] = {}
    for e in t.beta:
        fam.setdefault((e["center"], e["ligand"]), []).append(e)
    for (center, ligand), members in sorted(fam.items()):
        if (len(members) >= 2
                and all(float(m["logb"]) < _POOL_LOG_BETA for m in members)):
            t.pool_members.setdefault(center, [center])
            t.pool_ligands.setdefault(center, set()).add(ligand)
            for m in members:
                if m["complex"] not in t.pool_members[center]:
                    t.pool_members[center].append(m["complex"])
    for anchor, members in t.pool_members.items():
        t.pool_members[anchor] = sorted(set(members))
        for sp in t.pool_members[anchor]:
            t.pools[sp] = anchor
    for e in t.pka:
        t.pka_acid.setdefault(e["acid"], []).append(e)
        t.pka_base.setdefault(e["base"], []).append(e)
        for sp in (e["acid"], e["base"]):
            q = charge_of(sp)
            (t.cations if q > 0 else t.anions if q < 0 else set()).add(sp)
    # 晶格阴离子电对派生（热力学循环，非新数据）：电对 (ox, red=A) 的
    # red 为游离阴离子且 A 参与 Ksp 固相 M_xA_y 时，派生 (ox, red=solid)：
    # 氧化晶格中的 A 须先拆格子，aA = (Ksp/aM^x)^(1/y) 代入 Nernst 得
    # E0' = E0 + k·ν_A·pKsp/(n·y)（对离子活度 1；ν_A 为半反应中阴离子
    # 系数，如 I2/I- 为 2）。例：S/PbS（油画变黑 PbS + H2O2 -> 白色
    # PbSO4 的修复化学）、S/CuS（CuS 溶于热硝酸析 S）、S/HgS（E0'≈1.06，
    # HNO3/NO 0.96 不够——HgS 不溶于硝酸只溶王水的热力学根源）。
    # 只派生负氧族晶格（S2-）：硫化物矿石的氧化溶解是独立的化学现象类
    # （湿法冶金、定性分析）；卤素晶格电对（E0' 普遍 >1.2V）只会在卤素
    # 歧化/归中网络上叠加边际幻影通道（Cl2+Pb2+->PbCl2+ClO3- 之类），
    # 无对应实验化学，不派生。
    _an_solids: dict = {}
    for _e in t.ksp:
        _cat, _an = _e["pair"]
        if _cat != "H^+" and _an == "S^{2-}":
            _an_solids.setdefault(_an, []).append(_e)
    _existing = {(c["ox"], c["red"]) for c in t.couples}
    _extra = []
    for _c in t.couples:
        _red = _c["red"]
        if charge_of(_red) >= 0 or _red not in _an_solids:
            continue
        _el = elements_of(_red)
        if len(_el) != 1:
            continue   # 只处理单原子阴离子（S2-/卤离子等）
        _ela = next(iter(_el))
        _nu_a = elements_of(_c["ox"]).get(_ela, 0) / _el[_ela]
        if _nu_a <= 0:
            continue
        for _e in _an_solids[_red]:
            _solid = _e["solid"]
            if (_c["ox"], _solid) in _existing:
                continue
            _y = elements_of(_solid).get(_ela, 0)
            if _y <= 0:
                continue
            _extra.append({
                "ox": _c["ox"], "red": _solid, "n": _c["n"],
                "E0": round(_c["E0"] + K_NERNST_298 * _nu_a * _e["pKsp"]
                            / (_c["n"] * _y), 4),
                "note": f"晶格阴离子电对派生：{_c['ox']}/{_red} "
                        f"E0 + k·ν·pKsp({_solid})/(n·y)",
                "derived_from": _red})
            _existing.add((_c["ox"], _solid))
    t.couples.extend(_extra)
    for e in t.couples:
        for sp in (e["ox"], e["red"]):
            t.redox_species.add(sp)
            q = charge_of(sp)
            (t.cations if q > 0 else t.anions if q < 0 else set()).add(sp)
        # 物种作为氧化剂/还原剂参与的电势范围（惰性实现增益判定的近似配对用）
        t.redox_ox_E[e["ox"]] = max(t.redox_ox_E.get(e["ox"], -1e9), e["E0"])
        t.redox_red_E[e["red"]] = min(t.redox_red_E.get(e["red"], 1e9), e["E0"])
    # 补充：在 KSP 中无格子但参与拆解的离子
    t.anions.update(["NO_3^-", "ClO_4^-", "MnO_4^-", "ClO^-", "SO_3^{2-}", "HCO_3^-",
                     "HSO_3^-", "HS^-", "NO_2^-", "CN^-", "F^-", "Br^-", "I^-", "OH^-",
                     "BrO^-",
                     "CO_3^{2-}", "S^{2-}", "PO_4^{3-}", "HPO_4^{2-}", "H_2PO_4^-",
                     "SiO_3^{2-}", "H_2PO_2^-", "MnO_4^{2-}", "C_2^{2-}", "N^{3-}",
                     "CH_3COO^-", "H^-"])
    t.cations.update(["Na^+", "K^+", "NH_4^+", "H^+", "Li^+", "Sr^{2+}", "Ni^{2+}",
                      "Mn^{2+}", "Cr^{3+}", "Sn^{2+}", "Sc^{3+}", "Ti^{3+}", "V^{3+}",
                      "Au^{3+}"])
    for name, e in t.ex.items():
        if e.get("form") == "solid":
            t.solids.add(name)
        if e.get("form") == "gas":
            t.gases.add(name)
    # 数据源目录快照（DATA_DIR 全局可被后续 load_tables(data_dir=...)
    # 改写，表实例须记住自己的来源）
    t.src_dir = os.path.abspath(DATA_DIR)
    _compute_dH(t)
    _TABLES_CACHE[cache_key] = t
    return t


# ---------- van't Hoff 反应焓的运行时派生（Hess 定律，dH 单位 kJ/mol） ----------

def _half_balance(ox: str, red: str):
    """通用还原半反应配平：k·ox + h·H+ + ne·e- -> m·red + w·H2O。
    骨架元素依次尝试非 H/O 元素、O、H（O2/H2O 等全 H/O 体系靠后两者）。
    返回 (k, h, m, w, ne) 或 None（元素单边出现的非标准形，走配体释放路径）。
    w < 0 表示 H2O 在反应物侧（如 Fe3O4 + 2H2O + 2H+ + 2e- -> 3Fe(OH)2），
    化学上仍为有效半反应，正常返回。"""
    eo, er = elements_of(ox), elements_of(red)
    others = [X for X in set(eo) | set(er) if X not in ("H", "O")]
    for sk in others + ["O", "H"]:
        co0, cr0 = eo.get(sk, 0), er.get(sk, 0)
        if co0 == 0 or cr0 == 0:
            continue
        for k in (1, 2, 3, 4):
            if (k * co0) % cr0:
                continue
            m = (k * co0) // cr0
            if any(k * eo.get(X, 0) != m * er.get(X, 0) for X in others):
                break          # 骨架选择不当（元素单边出现），换下一骨架
            w = k * eo.get("O", 0) - m * er.get("O", 0)
            h = m * er.get("H", 0) + 2 * w - k * eo.get("H", 0)
            ne = k * charge_of(ox) + h - m * charge_of(red)
            if h < 0 or ne <= 0:
                break
            return k, h, m, w, ne
    return None


def _couple_dH(t, e):
    """电对还原半反应 ΔH 与电子数：标准半反应 -> 配体释放 -> SHE 约定。"""
    ox, red = e["ox"], e["red"]
    if (ox, red) == ("H^+", "H_2"):
        return 0.0, 2          # SHE：ΔG/ΔH/ΔS 全温度按定义为零
    bal = _half_balance(ox, red)
    if bal is not None:
        k, _h, m, w, ne = bal
        vals = [_dhf(t, red), _dhf(t, ox)]
        if all(v is not None for v in vals):
            return m * vals[0] + w * t.thermo["H_2O"] - k * vals[1], ne
        return None, None
    # 配体释放形：complex + ne·e- -> metal + nu·ligand（[Ag(NH3)2]+/Ag、
    # [AuCl4]-/Au 等）。配体数取自配合物化学式而非 beta 表的 nu（后者是
    # logβ 级数，可能与式中配体数不一致——[PtCl6]2- 的 beta nu=4 而式含
    # 6 Cl）；ne 由电荷平衡：ne = q(ox) − q(red) − nu·q(ligand)
    b = t.beta_by_complex.get(ox)
    if (b is not None and red in _ELEMENTS
            and elements_of(b["center"]).get(red)):
        lel = next((X for X in elements_of(b["ligand"]) if X != "H"), None)
        nu = (elements_of(ox).get(lel, 0) - elements_of(b["center"]).get(lel, 0)
              if lel else 0)
        vals = [_dhf(t, red), _dhf(t, b["ligand"]), _dhf(t, ox)]
        if nu > 0 and all(v is not None for v in vals):
            ne = charge_of(ox) - charge_of(red) - nu * charge_of(b["ligand"])
            if ne > 0:
                return vals[0] + nu * vals[1] - vals[2], ne
    return None, None


def _dhf(t, sp):
    return dhf_of(t, sp)


def dhf_of(t, sp: str) -> float | None:
    """物种标准生成焓 ΔHf°（kJ/mol）——**单一事实源**入口。

    数据面：thermo.json（裸名 = 引擎账本态：水溶/凝聚/永久气体；
    "X(g)" 后缀键 = 逸出气相态，仅供热模块的气体计账）。
    单质标准态（_ELEMENTS 集合）恒为 0；缺数据返回 None（调用方
    各自回退：van't Hoff 缺 dH 时用 Nernst ΔS≈0，热分析报
    宁缺毋假）。thermo 模块与本函数共用同一份数据，不再自带副本。"""
    v = t.thermo.get(sp)
    if v is not None:
        return v
    return 0.0 if sp in _ELEMENTS else None


def _compute_dH(t):
    """四表反应焓运行时派生（不写回数据库）：
    couples 还原半反应（dH/dH_n）；pka 酸式电离 acid + w·H2O -> base + n·H+
    （水合物种不可省略：CO2 + H2O -> HCO3- + H+ 与 H2CO3 差一个 H2O）；
    ksp 溶解；beta 配位。缺 ΔHf 的条目不设 dH，引擎回退旧行为。"""
    for e in t.couples:
        v, ne = _couple_dH(t, e)
        if v is not None:
            e["dH"], e["dH_n"] = round(v, 1), ne
    for e in t.pka:
        a, b = _dhf(t, e["acid"]), _dhf(t, e["base"])
        if a is None or b is None:
            continue
        n = e.get("n", 1)
        w = (elements_of(e["base"]).get("H", 0) + n
             - elements_of(e["acid"]).get("H", 0)) / 2
        if w == int(w):
            e["dH"] = round(b - a - int(w) * t.thermo["H_2O"], 1)
    for e in t.ksp:
        cat, an = e["pair"]
        vals = [_dhf(t, e["solid"]), _dhf(t, cat), _dhf(t, an)]
        if all(v is not None for v in vals):
            qc, qa = charge_of(cat), -charge_of(an)
            g = gcd(qc, qa)
            e["dH"] = round(qa // g * vals[1] + qc // g * vals[2] - vals[0], 1)
    for e in t.beta:
        vals = [_dhf(t, e["complex"]), _dhf(t, e["center"]), _dhf(t, e["ligand"])]
        if all(v is not None for v in vals):
            e["dH"] = round(vals[0] - vals[1] - e["nu"] * vals[2], 1)
