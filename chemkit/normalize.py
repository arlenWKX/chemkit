"""chemkit.normalize：投料规范化（投料 → 引擎账本）。

职能：
    · normalize：投料列表 → (ledger, H_excess, steps0, unknown)
      —— 可溶性强电解质电离（ex 表/Ksp 溶度分级）、浓酸分子态曲线
      （_mol_fraction）、强酸强碱中和入 steps0、未知物种报告
    · _split_salt/_split_acid/_salt_ksp_cell：电离拆解与溶度判级
    · _ionize_map：方程式呈现用的电离映射（Reaction.ionize 源）

依赖 candidates（常数）/core/data。
"""
from __future__ import annotations
from math import gcd
from .core import elements_of, charge_of, PKW_298
from .data import Tables
from .candidates import (WATER, H_ION, X_MIN, WATER_MOL_PER_L,
                         STRONG_MOLECULAR_ACIDS, _split_salt_general,
                         _ELEMENTS, _MONO_ANION_Q)

# ========================================================== 规范化

def _mol_fraction(name: str, c: float, T) -> float:
    """强酸分子分数（水活度驱动的协同缔合代理曲线）：f = c^p/(c_half^p + c^p)。
    单一表观 Ka 无法同时拟合稀区（≤4M 分子分数<1%）与浓区（≥10M 显著
    分子化）——浓区分子化本质是活度系数/水活度效应而非稀溶液 Ka。
    两个常数均有物理含义：c_half = 半分子化浓度（HNO3 12M、H2SO4 12M，
    与 Raman 形态数据量级一致），p = 协同度。替代原 conc_M 二元开关。
    无 spec 字段的酸（HCl 等——浓盐酸实际仍近全电离）恒为 0。"""
    spec = (T.ex.get(name) or {}).get("spec")
    if spec is None or c <= 0.0:
        return 0.0
    cp = c ** spec["p"]
    return cp / (spec["c_half"] ** spec["p"] + cp)


def normalize(substances: list[dict], cond: dict, T,
              origins: dict | None = None) -> tuple[dict, float, list, list[str]]:
    """origins（可选）：给定一个 dict 时，按投料逐项记录物种来源
    {物种: {投料名, ...}}（H_ION 记录向质子账本贡献过 H+/OH- 的投料），
    供 reacted（狭义化学反应）的"跨投料相互作用"判据使用。"""
    V = cond["V_L"]
    ledger: dict[str, float] = {}
    unknown: list[str] = []
    pos = 0.0   # 强酸 H+ 贡献
    neg = 0.0   # 强碱 OH- 贡献

    def add_species(name: str, mol: float):
        nonlocal pos, neg
        if name == H_ION:
            pos += mol
        elif name == "OH^-":
            neg += mol
        else:
            ledger[name] = ledger.get(name, 0.0) + mol

    for item in substances:
        name, mol = item["name"], float(item["mol"])
        if mol <= 0:
            continue
        _pre = (dict(ledger), pos, neg) if origins is not None else None
        entry = T.ex.get(name)
        if charge_of(name) != 0 and all(e in _ELEMENTS for e in elements_of(name)):
            # 任何元素合法的带电物种直接入账：离子的存在不依赖注册
            # （注册只是反应数据入口；未注册离子为旁观，如 Fe^{6+}）
            add_species(name, mol)
        elif entry is not None:
            if "conc_forms" in entry:
                # 连续形态分布：分子分数由表观电离平衡给出（不再按 conc_M
                # 一刀切）；电离部分仍按全电离约定拆为离子
                f_mol = _mol_fraction(name, mol / V, T)
                m_mol = mol * f_mol
                if m_mol > 0.0:
                    ledger[name] = ledger.get(name, 0.0) + m_mol
                if mol - m_mol > 0.0:
                    for sp, nu in _split_acid(name, T).items():
                        add_species(sp, nu * (mol - m_mol))
                # 影子库存：该酸（仅本酸，不含盐类同离子）的分子+离子总量，
                # 供逐轮再平衡求解目标分子池（输出时过滤 "__" 前缀）
                ledger[f"__tot_{name}"] = ledger.get(f"__tot_{name}", 0.0) + mol
            elif entry["form"] == "ions":
                for ion, nu in _ion_coeffs(name, entry["ions"], T).items():
                    add_species(ion, nu * mol)
            else:  # molecule / solid / gas 原样入账
                ledger[name] = ledger.get(name, 0.0) + mol
        elif charge_of(name) == 0:
            parts = _split_acid(name, T) or _split_salt(name, T)
            if parts is None:
                unknown.append(name)
            elif parts == "solid":
                ledger[name] = ledger.get(name, 0.0) + mol
            else:
                # 微溶（slight）盐的溶解度判定：全拆浓度下离子积已超 Ksp，
                # 说明投料主要以固相存在——保持固相入账，溶解度交给
                # dissolve 候选求平衡；否则"拆开即沉淀"会产生幻影循环
                # （ consumed 报告从未投入过的离子、reacted 误判）
                ks = _salt_ksp_cell(parts, T)
                if ks is not None:
                    cell, cat, an, m, n = ks
                    if cell.get("slight"):
                        c_cat, c_an = m * mol / V, n * mol / V
                        if c_cat ** m * c_an ** n > 10.0 ** (-cell["pKsp"]):
                            ledger[name] = ledger.get(name, 0.0) + mol
                            parts = "solid"
                if parts != "solid":
                    for sp, nu in parts.items():
                        add_species(sp, nu * mol)
        else:
            unknown.append(name)
        if origins is not None and _pre is not None:
            pre_led, pre_pos, pre_neg = _pre
            for sp, m_new in ledger.items():
                if sp.startswith("__") or sp == WATER:
                    continue
                if m_new - pre_led.get(sp, 0.0) > 0:
                    origins.setdefault(sp, set()).add(name)
            if pos - pre_pos > 1e-12 or neg - pre_neg > 1e-12:
                origins.setdefault(H_ION, set()).add(name)

    # 用户条件给定的初始酸碱性
    if cond.get("c_H"):
        pos += cond["c_H"] * V
    if cond.get("c_OH"):
        neg += cond["c_OH"] * V
    if cond.get("pH") is not None:
        pos += 10.0 ** (-cond["pH"]) * V

    steps0 = []
    neutralized = min(pos, neg)
    if neutralized > 1e-9:
        steps0.append({"kind": "neutralize", "equation": "H+ + OH- -> H2O",
                       "logK": PKW_298, "S": PKW_298, "extent": round(neutralized, 6),
                       "conversion": 1.0})
    H_excess = pos - neg
    ledger[WATER] = ledger.get(WATER, 0.0) + WATER_MOL_PER_L * V
    return ledger, H_excess, steps0, unknown


def _ion_coeffs(name: str, ions: list, T) -> dict:
    """`ex["ions"]`（**只列离子种类、不带计量**）→ {离子: 计量}。

    0.4.x 原实现按 1:1 记账（`add_species(ion, mol)`），对化学式非 1:1
    的盐**凭空破坏电荷与质量守恒**：La(NO₃)₃ → La³⁺ + 1NO₃⁻（净 +2）、
    Na₂S₂O₃ → Na⁺ + 1S₂O₃²⁻（净 −1）、K₂Cr₂O₇ → K⁺ + 1Cr₂O₇²⁻（净 −1）。
    `ex` 表里 34 条 form=ions 条目中 **17 条**如此（La/Ce/In/Zr/Pd 的
    硝酸盐与氯化物、Na₂MoO₄/Na₂WO₄/Na₂S₂O₃/K₂Cr₂O₇/Tl₂SO₄/H₂PtCl₆…）。

    修正口径（**以化学式为唯一依据**，与 `_split_salt` 同源）：
      1. 优先用 `_split_acid`/`_split_salt` 按化学式配平（它们保证原子
         与电荷双守恒，且已覆盖注册表 + 元素周期表级兜底）；
      2. 公式法失败时退回 `ions` 列表并按**整比缩放**使净电荷为零
         （`ions` 本身就是该盐的离子组成清单，缩放不改变组成）。
    """
    parts = _split_acid(name, T) or _split_salt(name, T)
    if isinstance(parts, dict) and parts:
        if sum(charge_of(sp) * nu for sp, nu in parts.items()) == 0:
            return parts
    # 退回 ex.ions + 整比缩放（正负电荷总量配平）
    pos: list = []
    neg: list = []
    for ion in ions:
        (pos if charge_of(ion) > 0 else neg).append(ion)
    q_pos = sum(charge_of(i) for i in pos)
    q_neg = -sum(charge_of(i) for i in neg)
    if not pos or not neg or q_pos <= 0 or q_neg <= 0:
        return {i: 1.0 for i in ions}
    g = gcd(int(q_pos), int(q_neg))
    k_pos, k_neg = q_neg // g, q_pos // g
    out: dict = {}
    for i in pos:
        out[i] = out.get(i, 0.0) + k_pos
    for i in neg:
        out[i] = out.get(i, 0.0) + k_neg
    return out


def _split_acid(name: str, T) -> dict | None:
    elems = elements_of(name)
    if "H" not in elems:
        return None
    n_h = elems["H"]
    for an in T.anions:
        q = -charge_of(an)
        if q <= 0 or n_h % q != 0:
            continue
        k = n_h // q
        target = dict(elements_of(an))
        target["H"] = target.get("H", 0) + k * q
        if target == elems:
            return {H_ION: n_h, an: k}
    return None


def _salt_ksp_cell(parts: dict, T):
    """识别盐拆解结果 {阳离子: m, 阴离子: n} 对应的 Ksp 条目，
    返回 (cell, cat, an, m, n)；非盐拆解（酸等）或无 Ksp 条目返回 None。"""
    if len(parts) != 2 or H_ION in parts:
        return None
    (s1, m), (s2, n) = parts.items()
    for cat, an, mm, nn in ((s1, s2, m, n), (s2, s1, n, m)):
        cell = T.ksp_by_pair.get((cat, an))
        if cell is not None:
            return cell, cat, an, mm, nn
    return None


def _ionize_map(species: list[str], T) -> dict[str, dict[str, float]]:
    """净方程式书写所需的"分子态 → 离子形"改写图（单一数据源）。

    哪些分子物种在净离子方程式中应改写为离子形（强电解质式书写）：
      - conc_forms 浓酸（HNO3/H2SO4 的分子形态）：HNO3 → H+ + NO3-；
      - ex form=="ions" 的物种（按 ex ions 清单，与 normalize 同口径）；
      - 高溶解度盐（Ksp pKsp ≤ 0.3，溶解度 ≳0.5M，如 NaHCO3 pKsp=-0.08）：
        化学习惯写自由离子（饱和析出仅是浓度记账）。真沉淀（pKsp 大）
        保持固相写法。改写让 Na+/K+ 等旁观离子自然抵消。
    金属 / 气体 / 弱酸分子（H2S、CH3COOH 等）与真固相不改写。
    返回 {分子名: {离子: 计量}}。"""
    out: dict[str, dict[str, float]] = {}
    for sp in species:
        if sp in out or charge_of(sp) != 0:
            continue
        entry = T.ex.get(sp)
        k = T.ksp_by_solid.get(sp)
        if entry is not None and "conc_forms" in entry:
            parts = _split_acid(sp, T)
            if parts:
                out[sp] = parts
            continue
        # ex ions 形态仅在无 Ksp 条目时才改写为离子形：有 Ksp 的盐是
        # 溶解度平衡体系（固相是真实形态，如 TlCl/K_2PtCl_4 的沉淀产物）
        if (entry is not None and entry.get("form") == "ions"
                and k is None):
            out[sp] = {ion: 1.0 for ion in entry["ions"]}
            continue
        if k is not None and k["pKsp"] <= 0.3:
            cat, an = k["pair"]
            qc, qa = charge_of(cat), -charge_of(an)
            g = gcd(qc, qa)
            out[sp] = {cat: qa // g, an: qc // g}
    return out


def _split_salt(name: str, T):
    elems = elements_of(name)
    for cat in sorted(T.cations, key=lambda c: -charge_of(c)):
        qc = charge_of(cat)
        if qc <= 0:
            continue
        ec = elements_of(cat)
        if any(elems.get(k, 0) < v for k, v in ec.items()):
            continue
        for an in T.anions:
            qa = -charge_of(an)
            if qa <= 0:
                continue
            ea = elements_of(an)
            for m in range(1, 5):
                if (m * qc) % qa != 0:
                    continue
                n = (m * qc) // qa
                tot: dict = {}
                for k, v in ec.items():
                    tot[k] = tot.get(k, 0) + m * v
                for k, v in ea.items():
                    tot[k] = tot.get(k, 0) + n * v
                if tot == elems:
                    cell = T.ksp_by_pair.get((cat, an))
                    if cell is not None and not cell.get("slight"):
                        return "solid"
                    return {cat: m, an: n}
    # 注册表未覆盖：元素周期表级兑底拆解（旁观离子，保留电荷/质量守恒）
    return _split_salt_general(name, T)

