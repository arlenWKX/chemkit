"""chemkit.equations：方程式装配（步骤聚合 → 净离子方程式组）。

职能：
    · 方程式解析（_parse_equation）与 H/O 原子配平（_balance_h_o）
    · H+ 正则形 → OH- 还原（_restore_oh，按终态 pH 酸碱性书写）
    · 系数有理化（_rationalize/_try_integer_snap：最简整数比 snap）
    · 步骤聚合与方程式构建（_collect_steps/_collapse_transient_
      intermediates/_build_equations：净程度对消、瞬态中间体折叠、
      trace 物种过滤）

被 system.Reaction 消费；依赖 core/data/candidates。
"""
from __future__ import annotations
import re
from fractions import Fraction
from math import gcd
from functools import reduce

from .core import elements_of, charge_of
from .data import Tables, load_tables
import os as _os

_TRACE_EQ = bool(_os.environ.get("CHEM_TRACE_EQ"))

from dataclasses import dataclass

from .candidates import WATER, H_ION, X_MIN, ANN_MIN_EXTENT

from .data import load_tables as _load_tables

_SOLIDS_C: frozenset | None = None


def _solid_set() -> frozenset:
    """固相物种全集（`Tables.solids` 的键），模块级缓存（痕量新相闸门用）。"""
    global _SOLIDS_C
    if _SOLIDS_C is None:
        _SOLIDS_C = frozenset(_load_tables().solids)
    return _SOLIDS_C

# 净方程守恒的绝对容差 = 报告量子：consumption/production 由 round(x, 6)
# 产出，每个系数误差 ≤5e-7；含物种数与电荷项余量后取 1e-5。
_EQ_TOL_ABS = 1e-5

# ---------- 数据表预加载（import chemkit 时执行，全程序仅一次磁盘 I/O）----------
# 方程式装配层的表依赖常量（固态酸 pKa 等静态提取）用默认表构建；
# system 层再导出同一单例（interfaces.TABLES）
TABLES: Tables = load_tables()


def default_tables() -> Tables:
    """返回模块级默认数据表（预加载单例）。"""
    return TABLES

# ============================================================ 方程式解析
# 引擎 step.equation 形如 '2NO_3^- + 3Cu + 8H^+ -> 2NO + 3Cu^{2+}'
# neutralize 步骤用旧式 'H+ + OH- -> H2O'，需归一化到引擎标准记法。

_SPECIES_NORMALIZE = {
    "H2O": "H_2O",
    "H+": "H^+",
    "OH-": "OH^-",
}

# 匹配 "系数+物种"，系数为可选的整数或小数；物种以非数字开头（字母/[/(
_TERM_RE = re.compile(r"^(\d+(?:\.\d+)?)(\D.*)$")
_SIDE_SEP = " + "

WATER = "H_2O"
H_ION = "H^+"
OH_ION = "OH^-"


def _parse_side(s: str) -> dict[str, float]:
    """解析 '2A + 3B' → {A: 2.0, B: 3.0}。"""
    result: dict[str, float] = {}
    for term in s.split(_SIDE_SEP):
        term = term.strip()
        if not term:
            continue
        m = _TERM_RE.match(term)
        if m:
            coef = float(m.group(1))
            species = m.group(2).strip()
        else:
            coef = 1.0
            species = term
        species = _SPECIES_NORMALIZE.get(species, species)
        result[species] = result.get(species, 0.0) + coef
    return result


def _parse_equation(eq) -> tuple[dict[str, float], dict[str, float]]:
    """解析 '2A + 3B -> 4C + D' → ({A:2, B:3}, {C:4, D:1})。

    接受 `Equation` 结构（v0.5.0：方程式先结构化后渲染——结构对象直接
    返回其系数，不走字符串回读）；接受任一箭头记号（`->`/`<=>`/`→`/`⇌`），
    因为箭头是**渲染细节**（结构里只有 `reversible` 布尔）。"""
    if isinstance(eq, Equation):
        return dict(eq.left), dict(eq.right)
    if not isinstance(eq, str):
        eq = str(eq)
    for a in (ARROW_FWD, ARROW_EQ, "\u2192", "\u21cc", "="):
        if a in eq:
            lhs, rhs = eq.split(a, 1)
            return _parse_side(lhs), _parse_side(rhs)
    return {}, {}


# ============================================================ H/O 原子配平
# 引擎 step.equation 过滤了 H2O（溶剂），离子方程式需要显式 H2O。

def _balance_h_o(consumed: dict[str, float], produced: dict[str, float]
                 ) -> tuple[dict[str, float], dict[str, float]]:
    """用 H2O（必要时 H+）配平 H、O 原子。

    引擎方程已用 H+ 正则化（OH- → H2O/H+），所以 H/O 不平衡仅由 H2O 被过滤
    导致。正常情况 h_diff = 2 × o_diff，纯加 H2O 即可；异常时退化为 H2O + H+。
    """
    h_lhs = sum(nu * elements_of(sp).get("H", 0) for sp, nu in consumed.items())
    h_rhs = sum(nu * elements_of(sp).get("H", 0) for sp, nu in produced.items())
    o_lhs = sum(nu * elements_of(sp).get("O", 0) for sp, nu in consumed.items())
    o_rhs = sum(nu * elements_of(sp).get("O", 0) for sp, nu in produced.items())

    o_diff = o_lhs - o_rhs      # 正：产物侧缺 O → 加 H2O 到产物
    h_diff = h_lhs - h_rhs      # 正：产物侧缺 H

    if abs(h_diff - 2 * o_diff) < 1e-6:
        # 纯 H2O 配平
        if o_diff > 1e-9:
            produced[WATER] = produced.get(WATER, 0.0) + o_diff
        elif o_diff < -1e-9:
            consumed[WATER] = consumed.get(WATER, 0.0) + (-o_diff)
    else:
        # 先用 H2O 配 O，再用 H+ 配 H（防御性兜底）
        if o_diff > 1e-9:
            produced[WATER] = produced.get(WATER, 0.0) + o_diff
            h_rhs += 2 * o_diff
        elif o_diff < -1e-9:
            consumed[WATER] = consumed.get(WATER, 0.0) + (-o_diff)
            h_lhs += 2 * (-o_diff)
        h_diff = h_lhs - h_rhs
        if h_diff > 1e-9:
            produced[H_ION] = produced.get(H_ION, 0.0) + h_diff
        elif h_diff < -1e-9:
            consumed[H_ION] = consumed.get(H_ION, 0.0) + (-h_diff)
    return consumed, produced


# ============================================================ H+ 正则形 → OH- 还原

def _restore_oh(consumed: dict[str, float], produced: dict[str, float],
                trace: float, allow: bool = True) -> tuple[dict[str, float], dict[str, float]]:
    """把引擎的 H+ 正则形还原为化学习惯的 OH- 写法，并抵消跨侧 H+/OH-。

    两步：
      1. OH- 还原：引擎把 OH- 记为 H2O(反应物) − H+(即产物 H+)，
         如 M + nH2O -> M(OH)n + nH+ 实为 M + nOH- -> M(OH)n。
         H2O 在反应物侧且 H+ 在产物侧且量匹配时，合并为 OH- 到反应物侧。
    只有"H2O 在反应物侧 + H+ 在产物侧"这一种挂侧可以还原为 OH-；
    跨侧 H+/OH- 抵消（H+ 左 OH- 右，或反之）在数学上不等价于生成
    H2O（净向量差 2H+ 或 2OH-），会破坏配平，故不做。

    v0.4.5 试行后撤回：曾按"终态 pH > 7 才算碱性"给还原加闸，想让
    FeCl₃+KSCN（pH 1.45）写成 `… + 3H2O -> … + Fe(OH)_3 + 3H^+`。实测
    **95 例翻红**——语料自身的呈现惯例是**一律**把 `H2O(左)+H^+(右)` 写成
    OH⁻（如 `CO2 + OH^- -> HCO3^-`、`Al^{3+} + 3OH^- -> Al(OH)_3`，即便
    这些体系并未投加强碱，OH⁻ 记为水自电离来源）。该惯例是既定语义，
    改动属于呈现层的全库重写，不是缺陷修复。故保留无条件还原。
    """
    if not allow:
        consumed = {sp: v for sp, v in consumed.items() if v > trace}
        produced = {sp: v for sp, v in produced.items() if v > trace}
        return consumed, produced
    w = consumed.get(WATER, 0.0)
    h = produced.get(H_ION, 0.0)
    if w > trace and h > trace:
        merge = min(w, h)
        consumed[WATER] = w - merge
        produced[H_ION] = h - merge
        consumed[OH_ION] = consumed.get(OH_ION, 0.0) + merge
        if consumed[WATER] <= trace:
            consumed.pop(WATER, None)
        if produced[H_ION] <= trace:
            produced.pop(H_ION, None)
    consumed = {sp: v for sp, v in consumed.items() if v > trace}
    produced = {sp: v for sp, v in produced.items() if v > trace}
    return consumed, produced


# ============================================================ 净方程守恒契约

# 报告量子（mol）：`_finalize_result` 用 round(x, 6) 产出 consumption/
# production，API 的呈现精度是 1e-6。
REPORT_QUANTUM = 1e-5

# 呈现分辨率（相对方程尺度）：项的系数被归一化后按 3 位小数打印，故
# 小于 5e-4 × 尺度 的项在**呈现上不可见**。守恒契约的容差必须取这个，
# 而不是 API 舍入量子——首版拿 1e-5 当契约容差，把"反应未走满"的
# 2.2e-5 残差（E45：账本 HCO₃⁻ 0.999978 而方程写 1）判成越界，闸门
# 因此拒绝一切合法改写 ⟹ 全库 111 例翻红。契约应约束"方程承诺的原子"，
# 容差则应与**它自己怎么显示**一致。
DISPLAY_QUANTUM = 5e-4

# 显示地板：系数按 3 位小数打印，小于 5e-4 的项打印出来就是 "0"——
# 浮点兜底路径必须把它们剔掉，否则出现 `1.001ClO^- + CO2 -> … + 0CO3^{2-}`
# 这种既有非整系数、又有零系数项的式子（T25/H77）。
_DISPLAY_FLOOR = 5e-4


# 迹量档的相对保真下限：绝对报告量子只许"占反应自身尺度的一小部分"。
# 判据 `max(DISPLAY_QUANTUM·scale, min(REPORT_QUANTUM, _REL_FLOOR·scale))`：
# 尺度 ≥ 1e-4 时与旧口径逐字相同（报告量子 1e-5 生效）；尺度更小时容差按
# 比例收紧。**为什么必须有这一条**：报告量子是**绝对**的（API 把 mol 舍入到
# 1e-6，容差留 10 倍余量），而净差可以是迹量级的——P10 的整条反应只有
# 1.17e-6 mol，1e-5 的容差比它大 10 倍，于是闸门允许把 Fe(OH)₃/Fe³⁺ 整块
# 删掉，剩下的 {H⁺ 4e-6} → {H₂O 3e-6} 归一化后印成 `1.333H^+ -> H_2O`
# ——一个连电荷都不平的"方程式"。相对保真闸把它挡在 10% 以内。
_REL_FLOOR = 0.1


def _contract_tol(scale: float) -> float:
    """守恒契约的容差 = max(呈现分辨率 × 尺度, min(报告量子, 10% × 尺度))。"""
    s = max(scale, 0.0)
    return max(DISPLAY_QUANTUM * s, min(REPORT_QUANTUM, _REL_FLOOR * s))


def _side_totals(side: dict, els: set) -> dict:
    return {el: sum(elements_of(s).get(el, 0) * v for s, v in side.items())
            for el in els}


def _drop_budget_ok(base_c: dict, base_p: dict,
                    c: dict, p: dict) -> bool:
    """**唯一的守恒闸门**：化简是否仍在守恒契约内。

    设计意图（v0.5.0 重构）：`_build_equations` 里有六个依次改数的启发式
    化简（质子噪声对消、水解对消、痕量过滤×2、共轭对对齐、H/O 补配平），
    此前**各自为政**——每个都能破坏守恒，没有任何一处对"真实净差"负责，
    容差又是 3% 相对值，于是 2.3% 的 Fe(OH)₃ 被删掉也照样通过并被固化成
    4 条测试标准。现在化简只是**候选**，唯一裁决在此。

    **判据（唯一一条）**：呈现出来的方程**自身必须配平**——非 H/O 元素与
    电荷两侧守恒到呈现分辨率 `_contract_tol`。H 与 O 不受约束（溶剂水的
    组成，由 `_balance_h_o` 用 H₂O/H⁺/OH⁻ 补齐，水活度 1、体相量级）。

    两次错误设计的记录（都实测过，别再走回去）：
      · 首版把 H/O 也纳入 → 把合法的补配平判成违规，`2H2S + SO2 -> 3S`
        连 `2H2O` 都丢了，近百例翻红；
      · 次版加"每侧元素总量相对 base（真值）不得漂移"的预算 → 但化简
        启发式（共轭对**对齐**、对消）**本来就会改量**，于是 T07
        `HCO3^- + OH^- -> CO3^{2-}` 这类完全正确的式子因"C 偏离真值
        0.014"被拒，浓硫酸/亚硝酸族同理，47 例翻红。
      **契约只约束"方程是不是一个配平的化学方程式"，不约束它等于账本的
      哪一部分**——后者是启发式（呈现选择）的职责，不是守恒的职责。
    真错误仍会被抓：Fe(OH)₃ 被水解对消删掉时，Fe 两侧差 0.528 ⟹ 拒绝。
    """
    els: set = set()
    for sp in list(c) + list(p):
        els |= set(elements_of(sp))
    els -= {"H", "O"}
    scale = max([abs(v) for v in list(base_c.values()) + list(base_p.values())]
                or [0.0])
    tol = _contract_tol(scale)
    for el in els:
        lhs = sum(elements_of(s).get(el, 0) * v for s, v in c.items())
        rhs = sum(elements_of(s).get(el, 0) * v for s, v in p.items())
        if abs(lhs - rhs) > tol:
            return False
    qc = sum(charge_of(s) * v for s, v in c.items())
    qp = sum(charge_of(s) * v for s, v in p.items())
    return abs(qc - qp) <= tol


# ============================================================ 系数有理化

def _fmt_term(nu: float, species: str) -> str:
    """格式化方程式一项：系数 1 省略，整数显示整数，浮点定点小数
    （禁用科学计数法，保持可解析）。"""
    if abs(nu - round(nu)) < 1e-6:
        nu = int(round(nu))
    if nu == 1:
        return species
    if isinstance(nu, int) or (isinstance(nu, float) and nu == int(nu)):
        return f"{int(nu)}{species}"
    return f"{nu:.3f}".rstrip("0").rstrip(".") + species


def _balanced_upto(consumed: dict, produced: dict) -> bool:
    """元素/电荷守恒判据，**恒用报告量子绝对容差**（不因系数恰为整数而改走
    精确档）。

    用途：校验"整数化会丢掉哪些痕量"。整数形式把 1e-6 级的痕量四舍五入为
    0 是**正当**的（该量低于 API 的 6 位小数报告精度），此时整数向量相对
    浮点向量天然差一个痕量；零容差会把这种正当舍入判成"不守恒"，于是所有
    含痕量的体系全退回浮点式（v0.4.5 实测：D11 出现 `7.999H^+`、
    D12 出现 `8.001H^+`，断言 2→18）。
    真正的错解（如 99:52 缺 1 个 Fe、电荷差 3）量级在 1 mol 尺度上，远大于
    1e-5，仍被拦下。
    """
    species = list(consumed) + list(produced)
    if not species:
        return True
    els: set = set()
    for sp in species:
        els |= set(elements_of(sp))
    for el in els:
        lhs = sum(elements_of(s).get(el, 0) * v for s, v in consumed.items())
        rhs = sum(elements_of(s).get(el, 0) * v for s, v in produced.items())
        if abs(lhs - rhs) > _EQ_TOL_ABS:
            return False
    return abs(sum(charge_of(s) * v for s, v in consumed.items())
               - sum(charge_of(s) * v for s, v in produced.items())) \
        <= _EQ_TOL_ABS


def _ints_balanced(names_c: list, ints_c: list,
                   names_p: list, ints_p: list) -> bool:
    """整数系数向量的元素与电荷守恒检查（零容差，精确整数算术）。

    用于**整数 snap** 路径——那里的系数本身就是化学计量整数，没有"近似
    配平"的余地（99:52 缺 1Fe、电荷 −3 即此档漏网并被固化成 4 条标准）。
    但 `_rationalize` 的整数化不同：它会把低于报告精度的痕量舍入为 0，
    那种情况用 `_balanced_upto`（见其文档）。
    """
    els: set = set()
    for s in names_c + names_p:
        els |= set(elements_of(s))
    for el in els:
        d = (sum(k * elements_of(s).get(el, 0) for s, k in zip(names_p, ints_p))
             - sum(k * elements_of(s).get(el, 0)
                   for s, k in zip(names_c, ints_c)))
        if d:
            return False
    return (sum(k * charge_of(s) for s, k in zip(names_p, ints_p))
            == sum(k * charge_of(s) for s, k in zip(names_c, ints_c)))


def _rationalize(vals: list[float], validate=None) -> list[int] | None:
    """浮点系数列表 → 最简整数比列表。

    两档精度逐级尝试（按最小值归一化后）：
      1. 细档：limit_denominator(8)（允许 1/8 级分数），相对误差 ≤2%——
         保留 25:17 这类真实非整数化学计量（CaCO3 溶解的 HCO3-/CO3^2- 分配），
         同时把 17.48 这类数值噪声收进 17.5；
      2. 粗档：limit_denominator(2)（整数/半整数），相对误差 ≤5%——
         用于噪声更大的近边界体系。
    两档都失败返回 None（调用方退化为浮点系数）。

    v0.4.5：新增 `validate`（整数向量 → 是否接受）。独立逐项近似会破坏
    浮点向量的守恒性（见 `_ints_balanced`），故由调用方传入守恒校验，
    **只接受守恒的整数向量**；两档都不守恒则返回 None，调用方退化为
    浮点系数（诚实呈现优于不守恒的"漂亮"整数）。
    """
    pos_vals = [v for v in vals if v > 0]
    if not pos_vals:
        return None
    mx_v = max(pos_vals)
    # 归一化基准取**非痕量**的最小值：若按全局最小值归一化，而最小值恰好是
    # 一个痕量（D11 的 Cr(OH)₃ 3.3e-4 对 H⁺ 7.999），主系数会被放大到数万，
    # `limit_denominator(8)` 只能给出 24021 这类巨型整数，调用方按
    # `max(int_vals) > 1000` 拒绝 ⟹ 全体系退回浮点式（`7.999H^+`）。
    _sig = [v for v in pos_vals if v >= 1e-3 * mx_v]
    min_v = min(_sig) if _sig else min(pos_vals)
    norm = [v / min_v for v in vals]
    for den_max, tol in ((8, 0.02), (2, 0.05)):
        fracs = []
        ok = True
        for v, raw in zip(norm, vals):
            if v <= 0:
                fracs.append(Fraction(0))
                continue
            # 痕量（< 0.1% 峰值）不参与"能否有理化"的判据：它们在任何整数
            # 形式里都必然被舍入为 0，而 `limit_denominator` 对 3e-6 这类
            # 相对量只能给出 0（相对误差 100%）⟹ 整档被判失败 ⟹ 全体系退回
            # 浮点式（实测 D11 出现 `7.999H^+`、D12 出现 `8.001H^+`）。
            # 舍入为 0 是否可接受，由 `validate`（守恒到报告量子）判定。
            if raw < 1e-3 * mx_v:
                fracs.append(Fraction(0))
                continue
            f = Fraction(v).limit_denominator(den_max)
            if abs(float(f) - v) / v > tol:
                ok = False
                break
            fracs.append(f)
        if not ok:
            continue
        pos = [f for f in fracs if f > 0]
        if not pos:
            return None
        min_frac = min(pos)
        normed = [f / min_frac for f in fracs]
        common = reduce(lambda a, b: a * b // gcd(a, b),
                        (f.denominator for f in normed), 1)
        int_vals = [int(f * common) for f in normed]
        g = reduce(gcd, [v for v in int_vals if v > 0], 0)
        if g == 0:
            return None
        int_vals = [v // g for v in int_vals]
        if validate is not None and not validate(int_vals):
            continue          # 本档整数比不守恒 → 试下一档
        return int_vals
    return None


def _try_integer_snap(consumed: dict[str, float], produced: dict[str, float],
                      max_coef_limit: int = 24) -> list[int] | None:
    """尝试将各系数四舍五入到最近的整数，再验证元素/电荷平衡。

    化学事实：净离子方程式系数应为整数比。引擎的 equilibrium 求解常给出
    0.986/1.014/0.014 这类"近整数但带 1-5% 残差"的值（弱酸/弱碱电离、
    配位平衡、同离子效应下的微小未转化部分）。常规 _rationalize 用
    limit_denominator(8) 把这些残差放大成 737:565:8 这种巨型系数。

    本函数策略：
      1. 多 scale 尝试（max/k for k=1..6）——适应不同化学计量比（1:1、2:1、
         3:2、5:1 等天然尺度）；
      2. 多 tol 尝试（2%/4%/6%/8%/10%/12%/15%）——从严格到宽松，取第一个
         能通过配平验证的；
      3. 每次都验证元素/电荷平衡（tol=0.5%）——防止错误 snap。
    返回与 consumed.values()+produced.values() 对齐的整数列表（含 0 表示
    被剔除的痕量物种），或 None。
    """
    all_vals = list(consumed.values()) + list(produced.values())
    pos = [v for v in all_vals if v > 0]
    if not pos:
        return None
    mx = max(pos)

    best: tuple[int, list[int]] | None = None

    for k in (1, 2, 3, 4, 5, 6):
        scale = mx / k
        if scale < 1e-9:
            continue
        for tol_rel in (0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15):
            snap_ints: list[int] = []
            cons_snap: dict[str, int] = {}
            prod_snap: dict[str, int] = {}
            ok = True
            for sp, v in consumed.items():
                if v <= 0:
                    snap_ints.append(0)
                    continue
                nv = v / scale
                r = round(nv)
                if r > 0 and abs(nv - r) <= tol_rel * max(r, 1):
                    snap_ints.append(r)
                    cons_snap[sp] = r
                elif r == 0 and nv <= 2 * tol_rel:
                    snap_ints.append(0)
                else:
                    ok = False
                    break
            if not ok:
                continue
            for sp, v in produced.items():
                if v <= 0:
                    snap_ints.append(0)
                    continue
                nv = v / scale
                r = round(nv)
                if r > 0 and abs(nv - r) <= tol_rel * max(r, 1):
                    snap_ints.append(r)
                    prod_snap[sp] = r
                elif r == 0 and nv <= 2 * tol_rel:
                    snap_ints.append(0)
                else:
                    ok = False
                    break
            if not ok or not cons_snap or not prod_snap:
                continue
            # 严格配平验证（元素 + 电荷）
            if not _balanced_quick(cons_snap, prod_snap, tol=0.005):
                continue
            mc = max(snap_ints)
            if mc > max_coef_limit:
                continue
            # 取 max coef 最小的方案
            if best is None or mc < best[0]:
                best = (mc, snap_ints)
            break  # 同 scale 下取第一个能配平的 tol（更低 tol 更安全）

    if best is None:
        return None
    snap_ints = best[1]
    g = reduce(gcd, [v for v in snap_ints if v > 0], 0)
    if g == 0:
        return None
    return [v // g for v in snap_ints]


def _project_balanced(species: list[str], signs: list[int],
                      amounts: list[float]) -> list[float] | None:
    """把量投影到元素/电荷平衡子空间（约束最小二乘 ν* = a − Aᵀ(AAᵀ)⁻¹Aa）。

    species/signs/amounts 平行列表（sign +1=消耗侧、−1=生产侧）。
    返回投影后的 ν*（平衡子空间内、与 a 欧氏最近）；子空间退化（仅零点
    ——剩余物种凑不出任何配平方程）返回 None。"""
    els: set[str] = set()
    for sp in species:
        els |= set(elements_of(sp))
    rows = []
    for el in sorted(els):
        rows.append([signs[i] * elements_of(sp).get(el, 0)
                     for i, sp in enumerate(species)])
    rows.append([signs[i] * charge_of(sp) for i, sp in enumerate(species)])
    n = len(species)
    m = len(rows)
    # AAᵀ（m×m）· y = A·a → ν* = a − Aᵀy（法方程解投影；m ≤ ~8 小系统）
    AAt = [[sum(rows[r][k] * rows[s][k] for k in range(n))
            for s in range(m)] for r in range(m)]
    Aa = [sum(rows[r][k] * amounts[k] for k in range(n)) for r in range(m)]
    # 高斯消元解 AAt·y = Aa（奇异 → 加 Tikhonov 投影到最近的平衡点）
    for col in range(m):
        piv = max(range(col, m), key=lambda r: abs(AAt[r][col]))
        if abs(AAt[piv][col]) < 1e-12:
            AAt[col][col] += 1e-9
            piv = col
        AAt[col], AAt[piv] = AAt[piv], AAt[col]
        Aa[col], Aa[piv] = Aa[piv], Aa[col]
        inv = 1.0 / AAt[col][col]
        for r in range(m):
            if r == col:
                continue
            f = AAt[r][col] * inv
            if f == 0.0:
                continue
            for cc in range(col, m):
                AAt[r][cc] -= f * AAt[col][cc]
            Aa[r] -= f * Aa[col]
    y = [Aa[i] / AAt[i][i] for i in range(m)]
    nu = [amounts[i] - sum(rows[r][i] * y[r] for r in range(m))
          for i in range(n)]
    # 投影退化检查：全零或负分量占比过高（正锥外——无法作为方程系数）
    mx = max((abs(v) for v in nu), default=0.0)
    if mx < 1e-12:
        return None
    neg = sum(1 for v in nu if v < -0.05 * mx)
    if neg:
        return None
    return [max(v, 0.0) for v in nu]


def _beautify_big_coeff(consumed: dict[str, float],
                        produced: dict[str, float]) -> tuple[dict, dict] | None:
    """大系数美化（v0.3.8）：浓度依赖的物种分配使真实净比为无理数
    （N10 PbBr₂+[PbBr₃]⁻ 混合 296:133:105:28）。教科书口径：净离子方程
    呈现主通道的干净计量，混合细节由多步方程列表承载。

    方法：逐个（先小者）再成对移除次要物种（量 <30% 峰值），把剩余量
    投影到元素/电荷平衡子空间（_project_balanced），对投影方向做小分母
    有理重构（≤24），精确配平校验后取最小整数比。全程配平硬保证——
    化学计量关系不因美化失真。找不到干净形式返回 None（保持原样）。

    **0.4.4 加宽移除集**：原先只试单/双物种移除，容量不足以覆盖 3 个混合
    配位通道——Fe³⁺+SCN⁻ 的 1/2/3 配位共存体系里，单/双移除后仍剩两个
    通道，snap 到小整数时缺 1Fe、电荷 −3（`99:52 -> 21:18:12`），而
    `_balanced_quick` 的相对容差没拦住它。补三连移除后该族塌缩到干净的
    单通道式 `3SCN^- + Fe^{3+} -> [Fe(SCN)_3]`（精确守恒）。
    """
    all_sp = list(consumed) + list(produced)
    a_all = [consumed.get(s, 0.0) + produced.get(s, 0.0) for s in all_sp]
    mx = max(a_all, default=0.0)
    if mx <= 0:
        return None
    minor = [s for s, v in zip(all_sp, a_all) if v < 0.3 * mx]
    if not minor:
        return None
    minor.sort(key=lambda s: consumed.get(s, 0.0) + produced.get(s, 0.0))
    removal_sets: list[tuple[str, ...]] = [(s,) for s in minor]
    removal_sets += [(x, y) for i, x in enumerate(minor)
                     for y in minor[i + 1:]]
    triples = sorted(((x, y, z) for i, x in enumerate(minor)
                      for j, y in enumerate(minor[i + 1:], i + 1)
                      for z in minor[j + 1:]),
                     key=lambda t: sum(consumed.get(s, 0.0)
                                       + produced.get(s, 0.0) for s in t))
    removal_sets += triples[:60]
    for rm in removal_sets:
        rms = set(rm)
        c2 = {s: v for s, v in consumed.items() if s not in rms}
        p2 = {s: v for s, v in produced.items() if s not in rms}
        if not c2 or not p2:
            continue
        sp2 = list(c2) + list(p2)
        sg2 = [1] * len(c2) + [-1] * len(p2)
        am2 = [c2.get(s, 0.0) + p2.get(s, 0.0) for s in sp2]
        nu = _project_balanced(sp2, sg2, am2)
        if nu is None:
            continue
        # 投影方向的小分母有理重构：逐分母尝试，精确配平 + 方向贴近
        mxp = max(nu)
        if mxp <= 0:
            continue
        unit = [v / mxp for v in nu]
        for den in range(1, 25):
            frs = []
            ok = True
            for v in unit:
                fr = Fraction(v).limit_denominator(den)
                if fr <= 0 or fr > 25:
                    ok = False
                    break
                frs.append(fr)
            if not ok:
                continue
            # LCM 分母整数化 → gcd 约减 → 最小整数比
            L = 1
            for f in frs:
                L = L * f.denominator // gcd(L, f.denominator)
            ints_i = [int(f * L) for f in frs]
            g = reduce(gcd, ints_i)
            if g == 0:
                continue
            ints_i = [v // g for v in ints_i]
            if max(ints_i) > 20 or min(ints_i) < 1:
                continue
            imx = max(ints_i)
            if any(abs(v / imx - u) > 0.1 for v, u in zip(ints_i, unit)):
                continue   # 方向偏离投影 >10%：重构失真
            scale = [float(v) for v in ints_i]
            nc = {s: k for s, k in zip(sp2[:len(c2)], scale[:len(c2)])}
            np_ = {s: k for s, k in zip(sp2[len(c2):], scale[len(c2):])}
            if not _balanced_quick(nc, np_, tol=0.001):
                continue
            # v0.5.0：整数化的美化产物走**报告量子绝对容差**——0.1% 相对
            # 容差正是 `99SCN^- + 52Fe^{3+}`（缺 1Fe、电荷 −3）当年得以
            # 通过的漏洞（见 _balanced_quick 文档）。
            if not _ints_balanced(list(nc), [int(v) for v in nc.values()],
                                  list(np_), [int(v) for v in np_.values()]):
                continue
            return nc, np_
    return None


# 可逆箭头（v0.5.0 用户口径）：**平衡过程用可逆符号**，完全反应用单向箭头。
# 用 ASCII 的 `<=>` 而不是 U+21CC：渲染结果要进日志/断言/GBK 控制台，
# ASCII 不会在打印时炸掉（语义相同）。
ARROW_FWD = "->"
ARROW_EQ = "<=>"


def species_plain(sp: str) -> str:
    """物种名 → **纯字符串**（连上下标标注都不用，最大兼容性）。

    只去掉**标记字符** `_ ^ { }`，内容一律保留：
    `H^+` → `H+`、`SO_4^{2-}` → `SO42-`、`CH_3COOH` → `CH3COOH`、
    `[Fe(SCN)_3]` → `[Fe(SCN)3]`。这是**手动读取的 fallback**（默认渲染是
    TeX，见 `Equation.__str__`）：供不支持任何标记的场合（纯文本日志、
    文件名、旧终端）使用。"""
    for ch in "_^{}":
        sp = sp.replace(ch, "")
    return sp


def species_tex(sp: str) -> str:
    r"""物种名 → TeX（**TeX 模式**的语法细节集中在这里）。

    引擎物种记法（`H^+` / `SO_4^{2-}` / `[Fe(SCN)_3]` / `CH_3COOH`）→ TeX：
      · 配方部分整体进 `\mathrm{}`（正体，符合化学排版惯例）；
      · `_x` → `_{x}`、`^x` → `^{x}`（多字符必须加花括号，否则 TeX 只取首字符）；
      · 电荷挂在 `\mathrm{}` **外面**（`\mathrm{SO_4}^{2-}`），与教科书一致。
    纯字符串模式不走这里（用户口径：纯字符串不考虑上下标问题，求最大兼容性）。
    """
    body, charge = sp, ""
    i = sp.rfind("^")
    if i >= 0:
        body, charge = sp[:i], sp[i + 1:]
    # 下标补花括号：_4 → _{4}（已带花括号的原样保留）
    out = []
    k = 0
    while k < len(body):
        ch = body[k]
        if ch in "_^" and k + 1 < len(body):
            nxt = body[k + 1]
            if nxt == "{":
                j = body.find("}", k)
                j = len(body) if j < 0 else j + 1
                out.append(f"{ch}{body[k + 1:j]}")
                k = j
                continue
            out.append(f"{ch}{{{nxt}}}")
            k += 2
            continue
        out.append(ch)
        k += 1
    tex = "\\mathrm{" + "".join(out) + "}"
    if charge:
        c = charge[1:-1] if charge.startswith("{") and charge.endswith("}") \
            else charge
        tex += "^{" + c + "}"
    return tex


def render_equation(left: dict, right: dict, reversible: bool = False,
                    mode: str = "tex") -> str | None:
    """把两侧系数渲染成离子方程式字符串（**渲染与结构分离**，v0.5.0）。

    箭头由**布尔**决定（用户口径：结构里存 reversible，符号只在渲染时出现）：
    `reversible=True` → `<=>`（平衡过程）；`False` → `->`（完全反应）。
    用 ASCII `<=>` 而非 U+21CC：渲染结果会进日志/断言/GBK 控制台，
    ASCII 不会在打印时炸掉（语义相同）。

    物种排序：先按系数降序，再按名称字典序（确定性输出，便于断言）。"""
    if not left or not right:
        return None
    c_items = sorted(((sp, v) for sp, v in left.items() if v > 0),
                     key=lambda x: (-x[1], x[0]))
    p_items = sorted(((sp, v) for sp, v in right.items() if v > 0),
                     key=lambda x: (-x[1], x[0]))
    if not c_items or not p_items:
        return None
    if mode == "plain":
        # 纯字符串 fallback：无任何标记（不写 `_`/`^`，也不加 \mathrm）
        def _p(v, sp):
            coef = "" if abs(v - 1.0) < 1e-12 else f"{_fmt_term(v, '')}"
            return f"{coef}{species_plain(sp)}"
        arrow = ARROW_EQ if reversible else ARROW_FWD
        lhs = " + ".join(_p(v, sp) for sp, v in c_items)
        rhs = " + ".join(_p(v, sp) for sp, v in p_items)
        return f"{lhs} {arrow} {rhs}"
    # 默认 TeX：放进 Markdown 的 $$…$$ 块即可正常显示
    def _t(v, sp):
        coef = "" if abs(v - 1.0) < 1e-12 else f"{_fmt_term(v, '')}"
        return f"{coef}{species_tex(sp)}"
    arrow = "\\rightleftharpoons" if reversible else "\\rightarrow"
    lhs = " + ".join(_t(v, sp) for sp, v in c_items)
    rhs = " + ".join(_t(v, sp) for sp, v in p_items)
    return f"{lhs} {arrow} {rhs}"


def make_equation(consumed: dict, produced: dict, reversible: bool = False,
                  kind: str = "", extent: float | None = None) -> "Equation | None":
    """由净差字典构造 `Equation`（规整系数后结构化存储）；空净差返回 None。"""
    norm = _normalize_equation(consumed, produced)
    if norm is None:
        return None
    return Equation(norm[0], norm[1], reversible, kind, extent)


@dataclass(frozen=True)
class Equation:
    """**结构化**离子方程式（先存结构，后渲染——用户 v0.5.0 口径）。

    left/right  物种 → 系数（呈现系数：整数 snap / 有理化 / 浮点兜底）
    reversible  **布尔**：True = 平衡过程（渲染成 `<=>`）；False = 完全反应（`->`）
    kind        机制标签（precip/dissolve/proton/redox/complex/…，可为空）
    extent      该步程度（mol，路由步骤用；净方程为 None）

    结构里**不存符号**（用户口径：可逆是语义、符号是渲染细节）；字符串只作
    呈现视图 `str(eq)`，比较/断言走结构（见 testsuit）。
    """
    left: dict
    right: dict
    reversible: bool = False
    kind: str = ""
    extent: float | None = None

    def __str__(self) -> str:
        """**默认 = TeX 模式**：直接放进 Markdown 的 `$$…$$` 块即可显示。"""
        return self.tex()

    def tex(self) -> str:
        r"""TeX 渲染（`\mathrm{}` 正体 + 规范上下标 + `\rightarrow`/`\rightleftharpoons`）。"""
        return render_equation(self.left, self.right, self.reversible,
                               mode="tex") or ""

    def plain(self) -> str:
        """**纯字符串 fallback**（手动调用）：不带任何上下标标注，最大兼容性。"""
        return render_equation(self.left, self.right, self.reversible,
                               mode="plain") or ""

    def __repr__(self) -> str:
        return f"Equation({self.plain()!r}, reversible={self.reversible})"


def _normalize_equation(consumed: dict[str, float], produced: dict[str, float]
                        ) -> tuple[dict, dict] | None:
    """把净差规整为**呈现系数**（结构，不拼接字符串）。

    优先级：① 整数 snap；② `_rationalize` 分数有理化；③ 浮点系数兜底。
    返回 (left, right) 两个"物种→系数"字典；调用方交给 `render_equation`
    渲染（v0.5.0：**先结构化存储、后渲染**——用户口径）。

    排序与渲染集中在 `render_equation`（先按系数降序、再按名称字典序）。"""
    if not consumed or not produced:
        return None
    # 显示级痕量剔除：<0.1% 峰值的项是数值噪声（近中性点体系的 H+ 残差
    # 会把有理化放大成巨型系数）；剔除后必须仍配平，否则保留原样
    all_vals = list(consumed.values()) + list(produced.values())
    mx = max(all_vals)
    tiny = {sp for side in (consumed, produced)
            for sp, v in side.items() if v < 1e-3 * mx}
    if tiny:
        c2 = {s: v for s, v in consumed.items() if s not in tiny}
        p2 = {s: v for s, v in produced.items() if s not in tiny}
        if c2 and p2 and _balanced_quick(c2, p2):
            consumed, produced = c2, p2
            all_vals = list(consumed.values()) + list(produced.values())

    # ① 整数 snap 优先（最干净的化学计量比，处理 equilibrium 残差）
    #    v0.4.5：整数化产物必须**精确守恒**才接受（见 _ints_balanced）
    _cn, _pn = list(consumed), list(produced)
    _nc = len(_cn)
    snap_ints = _try_integer_snap(consumed, produced)
    if snap_ints is not None and not _ints_balanced(
            _cn, snap_ints[:len(_cn)], _pn, snap_ints[len(_cn):]):
        snap_ints = None
    if snap_ints is not None:
        n_cons = len(consumed)
        cons_ints = snap_ints[:n_cons]
        prod_ints = snap_ints[n_cons:]
        cons_items = sorted([(s, v) for s, v in zip(consumed.keys(), cons_ints) if v > 0],
                            key=lambda x: (-x[1], x[0]))
        prod_items = sorted([(s, v) for s, v in zip(produced.keys(), prod_ints) if v > 0],
                            key=lambda x: (-x[1], x[0]))
        return dict(cons_items), dict(prod_items)

    # ② _rationalize 分数有理化（真实非整数比，如 1:2:1.5）——
    #    只接受**守恒**的整数向量（独立逐项近似会破坏守恒，见 _ints_balanced）
    def _ok_ints(iv: list) -> bool:
        c = {s: v for s, v in zip(_cn, iv[:_nc]) if v > 0}
        p = {s: v for s, v in zip(_pn, iv[_nc:]) if v > 0}
        return _balanced_upto(c, p)

    int_vals = _rationalize(all_vals, validate=_ok_ints)
    if int_vals is None or max(int_vals) > 1000:
        # ③ 浮点系数兜底（真实非化学计量混合，如 Fe→Fe2+/Fe3+ 混合价态）
        # 归一化用**稳健尺度**：取"显示显著"物种（≥0.1% 峰值）的最小值，
        # 而不是全局最小值。全局 min 一旦落在阈值边界的痕量物种上
        # （如 round(x,6) 恰好留下的 1e-6），整式会被放大约 10⁶ 倍——
        # v0.4.5 首轮试修即此：D12 出现 `941116.647H^+`、T78 出现
        # `31123.319H^+`，断言 2→20 全红。系数比不变 ⟹ 守恒性不受影响。
        _mx = max(all_vals)
        _sig = [v for v in all_vals if v >= 1e-3 * _mx]
        scale = min(_sig) if _sig else min(all_vals)
        if scale <= 0:
            return None
        cons_items = sorted(((sp, v / scale) for sp, v in consumed.items()
                             if v / scale > _DISPLAY_FLOOR),
                            key=lambda x: (-x[1], x[0]))
        prod_items = sorted(((sp, v / scale) for sp, v in produced.items()
                             if v / scale > _DISPLAY_FLOOR),
                            key=lambda x: (-x[1], x[0]))
        if not cons_items or not prod_items:
            # 显示阈值把一侧清空（全痕量体系）：放弃阈值，保底显示（诚实优先）
            cons_items = sorted(((sp, v / scale) for sp, v in consumed.items()),
                                key=lambda x: (-x[1], x[0]))
            prod_items = sorted(((sp, v / scale) for sp, v in produced.items()),
                                key=lambda x: (-x[1], x[0]))
    else:
        n_cons = len(consumed)
        cons_ints = int_vals[:n_cons]
        prod_ints = int_vals[n_cons:]
        # 零系数不显示（整数化/抵消后会出现 0，如 D12 的 `0MnCO_3`、
        # T78 的 `0NO_2^-`）——snap 路径本来就有 `if v > 0`，此处补齐
        cons_items = sorted([(s, v) for s, v in zip(consumed.keys(), cons_ints)
                             if v > 0], key=lambda x: (-x[1], x[0]))
        prod_items = sorted([(s, v) for s, v in zip(produced.keys(), prod_ints)
                             if v > 0], key=lambda x: (-x[1], x[0]))
    return dict(cons_items), dict(prod_items)


# ============================================================ 步骤聚合与方程式构建

# 净差阈值（mol）：低于此值的物种视为数值噪声，不进入方程。
_TRACE = 0.01

# 步骤显著性阈值（mol）：extent 低于此值视为数值噪声。
STEP_MIN = 1e-3


def _step_key(reactants: dict, products: dict) -> tuple:
    """步骤的规范键（忽略 H2O/H+ 挂侧差异的净反应键）。"""
    r2 = tuple(sorted((s, round(n, 6)) for s, n in reactants.items()))
    p2 = tuple(sorted((s, round(n, 6)) for s, n in products.items()))
    return (r2, p2)


def _collect_steps(steps: list[dict]) -> list[tuple[dict, dict, float]]:
    """从原始 steps 收集显著步骤：(reactants, products, extent)。

    - extent < STEP_MIN：噪声，跳过；
    - kind == 'dissolve'：物理溶解（NaHCO3 → Na+ + HCO3-），离子方程式
      不体现——真正参与反应的是溶解后的离子，由后续步骤表达；
    - 同一净反应（含 H2O 挂侧差异）聚合：extent 求和；
    - 互为逆反应的对子在聚合时自然抵消（extent 相减）。
    - 瞬态中间体塌缩：复杂形成+解离对（如 Ag+ + 2I- → [AgI2]-,
      [AgI2]- + Ag+ → 2AgI）塌缩为直接沉淀步（Ag+ + I- → AgI），
      避免多步方程中暴露瞬态中间体（化学事实上中间体不出现在净方程中）。
    """
    agg: dict[tuple, float] = {}
    reps: dict[tuple, tuple[dict, dict]] = {}
    for st in steps:
        ext = st.get("extent", 0.0)
        if ext < STEP_MIN or st.get("kind") == "dissolve":
            continue
        eq = st.get("equation", "")
        if not eq:
            continue
        r, p = _parse_equation(eq)
        if not r and not p:
            continue
        k = _step_key(r, p)
        rev = (k[1], k[0])
        if rev in agg:
            agg[rev] -= ext
            if abs(agg[rev]) < STEP_MIN:
                agg.pop(rev)
                reps.pop(rev, None)
            continue
        agg[k] = agg.get(k, 0.0) + ext
        reps.setdefault(k, (r, p))
    collected = [(reps[k][0], reps[k][1], ext) for k, ext in agg.items() if ext >= STEP_MIN]
    return _collapse_transient_intermediates(collected)


def _collapse_transient_intermediates(
        collected: list[tuple[dict, dict, float]]
        ) -> list[tuple[dict, dict, float]]:
    """瞬态中间体塌缩：Hess 风格的中间体消去。

    化学事实：如 Ag+ + 2I- → [AgI2]-（β2 形成）后立即 [AgI2]- + Ag+ → 2AgI
    （归中沉淀），中间体 [AgI2]- 仅是数值求解的桥接，不出现在净方程。
    本函数对每个 complex 物种 C：
      1. 收集所有涉及 C 的步骤（C 在反应物或产物侧）
      2. 计算 C 的形成总量（产 C 的 Σ nu×ext）与消耗总量（耗 C 的 Σ nu×ext）
      3. 若两者差额 < 5% max（C 是瞬态），合并所有涉及 C 的步骤为一个净步
         （消去 C 自身，保留其他物种的净消耗/生成）
      4. 若差额 > 5% max（C 是真实中间体或终产物），保留原步不动
    """
    if not collected:
        return collected

    def _is_complex(sp: str) -> bool:
        # 形如 [Ag(NH_3)_2]^+ / [AgI_2]^- / [CuCl_4]^{2-} 等配离子
        return sp.startswith("[")

    # 找所有 complex 物种
    all_complexes: set[str] = set()
    for r, p, e in collected:
        for sp in list(r) + list(p):
            if _is_complex(sp):
                all_complexes.add(sp)

    out = list(collected)
    for c in all_complexes:
        # 收集所有涉及 c 的步骤
        involving = [(r, p, e) for r, p, e in out if c in r or c in p]
        if not involving:
            continue
        # 计算 c 的形成和消耗总量
        total_formed = sum(p.get(c, 0.0) * e for r, p, e in involving)
        total_consumed = sum(r.get(c, 0.0) * e for r, p, e in involving)
        net = total_formed - total_consumed
        max_amt = max(total_formed, total_consumed, 1e-9)
        # 若 c 是瞬态（净 ≈ 0），合并所有涉及 c 的步骤为单一净步
        if abs(net) / max_amt < 0.05:
            # 计算净消耗/生成（去除 c 自身）
            net_r: dict[str, float] = {}
            net_p: dict[str, float] = {}
            for r, p, e in involving:
                for s, nu in r.items():
                    if s == c:
                        continue
                    net_r[s] = net_r.get(s, 0.0) + nu * e
                for s, nu in p.items():
                    if s == c:
                        continue
                    net_p[s] = net_p.get(s, 0.0) + nu * e
            # 抵消跨侧物种
            for s in set(net_r) & set(net_p):
                m = min(net_r[s], net_p[s])
                net_r[s] -= m
                net_p[s] -= m
            # 去除零项
            net_r = {s: v for s, v in net_r.items() if v > 1e-9}
            net_p = {s: v for s, v in net_p.items() if v > 1e-9}
            # 从 out 移除 involving，加入新的合并步
            out = [(r, p, e) for r, p, e in out
                   if not any(r is r2 and p is p2 and e == e2
                              for r2, p2, e2 in involving)]
            if net_r and net_p:
                # extent 取形成/消耗总量的最大值
                new_ext = max_amt
                out.append((net_r, net_p, new_ext))
    return out


def _step_to_ionic(reactants: dict, products: dict, extent: float
                   ) -> tuple[dict, dict, str | None]:
    """单个聚合步骤 → (consumed, produced, 纯字符串视图)。

    含 H2O 配平与 OH- 还原；第三个返回值只用于**步骤去重**（结构才是权威，
    见 `Equation`），按纯字符串（无标记）渲染。"""
    consumed = {sp: nu * extent for sp, nu in reactants.items()}
    produced = {sp: nu * extent for sp, nu in products.items()}
    consumed, produced = _balance_h_o(consumed, produced)
    consumed, produced = _restore_oh(consumed, produced, _TRACE)
    if not consumed or not produced:
        return {}, {}, None
    _n = _normalize_equation(consumed, produced)
    return consumed, produced, (None if _n is None
                                else render_equation(_n[0], _n[1], mode="plain"))


# 净方程式中"分子态 → 离子形"的改写不再维护本地表（旧 _NET_ACID_SPLIT /
# _NET_SALT_SPLIT 63 条与引擎拆盐知识重复，属双重记账）。单一数据源：
# 引擎 judge 结果的 "ionize" 图（浓酸分子形态 / ex ions 形态 / 高溶解度盐，
# 见 engine._ionize_map），本模块仅消费。


def _pool_ligand_flux(consumed: dict[str, float], produced: dict[str, float],
                      r: dict) -> dict[str, float]:
    """池内配合物的净配体摄入（带符号；正=游离配体被池吸入，负=释放）。

    元素守恒配套：[ZnCl]⁺ 净增 Δm ⟹ 吸入 Δm×ν(Cl⁻)——折叠把 [ZnCl]⁺
    并入 Zn²⁺ 名义时，这部分 Cl 原子的去向从 Cl⁻ 侧同步扣除，方程才
    平衡（T105 的 10Cl⁻ 旁观者消除、NR19 的配体消耗归零）。"""
    pool_nu = r.get("pool_nu") or {}
    if not pool_nu:
        return {}
    flux: dict[str, float] = {}
    for comp, (lig, nu) in pool_nu.items():
        d = produced.get(comp, 0.0) - consumed.get(comp, 0.0)
        if abs(d) > 1e-12:
            flux[lig] = flux.get(lig, 0.0) + nu * d
    return flux


def _fold_pools(consumed: dict[str, float], produced: dict[str, float],
                r: dict) -> tuple[dict[str, float], dict[str, float]]:
    """弱池折叠（v0.4.0 形态池语义，两侧 + 配体守恒配套）。

    池成员净差并入 anchor（center）名下；池内配合物的配体摄入/释放
    同步从游离配体侧扣除（元素守恒）。池内再分布在后续跨侧抵消中
    归零（NR19 无净方程）；池总量变化以 anchor 名义进入事件口径
    （O02 金属→池 0.95 ⟹ 2H⁺+Zn→H₂+Zn²⁺ 教科书形式回归）。"""
    pools = r.get("pool_map") or {}
    if not pools:
        return consumed, produced

    def _fold(side: dict[str, float]) -> dict[str, float]:
        out: dict[str, float] = {}
        for sp, v in side.items():
            a = pools.get(sp, sp)
            out[a] = out.get(a, 0.0) + v
        return out

    flux = _pool_ligand_flux(consumed, produced, r)
    consumed = _fold(consumed)
    produced = _fold(produced)
    for lig, x in flux.items():
        if abs(x) <= 1e-9:
            continue
        if x > 0:      # 池吸入：游离配体的事件消耗扣除吸入量
            consumed[lig] = consumed.get(lig, 0.0) - x
            if consumed[lig] <= 1e-9:
                consumed.pop(lig, None)
        else:          # 池释放：配体回归不计事件生成
            produced[lig] = produced.get(lig, 0.0) + x   # x<0，扣除
            if produced[lig] <= 1e-9:
                produced.pop(lig, None)
    return consumed, produced


def _pool_step_key(sp_r: set[str], sp_p: set[str],
                   pools: dict[str, str], pool_members: dict[str, list[str]],
                   pool_ligands: dict[str, set[str]]) -> bool:
    """步骤是否为弱池内再分布步骤（事件口径下不存在）。

    判据：物种集（除 H₂O/H⁺/OH⁻）全部落在单一池成员集 ∪ 该池族的游离
    配体集合内，且含至少一个池 complex（≠ anchor 的成员）——保证这是
    配位/解配/配体交换步而非跨池反应（Zn²⁺+C₂O₄²⁻→ZnC₂O₄ 沉淀保留、
    AgCl+2NH₃→[Ag(NH₃)₂]⁺ 强族事件保留）。"""
    sps = (sp_r | sp_p) - {WATER, H_ION, OH_ION}
    comp = {sp for sp in sps if sp in pools and pools[sp] != sp}
    if not comp:
        return False
    anchors = {pools[sp] for sp in comp}
    if len(anchors) != 1:
        return False
    a = next(iter(anchors))
    allowed = set(pool_members.get(a, ())) | set(pool_ligands.get(a, ()))
    return sps <= allowed


def _build_equations(steps: list[dict], r: dict
                     ) -> tuple[list, dict, dict, "Equation | None",
                                "Equation | None"]:
    """构建多步离子方程式列表 + 净离子方程式（**两版本**，v0.5.0）。

    返回 (equations, consumption, production, net_equation, net_equation_raw)：
      equations         各显著步骤的离子方程式（聚合、对消、按 extent 降序）；
      consumption/production  净消耗/生成（mol，含 H+/OH-/H2O）——来自账本净差
                        （consumption_raw/production_raw + H_excess 变化 + 中和步），
                        是真实摩尔量，振荡与中间体天然抵消；
      net_equation      净离子方程式（**精编版**）：按"有意义的平衡过程"口径
                        呈现——实质反应、解离/水解/配位、投料自带固相的溶解度
                        表达都给；只有**痕量副过程**（生成投料里没有的新相/
                        痕量配合物，且非该过程的全部）置 None；
      net_equation_raw  **原始版**：不做上述口径裁剪，永远给出账本净差的格式化
                        （用户要求："不带 raw 是当前行为，带 _raw 则三种情况
                        都考虑"）。两者不一致时正是"该体系只发生了痕量副过程"。
    """
    # ---- 多步方程式（步骤聚合路径）----
    collected = _collect_steps(steps)
    # v0.4.0 弱池内步骤过滤（事件口径）：Zn²⁺+Cl⁻⇌[ZnCl]⁺ 类再分布步
    # 不进多步方程（与 changed=False 同口径呈现）；跨池/强族步骤保留
    _pools = r.get("pool_map") or {}
    if _pools and collected:
        _pm = r.get("pool_members") or {}
        _pl = r.get("pool_ligands") or {}
        collected = [(rr, p, ext) for rr, p, ext in collected
                     if not _pool_step_key(set(rr), set(p), _pools, _pm, _pl)]
    # 箭头按**体系**判定（v0.5.0 用户口径）：完全反应 `->`，残存平衡 `<=>`。
    # 逐步判定需要每步的 conversion，而聚合后的步骤已混入多个来源；体系级
    # degree 是引擎自己的权威判定（0/1=未完成 ⟹ 平衡；2=完全），用它统一
    # 且自洽：净方程与多步叙述的箭头不会互相矛盾。
    # degree 是引擎给的三档**码**：2=完全（conv ≥ DEGREE_COMPLETE）、
    # 1=部分、0=无显著变化（见 engine 的 degree 计算）。注意别把它与
    # candidates 的转换率阈值 DEGREE_COMPLETE=0.99 混用。
    _rev = int(r.get("degree") or 0) < 2
    equations: list[Equation] = []
    if collected:
        # 步骤取舍用**引擎自己的显著程度判据**（`ANN_MIN_EXTENT`），不用
        # "主步骤的 5%"这个第二套阈值。两套阈值会互相打架：真实执行过、
        # 且远超显著线（0.04 mol 对 1e-3 的线 = 40 倍）的步骤，只因占主步骤
        # 比例小（H89 的 `Ca(OH)₂ + CO₃²⁻ + 2H⁺ → CaCO₃` 占 4.2%、
        # P11 的 `CO₃²⁻ + H⁺ → HCO₃⁻` 占 3.5%）就被剔除，而 `eq_has` 断言的
        # 正是"这条真实步骤必须在多步叙述里"——叙述不该吞掉自己执行过的步。
        # 痕量步骤（H89 尾部的 2.6e-5 级振荡步）仍被 ANN_MIN_EXTENT 挡掉 ✓。
        kept = [c for c in collected if c[2] >= ANN_MIN_EXTENT]
        kept.sort(key=lambda c: -c[2])
        seen_eq: set[str] = set()
        for rr, p, ext in kept:
            sc, sp2, eq = _step_to_ionic(rr, p, ext)
            if eq and sc and sp2 and eq not in seen_eq:
                seen_eq.add(eq)
                # **结构化存储**：存**规整系数**（不是 mol 量），渲染交给
                # Equation.__str__——结构才是权威，字符串只是视图。
                _sn = _normalize_equation(sc, sp2)
                equations.append(Equation(*(_sn if _sn else (sc, sp2)),
                                          reversible=_rev, extent=ext))

    # ---- 净方程（账本净差路径）----
    # 优先吃**精确净差**（`net_exact`，不 round/不设阈）：迹量反应在 1e-6
    # 报告口径下会丢真实项（P10 的 NH₃ 4e-7），净差因此电荷不平，任何呈现
    # 都配不平。没有该字段时（外部构造的 r）退回报告口径。
    _nx = r.get("net_exact")
    if _nx:
        consumed: dict[str, float] = {k: v for k, v in _nx["c"].items() if v > 0.0}
        produced: dict[str, float] = {k: v for k, v in _nx["p"].items() if v > 0.0}
    else:
        consumed = {e["name"]: e["mol"] for e in r.get("consumption", [])}
        produced = {e["name"]: e["mol"] for e in r.get("production", [])}
    # 分子态 → 离子形（单一数据源：engine ionize 图——浓酸分子形态 HNO3 →
    # H+ + NO3-、高溶解度盐 NaHCO3 → Na+ + HCO3- 等），让旁观离子自然抵消
    for sp, ions in r.get("ionize", {}).items():
        for side in (consumed, produced):
            x = side.pop(sp, 0.0)
            if x > 0.0:
                for ion, nu in ions.items():
                    side[ion] = side.get(ion, 0.0) + nu * x
    # v0.4.0 弱池折叠（净方程事件口径）：池成员净差并入 anchor，配体摄入
    # 同步抵扣（元素守恒配套）——NR19 无净方程；T105 的 10Cl⁻ 旁观者
    # 消除，教科书形式 2H⁺+Zn→H₂+Zn²⁺ 回归
    consumed, produced = _fold_pools(consumed, produced, r)
    # 质子账本净变化：He 是带符号账本（正=游离强酸），|dHe| 的离子形态
    # 由体系酸碱性 + 是否含强酸强碱中和步共同决定：
    #   - 纯碱性体系（无中和步）：dHe 表达 H+/OH- 在弱酸弱碱间的迁移，
    #     以 OH- 形态书写更自然（H+ 源自水自发电离，OH- 是真实形态）；
    #   - 含中和步（即强酸强碱滴定）：H+ 是真实反应物（来自 HCl/H2SO4 等），
    #     保持 H+ 形态——化学习惯酸碱滴定方程用 H+ 不用 OH-。
    He_i = _nx["He_i"] if _nx else r.get("H_excess_initial", 0.0)
    He_f = _nx["He_f"] if _nx else r.get("H_excess", 0.0)
    pH_f = r.get("final_pH")
    basic = pH_f is not None and pH_f > 7.0
    has_neutralize = any(st.get("kind") == "neutralize" for st in steps)
    basic_free = basic and not has_neutralize
    dHe = He_i - He_f          # >0: 体系酸性减弱；<0: 酸性增强
    if dHe > 1e-9:
        if basic_free:
            produced[OH_ION] = produced.get(OH_ION, 0.0) + dHe
        else:
            consumed[H_ION] = consumed.get(H_ION, 0.0) + dHe
    elif dHe < -1e-9:
        if basic_free:
            consumed[OH_ION] = consumed.get(OH_ION, 0.0) + (-dHe)
        else:
            produced[H_ION] = produced.get(H_ION, 0.0) + (-dHe)
    # 初始中和步（H+ + OH- -> H2O，normalize 阶段记账）
    for st in steps:
        if st.get("kind") == "neutralize":
            e = st.get("extent", 0.0)
            consumed[H_ION] = consumed.get(H_ION, 0.0) + e
            consumed[OH_ION] = consumed.get(OH_ION, 0.0) + e
            produced[WATER] = produced.get(WATER, 0.0) + e
    consumed = {s: v for s, v in consumed.items() if v > 1e-9}
    produced = {s: v for s, v in produced.items() if v > 1e-9}
    if not consumed or not produced:
        return equations, {}, {}, None, None
    # ---- 痕量净反应过滤：仅过滤形成配合物的痕量反应（如同离子效应下
    # 的 [AgCl2]- / [Ag(NH3)2]+ 形成痕迹量），不过滤纯溶解的 Ksp 表达
    # （如 BaSO4 → Ba2+ + SO4^2- 是 Ksp 平衡表达，应保留）。
    # 化学事实：B03 AgCl+0.1M NaCl 仅形成 4e-6 mol [AgCl2]- (0.04% 投料)，
    # 不构成化学反应方程式（同离子抑制）。但 BaSO4 纯水中溶解 1e-5 M
    # 虽痕量，Ksp 表达"BaSO4 ⇌ Ba2+ + SO4^2-"是化学标准写法。
    total_feed = 0.0
    for e in r.get("consumption", []):
        total_feed += e["mol"]
    for e in r.get("final", []):
        total_feed += e["mol"]
    feed_total = max(total_feed, 0.01)
    net_total = sum(consumed.values())
    # 仅当产物中含配合物离子且总消耗 <0.5% 总投料 且 < 5e-3 mol 时过滤
    has_complex_product = any(sp.startswith("[") for sp in produced)
    # 痕量副过程不再早退：置旗标，末尾统一裁决——这样**原始净方程**
    # （net_equation_raw，用户要求的"两版本"之一）仍可给出，而精编版
    # （net_equation）按口径置 None。consumption/production 保持真实净差
    # （呈现层仍能看到痕量量值），只是不写成方程式。
    _trace_net = bool(has_complex_product
                      and net_total < 0.005 * feed_total and net_total < 5e-3)

    # ============================================================
    # 守恒契约（v0.5.0 重构）
    # ------------------------------------------------------------
    # 到此处为止的 `consumed/produced`（账本净差 + 质子账本 + 中和步记账）
    # 就是**真实净差**——它的元素与电荷守恒由引擎账本保证。下面所有"化简"
    # 都只是**候选**：每一步算出的结果都要过唯一闸门 `_gate`，闸门以本快照
    # 为基准算**累计**偏差（不是逐步偏差），超过报告量子即整步回退。
    # 于是"丢项"只会发生在呈现精度之内，方程永远不会少一个原子。
    # ============================================================
    base_c, base_p = dict(consumed), dict(produced)

    def _gate(cur_c: dict, cur_p: dict, prev_c: dict, prev_p: dict):
        """唯一裁决：候选相对真值的累计偏差在报告量子内则采纳，否则回退。"""
        if cur_c and cur_p and _drop_budget_ok(base_c, base_p, cur_c, cur_p):
            return cur_c, cur_p
        return prev_c, prev_p

    def _step(fn, cur_c: dict, cur_p: dict, *a, **kw):
        tc, tp = fn(dict(cur_c), dict(cur_p), *a, **kw)
        tc = {s: v for s, v in tc.items() if v > 1e-9}
        tp = {s: v for s, v in tp.items() if v > 1e-9}
        return _gate(tc, tp, cur_c, cur_p)

    consumed, produced = _cancel_minor_protonation(consumed, produced)
    consumed, produced = _cancel_minor_hydrolysis(consumed, produced)
    # 同物种跨侧抵消（浓酸分子化/再电离等记账形态转换会在两侧各留一份）
    for sp in set(consumed) & set(produced):
        x = min(consumed[sp], produced[sp])
        consumed[sp] -= x
        produced[sp] -= x
    consumed = {s: v for s, v in consumed.items() if v > 1e-9}
    produced = {s: v for s, v in produced.items() if v > 1e-9}
    if not consumed or not produced:
        return equations, {}, {}, None, None
    # 痕量过滤 → 共轭酸碱对质量平衡调整 → 配平 → OH- 还原 → 终过滤
    # 全部走同一个闸门 `_step`（相对 base 的累计偏差 ≤ REPORT_QUANTUM）。
    # 各化简的**内在启发式**（阈值、配对规则）保持不变——它们只是候选生成器；
    # 是否允许落地由闸门说了算，不再各自带一套私有守卫。
    consumed, produced = _step(
        lambda c, p: _filter_trace(c, p, balanced=True), consumed, produced)
    # 共轭酸碱对不平衡调整：过滤掉 [Ag(NH3)2]+ (痕量) 后，消耗侧 NH3 残留
    # 与产物侧 NH4+ 不等（差额来自被过滤的产物）。化学事实上被过滤的产物
    # 是由共轭碱+其他离子形成的（如 2NH3 + Ag+ → [Ag(NH3)2]+），差额
    # 残留在共轭碱侧。同步调整较小者到较大者以恢复化学计量比。
    consumed, produced = _step(_align_conjugate_pair, consumed, produced)
    c2, p2 = _balance_h_o(dict(consumed), dict(produced))
    consumed, produced = _gate({s: v for s, v in c2.items() if v > 1e-9},
                               {s: v for s, v in p2.items() if v > 1e-9},
                               consumed, produced)
    consumed, produced = _restore_oh(consumed, produced, 1e-9, allow=True)
    # 配平兜底可能在对侧补出 H+：再做一次跨侧抵消
    for sp in set(consumed) & set(produced):
        x = min(consumed[sp], produced[sp])
        consumed[sp] -= x
        produced[sp] -= x
    # 终过滤（5%）：同样过闸
    consumed = {s: v for s, v in consumed.items() if v > 1e-9}
    produced = {s: v for s, v in produced.items() if v > 1e-9}
    consumed, produced = _step(
        lambda c, p: _filter_trace(c, p, frac=0.05, balanced=True),
        consumed, produced)
    if not consumed or not produced:
        return equations, {}, {}, None, None
    # H/O 未配平的净差（隐式水档）：先整数化再补水的候选。浮点档补水必然
    # 引入 H+/OH- 而被闸门拒（见 _stoichiometric_net），呈现只能退化成
    # 既非整数、又缺水的浮点式——本步把它换成教科书整数式。
    if not _h_o_balanced(consumed, produced):
        consumed, produced = _step(_stoichiometric_net, consumed, produced)
    # ---- 痕量**新相**闸门（v0.5.0，用户确认的边界）----------------------
    # 位置关键：必须在化简链**之后**——此时 consumed/produced 才是最终叙述形，
    # 用它判痕量才准（早期净差还含池内周转：NR19 的 Zn-Cl 配位周转使净差达
    # 0.78 mol，早期判据会误判为"非痕量"）。
    # 规则：生成"投料里没有的固相"，且净差**相对总投料 <0.5% 且 <5e-3 mol**
    # ⟹ 不叙述。化学依据：这类体系**在饱和线上**
    # （1 M CuSO4 自身水解到 pH 4.15 时 Q = 10^-19.7 = Ksp(Cu(OH)2)，析出
    # 5e-4 mol 化学上正确），但写成净方程与 changed=False 自相矛盾。
    # 两类反例天然不受影响：① 投料自带固相的溶解（E39 BaSO4、J03/J04
    # Ca(OH)2）产物里没有"新固相"；② **痕量投料**下的沉淀（H03 AgCl 1e-4、
    # H05/H06）：绝对量小，但相对其自身投料规模并不痕量（feed_total 有
    # 0.01 mol 下限、判据取 0.5%）⟹ 保留（它就是该测试的主题）。
    # 注：不能用"步的 conversion"当判据——痕量步的 x_max 本身也是痕量，
    # 比值可以接近 1（T92 的逆析出步 conversion 0.487 即此），会把该抑制的
    # 情形放行。
    _sol = _solid_set()
    _feed_sp = {e["name"] for e in r.get("initial", [])}
    _new_solid = bool(_sol) and any(sp in _sol and sp not in _feed_sp
                                    for sp in produced)
    _conv = max([st.get("conversion", 0.0) for st in steps
                 if st.get("kind") in ("precip", "dissolve")] or [0.0])
    _net_tot = sum(consumed.values())
    if _TRACE_EQ:
        print(f"  [eq-gate] net={_net_tot:.3g} feed={feed_total:.3g} "
              f"new_solid={_new_solid} conv={_conv:.3g} steps={len(steps)}")
    # **不抑制痕量新相**（用户 v0.5.0 口径修正）：饱和线上的痕量析出是
    # **有意义的平衡过程**（1 M CuSO₄ 的 pH 4.15 正是 Cu(OH)₂ 的 Ksp 线），
    # 精编版应把它讲出来，只是讲成**可逆平衡**（`reversible=True` ⟹ `<=>`）
    # 而不是"反应发生了"。`net_equation_raw` 仍给同一净差的原始渲染。
    # 保留的只有"痕量配合物"一条（B03 同离子隐蔽形态，不构成化学方程式）。
    _norm = _normalize_equation(consumed, produced)
    _nc, _np = (_norm if _norm is not None else (dict(consumed), dict(produced)))
    net_str = render_equation(_nc, _np)
    # ---- 大系数美化（v0.3.8）：浮点系数或 max>20 的混合通道净差 →
    # 主通道呈现（移除次要物种 + 平衡子空间投影 + 小整数比；配平硬保证，
    # 见 _beautify_big_coeff）。找不到干净形式保持原样（诚实优先）。
    # **入口判据看结构，不看渲染串**（v0.5.0 修复）：`net_str` 自 v0.5.0 起是
    # TeX（`\mathrm{…} \rightarrow …`），而 `_max_coef` 走 `_parse_equation`
    # 的箭头表（`->`/`<=>`/`→`/`⇌`/`=`），TeX 箭头一律认不出 ⟹ 恒返回 0
    # ⟹ "系数 >20 也要美化"这条分支**静默失效**，只有含小数点的串才进得来。
    # 后果实测：21 AgBr+浓氨水 的整数化结果 `34NH₃+13AgBr+8H₂O ⇌ 13Br⁻+
    # 13[Ag(NH₃)₂]⁺+8NH₄⁺+8OH⁻`（max=34，无小数点）漏出到净方程，而同一体系
    # 在 pKw 未锚定时因系数是 4.278 带小数点而正常美化 ⟹ 同一条化学结论
    # （2NH₃+AgBr ⇌ [Ag(NH₃)₂]⁺+Br⁻）随阈值出现/消失。
    _mx0, _int0 = _coef_stats(_nc, _np)
    if net_str is not None and (not _int0 or _mx0 > 20):
        bp = _beautify_big_coeff(consumed, produced)
        if bp is not None:
            bc, bp_ = bp
            _bn2 = _normalize_equation(bc, bp_)
            pretty = (None if _bn2 is None
                      else render_equation(_bn2[0], _bn2[1], mode="plain"))
            # 接受条件（v0.5.0 放宽并改为"看产物质量"）：
            #   ① 干净——系数全整数
            #   ② 小系数——max ≤ 20
            #   ③ 覆盖——被保留的物种承载净差的主要部分（≥50%），
            #      否则"主通道"名不副实，不如呈现诚实混合比
            # 原判据是"≤3 物种"，本意是"只取教科书单通道形式"，但它把
            # **多通道但通道清晰**的式子挡在外面：N09 的理想式
            # `14H^+ + 6Fe^{2+} + Cr_2O_7^{2-} -> 7H_2O + 6Fe^{3+} + 2Cr^{3+}`
            # 有 7 个物种、系数全小且守恒，却被拒 ⟹ 退化成
            # `505.535H^+ + 217.944Fe^{2+} + 36.324Cr_2O_7^{2-} -> …`
            # （非整系数 + 大系数）。物种数不是"是否干净"的判据，
            # 系数是否整数、是否够小、是否覆盖主通道才是。
            _mxb, _intb = _coef_stats(bc, bp_)
            if (_intb and _mxb <= 20
                    and _beautify_covers(consumed, produced, bc, bp_)):
                net_str = pretty
                _bn = _normalize_equation(bc, bp_)
                _nc, _np = (_bn if _bn is not None else (dict(bc), dict(bp_)))
    net_obj = Equation(_nc, _np, _rev) if net_str else None
    net_raw = Equation(_nc, _np, _rev) if net_str else None
    return (equations, consumed, produced,
            None if _trace_net else net_obj, net_raw)


def _beautify_covers(consumed: dict, produced: dict,
                     bc: dict, bp_: dict, frac: float = 0.5) -> bool:
    """美化式是否覆盖了主通道：被保留物种承载的净差 ≥ frac（投向两侧的
    较小者，避免"只留反应物"也算覆盖）。"""
    tot = sum(consumed.values()) + sum(produced.values())
    keep = (sum(bc.values()) + sum(bp_.values())) * 0.5
    return tot <= 0 or keep >= frac * tot * 0.5


def _coef_stats(*sides: dict) -> tuple[float, bool]:
    """系数统计 → (最大系数, 是否全整数)。

    v0.5.0：判"要不要美化/美化是否干净"必须看**结构**。原实现用
    `_max_coef(渲染串)`，而 `net_str` 自 v0.5.0 起是 TeX（`\\mathrm{…}
    \\rightarrow …`）⟹ `_parse_equation` 的箭头表（`->`/`<=>`/`→`/`⇌`/`=`）
    认不出 TeX 箭头，恒返回空 ⟹ 恒得 0.0（详见 `_build_equations` 里的
    大系数美化入口注释）。结构判定无此问题，也不依赖渲染模式。
    """
    vals = [v for side in sides for v in side.values()]
    return (max(vals, default=0.0),
            all(abs(v - round(v)) < 1e-9 for v in vals))


def _balanced_quick(consumed: dict, produced: dict,
                    tol: float | None = None) -> bool:
    """配平验证：元素与电荷两侧相等。

    **整数系数走精确整数算术；浮点系数走"报告量子"绝对容差**。

    0.4.x 一律用**相对**容差（默认 3%、美化路径 0.1%），容差被最大系数
    放大：对整数化的美化产物（最大系数 ~500）可达 0.1–15 eq，于是
    `99SCN^- + 52Fe^{3+} -> 21[Fe(SCN)]^{2+} + 18[Fe(SCN)_3] +
    12[Fe(SCN)_2]^+`（缺 1Fe、电荷 −3）与银氨族（缺 6H/2N）照样通过，
    被固化成测试标准。

    但**只把整数档改成精确还不够**：本函数的浮点档还被用作"痕量过滤是否
    允许"的守卫（L489/L903），3% 相对容差会让 2.3% 的真实失衡通过——
    FeCl₃+KSCN 的 `_filter_trace(frac=0.05)` 删掉 `Fe(OH)₃ 0.0118 mol`
    （占 Fe 的 2.3%）后守卫放行，净方程就少了 1 个 Fe
    （`8.231SCN^- + 4.351Fe^{3+} -> …`，Fe 差 0.097）。
    非整数档的正确容差是**报告量子**：`consumption/production` 由
    `round(x, 6)` 产出，每个系数误差 ≤5e-7，故绝对容差 `_EQ_TOL_ABS`
    （1e-5，含物种数余量）才是"守恒到打印精度"的判据。
    显式传 `tol` 的调用点（整数 snap / 美化投影）保持原相对口径。
    """
    species = list(consumed) + list(produced)
    if not species:
        return True
    vals = list(consumed.values()) + list(produced.values())
    if all(float(v).is_integer() for v in vals):
        els = set()
        for sp in species:
            els |= set(elements_of(sp))
        for el in els:
            lhs = sum(elements_of(s).get(el, 0) * int(n)
                      for s, n in consumed.items())
            rhs = sum(elements_of(s).get(el, 0) * int(n)
                      for s, n in produced.items())
            if lhs != rhs:
                return False
        return (sum(charge_of(s) * int(n) for s, n in consumed.items())
                == sum(charge_of(s) * int(n) for s, n in produced.items()))
    mx = max(vals, default=1.0)
    lim = (tol * mx) if tol is not None else _EQ_TOL_ABS
    els = set()
    for sp in species:
        els |= set(elements_of(sp))
    for el in els:
        lhs = sum(elements_of(s).get(el, 0) * n for s, n in consumed.items())
        rhs = sum(elements_of(s).get(el, 0) * n for s, n in produced.items())
        if abs(lhs - rhs) > lim:
            return False
    lhs = sum(charge_of(s) * n for s, n in consumed.items())
    rhs = sum(charge_of(s) * n for s, n in produced.items())
    return abs(lhs - rhs) <= lim


def _h_o_balanced(c: dict, p: dict) -> bool:
    """浮点净差的 H/O 是否已配平（显式水档 vs 隐式水档的判据）。"""
    for el in ("H", "O"):
        lhs = sum(elements_of(s).get(el, 0) * v for s, v in c.items())
        rhs = sum(elements_of(s).get(el, 0) * v for s, v in p.items())
        if abs(lhs - rhs) > _EQ_TOL_ABS:
            return False
    return True


def _stoichiometric_net(c: dict, p: dict) -> tuple[dict, dict]:
    """净差 → 整数化学计量 + 隐式水补齐（H/O 未配平的净差专用候选）。

    **为什么必须"先整数化再补水"**：账本净差的 H/O 由溶剂水承担（守恒契约
    只管非 H/O 元素 + 电荷），所以浮点净差通常**无法**同时满足"H/O 配平 +
    电荷配平"——比例里带 6e-4 残差时，补 H₂O 必然还要补 H⁺，电荷随即破契约
    （N15 实测：补水候选要求 H⁺ 0.001215 > 契约容差 1e-3 ⟹ 闸门拒绝），
    呈现于是退化成 `1.999CO_2 + SiO_3^{2-} -> 1.999HCO_3^- + H_2SiO_3`：
    既不是整数、也没有水，还差 6e-4 个电荷——正是"呈现把误差藏起来"的样子。

    整数向量在非 H/O 元素 + 电荷上是**精确**守恒的（这是化学方程式的硬
    条件，比任何容差都强），此时补水恰好是若干个 H₂O（不引入 H⁺/OH⁻），
    于是得到教科书形式 `2CO_2 + 2H_2O + SiO_3^{2-} -> 2HCO_3^- + H_2SiO_3`。

    候选生成器的自我约束（闸门之外）：每个整数系数必须落在原浮点值的
    **呈现量子** `_contract_tol` 内（最小二乘拟合尺度），否则说明"这不是
    同一个方程"，原样返回。系数上限 40（超过即放弃整数化，诚实呈现浮点）。
    """
    names_c, names_p = list(c), list(p)
    vals = [c[s] for s in names_c] + [p[s] for s in names_p]
    pos = [v for v in vals if v > 0.0]
    if not pos:
        return c, p
    mx = max(pos)
    sig = [v for v in pos if v >= 1e-3 * mx]
    scale0 = min(sig) if sig else min(pos)
    tol = _contract_tol(mx)
    ints: list[int] = []
    for v in vals:
        k = int(round(v / scale0))
        if k > 40:
            return c, p
        ints.append(max(k, 0))
    if not any(ints[:len(names_c)]) or not any(ints[len(names_c):]):
        return c, p
    # 最小二乘尺度：s = Σ(v·k)/Σ(k²)；逐系数偏差必须落在呈现量子内
    num = sum(v * k for v, k in zip(vals, ints))
    den = sum(k * k for k in ints)
    if den == 0:
        return c, p
    s = num / den
    if s <= 0.0 or any(abs(v - s * k) > tol for v, k in zip(vals, ints)):
        return c, p
    # 精确整数守恒：非 H/O 元素 + 电荷（化学方程式的硬条件）
    els: set = set()
    for sp in names_c + names_p:
        els |= set(elements_of(sp))
    els -= {"H", "O"}
    for el in els:
        lhs = sum(k * elements_of(s).get(el, 0)
                  for s, k in zip(names_c, ints[:len(names_c)]))
        rhs = sum(k * elements_of(s).get(el, 0)
                  for s, k in zip(names_p, ints[len(names_c):]))
        if lhs != rhs:
            return c, p
    if (sum(k * charge_of(s) for s, k in zip(names_c, ints[:len(names_c)]))
            != sum(k * charge_of(s)
                   for s, k in zip(names_p, ints[len(names_c):]))):
        return c, p
    c2 = {s: float(k) for s, k in zip(names_c, ints[:len(names_c)]) if k > 0}
    p2 = {s: float(k) for s, k in zip(names_p, ints[len(names_c):]) if k > 0}
    # 隐式水补齐（整数向量已电荷守恒 ⟹ 只需 H₂O；仍走通用补配平并回验）
    c3, p3 = _balance_h_o(dict(c2), dict(p2))
    c3 = {s: v for s, v in c3.items() if abs(v) > 1e-9}
    p3 = {s: v for s, v in p3.items() if abs(v) > 1e-9}
    if not c3 or not p3:
        return c, p
    if any(abs(v - round(v)) > 1e-6 for v in list(c3.values()) + list(p3.values())):
        return c, p                      # 补配平引入了分数 ⟹ 不是干净计量式
    if not (_h_o_balanced(c3, p3) and _balanced_quick(c3, p3)):
        return c, p
    return c3, p3


def _filter_trace(consumed: dict, produced: dict,
                  frac: float = 0.02,
                  balanced: bool = False) -> tuple[dict, dict]:
    """过滤净差中的痕量物种：< frac × 最大物种量视为噪声。

    `frac` 是**相对**阈值（稀体系（1e-4 M 级别）的小量反应不被绝对地板
    误杀）。但"小于 2% 峰值"不等于"可删"——删项破坏元素守恒，而守恒是
    净方程的硬要求。实测 FeCl₃+KSCN 的 `Fe(OH)₃ 0.0118 mol` 占峰值 1.2%
    却在 Fe 上占 2.3%，整批删掉后净方程少 1 个 Fe
    （`8.231SCN^- + 4.351Fe^{3+} -> …`），并被固化成 4 条测试标准。

    `balanced=True` 时改为**逐项守恒感知删除**：把候选按量升序逐个试删，
    删后仍满足守恒（`_EQ_TOL_ABS`，即报告量子）才真删；否则保留。
    于是 1e-6 级的真噪声照删，而 0.0118 的 Fe(OH)₃ 被留下——既清噪又不破坏守恒。
    """
    mx = max(list(consumed.values()) + list(produced.values()), default=0.0)
    thr = max(1e-6, frac * mx)
    if not balanced:
        return ({s: v for s, v in consumed.items() if v > thr},
                {s: v for s, v in produced.items() if v > thr})
    cands: list[tuple[float, str, bool]] = []      # (量, 物种, 是否在消耗侧)
    for s, v in consumed.items():
        if v <= thr:
            cands.append((v, s, True))
    for s, v in produced.items():
        if v <= thr:
            cands.append((v, s, False))
    c2, p2 = dict(consumed), dict(produced)
    for _v, s, in_cons in sorted(cands):
        c3, p3 = dict(c2), dict(p2)
        side = c3 if in_cons else p3
        side.pop(s, None)
        if _balanced_quick(c3, p3):
            c2, p2 = c3, p3
    return c2, p2

def _cancel_minor_hydrolysis(consumed: dict, produced: dict,
                             frac: float = 0.05) -> tuple[dict, dict]:
    """消去微量水解/溶沉记账噪声：氢氧化物固体与其配比的 H+ 成对微量
    出现时整体消去（FeCl2+Cl2 中 1% 的 Fe3+ 水解副反应不属于净方程）。
    对消后 H/O 差额由 balance_h_o 的 H2O 补齐，恰好回到教科书形式。
    仅当两者都相对主物种微量（<5%）时才消——纯水解体系（AlCl3 溶液）
    的水解本身就是主反应，不动。

    v0.4.5 守恒守卫：对消掉 `x` mol 氢氧化物固体等于从方程里删掉 `x` mol
    金属，**只在 x 低于报告量子时才允许**（否则方程少一个金属）。
    实测 FeCl₃+KSCN 的 `Fe(OH)₃ 0.0118 mol`（占 Fe 的 2.3%）被这里对消，
    净方程少 1 个 Fe（`8.231SCN^- + 4.351Fe^{3+} -> …`），并被固化成 4 条
    标准；round-4 曾误判为 `_filter_trace` 所致（过滤器的**输入**已缺该项，
    真凶在这一层）。
    """
    mx = max(list(consumed.values()) + list(produced.values()), default=0.0)
    if mx <= 0:
        return consumed, produced
    limit = frac * mx
    for e in TABLES.ksp:
        if e["pair"][1] != OH_ION:
            continue
        solid = e["solid"]
        n = charge_of(e["pair"][0])
        for side_solid, side_ion, ion in ((produced, produced, H_ION),
                                          (consumed, consumed, H_ION),
                                          (produced, produced, OH_ION),
                                          (consumed, consumed, OH_ION)):
            x = side_solid.get(solid, 0.0)
            if not (1e-9 < x < limit):
                continue
            y = side_ion.get(ion, 0.0)
            if y < limit and abs(y - n * x) <= 0.25 * max(n * x, 1e-9):
                tc, tp = dict(consumed), dict(produced)
                t_solid = tc if side_solid is consumed else tp
                t_ion = tc if side_ion is consumed else tp
                t_solid.pop(solid, None)
                rest = y - n * x
                if rest > 1e-9:
                    t_ion[ion] = rest
                else:
                    t_ion.pop(ion, None)
                    if rest < -1e-9:      # H+ 不足配比：差额挂到对侧
                        other = tp if side_ion is consumed else tc
                        other[ion] = other.get(ion, 0.0) + (-rest)
                if not _balanced_quick(tc, tp):
                    continue              # 对消会破坏守恒 → 保留该固相
                consumed.clear()
                consumed.update(tc)
                produced.clear()
                produced.update(tp)
                break
    return consumed, produced


def _cancel_minor_protonation(consumed: dict, produced: dict,
                              frac: float = 0.05) -> tuple[dict, dict]:
    """消去微量质子转移记账噪声（净离子方程式不写酸碱形态微调）。

    数据表中的每个 pKa 共轭对 acid ⇌ base + H+ 给出四种噪声模式（净差中
    三者按 1:1:1 同量出现且相对主物种微量时整体消去）：
      base + H+ → acid：消耗 H+/base，生成 acid
      acid → base + H+：消耗 acid，生成 base/H+
      base + H2O → acid + OH-：消耗 base，生成 acid/OH-
      acid + OH- → base + H2O：消耗 acid/OH-，生成 base
    """
    mx = max(list(consumed.values()) + list(produced.values()), default=0.0)
    if mx <= 0:
        return consumed, produced
    limit = frac * mx

    def _take(x, *side_sp):
        if x <= 1e-9 or x >= limit:
            return
        for side, sp in side_sp:
            side[sp] -= x
            if side[sp] <= 1e-9:
                side.pop(sp, None)

    for e in TABLES.pka:
        if e.get("n", 1) != 1:
            continue
        acid, base = e["acid"], e["base"]
        _take(min(consumed.get(H_ION, 0.0), produced.get(acid, 0.0),
                  consumed.get(base, 0.0)),
              (consumed, H_ION), (produced, acid), (consumed, base))
        _take(min(produced.get(H_ION, 0.0), consumed.get(acid, 0.0),
                  produced.get(base, 0.0)),
              (produced, H_ION), (consumed, acid), (produced, base))
        _take(min(consumed.get(base, 0.0), produced.get(acid, 0.0),
                  produced.get(OH_ION, 0.0)),
              (consumed, base), (produced, acid), (produced, OH_ION))
        _take(min(consumed.get(acid, 0.0), consumed.get(OH_ION, 0.0),
                  produced.get(base, 0.0)),
              (consumed, acid), (consumed, OH_ION), (produced, base))
        # 共轭对直接互转（质子经第三方传递，无净 H+/OH-）：
        # acid → base（如 NH4+ 供质子给沉淀溶解）：消耗 acid、生成 base
        _take(min(consumed.get(acid, 0.0), produced.get(base, 0.0)),
              (consumed, acid), (produced, base))
        # base → acid：消耗 base、生成 acid
        _take(min(consumed.get(base, 0.0), produced.get(acid, 0.0)),
              (consumed, base), (produced, acid))
    return consumed, produced


def _align_conjugate_pair(consumed: dict, produced: dict) -> tuple[dict, dict]:
    """痕量过滤后的共轭酸碱对齐：过滤掉 [Ag(NH3)2]+ 等痕量产物后，
    消耗侧 NH3 与产物侧 NH4+ 出现差额（来自被过滤的产物）。

    化学事实：被过滤的 [Ag(NH3)2]+ 系由 2NH3 + Ag+ 形成（消耗侧 NH3），
    过滤后消耗侧 NH3 残留 = 真实电离量 + 2×[Ag(NH3)2]+ 量。同步调整
    较小者到较大者（基于差额 < 5% 主物种阈值），以恢复 1:1 化学计量比。
    """
    if not consumed or not produced:
        return consumed, produced
    mx = max(list(consumed.values()) + list(produced.values()), default=0.0)
    if mx <= 0:
        return consumed, produced
    limit = 0.10 * mx  # 10% 阈值（放宽以容纳过滤痕量后的不平衡）
    for e in TABLES.pka:
        if e.get("n", 1) != 1:
            continue
        acid, base = e["acid"], e["base"]
        # base 消耗、acid 产生（如 NH3 → NH4+）
        if base in consumed and acid in produced:
            b = consumed[base]
            a = produced[acid]
            diff = abs(b - a)
            if diff < limit and diff > 1e-9:
                # 取较小者；较大者的差额视为被过滤产物的"幽灵"消耗
                m = min(b, a)
                consumed[base] = m
                produced[acid] = m
                if m <= 1e-9:
                    consumed.pop(base, None)
                    produced.pop(acid, None)
        # acid 消耗、base 产生（如 HAc → Ac-）
        if acid in consumed and base in produced:
            a = consumed[acid]
            b = produced[base]
            diff = abs(b - a)
            if diff < limit and diff > 1e-9:
                m = min(a, b)
                consumed[acid] = m
                produced[base] = m
                if m <= 1e-9:
                    consumed.pop(acid, None)
                    produced.pop(base, None)
    consumed = {s: v for s, v in consumed.items() if v > 1e-9}
    produced = {s: v for s, v in produced.items() if v > 1e-9}
    return consumed, produced


# ============================================================ Reaction

