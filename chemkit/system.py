"""chemkit.system：高层 API —— Reaction（一次反应）与 System（反应体系）。

两个对外对象：
    Reaction  一次反应的结果。数据分两层：
        - 非 raw（化学习惯、人类可读）：consumed / produced / final 中含
          H+ / OH- / H2O，方程式为净离子方程式，系数配平到最简整数比。
        - xxx_raw（引擎原始记账）：H2O 是溶剂不入账，H+/OH- 合记为单一
          带符号账本 H_excess_raw（正 = 残余游离强酸，负 = 残余游离强碱），
          consumed_raw / produced_raw / final_raw 不含 H+/OH-/H2O。
    System    反应体系：固定体积/温度/气压，add() 连续投料，每次投料后
              按累计投入量重新平衡，返回本次 Reaction。

一步式便捷函数：
    chemkit.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0) -> Reaction

方程式：
    reaction.equation   净离子方程式（全部步骤的净和），无显著反应时为 None
    reaction.equations  多步离子方程式列表（按贡献降序）——许多反应用多步
                        概括更贴近化学书写习惯（如多元弱酸分步中和、
                        沉淀-解配耦合），净方程难以表达时使用。
"""
from __future__ import annotations

import re
from fractions import Fraction
from math import gcd
from functools import reduce

from .data import Tables, load_tables
from .engine import judge
from .core import elements_of, charge_of

# ---------- 数据表预加载（import chemkit 时执行，全程序仅一次磁盘 I/O）----------
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
_ARROW = " -> "

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


def _parse_equation(eq: str) -> tuple[dict[str, float], dict[str, float]]:
    """解析 '2A + 3B -> 4C + D' → ({A:2, B:3}, {C:4, D:1})。"""
    if _ARROW not in eq:
        return {}, {}
    lhs, rhs = eq.split(_ARROW, 1)
    return _parse_side(lhs), _parse_side(rhs)


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
                trace: float) -> tuple[dict[str, float], dict[str, float]]:
    """把引擎的 H+ 正则形还原为化学习惯的 OH- 写法，并抵消跨侧 H+/OH-。

    两步：
      1. OH- 还原：引擎把 OH- 记为 H2O(反应物) − H+(即产物 H+)，
         如 M + nH2O -> M(OH)n + nH+ 实为 M + nOH- -> M(OH)n。
         H2O 在反应物侧且 H+ 在产物侧且量匹配时，合并为 OH- 到反应物侧。
    只有"H2O 在反应物侧 + H+ 在产物侧"这一种挂侧可以还原为 OH-；
    跨侧 H+/OH- 抵消（H+ 左 OH- 右，或反之）在数学上不等价于生成
    H2O（净向量差 2H+ 或 2OH-），会破坏配平，故不做。
    """
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


def _rationalize(vals: list[float]) -> list[int] | None:
    """浮点系数列表 → 最简整数比列表。

    两档精度逐级尝试（按最小值归一化后）：
      1. 细档：limit_denominator(8)（允许 1/8 级分数），相对误差 ≤2%——
         保留 25:17 这类真实非整数化学计量（CaCO3 溶解的 HCO3-/CO3^2- 分配），
         同时把 17.48 这类数值噪声收进 17.5；
      2. 粗档：limit_denominator(2)（整数/半整数），相对误差 ≤5%——
         用于噪声更大的近边界体系。
    两档都失败返回 None（调用方退化为浮点系数）。
    """
    pos_vals = [v for v in vals if v > 0]
    if not pos_vals:
        return None
    min_v = min(pos_vals)
    norm = [v / min_v for v in vals]
    for den_max, tol in ((8, 0.02), (2, 0.05)):
        fracs = []
        ok = True
        for v in norm:
            if v <= 0:
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
        return [v // g for v in int_vals]
    return None


def _format_equation(consumed: dict[str, float], produced: dict[str, float]) -> str | None:
    """将 consumed/produced 格式化为最简整数比的离子方程式字符串。

    物种排序：先按系数降序，再按名称字典序（确定性输出，便于测试断言）。
    系数比例无法整除（最大系数 >1000）时退化为浮点系数。
    """
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
    int_vals = _rationalize(all_vals)
    if int_vals is None or max(int_vals) > 1000:
        # 比例不整除（真实的非化学计量混合，如 Fe→Fe2+/Fe3+ 混合价态）：
        # 退化为浮点系数（按最小值归一化）
        scale = min(all_vals)
        if scale <= 0:
            return None
        cons_items = sorted(((sp, v / scale) for sp, v in consumed.items()),
                            key=lambda x: (-x[1], x[0]))
        prod_items = sorted(((sp, v / scale) for sp, v in produced.items()),
                            key=lambda x: (-x[1], x[0]))
    else:
        n_cons = len(consumed)
        cons_ints = int_vals[:n_cons]
        prod_ints = int_vals[n_cons:]
        cons_items = sorted(zip(consumed.keys(), cons_ints),
                            key=lambda x: (-x[1], x[0]))
        prod_items = sorted(zip(produced.keys(), prod_ints),
                            key=lambda x: (-x[1], x[0]))
    lhs = _SIDE_SEP.join(_fmt_term(v, sp) for sp, v in cons_items)
    rhs = _SIDE_SEP.join(_fmt_term(v, sp) for sp, v in prod_items)
    return f"{lhs}{_ARROW}{rhs}"


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
    return [(reps[k][0], reps[k][1], ext) for k, ext in agg.items() if ext >= STEP_MIN]


def _step_to_ionic(reactants: dict, products: dict, extent: float
                   ) -> tuple[dict, dict, str | None]:
    """单个聚合步骤 → (consumed, produced, 方程式)。含 H2O 配平与 OH- 还原。"""
    consumed = {sp: nu * extent for sp, nu in reactants.items()}
    produced = {sp: nu * extent for sp, nu in products.items()}
    consumed, produced = _balance_h_o(consumed, produced)
    consumed, produced = _restore_oh(consumed, produced, _TRACE)
    if not consumed or not produced:
        return {}, {}, None
    return consumed, produced, _format_equation(consumed, produced)


# 净方程式中把分子态强酸改写为离子形（浓酸记账形态 → 化学习惯离子形）
_NET_ACID_SPLIT = {
    "HNO_3": {H_ION: 1.0, "NO_3^-": 1.0},
    "H_2SO_4": {H_ION: 2.0, "SO_4^{2-}": 1.0},
}


def _build_equations(steps: list[dict], r: dict) -> tuple[list[str], dict, dict, str | None]:
    """构建多步离子方程式列表 + 净离子方程式。

    返回 (equations, consumed, produced, net_equation)：
      equations     各显著步骤的离子方程式（聚合、对消、按 extent 降序）；
      consumed/produced  净消耗/生成（mol，含 H+/OH-/H2O）——来自账本净差
                    （consumed_raw/produced_raw + H_excess 变化 + 中和步），
                    是真实摩尔量，振荡与中间体天然抵消；
      net_equation  净离子方程式（上述净差的格式化）。
    """
    # ---- 多步方程式（步骤聚合路径）----
    collected = _collect_steps(steps)
    equations: list[str] = []
    if collected:
        main_ext = max(ext for _, _, ext in collected)
        kept = [c for c in collected if c[2] >= main_ext * 0.05]
        kept.sort(key=lambda c: -c[2])
        seen_eq: set[str] = set()
        for rr, p, ext in kept:
            _, _, eq = _step_to_ionic(rr, p, ext)
            if eq and eq not in seen_eq:
                seen_eq.add(eq)
                equations.append(eq)

    # ---- 净方程（账本净差路径）----
    consumed: dict[str, float] = {e["name"]: e["mol"] for e in r.get("consumed", [])}
    produced: dict[str, float] = {e["name"]: e["mol"] for e in r.get("produced", [])}
    # 分子态强酸 → 离子形（HNO3 ⇌ H+ + NO3- 等）
    for acid, ions in _NET_ACID_SPLIT.items():
        for side in (consumed, produced):
            x = side.pop(acid, 0.0)
            if x > 0.0:
                for sp, nu in ions.items():
                    side[sp] = side.get(sp, 0.0) + nu * x
    # 质子账本净变化：He 是带符号账本（正=游离强酸），|dHe| 的离子形态
    # 由终态 pH 决定——偏酸体系记 H+、偏碱体系记 OH-（He 微小残差是记账
    # 幻影，终态 pH 才是真实的酸碱面貌）
    He_i = r.get("H_excess_initial", 0.0)
    He_f = r.get("H_excess", 0.0)
    pH_f = r.get("final_pH")
    basic = pH_f is not None and pH_f > 7.0
    dHe = He_i - He_f          # >0: 体系酸性减弱；<0: 酸性增强
    if dHe > 1e-9:
        if basic:
            produced[OH_ION] = produced.get(OH_ION, 0.0) + dHe
        else:
            consumed[H_ION] = consumed.get(H_ION, 0.0) + dHe
    elif dHe < -1e-9:
        if basic:
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
        return equations, {}, {}, None

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
        return equations, {}, {}, None
    # 痕量过滤 → 配平 → OH- 还原 → 终过滤。配平若被迫动用 H+ 兜底，
    # 说明过滤丢了携 H 物种，退回未过滤净差重来。
    raw_c, raw_p = dict(consumed), dict(produced)
    fc, fp = _filter_trace(consumed, produced)
    if not fc or not fp:
        # 过滤把一侧清空（全小量体系，如痕量溶解）：放弃过滤
        fc, fp = dict(consumed), dict(produced)
    consumed, produced = fc, fp
    c2, p2 = _balance_h_o(dict(consumed), dict(produced))
    if (H_ION in c2 and H_ION not in consumed) or (H_ION in p2 and H_ION not in produced):
        consumed, produced = raw_c, raw_p
        c2, p2 = _balance_h_o(dict(consumed), dict(produced))
    consumed, produced = c2, p2
    consumed, produced = _restore_oh(consumed, produced, 1e-9)
    # 配平兜底可能在对侧补出 H+：再做一次跨侧抵消
    for sp in set(consumed) & set(produced):
        x = min(consumed[sp], produced[sp])
        consumed[sp] -= x
        produced[sp] -= x
    # 终过滤（5%）：过滤后必须仍配平，否则退回未过滤版本
    consumed = {s: v for s, v in consumed.items() if v > 1e-9}
    produced = {s: v for s, v in produced.items() if v > 1e-9}
    fc, fp = _filter_trace(consumed, produced, frac=0.05)
    if fc and fp and _balanced_quick(fc, fp):
        consumed, produced = fc, fp
    if not consumed or not produced:
        return equations, {}, {}, None
    return equations, consumed, produced, _format_equation(consumed, produced)


def _balanced_quick(consumed: dict, produced: dict, tol: float = 0.03) -> bool:
    """快速配平验证：元素与电荷两侧相等（相对容差 3%）。"""
    species = list(consumed) + list(produced)
    mx = max(list(consumed.values()) + list(produced.values()), default=1.0)
    els = set()
    for sp in species:
        els |= set(elements_of(sp))
    for el in els:
        lhs = sum(elements_of(s).get(el, 0) * n for s, n in consumed.items())
        rhs = sum(elements_of(s).get(el, 0) * n for s, n in produced.items())
        if abs(lhs - rhs) > tol * mx:
            return False
    lhs = sum(charge_of(s) * n for s, n in consumed.items())
    rhs = sum(charge_of(s) * n for s, n in produced.items())
    return abs(lhs - rhs) <= tol * mx


def _filter_trace(consumed: dict, produced: dict,
                  frac: float = 0.02) -> tuple[dict, dict]:
    """过滤净差中的痕量物种：< 2% × 最大物种量视为噪声（纯相对阈值，
    稀体系（1e-4 M 级别）的小量反应不被绝对地板误杀）。"""
    mx = max(list(consumed.values()) + list(produced.values()), default=0.0)
    thr = max(1e-6, frac * mx)
    return ({s: v for s, v in consumed.items() if v > thr},
            {s: v for s, v in produced.items() if v > thr})


def _cancel_minor_hydrolysis(consumed: dict, produced: dict,
                             frac: float = 0.05) -> tuple[dict, dict]:
    """消去微量水解/溶沉记账噪声：氢氧化物固体与其配比的 H+ 成对微量
    出现时整体消去（FeCl2+Cl2 中 1% 的 Fe3+ 水解副反应不属于净方程）。
    对消后 H/O 差额由 balance_h_o 的 H2O 补齐，恰好回到教科书形式。
    仅当两者都相对主物种微量（<5%）时才消——纯水解体系（AlCl3 溶液）
    的水解本身就是主反应，不动。"""
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
                side_solid.pop(solid, None)
                rest = y - n * x
                if rest > 1e-9:
                    side_ion[ion] = rest
                else:
                    side_ion.pop(ion, None)
                    if rest < -1e-9:      # H+ 不足配比：差额挂到对侧
                        other = produced if side_ion is consumed else consumed
                        other[ion] = other.get(ion, 0.0) + (-rest)
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


# ============================================================ Reaction

class Reaction:
    """一次反应的结果（一次投料平衡后的完整快照）。

    非 raw 属性（化学习惯、人类可读，含 H+/OH-/H2O）：
        reacted       体系是否发生显著净变化（bool，净账本判定）
        chemical_reaction  其中是否包含狭义化学反应（bool：氧化还原/
                      中和/跨投料的沉淀配位等；纯溶解、电离、水解等
                      单投料形态变化为 False）
        degree        程度：complete / incomplete / hardly / none
        consumed      净消耗 {化学式: mol}
        produced      净生成 {化学式: mol}
        final         终态组成 {化学式: mol}（H2O 为溶剂不入）
        pH            终态 pH（None 表示不适用，如 OVERRIDE 路径）
        equation      净离子方程式字符串，无显著反应时为 None
        equations     多步离子方程式列表（按贡献降序）
        annotations   标注列表（slow / blocked 等）
        override      命中的 OVERRIDE id（或 None）
        escaped       逸出气相 {化学式: mol}（泡点扫气：超过 H(T)·p_ext
                      溶解上限的自产气体，逸出即消失，不再参与反应；
                      final 中同种气体只剩溶解态，produced 含逸出部分）

    raw 属性（引擎原始记账，H2O 不入账、H+/OH- 合记为 H_excess）：
        consumed_raw  消耗 {化学式: mol}
        produced_raw  生成 {化学式: mol}
        final_raw     终态 {化学式: mol}
        H_excess_raw  终态带符号质子账本（正 = 残余游离强酸 mol，
                      负 = 残余游离强碱 mol）
        steps         逐步反应过程（kind/equation/logK/S/extent/conversion）
        raw           judge() 原始 dict（备用，不鼓励直接读取）
    """
    __slots__ = (
        # 非 raw
        "reacted", "chemical_reaction", "degree", "consumed", "produced",
        "final", "pH", "annotations", "override", "escaped",
        # raw
        "consumed_raw", "produced_raw", "final_raw", "H_excess_raw", "raw",
        # 缓存
        "_equation", "_equations",
    )

    def __init__(self, r: dict):
        self.raw = r
        self.reacted: bool = r["reacted"]
        self.chemical_reaction: bool = r.get("chemical_reaction", r["reacted"])
        self.degree: str = r["degree"]
        self.pH: float | None = r["final_pH"]
        self.annotations: list[str] = list(r["annotations"])
        self.override: str | None = r.get("override")
        self.escaped: dict[str, float] = {e["name"]: e["mol"]
                                          for e in r.get("escaped", [])}

        # ---- raw（引擎记账）----
        self.consumed_raw: dict[str, float] = {e["name"]: e["mol"] for e in r["consumed"]}
        self.produced_raw: dict[str, float] = {e["name"]: e["mol"] for e in r["produced"]}
        self.final_raw: dict[str, float] = {e["name"]: e["mol"] for e in r["final"]}
        self.H_excess_raw: float = r.get("H_excess", 0.0)

        # ---- 非 raw（化学习惯）----
        self._equations, self.consumed, self.produced, self._equation = \
            _build_equations(r["steps"], r)

        # OVERRIDE 路径无 steps，退化为 raw 层面的化学式方程式
        if (not self.consumed and not self.produced and self.consumed_raw
                and self.override):
            self.consumed = dict(self.consumed_raw)
            self.produced = dict(self.produced_raw)
            self._equation = _format_equation(self.consumed, self.produced)
            if self._equation and not self._equations:
                self._equations = [self._equation]

        # 终态组成：raw final + H+/OH-（H2O 为溶剂不入）
        self.final: dict[str, float] = dict(self.final_raw)
        if self.H_excess_raw > 1e-6:
            self.final[H_ION] = self.final.get(H_ION, 0.0) + self.H_excess_raw
        elif self.H_excess_raw < -1e-6:
            self.final[OH_ION] = self.final.get(OH_ION, 0.0) + (-self.H_excess_raw)

    # ---- 便捷转发 ----
    @property
    def steps(self) -> list[dict]:
        """逐步反应过程（引擎记账）。每项含 kind/equation/logK/S/extent/
        conversion 等字段。"""
        return self.raw["steps"]

    @property
    def equation(self) -> str | None:
        """净离子方程式（如 'Zn + 2H^+ -> H_2 + Zn^{2+}'）。

        全部显著步骤的净和：中间体自然抵消，H2O 显式配平，OH- 从引擎的
        H+ 正则形还原，系数有理化到最简整数比。无显著反应时返回 None。
        """
        return self._equation

    @property
    def equations(self) -> list[str]:
        """多步离子方程式列表（按贡献降序）。

        许多反应用多步概括更自然：如 Ca(OH)2+CO2 是
        'CO2 + 2OH^- -> CO_3^{2-} + H2O' 与 'Ca^{2+} + CO_3^{2-} -> CaCO3'
        两步，净方程也能写但不直观。无显著反应时为空列表。
        """
        return list(self._equations)

    # ---- 魔术方法 ----
    def __bool__(self) -> bool:
        return self.reacted

    def __repr__(self) -> str:
        pro = ", ".join(f"{k}×{v:.3g}" for k, v in self.produced.items())
        return (f"<Reaction {'反应' if self.reacted else '不反应'} "
                f"{self.degree} [{pro}] pH={self.pH}>")


# ============================================================ System

class System:
    """反应体系：固定体积/温度/气压，支持连续投料。

    参数：
        substances  初始投料 {化学式: mol}，默认空（纯水）
        V           溶液体积（L），默认 1.0
        T           温度（K），默认 298.15（25 °C）
        T_C         温度（°C），若给定则覆盖 T
        p           外界气压（kPa），默认 101.3；影响气体逸出阈值
        tables      自定义数据表（默认用模块级 TABLES）

    建立时（若给了 substances）与每次 add() 自动触发反应——按累计投入量
    重新平衡（化学上等价于连续投料的再平衡），返回本次 Reaction。
    全部历史结果保存在 history。

    示例：
        sys = chemkit.System(V=1.0)
        sys.add("NaOH", 0.1)              # 纯水 + 0.1 mol NaOH
        r = sys.add("HCl", 0.15)          # 再投入 HCl，体系重新平衡
        r.equation, r.pH
    """

    def __init__(self,
                 substances: dict[str, float] | None = None,
                 *,
                 V: float = 1.0,
                 T: float = 298.15,
                 T_C: float | None = None,
                 p: float = 101.3,
                 tables: Tables | None = None):
        self.V_L: float = float(V)
        self.T_K: float = float(T) if T_C is None else float(T_C) + 273.15
        self.p_kpa: float = float(p)
        self._tables: Tables = tables if tables is not None else TABLES
        self._feeds: dict[str, float] = {}
        self.history: list[Reaction] = []
        self.result: Reaction | None = None
        if substances:
            for name, mol in substances.items():
                self._feeds[name] = self._feeds.get(name, 0.0) + float(mol)
            self._react()

    def add(self, name: str, mol: float) -> Reaction:
        """加入物质（mol），自动触发反应，返回本次 Reaction。"""
        self._feeds[name] = self._feeds.get(name, 0.0) + float(mol)
        return self._react()

    def _react(self) -> Reaction:
        subs = [{"name": n, "mol": m} for n, m in self._feeds.items()]
        cond = {"V_L": self.V_L, "T_K": self.T_K, "p_kpa": self.p_kpa}
        r = judge(subs, cond, self._tables)
        self.result = Reaction(r)
        self.history.append(self.result)
        return self.result

    @property
    def feeds(self) -> dict[str, float]:
        """累计投料 {化学式: mol}（副本，外部修改不影响内部状态）。"""
        return dict(self._feeds)

    def __repr__(self) -> str:
        fd = ", ".join(f"{k}×{v:.3g}" for k, v in self._feeds.items())
        return f"<System [{fd}] V={self.V_L}L T={self.T_K}K p={self.p_kpa}kPa>"


# ============================================================ 一步式 API

def react(substances: dict[str, float],
          *,
          V: float = 1.0,
          T: float = 298.15,
          T_C: float | None = None,
          p: float = 101.3,
          tables: Tables | None = None) -> Reaction:
    """一步式反应：建立体系并立即反应，返回 Reaction。

    参数与 System 一致（substances 必填）。

    示例：
        r = chemkit.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0)
        print(r.equation)   # 'Zn + 2H^+ -> H_2 + Zn^{2+}'
    """
    return System(substances, V=V, T=T, T_C=T_C, p=p, tables=tables).result
