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
from .candidates import WATER, H_ION, X_MIN

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
            return nc, np_
    return None


def _format_equation(consumed: dict[str, float], produced: dict[str, float]) -> str | None:
    """将 consumed/produced 格式化为最简整数比的离子方程式字符串。

    物种排序：先按系数降序，再按名称字典序（确定性输出，便于测试断言）。
    优先级：① 整数 snap（处理 equilibrium 残差，干净小整数）；
            ② _rationalize 分数有理化（处理真实非整数比，如 1:2:1.5）；
            ③ 浮点系数兜底（真实非化学计量混合，如 Fe→Fe2+/Fe3+）。
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

    # ① 整数 snap 优先（最干净的化学计量比，处理 equilibrium 残差）
    snap_ints = _try_integer_snap(consumed, produced)
    if snap_ints is not None:
        n_cons = len(consumed)
        cons_ints = snap_ints[:n_cons]
        prod_ints = snap_ints[n_cons:]
        cons_items = sorted([(s, v) for s, v in zip(consumed.keys(), cons_ints) if v > 0],
                            key=lambda x: (-x[1], x[0]))
        prod_items = sorted([(s, v) for s, v in zip(produced.keys(), prod_ints) if v > 0],
                            key=lambda x: (-x[1], x[0]))
        lhs = _SIDE_SEP.join(_fmt_term(v, sp) for sp, v in cons_items)
        rhs = _SIDE_SEP.join(_fmt_term(v, sp) for sp, v in prod_items)
        return f"{lhs}{_ARROW}{rhs}"

    # ② _rationalize 分数有理化（真实非整数比，如 1:2:1.5）
    int_vals = _rationalize(all_vals)
    if int_vals is None or max(int_vals) > 1000:
        # ③ 浮点系数兜底（真实非化学计量混合，如 Fe→Fe2+/Fe3+ 混合价态）
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
    """单个聚合步骤 → (consumed, produced, 方程式)。含 H2O 配平与 OH- 还原。"""
    consumed = {sp: nu * extent for sp, nu in reactants.items()}
    produced = {sp: nu * extent for sp, nu in products.items()}
    consumed, produced = _balance_h_o(consumed, produced)
    consumed, produced = _restore_oh(consumed, produced, _TRACE)
    if not consumed or not produced:
        return {}, {}, None
    return consumed, produced, _format_equation(consumed, produced)


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


def _build_equations(steps: list[dict], r: dict) -> tuple[list[str], dict, dict, str | None]:
    """构建多步离子方程式列表 + 净离子方程式。

    返回 (equations, consumption, production, net_equation)：
      equations         各显著步骤的离子方程式（聚合、对消、按 extent 降序）；
      consumption/production  净消耗/生成（mol，含 H+/OH-/H2O）——来自账本净差
                        （consumption_raw/production_raw + H_excess 变化 + 中和步），
                        是真实摩尔量，振荡与中间体天然抵消；
      net_equation      净离子方程式（上述净差的格式化）。
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
    consumed: dict[str, float] = {e["name"]: e["mol"] for e in r.get("consumption", [])}
    produced: dict[str, float] = {e["name"]: e["mol"] for e in r.get("production", [])}
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
    He_i = r.get("H_excess_initial", 0.0)
    He_f = r.get("H_excess", 0.0)
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
        return equations, {}, {}, None
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
    if has_complex_product and net_total < 0.005 * feed_total and net_total < 5e-3:
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
    # 痕量过滤 → 共轭酸碱对质量平衡调整 → 配平 → OH- 还原 → 终过滤
    raw_c, raw_p = dict(consumed), dict(produced)
    fc, fp = _filter_trace(consumed, produced)
    if not fc or not fp:
        # 过滤把一侧清空（全小量体系，如痕量溶解）：放弃过滤
        fc, fp = dict(consumed), dict(produced)
    # 共轭酸碱对不平衡调整：过滤掉 [Ag(NH3)2]+ (痕量) 后，消耗侧 NH3 残留
    # 与产物侧 NH4+ 不等（差额来自被过滤的产物）。化学事实上被过滤的产物
    # 是由共轭碱+其他离子形成的（如 2NH3 + Ag+ → [Ag(NH3)2]+），差额
    # 残留在共轭碱侧。同步调整较小者到较大者以恢复化学计量比。
    fc, fp = _align_conjugate_pair(fc, fp)
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
    net_str = _format_equation(consumed, produced)
    # ---- 大系数美化（v0.3.8）：浮点系数或 max>20 的混合通道净差 →
    # 主通道呈现（移除次要物种 + 平衡子空间投影 + 小整数比；配平硬保证，
    # 见 _beautify_big_coeff）。找不到干净形式保持原样（诚实优先）。
    if net_str is not None and ("." in net_str
                                or _max_coef(net_str) > 20):
        bp = _beautify_big_coeff(consumed, produced)
        if bp is not None:
            bc, bp_ = bp
            # ≤3 物种才替换（教科书单通道形式：Pb²⁺+2Br⁻→PbBr₂ 类）；
            # ≥4 物种的"美化"仍是混合物（Fe/SCN 配位阶梯），不如保留
            # 诚实混合比（测试有意锁定定量热力学呈现，N30/W11 族）
            if len(bc) + len(bp_) <= 3:
                pretty = _format_equation(bc, bp_)
                if pretty is not None and _max_coef(pretty) <= 20:
                    net_str = pretty
    return equations, consumed, produced, net_str


def _max_coef(eq: str) -> float:
    """方程式字符串的最大系数（解析失败返回 0）。"""
    rr, pp = _parse_equation(eq)
    if not rr and not pp:
        return 0.0
    return max(list(rr.values()) + list(pp.values()))


def _balanced_quick(consumed: dict, produced: dict, tol: float = 0.03) -> bool:
    """配平验证：元素与电荷两侧相等。

    **整数系数走精确整数算术，非整数走相对容差**——这是"配平硬保证"的
    真正落点。0.4.x 一律用相对容差（默认 3%、美化路径 0.1%），容差被
    最大系数放大：对整数化的美化产物（最大系数 ~500）可达 0.1–15 eq，
    于是 `99SCN^- + 52Fe^{3+} -> 21[Fe(SCN)]^{2+} + 18[Fe(SCN)_3] +
    12[Fe(SCN)_2]^+`（缺 1Fe、电荷 −3）与银氨族（缺 6H/2N）照样通过，
    被固化成测试标准。整数候选没有"近似配平"的余地（系数本身就是
    化学计量），故零容差；浮点候选（真实非化学计量混合比）保留容差。
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

