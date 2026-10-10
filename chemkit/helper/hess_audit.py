"""Hess 一致性审计：候选（派生 + 静态模板）的 logK 是否与其组元的线性组合相符。

原理（这正是"派生行秩亏"的判据）：候选的净变换向量必然落在**基候选**
（beta/pKa/Ksp/水）张成的反应子空间里——它本来就是这样造出来的。于是对每条
候选 d，解线性方程组

    Σ_i a_i · vec(base_i) = vec(d)        （vec = 物种计量的向量）

若可解（基集确实张成 d），则必须同时满足

    logK(d, T) == Σ_i a_i · logK(base_i, T)

**行秩亏的含义**：派生行是基行的线性组合，所以任何"秩检查"都看不出问题；
但组合系数若错（例如少乘一个 pKw），logK 就与基行不自洽——方程组变得
矛盾，走步/联立在这个方向上永远无法收敛。H46 即此：beta_ksp 通道把
`n_OH` 取成中心离子电荷而非固相 OH 数，Ag2O（x=2,y=2）少减一个 pKw，
logK 偏 +14。

**温度维（§7 X-31）**：上面这组等式对**任意温度**都必须成立。只在 298.15 K
审计会漏掉整类缺陷——OH⁻→H⁺ 折算（pKw）若被同时烘进基值与 pkw_coeff，
298.15 K 恰好恒等（pKw(298.15) ≡ 14.0），偏差 ∝ (pKw(T)−14.0)，只有换温度
才现形。templates.py 的 OH⁻ 型沉淀/溶解/配位四条模板正是如此（363 K 偏
ν·1.573，K 差 10^ν·1.573）。所以 `--T` 是必查项，不是可选项。

温度项按**引擎自己的两个通道**复算，否则审计会与引擎的口径打架：
  · pKw 通道：水行的系数 a_w 贡献 a_w·(pKw(298) − pKw(T))；
  · dH 通道：组合的 dH = Σ a_i·dH_i，但**任一带系数组元缺 dH 则整条不修正**
    （= candidates._hess_dH 的"缺数据→None 整体回退"语义）。
数据表 dH 不全（ksp/beta 有缺项），若按逐行各自修正，缺数据行会与引擎的
"整体不修正"口径产生上千条假阳性——首版即栽在此。

用法：
    python chemkit/helper/hess_audit.py                       # 298.15 K，派生候选
    python chemkit/helper/hess_audit.py --T 363.15 --with-templates
    python chemkit/helper/hess_audit.py --T 363.15 --templates        # 只审静态模板族
    python chemkit/helper/hess_audit.py --T 363.15 --src Tl_2S --show # 逐条打印组合
    python chemkit/helper/hess_audit.py --exact                       # 有理精确复核（慢）
"""
from __future__ import annotations

import sys
from fractions import Fraction as F

sys.path.insert(0, ".")

from chemkit.candidates import (WATER, build_derived, _ksp_xy,  # noqa: E402
                                PKW_298, logK_T)
from chemkit.core import pKw_of, _vant                          # noqa: E402
from chemkit.data import load_tables                            # noqa: E402
from chemkit.templates import _build_static_cands               # noqa: E402

H_ION = "H^+"
OH_ION = "OH^-"


def _vec(d: dict) -> dict:
    out: dict = {}
    for s, v in d.items():
        f = F(v).limit_denominator(10000)
        if f:
            out[s] = out.get(s, F(0)) + f
    return {s: v for s, v in out.items() if v}


def _ksp_vec(solid: str, cat: str, an: str, x: int, y: int) -> dict:
    """Ksp 溶解行的净生成向量，**按 H/O 配平补上水**。

    表里的 pKsp 对应的是配平式 `solid + w·H2O → x·cat + y·an`：
    M(OH)_y 型 w=0，而 M₂O 型（Ag₂O、Cu₂O）必须补 w=1——漏掉水会让
    这一行自身不守恒，任何组合校验都失去意义（首版审计即栽在此）。
    """
    from chemkit.core import elements_of
    es, ec, ea = elements_of(solid), elements_of(cat), elements_of(an)
    # w 由 H 平衡定（O 平衡对氧化物/氢氧化物恒同解）
    w = (ec.get("H", 0) * x + ea.get("H", 0) * y - es.get("H", 0)) / 2.0
    v = {solid: -1.0, cat: float(x), an: float(y)}
    if w:
        v[WATER] = v.get(WATER, 0.0) - w
    return _vec(v)


def base_rows(T) -> list[tuple[str, dict, float, float | None]]:
    """基候选：(名字, 反应向量, logK(298), dH|None)。

    向量约定统一为 **净生成量** `pr − r`（与派生候选同口径）——两侧同号会把
    所有组合都判成不可解。dH 是**该行自身**的焓（水行恒 None：其温度项走
    pKw 通道，不进 dH 通道，见 candidates._hess_dH 的注释）。
    """
    rows: list[tuple[str, dict, float, float | None]] = []
    # 水自电离：H2O -> H+ + OH-，logK = -pKw
    rows.append(("water", _vec({WATER: -1, H_ION: 1, OH_ION: 1}),
                 -PKW_298, None))
    for e in T.pka:
        # acid -> base + n·H+，logK = -pKa
        rows.append((f"pka:{e['acid']}->{e['base']}",
                     _vec({e["acid"]: -1, e["base"]: 1,
                           H_ION: e.get("n", 1)}),
                     -e["pka"], e.get("dH")))
    for b in T.beta:
        rows.append((f"beta:{b['complex']}",
                     _vec({b["center"]: -1, b["ligand"]: -b["nu"],
                           b["complex"]: 1}),
                     b["logb"], b.get("dH")))
    for e in T.ksp:
        cat, an = e["pair"]
        x, y = _ksp_xy(e)
        rows.append((f"ksp:{e['solid']}",
                     _ksp_vec(e["solid"], cat, an, x, y),
                     -e["pKsp"], e.get("dH")))
    return rows


def _vec_of_cand(c) -> dict:
    """候选的净生成向量 pr − r。"""
    dv = _vec(dict(c.r))
    for s, v in _vec(dict(c.pr)).items():
        dv[s] = dv.get(s, F(0)) - v
    return {s: -v for s, v in dv.items() if v}


def template_rows(T, T_K: float) -> list[tuple[str, dict, float]]:
    """静态模板候选（§4.2–4.4：质子/沉淀/溶解/配位/解离）作为被检行。

    氧化物酸溶（`T.ex` 的 M_xO_y + 2yH⁺）**排除**：它的 logK 是
    `x·(n·14 − pKsp − 1.5)` 的经验估计（含 −1.5 手调项），本就不是基行的
    线性组合，纳入只会制造已知的假阳性。其余模板全部应逐位 Hess 自洽。
    """
    out: list[tuple[str, dict, float]] = []
    for cand, _req in _build_static_cands(T_K, T, pKw_of(T_K)):
        if cand.kind == "dissolve" and cand.meta.get("solid") in T.ex:
            continue
        tag = cand.meta.get("solid") or cand.meta.get("dir") or ""
        out.append((f"{cand.kind}{'/' + str(tag) if tag else ''}",
                    _vec_of_cand(cand), logK_T(cand, T_K)))
    return out


def rref_solve(A: list[list[F]], b: list[F]):
    """精确有理 RREF 解 A·a = b。返回 (解, 是否相容)。自由变量取 0。"""
    m = len(A)
    n = len(A[0]) if m else 0
    M = [A[i][:] + [b[i]] for i in range(m)]
    piv = []
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, m):
            if M[i][c] != 0:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = M[r][c]
        M[r] = [v / inv for v in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [M[i][j] - f * M[r][j] for j in range(n + 1)]
        piv.append(c)
        r += 1
        if r == m:
            break
    for i in range(r, m):
        if all(M[i][j] == 0 for j in range(n)) and M[i][n] != 0:
            return None, False
    a = [F(0)] * n
    for i, c in enumerate(piv):
        a[c] = M[i][n]
    return a, True


def rref_solve_f(A: list[list[float]], b: list[float], tol: float = 1e-9):
    """浮点 RREF 解 A·a = b（部分主元）。返回 (解, 是否相容)。

    精确有理版单行 ~150 ms（max_sub≈130 时），1.4 万行要跑半小时——守卫跑不动
    就等于没有守卫。判据本身是"是否张成 + 预测 logK 差 ≤ 0.05 个对数单位"，
    浮点足够；病态情形只会被判成"不可解"（计入 unsolved，不产生假阳性）。
    需要逐位精确复核时用 `--exact`。
    """
    m = len(A)
    n = len(A[0]) if m else 0
    M = [list(A[i]) + [b[i]] for i in range(m)]
    piv = []
    r = 0
    for c in range(n):
        pr, best = None, tol
        for i in range(r, m):
            if abs(M[i][c]) > best:
                pr, best = i, abs(M[i][c])
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = 1.0 / M[r][c]
        M[r] = [v * inv for v in M[r]]
        for i in range(m):
            f = M[i][c]
            if i != r and f != 0.0:
                M[i] = [M[i][j] - f * M[r][j] for j in range(n + 1)]
        piv.append(c)
        r += 1
        if r == m:
            break
    for i in range(r, m):
        if all(abs(M[i][j]) <= tol for j in range(n)) and abs(M[i][n]) > tol:
            return None, False
    a = [0.0] * n
    for i, c in enumerate(piv):
        a[c] = M[i][n]
    return a, True


def main(argv: list[str]) -> None:
    tol = 0.05
    top = 25
    src_filter = None
    T_K = 298.15
    show = "--show" in argv
    exact = "--exact" in argv
    only_templates = "--templates" in argv
    with_templates = only_templates or "--with-templates" in argv
    if "--tol" in argv:
        tol = float(argv[argv.index("--tol") + 1])
    if "--top" in argv:
        top = int(argv[argv.index("--top") + 1])
    if "--src" in argv:
        src_filter = argv[argv.index("--src") + 1]
    if "--T" in argv:
        T_K = float(argv[argv.index("--T") + 1])

    T = load_tables()
    bases = base_rows(T)
    derived = build_derived(T)
    tmpls = template_rows(T, T_K)
    print(f"T = {T_K} K（pKw = {pKw_of(T_K):.4f}）；基候选 {len(bases)} 条；"
          f"派生候选 {len(derived)} 条；静态模板 {len(tmpls)} 条")

    pool: list[tuple[str, dict, float]] = []
    if not only_templates:
        pool += [(d.meta.get("src", ""), _vec_of_cand(d), logK_T(d, T_K))
                 for d in derived]
    if with_templates:
        pool += tmpls

    pk_t = pKw_of(T_K)
    bad: list[tuple[float, str, float, float, list]] = []
    unsolved = 0
    checked = 0
    for src, dv, lgk in pool:
        if src_filter and src_filter not in src:
            continue
        sp = set(dv)
        # 只用与 d 物种有交集的基行（其余系数必为 0）
        sub = [i for i, r in enumerate(bases) if sp & set(r[1])]
        if not sub:
            continue
        species = sorted(set().union(*[set(bases[i][1]) for i in sub]))
        if exact:
            A = [[bases[i][1].get(s, F(0)) for i in sub] for s in species]
            b = [dv.get(s, F(0)) for s in species]
            a, ok = rref_solve(A, b)
        else:
            A = [[float(bases[i][1].get(s, 0)) for i in sub] for s in species]
            b = [float(dv.get(s, 0)) for s in species]
            a, ok = rref_solve_f(A, b)
        if not ok or a is None:
            unsolved += 1
            continue
        resid = 0.0
        for k, s in enumerate(species):
            lhs = sum(float(a[j]) * float(bases[sub[j]][1].get(s, 0))
                      for j in range(len(sub)))
            resid = max(resid, abs(lhs - float(b[k])))
        if resid > 1e-6:
            unsolved += 1
            continue
        checked += 1
        terms = [(bases[sub[j]][0], float(a[j]))
                 for j in range(len(sub)) if abs(float(a[j])) > 1e-9]
        # 预测 logK(T)：298 基值线性组合 + pKw 通道（水行）+ dH 通道
        pred = sum(float(a[j]) * bases[sub[j]][2] for j in range(len(sub)))
        a_w = sum(float(a[j]) for j in range(len(sub))
                  if bases[sub[j]][0] == "water")
        pred += a_w * (PKW_298 - pk_t)
        dhs = [None if bases[sub[j]][3] is None
               else float(a[j]) * bases[sub[j]][3]
               for j in range(len(sub))
               if abs(float(a[j])) > 1e-9 and bases[sub[j]][0] != "water"]
        if all(d is not None for d in dhs):
            pred += _vant(sum(dhs), T_K)
        if show:
            print(f"  {src}\n     logK={lgk:+9.3f} 组合={pred:+9.3f}  "
                  + "  ".join(f"{v:+g}×[{n}]" for n, v in terms))
        if abs(pred - lgk) > tol:
            bad.append((abs(pred - lgk), src, lgk, pred, terms))

    bad.sort(key=lambda t: -t[0])
    # 已知的有意偏离：H_2SiO_3 族用**手调表观 pKa**（templates._solid_acid_pka：
    # pKa_eff = pKsp − 下游 pKa 之和 = 24.5 − 12.0 = 12.5），把"固相酸脱质子"
    # 与溶解度耦合起来，本就不是基行的线性组合。只有这 3 条，单独标注，
    # 免得后来者反复追它（其余任何偏离都必须给出化学理由）。
    known = [t for t in bad if any("SiO_3" in n for n, _ in t[4])]
    bad = [t for t in bad if t not in known]
    fam: dict[str, list[float]] = {}
    for delta, src, *_ in bad:
        key = src.split(":")[0] or src
        fam.setdefault(key, []).append(delta)
    print(f"可解并校验 {checked} 条；不可解（基集未张成）{unsolved} 条；"
          f"**logK 不自洽 {len(bad)} 条**（阈 {tol}）"
          f"；已知手调族 {len(known)} 条（H_2SiO_3 表观 pKa，不计入）")
    if fam:
        print("  按族：")
        for k in sorted(fam, key=lambda k: -max(fam[k])):
            print(f"    {k:<22} {len(fam[k]):>5} 条  maxΔ={max(fam[k]):.3f}")
    for delta, src, got, pred, terms in bad[:top]:
        print(f"\n  Δ{delta:8.3f}  {src}")
        print(f"     引擎 logK={got:+9.3f}   Hess 组合={pred:+9.3f}")
        print("     " + "  ".join(f"{a:+g}×[{n}]" for n, a in terms))


if __name__ == "__main__":
    main(sys.argv[1:])
