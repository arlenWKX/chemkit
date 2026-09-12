"""Hess 一致性审计：派生候选的 logK 是否与其组元的线性组合相符。

原理（这正是"派生行秩亏"的判据）：派生候选的净变换向量必然落在**基候选**
（beta/pKa/Ksp/水）张成的反应子空间里——它本来就是这样造出来的。于是对每条
派生候选 d，解线性方程组

    Σ_i a_i · vec(base_i) = vec(d)        （vec = 物种计量的精确有理向量）

若可解（基集确实张成 d），则必须同时满足

    logK(d) == Σ_i a_i · logK(base_i)

**行秩亏的含义**：派生行是基行的线性组合，所以任何"秩检查"都看不出问题；
但组合系数若错（例如少乘一个 pKw），logK 就与基行不自洽——方程组变得
矛盾，走步/联立在这个方向上永远无法收敛。H46 即此：beta_ksp 通道把
`n_OH` 取成中心离子电荷而非固相 OH 数，Ag2O（x=2,y=2）少减一个 pKw，
logK 偏 +14。

用法：
    python tools/hess_audit.py [--tol 0.05] [--top 30] [--src beta_ksp]
"""
from __future__ import annotations

import sys
from fractions import Fraction as F

sys.path.insert(0, ".")

from chemkit.candidates import (WATER, build_derived, _ksp_xy,  # noqa: E402
                                PKW_298)
from chemkit.core import charge_of                              # noqa: E402
from chemkit.data import load_tables                            # noqa: E402

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


def base_rows(T) -> list[tuple[str, dict, float]]:
    """基候选：(名字, 反应向量, logK)。

    向量约定统一为 **净生成量** `pr − r`（与派生候选的 dv 同口径）——
    两侧同号会把所有组合都判成不可解。
    """
    rows: list[tuple[str, dict, float]] = []
    # 水自电离：H2O -> H+ + OH-，logK = -pKw
    rows.append(("water", _vec({WATER: -1, H_ION: 1, OH_ION: 1}), -PKW_298))
    for e in T.pka:
        # acid -> base + n·H+，logK = -pKa
        rows.append((f"pka:{e['acid']}->{e['base']}",
                     _vec({e["acid"]: -1, e["base"]: 1,
                           H_ION: e.get("n", 1)}), -e["pka"]))
    for b in T.beta:
        rows.append((f"beta:{b['complex']}",
                     _vec({b["center"]: -1, b["ligand"]: -b["nu"],
                           b["complex"]: 1}), b["logb"]))
    for e in T.ksp:
        cat, an = e["pair"]
        x, y = _ksp_xy(e)
        rows.append((f"ksp:{e['solid']}",
                     _ksp_vec(e["solid"], cat, an, x, y), -e["pKsp"]))
    return rows


def rref_solve(A: list[list[F]], b: list[F]):
    """RREF 解 A·a = b（精确有理）。返回 (解, 是否相容)。自由变量取 0。"""
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


def main(argv: list[str]) -> None:
    tol = 0.05
    top = 30
    src_filter = None
    if "--tol" in argv:
        tol = float(argv[argv.index("--tol") + 1])
    if "--top" in argv:
        top = int(argv[argv.index("--top") + 1])
    if "--src" in argv:
        src_filter = argv[argv.index("--src") + 1]

    T = load_tables()
    bases = base_rows(T)
    print(f"基候选 {len(bases)} 条；派生候选 {len(build_derived(T))} 条")

    if src_filter:
        for e in T.ksp:
            if "Ag_2O" in e["solid"]:
                print("  Ag2O ksp 条目:", e, "_ksp_xy=", _ksp_xy(e))

    bad: list[tuple[float, str, float, float, list]] = []
    unsolved = 0
    checked = 0
    for d in build_derived(T):
        src = d.meta.get("src", "")
        if src_filter and src_filter not in src:
            continue
        dv = _vec(dict(d.r))
        for s, v in _vec(dict(d.pr)).items():
            dv[s] = dv.get(s, F(0)) - v
        dv = {s: -v for s, v in dv.items() if v}
        # 只用与 d 物种有交集的基行（其余系数必为 0）
        sp = set(dv)
        sub = [i for i, (_n, bv, _k) in enumerate(bases) if sp & set(bv)]
        if not sub:
            continue
        species = sorted(set().union(*[set(bases[i][1]) for i in sub]))
        A = [[bases[i][1].get(s, F(0)) for i in sub] for s in species]
        b = [dv.get(s, F(0)) for s in species]
        a, ok = rref_solve(A, b)
        if not ok or a is None:
            unsolved += 1
            continue
        # 回代校验（自由变量取 0 可能不张成）
        resid = 0.0
        for k, s in enumerate(species):
            lhs = sum(a[j] * bases[sub[j]][1].get(s, F(0))
                      for j in range(len(sub)))
            resid = max(resid, abs(float(lhs - b[k])))
        if resid > 1e-9:
            unsolved += 1
            if src_filter:
                print(f"  未张成 {src} resid={resid:g}")
            continue
        checked += 1
        pred = sum(float(a[j]) * bases[sub[j]][2] for j in range(len(sub)))
        if src_filter:
            terms = [(bases[sub[j]][0], float(a[j]))
                     for j in range(len(sub)) if a[j] != 0]
            print(f"  {src}\n     logK={d.logK:+9.3f} 组合={pred:+9.3f}  "
                  + "  ".join(f"{v:+g}×[{n}]" for n, v in terms))
        if abs(pred - d.logK) > tol:
            terms = [(bases[sub[j]][0], float(a[j]))
                     for j in range(len(sub)) if a[j] != 0]
            bad.append((abs(pred - d.logK), src, d.logK, pred, terms))

    bad.sort(key=lambda t: -t[0])
    print(f"可解并校验 {checked} 条；不可解（基集未张成）{unsolved} 条；"
          f"**logK 不自洽 {len(bad)} 条**（阈 {tol}）")
    for delta, src, got, pred, terms in bad[:top]:
        print(f"\n  Δ{delta:8.3f}  {src}")
        print(f"     引擎 logK={got:+9.3f}   Hess 组合={pred:+9.3f}")
        print("     " + "  ".join(f"{a:+g}×[{n}]" for n, a in terms))


if __name__ == "__main__":
    main([x for x in sys.argv[1:] if not x.startswith("-") or x.startswith("--")]
         or sys.argv[1:])
