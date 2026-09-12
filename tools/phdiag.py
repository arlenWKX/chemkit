"""pH 机器（启发式分支机）vs 电荷平衡精确解（acidbase.charge_pH）全量对照。

目的（0.5.x 施工顺序 ③ 的量化前置）：在**引擎自己的退出账本**上比较两条
pH 求法，回答两个问题：

  1. 分歧有多大、集中在哪些用例（决定替换风险面）；
  2. 分歧例上谁是化学正确的一方（决定替换方向）。

对照口径：probe 记录的 `ledger` 已经过 `_respeciate_strong_acids` +
`_presentation_He`，`probe["pH"]` 即机器在**同一账本**上的输出——两者
输入完全相同，唯一差别是算法，故差值即纯算法差。

用法：
    python tools/phdiag.py [-n 20] [--cases 路径] [--only 前缀,前缀]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter

sys.path.insert(0, ".")

from chemkit.acidbase import charge_pH, build_families          # noqa: E402
from chemkit.core import charge_of                              # noqa: E402
from chemkit.data import load_tables                            # noqa: E402
from chemkit.engine import judge                                # noqa: E402
from chemkit.testsuit import load_cases                         # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=20)
    ap.add_argument("--cases", default=None)
    ap.add_argument("--only", default=None, help="逗号分隔的用例名前缀")
    ap.add_argument("--tol", type=float, default=0.05)
    args = ap.parse_args()

    T = load_tables()
    fams = build_families(T)
    cases = load_cases(args.cases)
    if args.only:
        pre = tuple(x.strip() for x in args.only.split(",") if x.strip())
        cases = [c for c in cases if c["name"].startswith(pre)]

    rows = []
    none_n = 0
    for c in cases:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        probe: dict = {}
        judge(subs, cond, T, _probe=probe)
        if "ledger" not in probe:
            continue                        # OVERRIDE 直出，无账本画像
        V = float(cond.get("V_L", 1.0) or 1.0)
        T_K = float(cond.get("T_K", 298.15) or 298.15)
        led = probe["ledger"]
        He = probe["H_excess"]
        p_m = probe["pH"]
        c_H = float(cond.get("c_H", 0.0) or 0.0)
        # 正确口径：账本自身电荷平衡定 pH（H_excess 是待求量的另一半，
        # 传进去即双重记账）；c_H 是无阴离子的强酸条件量，须补回。
        p_c = charge_pH(led, V, T, T_K, c_H=c_H)
        if p_c is None:
            none_n += 1
        net = sum(charge_of(s) * m for s, m in led.items()
                  if not s.startswith("__") and s not in ("H^+", "OH^-"))
        fam_hit = sum(1 for s in led if s in fams)
        rows.append({
            "name": c["name"], "V": V, "T_K": T_K, "He": He, "c_H": c_H,
            "machine": p_m, "charge": p_c,
            "d": None if p_c is None else round(p_c - p_m, 3),
            "families": fam_hit, "nsp": len(led),
            "cons": abs(net + He + c_H),
            "note": (c.get("note") or "")[:70],
        })

    got = [r for r in rows if r["d"] is not None]
    big = sorted(got, key=lambda r: -abs(r["d"]))
    print(f"用例 {len(rows)}；charge_pH 无括号 {none_n}；可对照 {len(got)}")
    bad = [r for r in rows if r["cons"] > 1e-6]
    print(f"账本守恒违规 |Σz·n+He|>1e-6 : {len(bad)} 例"
          f"（最大 {max((r['cons'] for r in rows), default=0):.3g}）")
    hist = Counter()
    for r in got:
        a = abs(r["d"])
        k = ("<0.01" if a < 0.01 else "<0.05" if a < 0.05 else "<0.2" if a < 0.2
             else "<0.5" if a < 0.5 else "<1" if a < 1 else ">=1")
        hist[k] += 1
    for k in ("<0.01", "<0.05", "<0.2", "<0.5", "<1", ">=1"):
        print(f"  |ΔpH| {k:6s} : {hist[k]:5d}")
    print(f"  分歧 >{args.tol} : "
          f"{sum(1 for r in got if abs(r['d']) > args.tol)}")
    print(f"\n== 分歧最大 {args.n} 例（Δ = charge − machine）==")
    for r in big[:args.n]:
        print(f"  Δ{r['d']:+8.3f}  机器 {r['machine']:6.2f} → 电荷 {r['charge']:6.2f} "
              f" He={r['He']:+10.6f} 族{r['families']}/物种{r['nsp']:3d}  {r['name'][:52]}")
    print(f"\n== 一致例（|Δ|<0.01）样例 ==")
    for r in [x for x in got if abs(x["d"]) < 0.01][:args.n]:
        print(f"  pH {r['machine']:6.2f}  He={r['He']:+10.6f}  {r['name'][:56]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
