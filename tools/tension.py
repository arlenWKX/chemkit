"""张力例归因：`resid_live` 高的用例，究竟是**哪一条平衡**在驱动、为什么没走。

**为什么需要**（architecture §7 T-4）：`resid_live`（converg 的质量口径）只说
"退出账本上仍有值得解、walk 也会解的平衡"，不说**为什么**它没被解。可能的
原因完全不同，修法也完全不同：

  · `disabled`（单向或双向）——微步/振荡禁用把它锁死（J06 型真实病灶）；
  · `frozen` / `slow` / `blocked`——引擎既定语义（不计入 live，但要看得到）；
  · `ext_max < ANN_MIN_EXTENT`——痕量方向，无关的驱动（D38 型）；
  · 候选根本没被枚举出来（模板/在场/膜闸门）——枚举缺口。

本工具按 |S| 排序打印**退出账本**上的活跃平衡及其全部旗标，再对最强者直接
调 `solve_extent` 看"走步眼里的程度"（复用 tools/extent.py 的判据：x* 是否
被微步阈值挡下）。

用法：
    python tools/tension.py                 # 全量，取前 20
    python tools/tension.py -n 40
    python tools/tension.py --case RX12     # 指定用例前缀
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from chemkit.candidates import ANN_MIN_EXTENT, WATER, X_MIN   # noqa: E402
from chemkit.data import load_tables                          # noqa: E402
from chemkit.engine import (judge, solve_extent, S_of,        # noqa: E402
                            _respeciate_strong_acids, _presentation_He,
                            _fmt)
from chemkit.speciation import estimate_pH                    # noqa: E402
from chemkit.templates import enumerate_candidates            # noqa: E402
from chemkit.testsuit import load_cases                       # noqa: E402


def _why(a: dict) -> str:
    if a["frozen"]:
        return "frozen（引擎宣告平衡止震）"
    if a.get("slow"):
        return "slow（动力学标注）"
    if a.get("blocked"):
        return "blocked（膜封锁溶剂通道）"
    if a["ext_max"] < ANN_MIN_EXTENT:
        return f"痕量方向（ext_max {a['ext_max']:.2g} < {ANN_MIN_EXTENT:g}）"
    if a.get("dis_fwd") and a.get("dis_rev"):
        return "**双向禁用**"
    if a.get("dis_fwd") or a.get("dis_rev"):
        return f"**单向禁用**（fwd={a.get('dis_fwd')} rev={a.get('dis_rev')}）"
    return "**无旗标却没走**（枚举/仲裁缺口）"


def main() -> None:
    n = 20
    args = sys.argv[1:]
    if "-n" in args:
        n = int(args[args.index("-n") + 1])
    want: tuple[str, ...] = ()
    if "--case" in args:
        want = tuple(args[args.index("--case") + 1:])

    T = load_tables()
    cases = load_cases(None)
    if want:
        cases = [c for c in cases if c["name"].startswith(want)]
    rows = []
    for c in cases:
        subs = [{"name": nm, "mol": m} for nm, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        probe: dict = {}
        judge(subs, cond, T, _probe=probe)
        if not probe or not probe.get("active"):
            continue
        act = [a for a in probe["active"] if a["two_sided"]
               and not a["frozen"] and not a.get("slow") and not a.get("blocked")
               and a.get("ext_max", 9e9) >= ANN_MIN_EXTENT]
        if not act:
            continue
        top = max(act, key=lambda a: abs(a["S"]))
        rows.append((abs(top["S"]), c["name"], probe, top, len(act)))
    rows.sort(key=lambda r: -r[0])
    print(f"张力例（live > 0，按 |S| 排序）：共 {len(rows)} 例，列前 {n}\n")
    for val, name, probe, top, nact in rows[:n]:
        print(f"== {name[:46]}  live={val:.3f}  pH={probe['pH']} "
              f"exit={probe['exit']} iters={probe['iters']} 活跃={nact}")
        print(f"   |S|={val:.3f} kind={top['kind']:9s} ext_max={top['ext_max']:.3g}"
              f"  {top['eq'][:96]}")
        print(f"   → {_why(top)}")
        # 走步眼里的程度（与 tools/extent.py 同判据）
        if want:
            led = dict(probe["ledger"])
            V = float(cond.get("V_L", 1.0) or 1.0)
            T_K = float(cond.get("T_K", 298.15) or 298.15)
            He = _respeciate_strong_acids(led, probe["H_excess"], V, T)
            He = _presentation_He(led, He, V, T, T_K)
            pH = estimate_pH(led, He, V, T, T_K)
            for cd in enumerate_candidates(led, He, pH, V, T_K, T, True):
                if _fmt_eq(cd) != top["eq"]:
                    continue
                pres = cd.pres_specs
                if not all(led.get(s, 0.0) > X_MIN for s in pres[0] + pres[1]):
                    continue
                S = S_of(cd, led, V, pH, T_K, T, (), 101.325, True)
                d = 1 if S > 0 else -1
                x, xm = solve_extent(cd, d, led, He, V, T_K, T, frozenset(),
                                     iters=60, p_ext_kpa=101.325)
                thr = max(X_MIN, 1e-4 * xm)
                print(f"     走步判据：x*={x:.4g} x_max={xm:.4g} 微步阈={thr:.3g}"
                      f" → {'微步挡下' if x <= thr else '会执行'}")
                break


def _fmt_eq(c) -> str:
    r = " + ".join(f"{_fmt(nu)}{s}" for s, nu in c.r.items() if s != WATER)
    p = " + ".join(f"{_fmt(nu)}{s}" for s, nu in c.pr.items() if s != WATER)
    return f"{r} -> {p}"


if __name__ == "__main__":
    main()
