"""对退出账本上"仍有驱动"的候选直接调 solve_extent，看走步眼里的程度。

用途：区分"平衡已到位"与"微步阈值把它挡下"。
用法：python tools/extent.py N23 [N23 ...]
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from chemkit.candidates import WATER, X_MIN            # noqa: E402
from chemkit.data import load_tables                    # noqa: E402
from chemkit.engine import (judge, solve_extent, S_of,  # noqa: E402
                            _respeciate_strong_acids, _presentation_He)
from chemkit.speciation import estimate_pH              # noqa: E402
from chemkit.templates import enumerate_candidates      # noqa: E402
from chemkit.testsuit import load_cases                 # noqa: E402


def main(prefixes: list[str]) -> None:
    T = load_tables()
    for case in load_cases(None):
        if not case["name"].startswith(tuple(prefixes)):
            continue
        subs = [{"name": n, "mol": m} for n, m in case["subs"]]
        cond = case.get("cond") or {"V_L": 1.0}
        V = float(cond.get("V_L", 1.0) or 1.0)
        T_K = float(cond.get("T_K", 298.15) or 298.15)
        probe: dict = {}
        judge(subs, cond, T, _probe=probe)
        led = dict(probe["ledger"])
        He = _respeciate_strong_acids(led, probe["H_excess"], V, T)
        He = _presentation_He(led, He, V, T, T_K)
        pH = estimate_pH(led, He, V, T, T_K)
        print(f"\n== {case['name'][:44]}  pH={pH:.3f} He={He:.3g}")
        for c in enumerate_candidates(led, He, pH, V, T_K, T, True):
            pr = all(led.get(s, 0.0) > X_MIN for s in c.pres_specs[0])
            pp = all(led.get(s, 0.0) > X_MIN for s in c.pres_specs[1])
            if not (pr and pp):
                continue
            S = S_of(c, led, V, pH, T_K, T, (), 101.325, True)
            if abs(S) < 0.5:
                continue
            d = 1 if S > 0 else -1
            x, xm = solve_extent(c, d, led, He, V, T_K, T, frozenset(),
                                 iters=60, p_ext_kpa=101.325,
                                 gas_escape=True)
            solid = [s for s in list(c.r) + list(c.pr) if s in T.solids]
            thr = max(X_MIN, 1e-4 * xm)
            print(f"   S={S:+8.3f} kind={c.kind:8s} solid={solid} "
                  f"x*={x:.4g} x_max={xm:.4g} 阈值={thr:.4g} "
                  f"{'→ 微步挡下' if x <= thr else '→ 会执行'}")


if __name__ == "__main__":
    main(sys.argv[1:] or ["N23"])
