"""候选全量深探：给一个用例，列出退出点上"两侧在场"的候选及其完整
计量/ logK / S / 冻结状态，并检查每条净变换的元素与电荷守恒。

用法：python tools/cands.py H46 [--all]
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from chemkit.candidates import WATER, logK_T                # noqa: E402
from chemkit.core import charge_of, elements_of             # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.engine import (judge, S_of, X_MIN,             # noqa: E402
                            _respeciate_strong_acids, _presentation_He)
from chemkit.speciation import estimate_pH                  # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def fmt(d: dict) -> str:
    return " + ".join(f"{v:g}{s}" for s, v in sorted(d.items()) if s != WATER)


def net_ok(c) -> str:
    """净变换的元素与电荷守恒检查（r -> pr 的净变化应为零）。"""
    bad = []
    for el in set().union(*[set(elements_of(s)) for s in list(c.r) + list(c.pr)]):
        d = sum(nu * elements_of(s).get(el, 0) for s, nu in c.pr.items()) - \
            sum(nu * elements_of(s).get(el, 0) for s, nu in c.r.items())
        if abs(d) > 1e-9:
            bad.append(f"{el}{d:+g}")
    dq = sum(nu * charge_of(s) for s, nu in c.pr.items()) - \
        sum(nu * charge_of(s) for s, nu in c.r.items())
    if abs(dq) > 1e-9:
        bad.append(f"q{dq:+g}")
    return ",".join(bad)


def main(prefix: str, show_all: bool) -> None:
    T = load_tables()
    for case in load_cases(None):
        if not case["name"].startswith(prefix):
            continue
        subs = [{"name": n, "mol": m} for n, m in case["subs"]]
        cond = case.get("cond") or {"V_L": 1.0}
        V = float(cond.get("V_L", 1.0) or 1.0)
        T_K = float(cond.get("T_K", 298.15) or 298.15)
        probe: dict = {}
        r = judge(subs, cond, T, _probe=probe)
        led = dict(probe["ledger"])
        He = _respeciate_strong_acids(led, probe["H_excess"], V, T)
        He = _presentation_He(led, He, V, T, T_K)
        pH = estimate_pH(led, He, V, T, T_K)
        print(f"== {case['name']}  pH={pH:.3f} He={He!r}  "
              f"steps={len(r.get('steps', []))} frozen={probe['frozen_n']}")
        cands = enumerate_candidates(led, He, pH, V, T_K, T, True)
        rows = []
        for c in cands:
            pr = all(led.get(s, 0.0) > X_MIN for s in c.pres_specs[0])
            pp = all(led.get(s, 0.0) > X_MIN for s in c.pres_specs[1])
            if not (pr or pp):
                continue
            S = S_of(c, led, V, pH, T_K, T, (), 101.325, True)
            rows.append((S, c, pr and pp))
        rows.sort(key=lambda t: -abs(t[0]))
        for S, c, two in rows:
            if not show_all and not two:
                continue
            bad = net_ok(c)
            k = logK_T(c, T_K)
            print(f"  S={S:+9.3f} logK={k:+8.3f} {c.kind:7s} "
                  f"{'two' if two else 'one'} {'BAD:' + bad if bad else ''}")
            print(f"      {fmt(c.r)}  ->  {fmt(c.pr)}")
        break


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("-")]
    main(a[0] if a else "H46", "--all" in sys.argv)
