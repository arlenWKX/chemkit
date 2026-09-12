"""单例账本/候选深探：把退出账本、机器 pH、电荷平衡 pH、候选 S 一次打全。

用法：python tools/case.py H46 [J06 ...]
"""
from __future__ import annotations

import json
import sys

sys.path.insert(0, ".")

from chemkit.acidbase import charge_pH, build_families, eff_pka   # noqa: E402
from chemkit.core import charge_of, pKw_of                        # noqa: E402
from chemkit.data import load_tables                              # noqa: E402
from chemkit.engine import judge                                  # noqa: E402
from chemkit.testsuit import load_cases                           # noqa: E402


def show(prefixes: list[str], verbose: bool = False) -> None:
    T = load_tables()
    fams = build_families(T)
    for c in load_cases(None):
        if not c["name"].startswith(tuple(prefixes)):
            continue
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        probe: dict = {}
        r = judge(subs, cond, T, _probe=probe)
        V = float(cond.get("V_L", 1.0) or 1.0)
        T_K = float(cond.get("T_K", 298.15) or 298.15)
        print(f"\n{'='*78}\n{c['name']}   cond={cond}")
        print(f"  note: {(c.get('note') or '')[:200]}")
        print(f"  degree={r['degree']} changed={r['changed']} "
              f"final_pH={r.get('final_pH')} steps={len(r.get('steps', []))}")
        print(f"  net_equation: {r.get('net_equation')}")
        print(f"  exit={probe.get('exit')} iters={probe.get('iters')} "
              f"resid={probe.get('max_abs_S')} frozen={probe.get('frozen_n')} "
              f"disabled={probe.get('disabled_n')}")
        if "ledger" not in probe:
            continue
        led = probe["ledger"]
        He = probe["H_excess"]
        p_m = probe["pH"]
        p_c = charge_pH(led, V, T, T_K, c_H=float(cond.get("c_H", 0.0) or 0.0))
        p_c2 = charge_pH(led, V, T, T_K)
        print(f"  机器 pH={p_m}  charge_pH(He)={p_c}  charge_pH(no He)={p_c2}"
              f"  pKw={pKw_of(T_K):.3f}")
        net = sum(charge_of(s) * m for s, m in led.items())
        print(f"  He={He!r}  Σz·n={net:.6f}  Σz·n+He={net + He:.2e}")
        print(f"  --- 账本 ({len(led)}) ---")
        for s, m in sorted(led.items(), key=lambda kv: -kv[1]):
            fi = fams.get(s)
            tag = ""
            if fi is not None:
                fid, order, idx, pka, dh = fi[:5]
                ek = eff_pka(pka, dh, T_K)
                tag = (f" [族 {fid}: {list(order)} pKa="
                       f"{[round(x, 2) for x in ek]} q="
                       f"{[charge_of(x) for x in order]}"
                       f" sgn={list(fi[5])}]")
            print(f"    {s:24s} {m:12.6g}  z={charge_of(s):+d}{tag}")
        if verbose:
            print("  --- 候选 S（两侧在场）---")
            for a in probe["active"]:
                if a["two_sided"]:
                    print(f"    S={a['S']:+8.3f} froz={int(a['frozen'])} "
                          f"dis={int(a['disabled'])} {a['eq'][:90]}")
        print("  has:", json.dumps(c.get("has") or c.get("has_range") or {},
                                  ensure_ascii=False))
        print("  has_not:", json.dumps(c.get("has_not") or {}, ensure_ascii=False))
        if c.get("ph"):
            lo, hi = c["ph"]
            ok_m = lo <= r.get("final_pH", -99) <= hi
            ok_c = p_c is not None and lo <= p_c <= hi
            print(f"  ph 期望 [{lo}, {hi}]：机器 {r.get('final_pH')} "
                  f"{'✓' if ok_m else '✗'}  电荷 {None if p_c is None else round(p_c, 3)} "
                  f"{'✓' if ok_c else '✗'}")


if __name__ == "__main__":
    a = sys.argv[1:]
    v = "--v" in a
    a = [x for x in a if not x.startswith("-")]
    show(a or ["H46"], v)
