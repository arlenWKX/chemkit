# -*- coding: utf-8 -*-
"""第 263 轮 · N15 回归溯源：新路径（两性支路精确解）在哪些状态上触发？

记录每一次 `_exact_ok(..., min_fams=1)` 的判定与随后的精确解，并标出
**该账本里有没有固相**（以及最大固相量）——用于回答"N15 的极限环是不是
被新路径带进来的、且是否越过了 `_exact_ok` 的固相闸"。

用法：python tools/n15_trace.py N15
"""
from __future__ import annotations

import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit import speciation as _sp            # noqa: E402
from chemkit.candidates import WATER, X_MIN      # noqa: E402
from chemkit.data import load_tables             # noqa: E402

_ORIG_OK = _sp._exact_ok
_ORIG_CHG = _sp.charge_pH
LOG: list = []


def main(argv: list[str]) -> int:
    name = argv[0] if argv else "N15"
    T = load_tables()
    solids = frozenset(T.solids)

    def ok(ledger, He_res, Tt, V, min_fams=2):
        r = _ORIG_OK(ledger, He_res, Tt, V, min_fams)
        if min_fams == 1:
            smax = 0.0
            snames = []
            for sp, m in ledger.items():
                if sp in solids and m > 0.0:
                    snames.append((sp, m))
                    smax = max(smax, m)
            LOG.append({"min_fams": 1, "ok": r, "solid_max": smax,
                        "solids": sorted(snames, key=lambda t: -t[1])[:3],
                        "He_res": He_res,
                        "n": len([1 for s, m in ledger.items()
                                  if m > X_MIN and s != WATER])})
        return r

    def chg(ledger, V, Tt, T_K, **kw):
        r = _ORIG_CHG(ledger, V, Tt, T_K, **kw)
        if kw.get("fast"):
            LOG.append({"chg": True, "pH": r})
        return r

    _sp._exact_ok = ok
    _sp.charge_pH = chg
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = {c["name"]: c for c in json.load(fh)}
    hit = [c for k, c in cases.items() if k.startswith(name)][0]
    from chemkit.engine import judge
    pr = {}
    subs = [{"name": s[0], "mol": float(s[1])} for s in hit["subs"]]
    r = judge(subs, hit.get("cond") or {}, T, _probe=pr)
    _sp._exact_ok = _ORIG_OK
    _sp.charge_pH = _ORIG_CHG

    calls = [e for e in LOG if e.get("min_fams") == 1]
    fires = [e for e in calls if e["ok"]]
    withsolid = [e for e in fires if e["solid_max"] > 0.0]
    print(f"=== {hit['name']} ===")
    print(f"pH={r.get('final_pH')} iters={pr.get('iters')} exit={pr.get('exit')}")
    print(f"两性支路被问次数 {len(calls)}；判为可信 {len(fires)}；"
          f"其中**账本含固相**的 {len(withsolid)}；X_MIN={X_MIN:g}")
    mx = max((e["solid_max"] for e in fires), default=0.0)
    print(f"可信调用里最大固相量 {mx:.6g}")
    print("\n前 12 次可信调用：")
    for e in fires[:12]:
        print(f"   solid_max={e['solid_max']:<12.4g} He_res={e['He_res']:+.3e} "
              f"n_sp={e['n']:<4} {e['solids']}")
    print("\n最后 12 次可信调用：")
    for e in fires[-12:]:
        print(f"   solid_max={e['solid_max']:<12.4g} He_res={e['He_res']:+.3e} "
              f"n_sp={e['n']:<4} {e['solids']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
