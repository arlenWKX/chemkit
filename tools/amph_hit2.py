# -*- coding: utf-8 -*-
"""第 263 轮 · 抓"两性支路精确解"接管时**真正喂给 `charge_pH` 的账本**。

`estimate_state` 内部先跑 `_buffer_titration` 得到虚拟账本，分支 4 用的是
**那个**账本（不是调用方传进来的）。故包 `charge_pH` 才能看到真身。

用法：python tools/amph_hit2.py N15
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

LOG: list = []
_ORIG = _sp.charge_pH


def main(argv: list[str]) -> int:
    name = argv[0] if argv else "N15"
    T = load_tables()
    fams = _sp.build_families(T)
    amph = set(T.pka_acid) & set(T.pka_base)
    solids = frozenset(T.solids)

    def wrap(ledger, V, Tt, T_K, **kw):
        r = _ORIG(ledger, V, Tt, T_K, **kw)
        if kw.get("fast") and r is not None:
            present_amph = [s for s, m in ledger.items()
                            if m > X_MIN and s in amph]
            if present_amph:
                LOG.append({
                    "pH": r, "He_arg": kw.get("c_H"),
                    "amph_present": present_amph,
                    "solids": {s: m for s, m in ledger.items()
                               if s in solids and m > 0.0},
                    "n_fam": sum(1 for s, m in ledger.items()
                                 if m > 0.0 and s in fams),
                    "ledger": {s: round(m, 8) for s, m in
                               sorted(ledger.items())
                               if m > 0.0 and s != WATER},
                })
        return r

    _sp.charge_pH = wrap
    try:
        with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                     encoding="utf-8") as fh:
            cases = {c["name"]: c for c in json.load(fh)}
        hit = [c for k, c in cases.items() if k.startswith(name)][0]
        from chemkit.engine import judge
        pr = {}
        subs = [{"name": s[0], "mol": float(s[1])} for s in hit["subs"]]
        res = judge(subs, hit.get("cond") or {}, T, _probe=pr)
    finally:
        _sp.charge_pH = _ORIG
    print(f"=== {hit['name']} ===  pH={res.get('final_pH')} "
          f"iters={pr.get('iters')} exit={pr.get('exit')}")
    print(f"含两性物种的 charge_pH 调用：{len(LOG)}")
    for i, e in enumerate(LOG[:10], 1):
        led = "  ".join(f"{k}={v:g}" for k, v in e["ledger"].items())
        print(f"\n[{i}] charge_pH -> {e['pH']:.6f}   n_fam={e['n_fam']}  "
              f"amph={e['amph_present']}  solids={e['solids']}")
        print(f"    {led}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
