# -*- coding: utf-8 -*-
"""第 263 轮 · N15 极限环是不是**预先存在**的走步脆弱性？

做法：把 N15 的 CO₂ 投料 ±1% 微扰，在**新路径关闭**（= 改前行为）下扫一遍。
若某个微扰点也落进同一个 `2H^+ + SiO_3^{2-} -> H_2SiO_3` /
`HSiO_3^- -> SiO_3^{2-} + H^+` 等量交替环（exit=max-iter、iters≈3000），
则"极限环"是**走步层既有的**，第 263 轮的精确解只是把轨迹推到了它的吸引域里
——而不是新路径本身算错了 pH。

用法：python tools/n15_cycle.py
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
from chemkit.data import load_tables             # noqa: E402

_ORIG = _sp._exact_ok


def _off(ledger, He_res, T, V, min_fams=2):
    if min_fams == 1:
        return False
    return _ORIG(ledger, He_res, T, V, min_fams)


def run(T, subs, cond, resid_live_from_probe=True):
    from chemkit.converg import _live
    from chemkit.engine import judge
    pr = {}
    r = judge(subs, cond, T, _probe=pr)
    steps = [(s["equation"], s["extent"]) for s in (r.get("steps") or [])]
    tail = steps[-6:]
    cyc = (len(tail) == 6
           and abs(tail[-1][1] - tail[-3][1]) < 1e-9
           and abs(tail[-2][1] - tail[-4][1]) < 1e-9)
    return {
        "pH": r.get("final_pH"), "resid": abs(_live(pr.get("active") or [])),
        "iters": pr.get("iters"), "exit": pr.get("exit"),
        "n_steps": len(steps), "cycle": cyc,
        "h2sio3": next((e["mol"] for e in (r.get("final") or [])
                        if e["name"] == "H_2SiO_3"), 0.0),
    }


def main() -> int:
    T = load_tables()
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = {c["name"]: c for c in json.load(fh)}
    hit = [c for k, c in cases.items() if k.startswith("N15")][0]
    cond = hit.get("cond") or {}
    base = [{"name": s[0], "mol": float(s[1])} for s in hit["subs"]]
    print(f"=== {hit['name']} ===")
    print("投料:", [(s['name'], s['mol']) for s in base],
          " cond:", cond)
    print("\n【新路径关闭 = 改前行为】CO₂ 量微扰扫描")
    print(f"{'CO2/mol':>12} {'pH':>7} {'resid':>8} {'iters':>6} {'exit':>10} "
          f"{'n_steps':>8} {'环?':>4} {'H2SiO3':>8}")
    _sp._exact_ok = _off
    try:
        for f in (0.990, 0.995, 0.999, 1.0, 1.001, 1.005, 1.010,
                  0.97, 0.95, 1.03, 1.05):
            subs = [dict(s) for s in base]
            for s in subs:
                if s["name"].startswith("CO_2"):
                    s["mol"] = round(s["mol"] * f, 6)
            o = run(T, subs, cond)
            print(f"{subs[-1]['mol']:>12.6g} {str(o['pH']):>7} "
                  f"{o['resid']:>8.4f} {str(o['iters']):>6} "
                  f"{str(o['exit'])[:10]:>10} {o['n_steps']:>8} "
                  f"{'是' if o['cycle'] else '否':>4} {o['h2sio3']:>8.4f}")
    finally:
        _sp._exact_ok = _ORIG
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
