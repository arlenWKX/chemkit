"""跨版本逐例对照：按 He 量级分层统计 pH 位移，并给指定用例的明细。

用法：python tools/xver.py [旧快照.json 新快照.json] [用例名前缀 ...]
"""
from __future__ import annotations

import json
import sys
from collections import Counter

_DEF = ("_snap_before.json", "_snap_after2.json")


def band(h):
    if h is None:
        return "无"
    a = abs(h)
    if a == 0.0:
        return "He=0"
    if a < 1e-6:
        return "<1e-6"
    if a < 1e-4:
        return "1e-4..1e-6"
    if a < 1e-3:
        return "1e-3..1e-4"
    if a < 1e-2:
        return "1e-2..1e-3"
    return ">=1e-2"


def main(argv):
    js = [a for a in argv if a.endswith(".json")]
    fa, fb = (js + list(_DEF))[:2]
    prefixes = [a for a in argv if not a.endswith(".json")]
    A = json.load(open(fa, encoding="utf-8"))
    B = json.load(open(fb, encoding="utf-8"))
    print(f"A={fa}  B={fb}")

    hb, shift, dmax = Counter(), Counter(), Counter()
    for n, a in A.items():
        b = B.get(n)
        if b is None or a["He"] is None:
            continue
        k = band(a["He"])
        hb[k] += 1
        d = (b["pH"] or 0) - (a["pH"] or 0)
        if abs(d) > 0.005:
            shift[k] += 1
            dmax[k] = max(dmax[k], abs(d))
    print("\n按**旧版** He 量级分层的 pH 位移分布：")
    print(f"{'band':14s} {'例数':>6s} {'pH变':>6s} {'max|dPH|':>10s}")
    for k in ("He=0", "<1e-6", "1e-4..1e-6", "1e-3..1e-4",
              "1e-2..1e-3", ">=1e-2", "无"):
        if hb[k]:
            print(f"{k:14s} {hb[k]:6d} {shift[k]:6d} {dmax[k]:10.3f}")

    ra = [a["resid"] for a in A.values() if a["resid"] is not None]
    rb = [b["resid"] for b in B.values() if b["resid"] is not None]
    print(f"\n残差 max|S|: 旧 max={max(ra):.1f} >0.1={sum(1 for x in ra if x > 0.1)}"
          f" >1={sum(1 for x in ra if x > 1)}   "
          f"新 max={max(rb):.1f} >0.1={sum(1 for x in rb if x > 0.1)}"
          f" >1={sum(1 for x in rb if x > 1)}")
    ita = sum(a["iters"] or 0 for a in A.values())
    itb = sum(b["iters"] or 0 for b in B.values())
    print(f"iters 合计 {ita} → {itb} ({(itb - ita) / ita * 100:+.1f}%)")

    if prefixes:
        for n, a in A.items():
            if not n.startswith(tuple(prefixes)):
                continue
            b = B.get(n, {})
            print(f"\n== {n}")
            print(f"   pH  {a['pH']} → {b.get('pH')}   He  {a['He']} → {b.get('He')}"
                  f"   旧band {band(a['He'])}")
            print(f"   deg {a['deg']}→{b.get('deg')} changed {a['changed']}→"
                  f"{b.get('changed')} resid {a['resid']}→{b.get('resid')} "
                  f"iters {a['iters']}→{b.get('iters')} steps {a['steps']}→"
                  f"{b.get('steps')}")
            la, lb = a["led"], b.get("led") or {}
            for k in sorted(set(la) | set(lb)):
                va, vb = la.get(k, 0.0), lb.get(k, 0.0)
                if abs(va - vb) > 1e-9:
                    print(f"     {k:24s} {va:12.6g} → {vb:12.6g}")


if __name__ == "__main__":
    main([x for x in sys.argv[1:] if not x.startswith("-")])
