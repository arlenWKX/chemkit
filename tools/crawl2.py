# -*- coding: utf-8 -*-
"""第 264 轮 · N15 爬行态的**末态探针**：两个交替候选的 S / x_max / 账本。

`[pick]` 只给被选中者的 S；要判"为什么每步只有 2e-4"，必须同时拿到
**两个候选的 x_max**（`solve_extent` 的返回值）与末态账本。

用法：python tools/crawl2.py N15
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

from chemkit.data import load_tables             # noqa: E402
from chemkit.converg import _live                 # noqa: E402


def main(argv: list[str]) -> int:
    pre = argv[0] if argv else "N15"
    T = load_tables()
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = {c["name"]: c for c in json.load(fh)}
    c = [v for k, v in cases.items() if k.startswith(pre)][0]
    from chemkit.engine import judge
    pr = {}
    subs = [{"name": s[0], "mol": float(s[1])} for s in c["subs"]]
    r = judge(subs, c.get("cond") or {}, T, _probe=pr)
    print(f"=== {c['name']} ===  pH={r.get('final_pH')} "
          f"iters={pr.get('iters')} exit={pr.get('exit')} "
          f"resid={abs(_live(pr.get('active') or [])):.4f}")
    print("\n-- probe 字段 --")
    for k in sorted(pr):
        v = pr[k]
        if k == "active":
            continue
        print(f"   {k} = {str(v)[:160]}")
    act = pr.get("active") or []
    print(f"\n-- active 候选（{len(act)}）--")
    for a in act:
        if isinstance(a, dict):
            print("   " + json.dumps(a, ensure_ascii=False)[:300])
        else:
            print("   " + str(a)[:300])
    print("\n-- 末态账本（>1e-7）--")
    for e in sorted(r.get("final") or [], key=lambda e: -e["mol"]):
        if e["mol"] > 1e-7:
            print(f"   {e['name']:<22} {e['mol']:.8g}")
    print("\n-- 末态步表尾部 12 --")
    for s in (r.get("steps") or [])[-12:]:
        print(f"   {s['equation'][:74]:<74} x={s['extent']:.8g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
