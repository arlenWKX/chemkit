# -*- coding: utf-8 -*-
"""第 263 轮 · 抓出"电荷平衡精确解(两性)"接管时的**完整账本**。

`estimate_state` 被 `estimate_pH` 按模块全局名调用，故包一层即可零侵入
记录：命中新标签的那次调用，问的是什么账本、得了什么 pH。

用法：python tools/amph_hit.py N15
      python tools/amph_hit.py TC1
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
from chemkit.candidates import WATER             # noqa: E402
from chemkit.data import load_tables             # noqa: E402

TAG = "电荷平衡精确解(两性)"
LOG: list = []
_ORIG = _sp.estimate_state


def wrap(ledger, H_excess, V, T, T_K, *a, **k):
    _sp.PH_TAGS = []
    try:
        r = _ORIG(ledger, H_excess, V, T, T_K, *a, **k)
        tags = list(_sp.PH_TAGS)
    finally:
        _sp.PH_TAGS = None
    if TAG in tags:
        LOG.append({
            "He": H_excess,
            "ledger": {s: m for s, m in sorted(ledger.items())
                       if m > 0.0 and s != WATER},
            "pH_out": round(r[0], 6),
        })
    return r


def main(argv: list[str]) -> int:
    name = argv[0] if argv else "N15"
    T = load_tables()
    _sp.estimate_state = wrap
    try:
        with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                     encoding="utf-8") as fh:
            cases = {c["name"]: c for c in json.load(fh)}
        hit = [c for k, c in cases.items() if k.startswith(name)][0]
        from chemkit.engine import judge
        pr = {}
        subs = [{"name": s[0], "mol": float(s[1])} for s in hit["subs"]]
        r = judge(subs, hit.get("cond") or {}, T, _probe=pr)
    finally:
        _sp.estimate_state = _ORIG
    print(f"=== {hit['name']} ===  pH={r.get('final_pH')} "
          f"iters={pr.get('iters')} exit={pr.get('exit')}")
    print(f"新路径接管次数：{len(LOG)}")
    for i, e in enumerate(LOG[:20], 1):
        led = "  ".join(f"{k}={v:.6g}" for k, v in e["ledger"].items())
        print(f"\n[{i}] He={e['He']:+.6e}  pH_out={e['pH_out']}")
        print(f"    {led}")
    if len(LOG) > 20:
        print(f"\n… 另有 {len(LOG) - 20} 次")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
