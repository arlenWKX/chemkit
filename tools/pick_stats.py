# -*- coding: utf-8 -*-
"""第 264 轮 · 量「同产物可动量回退」的**触发面**（第 263 轮量分支占比的同款做法）。

`_pick_ev` 是 `engine.py` 的模块级函数、被 `engine.py` 内部按全局名调用，
故包一层即可零侵入统计：① 全库被问了多少次；② 回退真正**改写了拾取**多少次；
③ 落在多少个用例上。

⚠️ 串行跑（~3–4 min）。运行期间**冻结 chemkit/**。
输出 `logs/pick_stats.json`。

用法：python tools/pick_stats.py
"""
from __future__ import annotations

import collections
import io
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as _eng                     # noqa: E402
from chemkit.data import load_tables              # noqa: E402

_ORIG = _eng._pick_ev
STATS = {"calls": 0, "fallback": 0, "cand_sum": 0}
PER_CASE = collections.Counter()
CASE = [""]


def _wrap(evals, ledger):
    r = _ORIG(evals, ledger)
    STATS["calls"] += 1
    STATS["cand_sum"] += len(evals)
    top = max(evals, key=lambda e: e[2])
    if r is not top:
        STATS["fallback"] += 1
        PER_CASE[CASE[0]] += 1
    return r


def main() -> int:
    _eng._pick_ev = _wrap
    T = load_tables()
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    from chemkit.engine import judge
    t0 = time.time()
    rows = []
    for i, c in enumerate(cases, 1):
        CASE[0] = c["name"]
        try:
            judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                  c.get("cond") or {}, T)
            rows.append({"name": c["name"], "n": PER_CASE.get(c["name"], 0)})
        except Exception as exc:                                # noqa: BLE001
            rows.append({"name": c["name"], "err": repr(exc)})
        if i % 100 == 0:
            print(f"  {i}/{len(cases)} … {time.time() - t0:.0f}s", flush=True)
    _eng._pick_ev = _ORIG
    hot = sorted((r for r in rows if r.get("n")),
                 key=lambda r: -r["n"])
    out = {"calls": STATS["calls"], "fallback": STATS["fallback"],
           "cand_sum": STATS["cand_sum"],
           "rate": STATS["fallback"] / max(1, STATS["calls"]),
           "n_cases_fired": len(hot), "n_cases": len(cases),
           "hot": hot, "wall_s": round(time.time() - t0, 1)}
    with io.open(os.path.join(ROOT, "logs", "pick_stats.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(f"\n拾取调用 {out['calls']} 次（平均候选 {out['cand_sum'] / max(1, out['calls']):.2f} 个）")
    print(f"回退改写拾取 {out['fallback']} 次 = **{out['rate'] * 100:.4f}%**")
    print(f"落在 {out['n_cases_fired']} / {out['n_cases']} 个用例上；"
          f"串行墙钟 {out['wall_s']}s")
    print("\n触发最多的 20 例：")
    for r in hot[:20]:
        print(f"   {r['name'][:40]:<40} {r['n']}")
    print("\n-> logs/pick_stats.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
