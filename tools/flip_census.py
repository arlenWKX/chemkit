# -*- coding: utf-8 -*-
"""第 272 轮 · **通过性翻转普查**：两份 `logs/suite-*.json` 逐例对拍
`ok` / `errors` / 关键数值，把"改进换了多少断言"摆清楚。

`tools/suite_show.py` 的 diff 模式只吃并行口径的旧格式（`d["fails"]`），
`tools/acc_metrics.py` 只比残差 —— 缺的正是"谁翻红、红在哪一条断言"。

用法：python tools/flip_census.py logs/A.json logs/B.json
"""

from __future__ import annotations

import json
import os
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def _load(p: str) -> list:
    with open(p, encoding="utf-8") as fh:
        d = json.load(fh)
    return d if isinstance(d, list) else (d.get("rows") or d.get("cases"))


def _errs(r: dict) -> tuple:
    e = r.get("errors")
    if isinstance(e, list):
        return tuple(e)
    return (e,) if e else ()


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    A = {r["name"]: r for r in _load(argv[0])}
    B = {r["name"]: r for r in _load(argv[1])}
    only_a = sorted(set(A) - set(B))
    only_b = sorted(set(B) - set(A))
    if only_a or only_b:
        print(f"[集合不同] 仅 A {len(only_a)} 例、仅 B {len(only_b)} 例")
    new_fail = [k for k in B if k in A and A[k].get("ok") and not B[k].get("ok")]
    new_pass = [k for k in B if k in A and not A[k].get("ok") and B[k].get("ok")]
    print(f"A ok={sum(1 for r in A.values() if r.get('ok'))}/{len(A)}   "
          f"B ok={sum(1 for r in B.values() if r.get('ok'))}/{len(B)}")
    print(f"翻红 {len(new_fail)} 例   翻绿 {len(new_pass)} 例\n")
    print(f"=== 翻红 {len(new_fail)} 例（逐条列出新出现的断言失败）===")
    for k in new_fail:
        a, b = A[k], B[k]
        old = set(_errs(a))
        new = [e for e in _errs(b) if e not in old]
        dr = (b.get("resid_live") or 0.0) - (a.get("resid_live") or 0.0)
        print(f"\n-- {k}")
        print(f"   pH {a.get('pH')} -> {b.get('pH')}   "
              f"resid {(a.get('resid_live') or 0):.3f} -> "
              f"{(b.get('resid_live') or 0):.3f} ({dr:+.3f})   "
              f"iters {a.get('iters')} -> {b.get('iters')}")
        for e in new[:6]:
            print(f"   ✗ {e}")
        if len(new) > 6:
            print(f"   … 另 {len(new) - 6} 条")
    print(f"\n=== 翻绿 {len(new_pass)} 例 ===")
    for k in new_pass:
        a, b = A[k], B[k]
        print(f"   {k}   resid {(a.get('resid_live') or 0):.3f} -> "
              f"{(b.get('resid_live') or 0):.3f}   pH {a.get('pH')} -> {b.get('pH')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
