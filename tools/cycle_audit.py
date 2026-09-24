# -*- coding: utf-8 -*-
"""触发点④（净质子循环）在 EU01 空转中的处置普查（第 289 轮）。

回答两个问题：
  ① 循环被识别了多少次、s0 落在哪个分支（s0≈0 / s0<0 冻 / 跳）；
  ② s0≈0 时各腿的 |S| 有多大 —— 区分"真平衡"与"张力"（腿强而净零）。

用法：python tools/cycle_audit.py EU01
"""
from __future__ import annotations

import sys
from collections import Counter

sys.path.insert(0, ".")

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit import engine                     # noqa: E402
from chemkit.data import load_tables          # noqa: E402
from chemkit.testsuit import load_cases       # noqa: E402


def main(keys: list[str]) -> int:
    T = load_tables()
    cases = load_cases(None)
    picks = [c for c in cases if c["name"].startswith(tuple(keys))]
    for c in picks:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        engine.CYCLE_AUDIT = []
        pr: dict = {}
        engine.judge(subs, c.get("cond") or {}, T, _probe=pr)
        au = engine.CYCLE_AUDIT
        engine.CYCLE_AUDIT = None
        print("=" * 76)
        print(f"{c['name']}  iters={pr.get('iters')}  循环识别 {len(au)} 次")
        band = Counter()
        for a in au:
            legs = a["legs"]
            mx = max((abs(s) for _k, s in legs), default=0.0)
            if abs(a["s0"]) <= 0.05:
                band[f"s0≈0 腿max|s|={'>2' if mx > 2 else '<=2'}"] += 1
            else:
                band[f"|s0|>0.05 ({'+' if a['s0'] > 0 else '-'})"] += 1
        for k, v in band.most_common():
            print(f"  {k}: {v}")
        print("  -- 抽样（首 2 / 末 3 条）--")
        for a in (au[:2] + au[-3:]):
            print(f"  hist={a['hist']} s0={a['s0']:+.4f} out={a.get('out')}"
                  f" s_hi={a.get('s_hi')}")
            for k, s in a["legs"]:
                print(f"      S={s:+8.3f}  {k}")
        outc = Counter(a.get("out", "?").split(" ")[0] for a in au)
        print("  -- 出口分布 --")
        for k, v in outc.most_common():
            print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["EU01"]))
