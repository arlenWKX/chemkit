# -*- coding: utf-8 -*-
"""逐例跨档追踪：把若干个用例在**多份 `logs/suite-*.json` 留档**里的
关键指标并排打印，用于"某例在 A/B/C 三个版本里各是多少"的快速归因。

## 为什么需要

`tools/acc_metrics.py` 只给**全局**指标与"增减最多 8 例"；
`tools/flip_census.py` 只给通过性翻转；两者都回答不了
"**`H35` 在 r271 / r272 / r273 三版里各是多少**"。
而第 272–273 轮的归因恰恰全卡在这种逐例三档对照上
（`H35` 0.363 → 10.055 → 11.094，必须同时看到三档才能判断
"是重建引起的"还是"走步轨迹被动改变"）。

用法：
    python tools/case_track.py H35 Ni41 -- logs/suite-r271b-1331.json logs/suite-r272d-1337.json logs/suite-r273a-1337.json
    python tools/case_track.py H35 --glob r27
"""

from __future__ import annotations

import glob
import json
import os
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(p: str) -> dict:
    with open(p, encoding="utf-8") as fh:
        d = json.load(fh)
    rows = d if isinstance(d, list) else (d.get("rows") or d.get("cases"))
    return {r["name"]: r for r in rows}


def main(argv: list[str]) -> int:
    if "--" not in argv:
        print(__doc__)
        return 2
    i = argv.index("--")
    wants, paths = argv[:i], argv[i + 1:]
    if "--glob" in wants:
        j = wants.index("--glob")
        pat = wants[j + 1]
        wants = [w for k, w in enumerate(wants) if k not in (j, j + 1)]
        paths = sorted(glob.glob(os.path.join(ROOT, "logs", f"suite-*{pat}*.json")))
    if not paths:
        print("没给留档文件")
        return 2
    arch = [(os.path.basename(p).replace("suite-", "").replace(".json", ""),
             _load(p)) for p in paths]
    for w in wants:
        keys = [k for k in arch[0][1] if k.split()[0] == w] or \
               [k for k in arch[0][1] if k.startswith(w)]
        if not keys:
            # 该前缀可能在较新的档里才有 ⟹ 逐档找
            for _nm, d in arch:
                keys = [k for k in d if k.split()[0] == w]
                if keys:
                    break
        if not keys:
            print(f"[{w}] 各档都没有")
            continue
        key = keys[0]
        print(f"\n=== {key} ===")
        print(f"{'档':<20}{'ok':>5}{'resid':>10}{'pH':>9}{'iters':>8}"
              f"{'exit':>14}")
        for nm, d in arch:
            r = d.get(key)
            if r is None:
                print(f"{nm:<20}{'—':>5}{'（本档无此例）':>10}")
                continue
            print(f"{nm:<20}{('✓' if r.get('ok') else '✗'):>5}"
                  f"{(r.get('resid_live') or 0.0):>10.4f}"
                  f"{(r.get('pH') if r.get('pH') is not None else float('nan')):>9.3f}"
                  f"{(r.get('iters') or 0):>8}{str(r.get('exit')):>14}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
