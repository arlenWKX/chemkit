"""dump 计时/迭代对照：找出变慢的用例。

用法：python tools/perf.py 旧.json 新.json [-n 20]
"""
from __future__ import annotations

import json
import sys


def main(a: str, b: str, n: int = 20) -> None:
    A = {x["name"]: x for x in json.load(open(a, encoding="utf-8"))["cases"]}
    B = {x["name"]: x for x in json.load(open(b, encoding="utf-8"))["cases"]}
    rows = []
    for k, x in A.items():
        y = B.get(k)
        if y is None:
            continue
        rows.append((y["ms"] - x["ms"], k, x, y))
    rows.sort(key=lambda t: -t[0])
    print(f"== ms 增幅最大 {n} 例 ==")
    for d, k, x, y in rows[:n]:
        print(f"  {d:+9.1f}ms  {x['ms']:8.1f}→{y['ms']:8.1f}  "
              f"iters {x.get('iters')}→{y.get('iters')}  {k[:50]}")
    print(f"\n== ms 降幅最大 {n} 例 ==")
    for d, k, x, y in rows[-n:]:
        print(f"  {d:+9.1f}ms  {x['ms']:8.1f}→{y['ms']:8.1f}  "
              f"iters {x.get('iters')}→{y.get('iters')}  {k[:50]}")
    da = sum(x["ms"] for x in A.values())
    db = sum(y["ms"] for y in B.values())
    print(f"\n总时长 {da:.0f} → {db:.0f} ms ({(db - da) / da * 100:+.1f}%)")
    n_it = sum(1 for k, x in A.items()
               if B.get(k) and x.get("iters") != B[k].get("iters"))
    print(f"iters 有变化的用例 {n_it}/{len(A)}")
    ia = sum(x.get("iters", 0) for x in A.values())
    ib = sum(y.get("iters", 0) for y in B.values())
    print(f"iters 合计 {ia} → {ib} ({(ib - ia) / ia * 100:+.1f}%)")


if __name__ == "__main__":
    args = [x for x in sys.argv[1:] if not x.startswith("-")]
    nn = 20
    if "-n" in sys.argv:
        nn = int(sys.argv[sys.argv.index("-n") + 1])
    main(args[0], args[1], nn)
