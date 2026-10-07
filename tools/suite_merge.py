# -*- coding: utf-8 -*-
"""合并 `--shard=i/N` 的分片结果 → `logs/suite-latest.json`。

为什么需要它（第 296 轮）：套件写文件是共享资源，`dev.py suite` 每次都写
`logs/suite-latest.json`，于是**两个分片并行跑会互相覆盖** —— 这就是
"套件修改文件导致不能并行操作"的根因。分片各自写
`logs/suite-shard<i>-of<n>.json`，本脚本负责按 `index` 合并回单一全量留档。

用法：
    python tools/suite_merge.py 4            # 合并 shard0..3（各 of4）
    python tools/suite_merge.py 4 --check    # 顺带检查分片完整性（无缺失/重复）
"""
from __future__ import annotations

import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

LOGS = os.path.join(ROOT, "logs")


def main(argv: list[str]) -> int:
    n = int(argv[0]) if argv and argv[0].isdigit() else 4
    parts, missing = [], []
    for i in range(n):
        p = os.path.join(LOGS, f"suite-shard{i}-of{n}.json")
        if not os.path.exists(p):
            missing.append(p)
            continue
        with io.open(p, encoding="utf-8") as fh:
            parts.append(json.load(fh))
    if missing:
        print("!! 缺分片：" + ", ".join(os.path.basename(m) for m in missing))
        return 2

    rows: list[dict] = []
    for p in parts:
        rows.extend(p)
    rows.sort(key=lambda r: r.get("index", 0))
    seen: dict[int, int] = {}
    for r in rows:
        seen[r.get("index", 0)] = seen.get(r.get("index", 0), 0) + 1
    dup = {k: v for k, v in seen.items() if v > 1}
    n_ok = sum(1 for r in rows if r.get("ok"))
    print(f"[合并] {n} 个分片 → {len(rows)} 例；通过 {n_ok}/{len(rows)}")
    if dup:
        print(f"!! 重复 index {len(dup)} 个：{sorted(dup)[:10]}")
    idx = sorted(seen)
    gaps = [k for k in range(1, (max(idx) if idx else 0) + 1) if k not in seen]
    if gaps:
        print(f"!! 缺 index {len(gaps)} 个：{gaps[:10]}")
    if not dup and not gaps:
        print("[检查] index 连续且无重复 ✓")

    dst = os.path.join(LOGS, "suite-latest.json")
    if os.path.exists(dst) and "--force" not in argv:
        print(f"!! {os.path.relpath(dst, ROOT)} 已存在（可能是全量直跑的结果）。")
        print("   加 --force 覆盖；或先备份。")
        return 3
    with io.open(dst, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(rows, fh, ensure_ascii=False)
    print(f"[写出] {os.path.relpath(dst, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
