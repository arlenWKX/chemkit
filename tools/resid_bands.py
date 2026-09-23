# -*- coding: utf-8 -*-
"""第 230 轮 · `p90` 的真实驱动：0.1~0.5 带的归因构成。

第 229 轮算清：`resid_p90` 由排序后第 `int(n*0.9)` 个值决定（现 0.447），
要降到 0.1 需修好 **82 例 `>0.1`**；而"签名"只覆盖 32 例（且都是大幅度例）。
⟹ **决定 `p90` 的是长尾中段那一大批低幅度残差**，不是最大的几个。

本脚本（读档，不重跑）回答：**那一带是什么机制**——同一种（可一修），
还是散乱（说明是"普遍的小不收敛"，需另想办法）？

分桶：`[0.1,0.3)` / `[0.3,1)` / `[1,5)` / `[5,∞)`，每桶统计
用例数 + 归因反应（`resid_src_eq`）的**反应类型**分布 + 代表例。

用法： python tools/resid_bands.py
"""
import collections
import io
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

d = json.load(io.open(os.path.join(ROOT, "logs", "suite-latest.json"),
                      encoding="utf-8"))
try:
    cen = json.load(io.open(os.path.join(ROOT, "logs",
                                         "signature_census.json"),
                            encoding="utf-8"))
    hits = {h[0] for h in cen["hits"]}
except Exception:                                           # noqa: BLE001
    hits = set()

BANDS = [(0.1, 0.3), (0.3, 1.0), (1.0, 5.0), (5.0, 1e9)]
print(f"全库 {len(d)} 例；签名命中 {len(hits)} 例\n")

# 归因反应的"类型"粗分类
def rtype(eq: str) -> str:
    if not eq:
        return "(无归因)"
    if "-->" in eq or "->" in eq:
        pass
    if "H^+" in eq or "OH^-" in eq:
        pass
    for kw, name in (("(OH)", "羟合/水解"), ("Cl]", "氯合"), ("CO_3", "碳酸盐"),
                     ("SO_3", "亚硫酸盐"), ("SO_4", "硫酸盐"), ("SCN", "硫氰"),
                     ("NH_3", "氨"), ("HCO_3", "碳酸氢"), ("O_2", "氧/氧化"),
                     ("H_2", "氢"), ("S^2-", "硫"), ("PO_4", "磷酸")):
        if kw in eq:
            return name
    return "其它"


for lo, hi in BANDS:
    rows = [c for c in d
            if lo <= (c.get("resid_live") or 0.0) < hi]
    if not rows:
        continue
    nh = sum(1 for c in rows if c["name"] in hits)
    print("=" * 74)
    print(f"band [{lo}, {hi})：**{len(rows)} 例**（其中签名命中 {nh}）")
    print("=" * 74)
    tp = collections.Counter(rtype(c.get("resid_src_eq") or "") for c in rows)
    print("  归因类型分布:")
    for k, v in tp.most_common(8):
        print(f"    {k:12s} {v:3d}")
    print("  代表例（前 8，按 resid 降序）:")
    for c in sorted(rows, key=lambda c: -(c.get("resid_live") or 0))[:8]:
        mark = "★" if c["name"] in hits else " "
        print(f"   {mark}{c['name'][:34]:34s} "
              f"{(c.get('resid_live') or 0):7.3f}  "
              f"{str(c.get('resid_src_eq'))[:40]}")
    print()

print("=== 判读 ===")
print("  · 若 [0.1,0.3) 带（决定 p90 的那批）**类型集中** ⟹ 可能是同一机制，可一修。")
print("  · 若**类型散乱、且多为不同体系** ⟹ 是'普遍的小不收敛'，")
print("    修单个机制无法把 p90 压到 0.1，应把验收指标改为")
print("    `resid_max`/残差质量（那才是'正确性'），或分列两个指标。")
