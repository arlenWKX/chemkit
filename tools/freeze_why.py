# -*- coding: utf-8 -*-
"""第 274 轮 · **冻结守卫的"为什么没认作强驱动"普查**。

## 为什么需要

`tools/tension.py`（第 274 轮改成与 `engine.resid_live_ok` 同一判据后）实测
`H35 Ca(HCO3)2 加热 353K`：

    live=11.094  pH=11.101  exit=no-cands
      |S|=11.094 kind=derived  CaCO_3 + 2H^+ -> Ca^{2+} + CO_2
      → **冻结在强驱动上**（|S| > RESID_FROZEN_TOL=0.1，X-38 病根一）
        走步判据：x*=0.01021 = x_max ⟹ **已执行**

即：一条 |S| = 11.094、求解器判"跑到底"的通道**被冻结**了，
走步因此不执行它，而质量口径如实把它计入残差。

引擎自己的契约（X-38 病根二）说得很清楚：`_freeze` 必须跳过
`_strong_keys()` 里的键（"|S| > FREEZE_MAX_S 就是它没达成的明证"）。
所以问题只能是 **`_strong_keys()` 把它滤掉了**，而它有 7 条 `continue`。
本工具读 `engine.FREEZE_AUDIT`（生产路径恒 None），把每条候选的**滤除原因**
聚合出来 —— 直接回答"哪一条闸误伤了"。

用法：
    python tools/freeze_why.py H35
    python tools/freeze_why.py Ni41 --top 20
"""

from __future__ import annotations

import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    top = 14
    for a in argv:
        if a.startswith("--top="):
            top = int(a[6:])
    if not args:
        print(__doc__)
        return 2
    pre = args[0]
    T = load_tables()
    hit = [c for c in load_cases(None) if c["name"].split()[0] == pre]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    eng.FREEZE_AUDIT = []
    try:
        r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                      c.get("cond") or {}, T)
    finally:
        rec = eng.FREEZE_AUDIT
        eng.FREEZE_AUDIT = None
    print(f"=== {c['name']} ===  pH={r.get('final_pH')}  "
          f"守卫审计 {len(rec)} 条")

    by_why: dict = defaultdict(int)
    for x in rec:
        by_why[x["why"]] += 1
    print(f"\n{'滤除原因':<22}{'条数':>8}")
    for w, n in sorted(by_why.items(), key=lambda kv: -kv[1]):
        print(f"{w:<22}{n:>8}")

    # 指定候选的**逐次**轨迹：同一候选在每次 `_strong_keys()` 调用里
    # 是"被滤掉（哪条闸）"还是"**STRONG**"。用于判定"是闸误伤、还是
    # 冻结发生时该通道确实弱、之后才涨起来"——两种成因的修法完全不同。
    eqf = None
    for a in argv:
        if a.startswith("--eq="):
            eqf = a[5:]
    if eqf:
        print(f"\n候选轨迹（--eq={eqf!r}）：")
        hits = [x for x in rec if eqf in x["eq"]]
        if not hits:
            print("  （守卫审计里没有这条候选）")
        for x in hits:
            sv = f"{abs(x['S']):>8.3f}" if x["S"] is not None else "     —  "
            print(f"  {x['why']:<14} dir={x['dir']:+d} |S|={sv}  {x['eq']}")

    # 最刺眼的一档：S 很强却被滤掉（按 |S| 排序）
    print(f"\n被滤掉的候选里 |S| 最大的 {top} 条"
          f"（这些是\"本可解、却没被认作强驱动\"的通道）：")
    cand = [x for x in rec if x["why"] != "**STRONG**"
            and x["S"] is not None and x["why"] != "weak_S"]
    seen = set()
    cand.sort(key=lambda x: -abs(x["S"]))
    for x in cand:
        if x["eq"] in seen:
            continue
        seen.add(x["eq"])
        print(f"  |S|={abs(x['S']):>8.3f}  {x['why']:<14} "
              f"kind={x['kind']:<9} dir={x['dir']:+d}")
        print(f"      {x['eq']}")
        if len(seen) >= top:
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
