"""A/B：二分收敛容差到底省了多少？——用**确定性**指标裁决，不用墙钟。

背景（architecture §7 W-4）：共享沙箱里同一份代码的三次全量墙钟是
29.1 / 32.5 / 33.9 ms（均值逐轮漂 15%），墙钟不足以裁决"求根改动是否
真的省了"。`S_of` 调用次数（`engine.SOF_CALLS`）是确定性的：同一份代码
⇒ 同一计数，而且它正是走步每步的钱主要花在哪（§7 U-1 剖面）。

做法：同一进程内跑三遍全量——off(容差=0，跑满到浮点饱和) / on(1e-11)
/ off。第一遍 off 与第三遍 off 同时充当**过程内确定性校验**（同代码同
进程应给同一 digest 与同一 SOF 计数），并暴露"第二遍更热"的墙钟偏置
（这正是各版本墙钟数字不可比的原因）。

用法：python tools/sof_ab.py [--passes 3]
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from chemkit import engine                      # noqa: E402
from chemkit.converg import run_all             # noqa: E402

TOL_ON = 1e-11


def measure(tol: float) -> tuple[list[dict], float]:
    engine._EXTENT_TOL_REL = tol
    t0 = time.perf_counter()
    recs = run_all(count_sof=True)
    return recs, time.perf_counter() - t0


def main() -> None:
    n_passes = 3
    if "--passes" in sys.argv:
        n_passes = int(sys.argv[sys.argv.index("--passes") + 1])

    plan = [0.0, TOL_ON] + [0.0] * (n_passes - 2)
    runs = []
    for i, tol in enumerate(plan, 1):
        recs, dt = measure(tol)
        sof = sum(r["sof"] for r in recs)
        iters = sum(r.get("iters", 0) for r in recs)
        steps = sum(r["steps_n"] for r in recs)
        runs.append((tol, recs, dt, sof, iters, steps))
        tag = "OFF(跑满浮点饱和)" if tol == 0.0 else f"ON (tol={tol:g})"
        print(f"pass{i} {tag:22s} 墙钟 {dt:6.1f}s  SOF {sof:7d}  "
              f"iters {iters:6d}  steps {steps:6d}")
        sys.stdout.flush()

    # ---- 确定性校验：同代码同进程的两遍 off 必须逐例一致
    off1, off2 = runs[0][1], runs[-1][1]
    if runs[-1][0] == 0.0 and len(runs) > 2:
        d = [b["name"] for a, b in zip(off1, off2) if a["digest"] != b["digest"]]
        print(f"\n确定性校验（off 首遍 vs off 末遍）：digest 位移 {len(d)} 例"
              f"{'  ← 过程内非确定！' if d else '  ✓'}")
        print(f"  逐例 SOF 位移 {sum(1 for a, b in zip(off1, off2) if a['sof'] != b['sof'])} 例")

    # ---- 主裁决：off(首遍，最冷) vs on
    on = next(r for r in runs if r[0] != 0.0)
    a, b = runs[0][1], on[1]
    print(f"\n== 主裁决：off(首遍) → on ==")
    print(f"SOF 调用 {runs[0][3]} → {on[3]}  ({(on[3] - runs[0][3]) / runs[0][3] * 100:+.1f}%)")
    print(f"iters    {runs[0][4]} → {on[4]}  ({(on[4] - runs[0][4]) / runs[0][4] * 100:+.1f}%)")
    print(f"steps    {runs[0][5]} → {on[5]}  ({(on[5] - runs[0][5]) / runs[0][5] * 100:+.1f}%)")
    print(f"语义位移 {sum(1 for x, y in zip(a, b) if x['digest'] != y['digest'])} 例")
    sav = sorted(((x["sof"] - y["sof"], x["name"], x["sof"], y["sof"])
                  for x, y in zip(a, b)), reverse=True)
    print("\n省得最多的 15 例（off → on，SOF 调用数）：")
    for d, nm, so, sn in sav[:15]:
        print(f"  {d:7d}  {nm[:52]:54s} {so:7d} → {sn:7d}")
    print("\n反而变多的 10 例：")
    for d, nm, so, sn in sav[-10:]:
        print(f"  {d:7d}  {nm[:52]:54s} {so:7d} → {sn:7d}")


if __name__ == "__main__":
    main()
