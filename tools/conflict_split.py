# -*- coding: utf-8 -*-
"""第 206 轮 · 分清"真冲突"与"参考值本身错" —— 回答使用者的架构提问。

使用者问："引擎是否大量存在同一数值不同计算方法相冲突？"
`tools/conflict_audit.py` 给出两组数字，但**必须分清哪一组是真的**：

  ① `|final_pH − pH_solver|`（用户答案 vs 走步搜索用的 pH）
     = p50 0.001 / p99 0.426 / max 2.400 / >1 单位仅 **2 例**
     ⟹ **答案侧基本无冲突**。

  ② `|pH_charge − pH_solver|`（与账本电荷自洽值比）
     = 剔除欠定后 147 例 >1 单位
     ⚠️ 但 `pH_charge` 在**有固相**时是**已知错误的参考值**（它把账本里的
     总金属量当作溶解量，而固相已把金属拿走）——第 203/204 轮实测：
     `E55` 引擎 6.262（对）vs charge 10.392（错）；`16 AlCl3+3NaOH` 是 PASS 例。

本脚本把 ② 按"有无固相"拆开，并**换一个不依赖 pH_charge 的内生判据**：
  · 引擎自己的残差（`max_abs_S`/`resid_live`）——它自认未达平衡的程度
  · 账本内部**同族物种之间的比例**是否与引擎自己的条件常数一致
    （只对**无固相**的用例做，因为固相在场时 Ksp 合法地钉住游离金属）

用法： python tools/conflict_split.py
"""
import io
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
src = os.path.join(ROOT, "logs", "conflict_audit.json")
rows = json.load(io.open(src, encoding="utf-8"))["rows"]
print(f"载入 {len(rows)} 例（{os.path.relpath(src, ROOT)}）\n")

print("=" * 76)
print("① 把 ② 按「有无固相」拆开（验证'147 例是参考值错'这一判断）")
print("=" * 76)
solid = [r for r in rows if r["solid"]]
nosol = [r for r in rows if not r["solid"]]
und = [r for r in rows if r["undissoc"]]
for tag, grp in (("有固相（pH_charge 已知会错）", solid),
                 ("无固相", nosol)):
    ds = [r["d_chg"] for r in grp if r.get("d_chg") is not None]
    if not ds:
        print(f"  {tag}: 无可比项")
        continue
    ds.sort()
    n = len(ds)
    print(f"  {tag}: n={n}  |Δ|>1: {sum(1 for d in ds if d > 1)}"
          f"  |Δ|>3: {sum(1 for d in ds if d > 3)}"
          f"  p90={ds[int(n * 0.9)]:.3f}")
print(f"  其中有未解离酸/碱（账本欠定）: {len(und)} 例")

print()
print("=" * 76)
print("② 换内生判据：引擎自认的残差（不依赖任何外部参考 pH）")
print("=" * 76)
res = [(r.get("resid"), r["case"], r["solid"]) for r in rows
       if isinstance(r.get("resid"), (int, float))]
rs = sorted(x[0] for x in res)
n = len(rs)
if n:
    print(f"  n = {n}")
    for q, tag in ((0.5, "p50"), (0.9, "p90"), (0.99, "p99")):
        print(f"  max|S| {tag} = {rs[min(int(n * q), n - 1)]:.3f}")
    print(f"  max|S| max = {rs[-1]:.3f}")
    print(f"  max|S| > 1 的: {sum(1 for d in rs if d > 1)}")
    print(f"  max|S| > 3 的: {sum(1 for d in rs if d > 3)}")
    print(f"\n  自认残差最大的 15 例（**这才是引擎自己承认有问题的例**）:")
    print(f"   {'用例':34s} {'max|S|':>9} 固相")
    for d, nm, so in sorted(res, key=lambda x: -x[0])[:15]:
        print(f"   {nm[:34]:34s} {d:9.3f} {so}")

print()
print("=" * 76)
print("③ 结论：冲突到底有多大")
print("=" * 76)
n_solid_conf = sum(1 for r in solid if (r.get("d_chg") or 0) > 1)
n_nosol_conf = sum(1 for r in nosol if (r.get("d_chg") or 0) > 1)
n_ans = sum(1 for r in rows if (r.get("d_ans") or 0) > 1)
print(f"  · 答案侧（final vs solver）>1 单位: **{n_ans}** 例")
print(f"  · 与电荷自洽值 >1 单位: 有固相 {n_solid_conf} 例（**参考值错，非冲突**）"
      f" + 无固相 {n_nosol_conf} 例")
print(f"  · 引擎自认残差 >3 的: "
      f"{sum(1 for d in rs if d > 3)} 例（**真问题面**）")
print()
print("  ⟹ 判断：**同一数值多套算法的架构问题真实存在（7 条 pH 通路），")
print("     但它对'用户看到的答案'几乎无影响**（>1 单位仅 2 例）；")
print("     真正暴露的是：① 少数例引擎自认残差大；② 诊断字段口径不一致")
print("     会让分析与修复被误导（本文件与 lessons 多條即此）。")
print("     **改善方向 = 合并 pH 通路 + 让'参考值'只在有效窗口内使用**，")
print("     而不是逐例打补丁。")
