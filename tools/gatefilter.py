# -*- coding: utf-8 -*-
"""第 204 轮 · 合取条件过滤器设计（用已有的全量 sweep 数据，不重跑）。

背景（第 203 轮结论）：
  · 整体窄化 `exact_proton_pH` 的闸（L1191）已**实测否决**——子代理实测
    新放行 47 例中 24 例会变、13 例动 >1 pH 单位（EN02 −1.0→+6.7 等）。
  · F31 与它们的本质区别待确认：推测是"**账本里没有强酸储备**"。

本脚本在 `logs/phaudit_sweep_step1.json`（子代理产出的全量 1176 例逐例数据）
上直接套用**逐级收紧的合取条件**，统计每级剩下多少例、分别是哪些，
从而把"安全且充分的判据"找出来——**不重跑套件**。

条件（逐级叠加）：
  C1 现在被闸拒（`exact is None`）而**放宽版**给得出值（`exact_solid_only` 是 float）
  C2 无固相在场（`solids` 空）——C1 已隐含，但显式核对
  C3 **无强酸/强碱条件**（`c_H`/`c_OH` 均非正、`c_pH` 未设定）
  C4 账本电荷**自相矛盾**（|net + He| 超阈）—— 这才是要修的那种态
  C5 账本里**没有游离强酸/强碱储备**（用 `net` 与 He 的符号关系判）

用法： python tools/gatefilter.py [He 容差，默认 1e-6]
"""
import io
import json
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

TOL = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-6
src = os.path.join(ROOT, "logs", "phaudit_sweep_step1.json")
recs = json.load(io.open(src, encoding="utf-8"))["records"]
print(f"载入 {len(recs)} 例（{os.path.relpath(src, ROOT)}）\n")


def fnum(v):
    return v if isinstance(v, float) else None


def cond(name, rs):
    print(f"  {name:52s} {len(rs):4d} 例")
    return rs


def show(rs, tag, n=25):
    print(f"\n  --- {tag} 逐例（前 {n}）---")
    print(f"   {'用例':34s} {'现在pH':>8} {'放宽后':>9} {'Δ':>8} "
          f"{'net':>11} {'He':>11} 阳离子")
    for r in rs[:n]:
        e = fnum(r.get("est"))
        x = fnum(r.get("exact_solid_only"))
        d = f"{x - e:+8.3f}" if (e is not None and x is not None) else "-"
        print(f"   {str(r['name'])[:34]:34s} "
              f"{(f'{e:.4f}' if e is not None else '-'):>8} "
              f"{(f'{x:.4f}' if x is not None else '-'):>9} {d:>8} "
              f"{str(r.get('net'))[:11]:>11} "
              f"{str(r.get('probe_H_excess'))[:11]:>11} "
              f"{','.join(r.get('ksp_cats') or [])[:18]}")


S1 = [r for r in recs if r.get("exact") is None
      and isinstance(r.get("exact_solid_only"), float)]
print("=== 逐级收紧 ===")
cond("基线：全库", recs)
c1 = cond("C1 现闸拒 ∧ 放宽版给得出值（= narrowed_admit）", S1)
c2 = cond("C2 ∧ 无固相在场",
          [r for r in c1 if not (r.get("solids") or [])])
c3 = cond("C3 ∧ c_H/c_OH 均非正 ∧ c_pH 未设定",
          [r for r in c2
           if not (fnum(r.get("c_H")) or 0.0) > 0.0
           and not (fnum(r.get("c_OH")) or 0.0) > 0.0
           and r.get("c_pH") is None])
c4 = cond(f"C4 ∧ 账本电荷自相矛盾 |net+He| > {TOL:g}",
          [r for r in c3
           if abs((fnum(r.get("net")) or 0.0)
                  + (fnum(r.get("probe_H_excess")) or 0.0)) > TOL])
# C5：区分"账本自相矛盾来自被吃掉的碱"还是"来自游离强酸"。
# 判据：pH 侧读出的游离质子量 |He|/V 与账本 net 的**符号**是否一致。
# 账本 net>0 却 He≈0（或反向）⟹ 记账与组成打架 ⟹ 要修的那种态。
c5 = cond("C5 ∧ net 与 He 符号相反或 |He| 远小于 |net|",
          [r for r in c4
           if (fnum(r.get("net")) or 0.0) * (fnum(r.get("probe_H_excess")) or 0.0) < 0
           or abs(fnum(r.get("probe_H_excess")) or 0.0) < 0.01 * abs(fnum(r.get("net")) or 1.0)])

show(c1, "C1 narrowed_admit（整体窄化的风险面）")
show(c5, "C5 合取条件命中集（**本轮的候选修复面**）")

print("\n=== 判读 ===")
print(f"  整体窄化面 C1 = {len(c1)} 例；合取条件 C5 = {len(c5)} 例")
print("  目标：C5 足够小（≤5 例）且**含 F31** ⟹ 判据可用。")
print("  若 C5 仍很大 ⟹ 需要更强判据（或说明这类自相矛盾态在库里很常见，")
print("  修法应改为'不让走步产生这种态'而不是'事后放行精确解'）。")
f31 = [r for r in c5 if str(r["name"]).startswith("F31")]
print(f"\n  F31 是否在 C5 内: {'是 ✓' if f31 else '否 ✗（判据不充分）'}")

out = os.path.join(ROOT, "logs", "gatefilter.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"C1": [r["name"] for r in c1],
               "C5": [r["name"] for r in c5],
               "C5_detail": c5}, f, ensure_ascii=False, indent=1, default=str)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
