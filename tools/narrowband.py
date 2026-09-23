# -*- coding: utf-8 -*-
"""第 205 轮 · F31 修复的**窄带命中面测量**（改 estimate_pH 之前必做）。

第 204 轮定影的唯一缺陷：`estimate_pH(ledger, He)` 在 `He` 跨过**账本自身
电荷不平衡量** `−Σz·n` 时硬切换分支（F31：−0.024007→12.3803 ✓，
−0.020000→1.6197 ✗，中间是完全平台）。

修法要**窄**。本脚本按 gateface 的分类口径（无固相 ∧ 在场 Ksp-OH 阳离子
= 190 例）逐例测量：

  · `pH_engine`  = 引擎自报（走步口径 `pH_solver`）
  · `pH_charge`  = 账本电荷自洽值（唯一确定：净电荷 → [OH⁻] 或 [H⁺]）
  · `He`         = 探针导出的 He（呈现口径，**仅用于筛"无强酸储备"**）
  · `net`        = Σz·n（除 H⁺/OH⁻）
  · `solid`      = 终态是否有固相
  · `strong`     = 账本里是否有游离强酸/强碱储备（判据见下）
  · `gap`        = |pH_charge − pH_engine|

**关键判据（本轮要验证的合取条件）**：
  ① 无固相在场            ② 在场 Ksp-OH 阳离子
  ③ 无强酸/强碱储备        ④ 电荷不平衡（|net| 显著）
  ⑤ `gap` 大（引擎读错悬崖）
只有 ①∧②∧③∧④∧⑤ 全真的用例才是"该修的面"。数一数有多大。

用法： python tools/narrowband.py [gap 阈值，默认 2.0] [步长，默认 1]
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
from chemkit.testsuit import load_cases, run_case           # noqa: E402

GAP = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
STEP = int(sys.argv[2]) if len(sys.argv) > 2 else 1
T = load_tables()

SOLIDS = frozenset(T.solids)
KSP_CATS = frozenset(e["pair"][0] for e in T.ksp if e["pair"][1] == "OH^-")
# 游离强酸/强碱储备：账本里的分子态强酸（HNO3/H2SO4 等）或游离 H+/OH-
STRONG_MOL = frozenset(("HNO_3", "H_2SO_4", "HCl", "HBr", "HI", "HClO_4"))

cases = load_cases(None)
rows = []
scanned = 0
for i, c in enumerate(cases):
    if i % STEP:
        continue
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    T_K = float(cond.get("T_K", 298.15))
    # 有强酸/强碱条件（c_H/c_OH/c_pH）的直接排除（第 203 轮教训：那些例
    # 的电荷不平衡由 c_H 解释，判据应当拒绝）
    if cond.get("c_H") or cond.get("c_OH") or cond.get("pH") is not None:
        continue
    probe = {}
    try:
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  {"V_L": V, "T_K": T_K}, T, _probe=probe)
    except Exception:                                       # noqa: BLE001
        continue
    led = dict(probe.get("ledger") or {})
    if not led:
        continue
    scanned += 1
    if any(m > spec.X_MIN and sp in SOLIDS for sp, m in led.items()):
        continue                                    # ① 有固相 ⟹ 排除
    cats = [sp for sp, m in led.items()
            if m > spec.X_MIN and sp in KSP_CATS]
    if not cats:
        continue                                    # ② 无 Ksp-OH 阳离子 ⟹ 排除
    if any(led.get(sp, 0.0) > spec.X_MIN for sp in STRONG_MOL):
        continue                                    # ③ 有分子态强酸 ⟹ 排除
    net = 0.0
    for sp, m in led.items():
        # `__`-前缀是引擎内部记账伪物种（如 `__tot_HNO_3`），不是真实物种；
        # `core.parse_species` 会因 `_` 报 FormulaError ⟹ 必须跳过
        # （与 `speciation.exact_proton_pH` 的 `sp.startswith("__")` 同口径）。
        if sp == "H_2O" or m <= 0 or sp.startswith("__"):
            continue
        net += spec.charge_of(sp) * m
    pKw = spec.pKw_of(T_K)
    if abs(net) < 1e-6:
        pH_charge = None                            # ④ 电荷已平衡 ⟹ 无需切换
    elif net > 0:
        pH_charge = pKw + math.log10(net / V)
    else:
        pH_charge = -math.log10(-net / V)
    pH_eng = probe.get("pH_solver")
    if pH_charge is None or not isinstance(pH_eng, float):
        continue
    gap = abs(pH_charge - pH_eng)
    rows.append({"case": c["name"], "pH_eng": pH_eng, "pH_charge": pH_charge,
                 "gap": gap, "net": net, "cats": cats,
                 "He": probe.get("H_excess"),
                 "resid": probe.get("max_abs_S"),
                 "iters": probe.get("iters")})

rows.sort(key=lambda r: -r["gap"])
hit = [r for r in rows if r["gap"] > GAP]
print(f"扫描 {scanned} 例（已排除带 c_H/c_OH/c_pH 条件的用例）")
print(f"①无固相 ∧ ②在场 Ksp-OH 阳离子 ∧ ③无分子态强酸 ∧ ④电荷不平衡 的用例："
      f"**{len(rows)} 例**")
print(f"其中 ⑤gap > {GAP} 的：**{len(hit)} 例** ⟹ 这是候选修复面\n")

print(f"{'用例':34s} {'pH_eng':>8} {'pH_chg':>8} {'gap':>7} {'net':>12} "
      f"{'iters':>6} 阳离子")
for r in rows[:35]:
    mark = "★" if r["gap"] > GAP else " "
    print(f"{mark}{r['case'][:33]:33s} {r['pH_eng']:8.3f} {r['pH_charge']:8.3f} "
          f"{r['gap']:7.3f} {r['net']:12.3e} {str(r['iters']):>6} "
          f"{','.join(r['cats'])[:20]}")
if len(rows) > 35:
    print(f"  … 另 {len(rows) - 35} 例")

f31 = [r for r in rows if str(r["case"]).startswith("F31")]
print(f"\nF31: {'命中 ✓ gap=%.3f' % f31[0]['gap'] if f31 else '未命中 ✗'}")

print("\n=== 判读 ===")
print(f"  修复面 = {len(hit)} 例。**F31 必须在其中**，且总数应小到可逐例核对。")
print("  若远大于 20 例 ⟹ 判据不够窄，或说明'读错悬崖'在库里很常见，")
print("  则应改思路为'不让走步进入该态'而非'事后按合取条件切换分支'。")

out = os.path.join(ROOT, "logs", "narrowband.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"gap_thr": GAP, "scanned": scanned, "n_band": len(rows),
               "n_hit": len(hit), "rows": rows}, f, ensure_ascii=False,
              indent=1, default=str)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
