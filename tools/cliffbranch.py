# -*- coding: utf-8 -*-
"""第 207 轮 · 定位 F31 的 pH 悬崖发生在哪个分支（用 PH_TAGS 审计钩子）。

`tools/whystuck.py` 实测：走步过程中同一族账本，
  He = -0.031250 -> pH 12.4949   （正确，碱侧）
  He = -0.015625 -> pH  1.6193   （错，酸侧）★跳变发生在这一档
`estimate_state` 的分支结构（speciation.py L669-1020）：
  ① `_buffer_titration` 命中 -> 滴定/Henderson 直返
  ② `_role_free` 且 |He| ≥ 1e-3 -> **直读**（H⁺/OH⁻ 直读）
  ③ 否则落分支 4（各来源取最大）+ `_deg` 退化区内的精确解
本脚本逐 He 打印 (pH, 分支 tag)，看跳变处**换了哪条分支**。

`speciation.PH_TAGS` 是官方审计钩子（生产路径恒 None，零行为影响）。

用法： python tools/cliffbranch.py
"""
import io
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

import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15

# whystuck 实测的走步序列（Ga³⁺/Ga(OH)₄⁻ 互补，Cl⁻=Na⁺=3）
SEQ = [(1.000000, 0.000000), (0.500000, -1.000000), (0.375000, -0.500000),
       (0.312500, -0.250000), (0.281250, -0.125000), (0.265620, -0.062500),
       (0.257810, -0.031250), (0.253910, -0.015625), (0.255860, -0.023438),
       (0.256840, -0.027344), (0.253447, -0.024007), (0.250000, -0.010000),
       (0.260000, -0.040000), (0.240000, -0.005000)]

print(f"{'[Ga3+]':>9} {'He':>11} {'pH':>9}  {'分支 tag':28s} 账本")
for ga3, he in SEQ:
    ga4 = 1.0 - ga3
    led = {"Ga^{3+}": ga3, "Cl^-": 3.0, "Na^+": 3.0, "[Ga(OH)_4]^-": ga4,
           "H_2O": 55.6}
    spec.PH_TAGS = []
    try:
        ph = spec.estimate_pH(led, he, V, T, T_K)
    except Exception as exc:                                # noqa: BLE001
        ph = float("nan")
        print(f"  {ga3:9.5f} {he:+11.6f}   异常 {type(exc).__name__}: {exc}")
        spec.PH_TAGS = None
        continue
    tags = [t for t in (spec.PH_TAGS or []) if not t.startswith("__")]
    spec.PH_TAGS = None
    print(f"  {ga3:9.5f} {he:+11.6f} {ph:9.4f}  {str(tags):28s} "
          f"{{Ga3+ {ga3:.5f}, Ga4- {ga4:.5f}}}")

print("\n=== 关键：|He| = 1e-3 两侧的行为（阈值悬崖）===")
base = {"Cl^-": 3.0, "Na^+": 3.0, "H_2O": 55.6}
for he in (-3e-3, -1.1e-3, -1.0e-3, -0.9e-3, -5e-4, -1e-4, 0.0,
           1e-4, 5e-4, 9e-4, 1.0e-3, 1.1e-3, 3e-3):
    # 固定一个"含镓羟合物、无 pKa 族"的账本（复刻 F31 终态比例）
    led = dict(base)
    led["Ga^{3+}"] = 0.253447
    led["[Ga(OH)_4]^-"] = 0.743147
    spec.PH_TAGS = []
    try:
        ph = spec.estimate_pH(led, he, V, T, T_K)
    except Exception as exc:                                # noqa: BLE001
        ph = float("nan")
    tags = [t for t in (spec.PH_TAGS or []) if not t.startswith("__")]
    spec.PH_TAGS = None
    print(f"  He={he:+10.1e}  pH={ph:9.4f}  {tags}")

print("\n=== 判读 ===")
print("  若 He=-0.031250 走「OH⁻直读」而 He=-0.015625 落「酸侧max」⟹")
print("  **悬崖 = `_role_free` 早退分支与分支 4 之间的切换**：")
print("    · |He| ≥ 1e-3 且账本无 role 物种 ⟹ 直读（对）")
print("    · |He| < 1e-3 ⟹ 落分支 4，而分支 4 在此账本上给出酸侧 1.62（错）")
print("  注意 `_role_free` 由**账本物种**决定，与 He 无关 ⟹ 同一账本下")
print("  He 跨越 ±1e-3 就会换分支，而 F31 的平衡点恰好落在该阈值附近。")
