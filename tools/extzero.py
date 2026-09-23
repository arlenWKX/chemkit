# -*- coding: utf-8 -*-
"""第 209 轮 · 直接测 `solve_extent`：终态上沉淀腿为何给出 x*=0？

`tools/f31trace.py` 抓到的轨迹（关键三行）：
```
[pick] Ga^{3+} + 3H2O -> Ga(OH)_3 + 3H^+   1   28.32
  [ext] 0.0            pH 12.38   He -0.0241     <- 注意：是 12.38 不是 1.62
  [micro] **零推进（S 强而 x*≈0）** x_max=0.256 ext=0 S=+28.320
```
即：**在最终迭代里，候选被正确挑中、pH 也正确（12.38）、S 强正（+28.3），
但 `solve_extent` 给出 x* = 0** ⟹ 走步无步可走 ⟹ `no-cands`。
病根因此**不在候选筛选、也不在 pH 分支**，而在 `solve_extent` 的求根。

本脚本在终态账本上直接调 `solve_extent`，并打印 f(x) 在若干 x 上的取值，
看**括号内 S 的符号分布**（是否全程为正 → 无过零点；或某处跳变）。

用法： python tools/extzero.py
"""
import io
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
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15

# trace 里的终态账本（逐字抄自 [micro-led]）
LED = {"Cl^-": 3.0, "Ga^{3+}": 0.253447, "Na^+": 3.0,
       "[Ga(OH)]^{2+}": 0.00340625, "[Ga(OH)_4]^-": 0.743147, "H_2O": 55.6}
HE = -0.024
PH = 12.38

print("=== 候选枚举（终态）===")
cands = enumerate_candidates(LED, HE, PH, V, T_K, T, True)
target = None
for c in cands:
    if c.kind == "precip" and "Ga(OH)_3" in c.pr:
        target = c
        break
if target is None:
    print("  未找到 Ga(OH)_3 沉淀候选")
    sys.exit(1)
print(f"  kind={target.kind}  r={target.r}  pr={target.pr}")

print("\n=== S_of 随 x 的变化（手工按 _exec 语义推账本）===")


def step(led, c, d, x):
    out = dict(led)
    he = 0.0
    rr = c.r if d > 0 else c.pr
    pp = c.pr if d > 0 else c.r
    for s_, nu in rr.items():
        if s_ == eng.H_ION:
            he -= nu * x
        elif s_ != eng.WATER:
            out[s_] = out.get(s_, 0.0) - nu * x
    for s_, nu in pp.items():
        if s_ == eng.H_ION:
            he += nu * x
        elif s_ != eng.WATER:
            out[s_] = out.get(s_, 0.0) + nu * x
    return out, he


x_max = min(LED.get(s_, 0.0) / nu for s_, nu in target.r.items()
            if s_ not in (eng.WATER, eng.H_ION))
print(f"  计量上限 x_max = {x_max:.6g}")
print(f"  碱储备 = -He = {-HE:.6g} mol；每 mol 反应耗 {target.pr.get(eng.H_ION,0)} mol H⁺ "
      f"⟹ 预算上限 = {-HE / target.pr.get(eng.H_ION, 1):.6g}")
print(f"\n  {'x':>10} {'pH':>9} {'He':>11} {'S(x)':>12}")
for x in (0.0, 1e-6, 1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.256):
    led_x, dhe = step(LED, target, 1, x)
    He_x = HE + dhe
    try:
        ph_x = eng.estimate_pH(led_x, He_x, V, T, T_K)
        S_x = eng.S_of(target, led_x, V, ph_x, T_K, T, frozenset(),
                       eng.P_EXT_KPA, True, {})
        print(f"  {x:10.6g} {ph_x:9.3f} {He_x:+11.6f} {S_x:+12.3f}")
    except Exception as exc:                                # noqa: BLE001
        print(f"  {x:10.6g} 异常 {type(exc).__name__}: {exc}")

print("\n=== 引擎自己的 solve_extent ===")
for he_use in (HE, 0.0):
    for ph_use in (PH, 1.62):
        try:
            ext, xm = eng.solve_extent(target, 1, LED, he_use, V, T_K, T,
                                       frozenset())
            print(f"  He={he_use:+.4g} pH(引擎内部自算) -> "
                  f"ext={ext:.6g}  x_max={xm:.6g}")
            break
        except Exception as exc:                            # noqa: BLE001
            print(f"  He={he_use:+.4g} 异常 {type(exc).__name__}: {exc}")

print("\n=== 判读 ===")
print("  · 若 S(x) 在 (0, x_max] 全程为正 ⟹ `solve_extent` 无过零点、")
print("    按契约返回 '全程为正则取计量上限' —— 但实测返回 0 ⟹ **契约被违反**，")
print("    值得查 `solve_extent` 的早退分支（如 H⁺ 直读/微步止损）。")
print("  · 若 S(x) 穿零一次但根≈0 ⟹ 说明 S 对 x 的单调性被破坏。")
