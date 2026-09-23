# -*- coding: utf-8 -*-
"""第 212 轮 · 查 `estimate_state` 在 F31 终态上到底走了哪条分支。

第 209/211 轮事实：F31 终态账本上 `estimate_pH` 返回 **1.6197**（酸侧），
而 `acidbase.charge_pH` 返回 **12.3803**。
第 197 轮已在 `speciation.py` L905-1010 落地"退化区精确质子条件接管"：
```
_deg = (h_c > 0.0 and o_c > 0.0 and 0.1*o_c <= h_c <= 10.0*o_c)
if _deg or pin_mode:  ... charge_pH ...
```
若 F31 的 `h_c`/`o_c` 真的都在 0.024 量级（弱酸式给 0.0240、碱侧 −He=0.024），
则**比值恰为 1.0 ⟹ `_deg` 为真 ⟹ 本该用 charge_pH**。可是引擎给的是 1.6197。
**本脚本查清为什么没走那条路**：打印 h_c/o_c、各分支 tag、以及
`_deg` 三个子条件的真假。

用法： python tools/gatecheck.py
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

import chemkit.speciation as spec                           # noqa: E402
import chemkit.acidbase as ab                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
LED = {"Cl^-": 3.0, "Ga^{3+}": 0.253447, "Na^+": 3.0,
       "[Ga(OH)]^{2+}": 0.00340625, "[Ga(OH)_4]^-": 0.743147, "H_2O": 55.6}

print("=== 分支审计（PH_SRC / PH_TAGS 官方钩子）===")
for he in (-0.024, -0.024007, -0.03, -0.02, -0.01):
    spec.PH_TAGS = []
    spec.PH_SRC = []
    ph = spec.estimate_pH(dict(LED), he, V, T, T_K)
    tags = [t for t in (spec.PH_TAGS or []) if not t.startswith("__")]
    src = [(sp, k, val) for sp, k, val in (spec.PH_SRC or [])]
    spec.PH_TAGS = spec.PH_SRC = None
    print(f"\n  He={he:+.6f}  estimate_pH = {ph:.4f}  tags={tags}")
    for sp, k, val in src:
        print(f"      {str(sp):24s} {str(k):28s} {val}")

print("\n=== 逐项复算 h_c / o_c 的候选来源 ===")
# 酸侧：Ga3+ 弱酸式
Ka = 10.0 ** (11.4 - 14.0)
c = LED["Ga^{3+}"] / V
h_weak = (-Ka + math.sqrt(Ka * Ka + 4 * Ka * c)) / 2
print(f"  Ga³⁺ 弱酸式        : Ka={Ka:.6g}  c={c:.6g}  -> h_c={h_weak:.6g} "
      f"(pH {-math.log10(h_weak):.4f})")
# 共轭碱（若已登记）
cb = LED.get("[Ga(OH)]^{2+}", 0.0) / V
if cb > 0:
    h_pair = Ka * c / (Ka + cb)
    print(f"  含共轭碱式         : c_base={cb:.6g} -> h_c={h_pair:.6g} "
          f"(pH {-math.log10(h_pair):.4f})")
# 碱侧：残余强碱
print(f"  碱侧 −He/V         : o_c={-he:.6g} (pH {14.0 + math.log10(-he):.4f})")
print(f"\n  ⟹ h_c(弱酸式)={h_weak:.6g} 与 o_c={-he:.6g} 的比值 = "
      f"{h_weak / (-he):.4f}")
print(f"     `_deg` 判据 0.1*o_c <= h_c <= 10*o_c : "
      f"{0.1 * (-he) <= h_weak <= 10 * (-he)}")

print("\n=== charge_pH 直接调用 ===")
for he in (-0.024, -0.024007):
    try:
        print(f"  He={he:+.6f} -> charge_pH = "
              f"{ab.charge_pH(dict(LED), V, T, T_K)}")
    except Exception as exc:                                # noqa: BLE001
        print(f"  He={he:+.6f} -> 异常 {type(exc).__name__}: {exc}")

print("\n=== 判读 ===")
print("  · 若 tags 显示走了 '电荷平衡精确解' 却仍得 1.62 ⟹ 那条路的入口")
print("    条件与我的复算不一致（例如 h_c/o_c 的来源不是我以为的那两个）。")
print("  · 若 tags 显示 '酸侧max' ⟹ `_deg` 为假 ⟹ 打印的 h_c/o_c 与真实不符，")
print("    需要从 PH_SRC 读出**真实**的 h_c/o_c 来源（上面已打印）。")
