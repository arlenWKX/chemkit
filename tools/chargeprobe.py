# -*- coding: utf-8 -*-
"""第 203 轮 · 精确质子条件解（`acidbase.charge_pH`）在 F31 终态的表现。

为什么关键：PHREEQC 的做法（model.cpp L4602-4667 + `ph_unknown ==
charge_balance_unknown`）是**把 pH 当作未知量、由电荷平衡方程解出来**
（Newton-Raphson，**没有任何酸/碱分支**）。chemkit 的对等物就是
`acidbase.charge_pH`（精确质子条件）。若它在 F31 终态**能解**，则修法
很轻：把"启发式分支"换成/回退到精确解即可。

本脚本直接调 `charge_pH`（含 pinned 6 元组两种形态）测其在终态/起始态
/中间态的可解性与返回值。

用法： python tools/chargeprobe.py
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

import chemkit.engine as eng                                # noqa: E402
import chemkit.acidbase as ab                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15

STATES = {
    "起始态(3eq 碱,全 Ga3+)": {
        "Ga^{3+}": 1.0, "Cl^-": 3.0, "Na^+": 3.0, "H_2O": 55.6},
    "中间(0.5 镓酸根)": {
        "Ga^{3+}": 0.5, "Cl^-": 3.0, "Na^+": 3.0, "[Ga(OH)_4]^-": 0.5,
        "H_2O": 55.6},
    "引擎终态": {
        "Ga^{3+}": 0.253447023, "Cl^-": 3.0, "Na^+": 3.0,
        "[Ga(OH)_4]^-": 0.743146732, "[Ga(OH)]^{2+}": 0.003406245,
        "H_2O": 55.6},
}

print("=== `acidbase.charge_pH` 可用性 ===")
print(f"  签名: {ab.charge_pH.__doc__.splitlines()[0] if ab.charge_pH.__doc__ else '?'}")
import inspect                                              # noqa: E402
print(f"  inspect: {inspect.signature(ab.charge_pH)}")

print("\n=== 逐态调用 ===")
for tag, led in STATES.items():
    print(f"\n--- {tag} ---")
    print(f"    ledger = { {k: round(v, 6) for k, v in led.items() if k != 'H_2O'} }")
    # ① 位置参数（无 pinned）
    try:
        r = ab.charge_pH(led, V, T, T_K)
        print(f"    charge_pH(led,V,T,T_K)              -> {r}")
    except Exception as exc:                                # noqa: BLE001
        print(f"    charge_pH(led,V,T,T_K)              -> {type(exc).__name__}: {exc}")
    # ② pinned=None 关键字
    try:
        r = ab.charge_pH(led, V, T, T_K, pinned=None)
        print(f"    charge_pH(..., pinned=None)         -> {r}")
    except Exception as exc:                                # noqa: BLE001
        print(f"    charge_pH(..., pinned=None)         -> {type(exc).__name__}: {exc}")
    # ③ engine 的封装（走步实际调用的那个）
    for fn in ("closed_pH", "exact_proton_pH"):
        f = getattr(eng, fn, None)
        if f is None:
            print(f"    engine.{fn} [缺]")
            continue
        for he in (0.0, -0.024007, -3.0):
            try:
                r = f(led, he, V, T, T_K)
                print(f"    engine.{fn}(He={he:+.6f}) -> {r}")
            except Exception as exc:                        # noqa: BLE001
                print(f"    engine.{fn}(He={he:+.6f}) -> "
                      f"{type(exc).__name__}: {exc}")

print("\n=== 判读 ===")
print("  · 若 `charge_pH` 在终态返回一个**中性附近**的 pH（≈5.8~12.4）而")
print("    `closed_pH` 返回 None ⟹ 精确解**存在但没被走步用上**（接线问题，")
print("    修法轻：让 estimate_pH 在分支不可信时回退到 charge_pH）。")
print("  · 若 `charge_pH` 也返回 None ⟹ 精确解本身解不出来（缺能力，")
print("    修法重：要先让精确质子条件能处理'强碱+两性金属+羟合梯'这一档）。")
print("  · 对照 PHREEQC：它把 pH 当未知量由电荷平衡方程解（Newton），")
print("    **没有分支** ⟹ 我们的目标是让 estimate_pH 具备同样的'无分支'性质。")
