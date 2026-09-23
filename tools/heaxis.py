# -*- coding: utf-8 -*-
"""第 203 轮 · `estimate_pH` 的 He 轴连续性普查（修法可行性判定）。

已定影：F31 首步的二分根落在 `He = 0.0` 的**阶跃**上（f31_scan.py）。
修法首选是"让 estimate_pH 在 He=0 处连续"。本脚本回答两个问题：
  ① He=0 是**唯一**的悬崖，还是 He 轴上有一串？（决定修法范围）
  ② 该账本下 `He -> pH` 是否单调？非单调 ⟹ 二分/走步全建立在流沙上。

同时打印每一档的 `closed_pH`（精确解）与 `estimate_pH`（启发式）之差——
差大的档就是启发式不可信、应交给精确解的档。

用法： python tools/heaxis.py [用例前缀，默认 F31]
"""
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
import chemkit.normalize as _norm                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith(PRE + " ")][0]
V, T_K = 1.0, 298.15
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
led0, _he0, _st, _un = _norm.normalize(
    subs, {"V_L": V, "T_K": T_K, "c_H": None, "c_OH": None, "pH": None,
           "p_kpa": eng.P_EXT_KPA}, T, {})

print(f"=== 起始账本上的 He -> pH 曲线（{c['name']}）===")
print(f"  ledger = { {k: round(v, 6) for k, v in led0.items() if k != 'H_2O'} }")
print(f"\n  {'He':>12} {'estimate_pH':>12} {'closed_pH':>11} {'presentation':>12} "
      f"{'差(est-cls)':>12}")
prev = None
prev_he = None
bad = []
for he in (-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.2, -0.1, -0.05, -0.02,
           -0.01, -0.005, -0.002, -0.001, -5e-4, -2e-4, -1e-4, -1e-5, -1e-6,
           0.0,
           +1e-6, +1e-5, +1e-4, +2e-4, +5e-4, +1e-3, +0.002, +0.005,
           +0.01, +0.05, +0.1, +0.5, +1.0):
    try:
        e = eng.estimate_pH(led0, he, V, T, T_K)
    except Exception as exc:                                # noqa: BLE001
        e = f"异常{type(exc).__name__}"
    try:
        cl = eng.closed_pH(led0, he, V, T, T_K)
        cl = "None" if cl is None else cl[0]
    except Exception as exc:                                # noqa: BLE001
        cl = f"异常{type(exc).__name__}"
    try:
        pr = eng.presentation_pH(led0, he, V, T, T_K)
    except Exception as exc:                                # noqa: BLE001
        pr = f"异常{type(exc).__name__}"
    _d = (f"{e - cl:+.3f}" if isinstance(e, float) and isinstance(cl, float)
          else "-")
    print(f"  {he:+12.6g} {e if not isinstance(e, float) else round(e, 4):>12} "
          f"{cl if not isinstance(cl, float) else round(cl, 4):>11} "
          f"{pr if not isinstance(pr, float) else round(pr, 4):>12} {_d:>12}")
    if isinstance(e, float):
        if prev is not None and e < prev - 1.0:
            bad.append((prev_he, he, prev, e))
        prev, prev_he = e, he

print("\n=== ① 非单调下降（pH 应随 He 增大而单调下降）===")
if bad:
    for a, b, pa, pb in bad:
        print(f"  ⚠ He {a:+.6g} -> {b:+.6g} : pH {pa:.4f} -> {pb:.4f} "
              f"(跳 {pb - pa:+.3f})")
else:
    print("  未发现 >1 单位的反常下降")

print("\n=== ② He=0 邻域细扫（确认是阶跃还是尖点）===")
for he in (-1e-3, -1e-4, -1e-5, -1e-6, -1e-9, 0.0, 1e-9, 1e-6, 1e-5, 1e-4):
    try:
        e = eng.estimate_pH(led0, he, V, T, T_K)
        print(f"  He={he:+.1e}  estimate_pH = {e:.6f}")
    except Exception as exc:                                # noqa: BLE001
        print(f"  He={he:+.1e}  异常 {type(exc).__name__}: {exc}")

print("\n=== 判读 ===")
print("  He=0 的 estimate_pH 若与 He→0⁻ 的极限差 >3 个 pH 单位 ⟹ 阶跃坐实，")
print("  且阶跃宽度 <1e-9 ⟹ 不是数值噪声，是**分支判据**的硬切换")
print("  （`H_excess <= 0` 走分支 3 强碱，`> 0` 落分支 4）⟹ 修法必须改判据，")
print("  不能靠容差调参。")
