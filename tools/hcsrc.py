# -*- coding: utf-8 -*-
"""第 207 轮 · 抓 F31 悬崖处 h_c/o_c 的**来源**（PH_SRC 审计钩子）。

`tools/cliffbranch.py` 实测（同一族账本）：
  He = -0.031250 -> pH 12.4949 = pKw + log10(0.031250)  ✓ 碱侧
  He = -0.015625 -> pH  1.6192                          ✗ 酸侧
  He = -0.027344 -> pH 12.4369 = pKw + log10(0.027344)  ✓
  He = -0.010000 -> pH  1.6228                          ✗
⟹ 分界在 |He| ≈ 0.031（不是常数阈，像是两来源比大小）。
且 `pH 1.62 ⟺ h_c = 0.024`，而 0.024 恰好是 `|He|` 在 Walk 终态的量级
（**疑似"同一量的两种记法"竞争**，D13；或弱酸路线的值）。

`speciation.PH_SRC` 是官方审计钩子（生产恒 None）：记录 (物种, 贡献类型, 值)。

用法： python tools/hcsrc.py
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

print("=== 单看 Ga3+ 的 Ka（beta 路线的 h_c 潜力）===")
from chemkit.data import load_tables as _lt                  # noqa: E402
_first_k = {b["center"]: 10.0 ** (b["logb"] - 14.0) for b in T.beta
            if b["ligand"] == "OH^-" and b.get("nu") == 1 and b.get("m", 1) == 1}
print(f"  _first_k = {_first_k}")
for cat, ka in _first_k.items():
    for c in (0.2534, 0.24):
        h = (-ka + (ka * ka + 4 * ka * c) ** 0.5) / 2
        print(f"  弱酸二次式 {cat} c={c}: Ka={ka:.4g} -> h={h:.6g} "
              f"(pH {(-__import__('math').log10(h)):.4f})")
print("  ⟹ 若 h_c 由该式给出，pH 应为 ~2.1，**不是 1.62** ⟹ 另有来源。")

print("\n=== 逐档抓 PH_SRC（看 h_c/o_c 由谁贡献）===")
for ga3, he in ((0.25781, -0.031250), (0.25684, -0.027344),
                (0.25391, -0.015625), (0.25000, -0.010000),
                (0.25345, -0.024007)):
    led = {"Ga^{3+}": ga3, "Cl^-": 3.0, "Na^+": 3.0,
           "[Ga(OH)_4]^-": 1.0 - ga3, "H_2O": 55.6}
    spec.PH_TAGS = []
    spec.PH_SRC = []
    try:
        ph = spec.estimate_pH(led, he, V, T, T_K)
    except Exception as exc:                                # noqa: BLE001
        ph = float("nan")
        print(f"  He={he:+.6f} 异常 {type(exc).__name__}: {exc}")
        spec.PH_TAGS = spec.PH_SRC = None
        continue
    tags = [t for t in (spec.PH_TAGS or []) if not t.startswith("__")]
    src = list(spec.PH_SRC or [])
    spec.PH_TAGS = spec.PH_SRC = None
    print(f"\n  He={he:+.6f}  pH={ph:.4f}  tag={tags}")
    for sp, kind, val in src:
        if str(sp).startswith("__"):
            print(f"      [自由侧] {sp} = {val}")
        else:
            print(f"      {sp:22s} {kind:10s} -> {val:.6g}")
    if not src:
        print("      （无来源记录）")
