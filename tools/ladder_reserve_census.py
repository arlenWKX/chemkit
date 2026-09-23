# -*- coding: utf-8 -*-
"""第 276 轮 · **OH 梯配离子的滴定储备普查**：现行储备是"一步脱光"还是"逐级"。

## 背景

`_buffer_titration` 的**配离子碱储备**只读 `beta_pka:` 派生条目，而那一族是
**一步脱光到自由中心**：

    [M(OH)_k] + k·H⁺ -> M^{z+} + k·H₂O     （k = 梯级、z = 中心电荷）

于是滴定器**无法在中途停在氢氧化物固相上**。实测
`H45 Na[Al(OH)4]+HCl 半量`（`tools/ph_path.py`）：

    传入 He = +0.37472143（0.5 mol HCl 减去被吸收的）
    _buffer_titration: tit=None  He_res=0
    虚拟账本变化: {'[Al(OH)_4]^-': -0.09368, 'Al^{3+}': +0.09368}
    ⟹ 0.3747 / 0.09368 = 3.999... ⟹ **按 4 个 H⁺/铝酸根吸收**

而化学上第一步是 **1 个 H⁺**：

    [Al(OH)₄]⁻ + H⁺ -> Al(OH)₃(s) + H₂O
    ⟹ 0.3747 mol H⁺ 把 0.3747 mol 铝酸根变固相 ⟹ 铝酸根 0.5 / Al(OH)₃ 0.5
    ⟹ 与手算（电荷平衡 + K = β₄·Ksp）**逐位一致**

本工具列出每个 OH 梯配离子的：现行储备的 νH⁺ 与"终止于固相的逐级"νH⁺、
两种 logK，以及库内是否已有对应的 `ksp_beta`（溶回）条目。

**数据来源（全部库内，无新数据）**：
    逐级 `[M(OH)_k] + (k−z)·H⁺ -> M(OH)_z(s) + (k−z)·H₂O`
    logK = (k−z)·pKw + pKsp − logβ_k
校验 `[Al(OH)_4]^-`：1×14 + 33 − 34.5 = **12.5**，
即 `Al(OH)₃(s) + OH⁻ -> [Al(OH)₄]⁻` 的 logK = 14 − 12.5 = **+1.5** ✓
（与该 beta 条目自己的 `calibrated` 注记一致）。

用法：python tools/ladder_reserve_census.py
"""

from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.candidates import PKW_298                       # noqa: E402
from chemkit.core import charge_of                           # noqa: E402
from chemkit.data import load_tables                         # noqa: E402


def main() -> int:
    T = load_tables()
    ksp = {}
    for e in T.ksp:
        c, a = e["pair"]
        if a == "OH^-" and c not in ksp:
            ksp[c] = e
    print(f"{'配离子':<16}{'k':>3}{'z':>3}{'现行 νH+':>10}{'逐级 νH+':>10}"
          f"{'逐级 logK':>11}{'溶回 logK':>11}  固相")
    n = 0
    for b in T.beta:
        if b["ligand"] != "OH^-" or b.get("m", 1) != 1:
            continue
        k = b.get("nu", 1)
        z = charge_of(b["center"])
        if k <= z:
            continue
        e = ksp.get(b["center"])
        if e is None:
            continue
        nd = k - z
        lg = nd * PKW_298 + e["pKsp"] - b["logb"]
        n += 1
        print(f"{b['complex']:<16}{k:>3}{z:>3}{k:>10}{nd:>10}"
              f"{lg:>11.3f}{PKW_298 - lg:>11.3f}  {e['solid']}")
    print(f"\n共 {n} 个 OH 梯配离子的储备应当用**逐级（终止于固相）**"
          f"而不是\"一步脱光\"。")
    print("现行 `beta_pka:` 储备一律取 νH⁺ = k（见 speciation._buffer_titration）。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
