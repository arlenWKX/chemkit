# -*- coding: utf-8 -*-
"""第 275 轮 · **含氧酸根碱侧角色普查**：`bases_map` 新收了哪些物种、Kb 从哪来。

## 背景

`estimate_pH` 分支 4 的 `o_c` 只扫 `rolemap`（= `acids_map ∪ bases_map ∪
hyd_map ∪ amph_pH`）。两性金属的**含氧酸根**（`[Al(OH)₄]⁻`、`[Zn(OH)₄]²⁻`
等）三者都不在 ⟹ **没有角色** ⟹ `o_c` 看不见账本里那一大池碱。

第 275 轮把它们的碱式补上，Kb **完全由库内数据推出**（无新数据、无阈值）：

    [M(OH)_k] ⇌ [M(OH)_{k-1}] + OH⁻   ⟹ Kb = β_{k-1}/β_k   （β₀ ≡ 1）
    中间级缺数据（Al³⁺ 恰缺 ν=3）时退回固相：
    [M(OH)_k] ⇌ M(OH)_z(s) + OH⁻      ⟹ Kb = 1/(β_k·Ksp) = 10^(pKsp − logβ_k)

校验：`[Al(OH)_4]^-` 的 β₄ = 10^34.5、pKsp(Al(OH)₃) = 33 ⟹ Kb = 10^(33−34.5)
= 0.0316 ⟹ `Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻` 的 logK = **+1.5**，与该 `beta`
条目自己的 `calibrated` 注记（"整体 logK ≈ +1.5"）**逐位一致**。

用法：python tools/ladder_base_census.py
"""

from __future__ import annotations

import math
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

from chemkit.core import charge_of                        # noqa: E402
from chemkit.data import load_tables                      # noqa: E402
from chemkit.speciation import _pksp                      # noqa: E402

T_K = 298.15


def main() -> int:
    T = load_tables()
    ksp: dict = {}
    for e in T.ksp:
        c, a = e["pair"]
        if a == "OH^-" and c not in ksp:
            ksp[c] = e
    lad: dict = {}
    skipped = []
    for b in T.beta:
        if b["ligand"] != "OH^-" or b.get("m", 1) != 1:
            continue
        cx = b.get("complex")
        if not cx:
            continue
        if charge_of(cx) >= 0:
            continue
        lad[(b["center"], b.get("nu", 1))] = (cx, float(b["logb"]))
    print(f"{'含氧酸根':<16}{'中心':<12}{'Kb':>12}{'logK(溶回)':>12}  Kb 来源")
    n = 0
    for (ctr, k), (cx, lb) in sorted(lad.items()):
        prev = lad.get((ctr, k - 1))
        if prev is not None:
            kb = 10.0 ** (prev[1] - lb)
            src = f"beta 相邻级 nu={k - 1}"
        else:
            e = ksp.get(ctr)
            if e is None:
                skipped.append((cx, "既无相邻级、中心也无 Ksp(OH)"))
                continue
            kb = 10.0 ** (_pksp(e, T_K) - lb)
            src = f"Ksp {e['solid']}"
        if kb <= 0.0:
            skipped.append((cx, "Kb ≤ 0"))
            continue
        n += 1
        print(f"{cx:<16}{ctr:<12}{kb:>12.4g}{-math.log10(kb):>12.3f}  {src}")
    print(f"\n共 {n} 个含氧酸根获得碱侧角色")
    for cx, why in skipped:
        print(f"  [跳过] {cx}：{why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
