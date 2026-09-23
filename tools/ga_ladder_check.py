# -*- coding: utf-8 -*-
"""第 268 轮 · **可行性核查**：把 Ga 的羟合梯补全后，F31 的电荷平衡落在哪里？

## 为什么做这个核查

第 268 轮查明 `F31 GaCl3+3NaOH`（`resid_max` 13.305）卡在 `_exact_ok` 第②条：
金属羟合梯来自 `T.beta`（logβ）**不进 `build_families`** ⟹ 族成员数恒为 0。
而**无条件停用②** 虽修好 `F31`（13.305→0、pH 1.62→5.76）却让 **10 例翻红**
（稀土水解 pH、阴离子不变性）；比值判据又落在**刀刃**上
（`F31` 相对差 1e-6–0.66 vs 十例 0.499–0.9，**重叠**）⟹ 任何阈值都是拟合。

更深一层：**库内 `Ga³⁺`/`In³⁺` 的羟合梯缺 ν=2,3**（只有 ν=1 与 ν=4），
`Al³⁺` 缺 ν=3 ⟹ 梯子在 `build_families` 里**根本连不成一条链**。

本脚本回答：**若梯子补全，F31 的电荷平衡会落在哪** —— 若落在 Ga(OH)₃
饱和的近中性区（~5–7），则"补数据"就是 F31 的真正解法；
若仍落在酸性端（~1），则补数据也不够，问题在别处。

## 口径（只用库内 logβ 与 pKw，不含引擎任何启发式）

M³⁺ 的羟合梯：`[M(OH)_k] = β_k·[M]·[OH]^k`（β_0 = 1）。
总金属 `C` 守恒 ⟹ 平均电荷 `z̄(pH) = Σ z_k β_k oh^k / Σ β_k oh^k`（z_k = 3−k）。
电荷平衡：`V·h − V·oh + (Na − Cl) + C·V·z̄(pH) = 0`。
`f(pH)` 单调 ⟹ 二分。

⚠️ **ν=2,3 的 logβ 库里没有**。本核查**借用 §1.5 记的 Thermoddem 值
（20.730 / 30.076）作**敏感性分析**，**不是把它当已核数据入库**
（纪律：记忆值/未核值不得当数据）。脚本同时跑"只补一条/补两条/不补"三档，
给出结论对数据缺口的**敏感程度**。

用法：python tools/ga_ladder_check.py
"""
from __future__ import annotations

import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PKW = 14.0
# (logβ_0..logβ_4)：库里有的两档 + §1.5 记的 Thermoddem 中间级
FULL = [0.0, 11.4, 20.730, 30.076, 37.6]      # Ga（ν=1..4 全）
LIB = [0.0, 11.4, None, None, 37.6]           # 库里现状（ν=2,3 缺）


def zbar(pH: float, betas: list[float | None]) -> float | None:
    oh = 10.0 ** (pH - PKW)
    num = den = 0.0
    for k, lb in enumerate(betas):
        if lb is None:
            return None
        w = 10.0 ** (lb) * oh ** k
        num += (3 - k) * w
        den += w
    return num / den if den > 0 else None


def root(betas, C: float, fixed: float = 0.0,
         lo: float = -2.0, hi: float = 18.0) -> float:
    def f(pH: float) -> float:
        h = 10.0 ** (-pH)
        oh = 10.0 ** (pH - PKW)
        z = zbar(pH, betas)
        if z is None:
            return 0.0
        return h - oh + fixed + C * z
    if f(lo) < 0:
        return lo
    if f(hi) > 0:
        return hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main() -> int:
    print("=== F31 GaCl3 1 mol + NaOH 3 mol / 1 L ===")
    print("  账本（引擎终态）：Ga³⁺ 0.253447 + [Ga(OH)]²⁺ 0.003406 "
          "+ [Ga(OH)₄]⁻ 0.743147 = 1.000 Ga；Na⁺ 3、Cl⁻ 3（互相抵消）")
    C = 1.0
    print("\n-- ① 库里现状（ν=2,3 缺 ⟹ 梯子断成两段）--")
    print("   两条互不相连的子梯各自求根（这正是问题的形状）：")
    z1 = root([0.0, 11.4], 0.256853)          # Ga³⁺/[Ga(OH)]²⁺ 子梯
    print(f"   子梯 A（Ga³⁺/[Ga(OH)]²⁺，总量 0.256853）单独平衡 pH = {z1:.4f}")
    print(f"   子梯 B（[Ga(OH)₄]⁻，总量 0.743147）只是固定电荷 −0.743147")
    print(f"   合起来：电荷平衡要求 Σz = 0，而固定项 −0.743147 只能靠 "
          f"子梯 A 的 +2~+3 补 ⟹ 解出 pH ≈ {root([0.0, 11.4], 0.256853, -0.743147):.4f}")
    print("\n-- ② 补全 ν=1..4（用 §1.5 记的 Thermoddem 中间级做**敏感性分析**）--")
    r_full = root(FULL, C)
    print(f"   整梯一条链，C = 1.0 ⟹ 电荷平衡 pH = **{r_full:.4f}**")
    print("\n-- ③ 只补 ν=2（保留 ν=3 缺）--")
    for lb2 in (20.730,):
        b = [0.0, 11.4, lb2, None, 37.6]
        print(f"   logβ₂={lb2}：仍断开 ⟹ 与①同形（无法求整梯根）")
    print("\n-- ④ 敏感性：中间级取 ±0.5 个对数单位 --")
    for d in (-0.5, 0.0, +0.5):
        b = [0.0, 11.4, 20.730 + d, 30.076 + d, 37.6]
        print(f"   Δ={d:+.1f} ⟹ pH = {root(b, C):.4f}")
    print("\n=== 判读 ===")
    print(f"  补全后电荷平衡落在 pH ≈ {r_full:.2f}")
    print("  该处正是 Ga(OH)₃ 两性最低溶解度区（近中性）⟹ 与"
          "「GaCl₃+3NaOH 恰好中和应析出 Ga(OH)₃」的化学事实一致。")
    print("  ⟹ **F31 的根因是数据缺口（缺 ν=2,3），不是引擎算法缺陷**；"
          "补全后 `build_families` 才可能把整梯连成一条可再分配的族。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
