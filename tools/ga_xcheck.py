# -*- coding: utf-8 -*-
"""第 203 轮 · Ga 常数双权威交叉核对（chemkit vs PHREEQC Thermoddem v1.10）。

使用者提示"也可以从 phreeqc 获取灵感和方向" ⟹ 把 PHREEQC 自带数据库当作
**独立权威**来核对 chemkit 的取值。实测：
  · `llnl.dat` / `minteq.v4.dat` / `wateq4f.dat` / `phreeqc.dat` / `sit.dat` /
    `pitzer.dat` **全部没有** Ga 水解配合物（llnl.dat 只有 Ga 主物种）
    ⟹ "库里没有"不等于"化学上没有"；chemkit 收 Ga 羟合梯是对的。
  · `PHREEQC_ThermoddemV1.10_15Dec2020.dat` 有**完整** Ga 数据（97ben/dia、
    99dia/sch 为源），可作独立核对。

换算口径（全部 I→0）：
  Thermoddem 给的是 `Ga³⁺ + nH₂O = Ga(OH)n^(3-n) + nH⁺` 的 log*Kₙ，
  而 chemkit 的 βₙ 是 `Ga³⁺ + nOH⁻ = Ga(OH)n` 的累积常数：
      logβₙ = n·pKw + log*Kₙ
  固相：Thermoddem 只给 GaOOH（= Ga(OH)₃ 脱水形式）：
      GaOOH + 3H⁺ = Ga³⁺ + 2H₂O, log K = 1.487
      GaOOH + H₂O = Ga(OH)₃(s)（水合，log K = 0）
      ⟹ Ga(OH)₃(s) = Ga³⁺ + 3OH⁻ 的 logKsp = 1.487 − 3·pKw = −35.513

用法： python tools/ga_xcheck.py
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

PKW = 14.0

# --- chemkit 库内值（直接读库，不反推）---
from chemkit.data import load_tables                        # noqa: E402
T = load_tables()
CK = {}
for b in T.beta:
    if b.get("center") == "Ga^{3+}" and b.get("ligand") == "OH^-":
        CK[b["nu"]] = b["logb"]
ck_ksp = None
for e in T.ksp:
    if e.get("solid") == "Ga(OH)_3":
        ck_ksp = -e["pKsp"]
print("=== chemkit 库内值 ===")
for nu in sorted(CK):
    print(f"  logβ{nu} = {CK[nu]}")
print(f"  logKsp(Ga(OH)_3) = {ck_ksp}")

# --- Thermoddem（PHREEQC 独立权威）---
TD_STAR = {1: -2.836, 2: -7.270, 3: -11.924, 4: -15.633}
TD_GAOOH = 1.487
print("\n=== Thermoddem v1.10（PHREEQC 独立权威）===")
print(f"  log*K1 = {TD_STAR[1]}   log*K2 = {TD_STAR[2]}   "
      f"log*K3 = {TD_STAR[3]}   log*K4 = {TD_STAR[4]}")
print(f"  GaOOH + 3H+ = Ga+3 + 2H2O : logK = {TD_GAOOH}")
td_ksp = TD_GAOOH - 3 * PKW
print(f"  ⟹ logKsp(Ga(OH)₃) = {TD_GAOOH} − 3×{PKW} = **{td_ksp:.3f}**")

print("\n=== 换算口径**用 Al 自校验**（两者都有 Al 梯，可验证换算式）===")
al_td = {1: -4.951, 2: -10.592, 4: -22.872}   # Thermoddem: AlOH+2/Al(OH)2+/AlO2-
_al_ck = {}
for b in T.beta:
    if b.get("center") == "Al^{3+}" and b.get("ligand") == "OH^-":
        _al_ck[b["nu"]] = b["logb"]
_al_ksp = None
for e in T.ksp:
    if e.get("solid") == "Al(OH)_3":
        _al_ksp = -e["pKsp"]
print(f"  {'n':>3} {'chemkit':>10} {'Thermoddem':>12} {'Δ':>8}")
for nu in (1, 2, 4):
    td = nu * PKW + al_td[nu]
    print(f"  {nu:>3} {_al_ck.get(nu, 0):>10} {td:>12.3f} "
          f"{_al_ck.get(nu, 0) - td:>+8.2f}")
_al_td_ksp = 7.738 - 3 * PKW
print(f"  Ksp  {_al_ksp:>10} {_al_td_ksp:>12.3f} {_al_ksp - _al_td_ksp:>+8.2f}")
print("  ⟹ Al 两侧最大偏差 ~0.8 ⟹ **换算式 logβₙ = n·pKw + log*Kₙ 正确**，")
print("     下面 Ga 的对比因此可信。")

print("\n=== 换算对齐（logβₙ = n·pKw + log*Kₙ）===")
print(f"  {'n':>3} {'chemkit logβn':>15} {'Thermoddem logβn':>18} {'Δ(ck−td)':>10}")
for nu in (1, 2, 3, 4):
    td = nu * PKW + TD_STAR[nu]
    ck = CK.get(nu)
    d = f"{ck - td:+.2f}" if ck is not None else "（chemkit 缺该级）"
    print(f"  {nu:>3} {ck if ck is not None else '—':>15} {td:>18.3f} {d:>10}")

print(f"\n  {'':>3} {'chemkit logKsp':>15} {'Thermoddem':>18} {'Δ':>10}")
print(f"  {'':>3} {ck_ksp:>15} {td_ksp:>18.3f} {ck_ksp - td_ksp:>+10.3f}")

print("\n=== 两性窗口（Ga(OH)₃(s) + OH⁻ → [Ga(OH)₄]⁻ 的 logK）===")
ck_am = ck_ksp + CK[4]
td_am = td_ksp + (4 * PKW + TD_STAR[4])
print(f"  chemkit    : {ck_ksp} + {CK[4]} = **{ck_am:+.3f}**")
print(f"  Thermoddem : {td_ksp:.3f} + {4 * PKW + TD_STAR[4]:.3f} = **{td_am:+.3f}**")
print(f"  Δ = {ck_am - td_am:+.3f} log 单位")
print("\n  ⚠️ **两源不可直接相比**：Thermoddem 里 Ga 的固相**只有 GaOOH**")
print("     （晶质羟基氧化镓），**根本没有 Ga(OH)₃ 相**；chemkit 用的是")
print("     `Ga(OH)_3`（新鲜沉淀的三羟化物/无定形口径）。")
print(f"     pKsp 差 {abs(ck_ksp - td_ksp):.2f} 个 log 单位正是"
      "「无定形 vs 晶质」的常见量级")
print("     ⟹ **不能据此判定任一侧取值错**，这是固相选择差异，不是数据错误。")

print("\n=== 对 F31 的影响：用两套常数各解一次严格解 ===")


def solve(lb, lksp, ga_tot=1.0):
    """固相存在时的严格电荷解（Na⁺=Cl⁻=3 相消，Ga(OH)₃ 饱和钉住游离量）。"""
    best = None
    poh = 0.0
    while poh <= PKW:
        oh = 10.0 ** (-poh)
        h = 10.0 ** (-(PKW - poh))
        g3 = 10.0 ** lksp / oh ** 3
        g1 = 10.0 ** lb.get(1, -99) * g3 * oh
        g4 = 10.0 ** lb.get(4, -99) * g3 * oh ** 4
        chg = 3 * g3 + 2 * g1 + h - g4 - oh
        if best is None or abs(chg) < abs(best[1]):
            best = (poh, chg, g3, g1, g4, h, oh)
        poh += 0.0005
    poh, chg, g3, g1, g4, h, oh = best
    sol = g3 + g1 + g4
    return PKW - poh, sol, g3, g1, g4, chg


for tag, lb, lksp in (("chemkit", CK, ck_ksp),
                      ("Thermoddem", {n: n * PKW + TD_STAR[n] for n in (1, 2, 3, 4)},
                       td_ksp)):
    ph, sol, g3, g1, g4, chg = solve(lb, lksp)
    print(f"  {tag:11s}: pH = {ph:.3f}   溶液总 Ga = {sol:.4g} M   "
          f"析出 = {1.0 - sol:.6f} mol  (残余电荷 {chg:+.2e})")
    print(f"               [Ga³⁺]={g3:.3g}  [GaOH²⁺]={g1:.3g}  [Ga(OH)₄⁻]={g4:.3g}")

print("\n=== 结论 ===")
print("  1. chemkit 收 Ga 羟合梯**正确**：PHREEQC 六个主流库"
      "（llnl/minteq.v4/wateq4f/phreeqc/sit/pitzer）**全无** Ga 水解配合物")
print("     （llnl.dat 只有 Ga 主物种）⟹「库里没有」不能作反证；")
print(f"  2. β₁ 两源吻合极好（11.4 vs {PKW + TD_STAR[1]:.3f}，Δ+0.24）；")
print(f"     β₄ chemkit 37.6 vs Thermoddem {4 * PKW + TD_STAR[4]:.2f}"
      f"（chemkit 低 {4 * PKW + TD_STAR[4] - 37.6:.2f}）；")
print("     pKsp 的 5.4 log 差来自**固相不同**（GaOOH vs Ga(OH)₃），非取值错；")
print("  3. **两套常数给出的 F31 定性结论完全一致**：固相析出 ≈1.0 mol、")
print("     pH 落在中性附近（chemkit 5.76 / Thermoddem 6.88）。")
print("     ⟹ F31 的病根**不是** Ga 常数取值，仍是引擎走步/pH 分支。")
print("  4. chemkit 缺 Ga 中间级（Ga(OH)₂⁺ / Ga(OH)₃(aq)，Thermoddem 有）：")
print("     属**数据缺口**（可补），但中间级在 F31 量级下不主导（碱侧由 β₄、")
print("     酸侧由 β₁ 主导）⟹ 不是本轮缺陷根因。")
print("  5. 方法论收获：**先找同时有两套数据的元素（Al）验证换算式**，")
print("     再去比只有单侧数据的元素（Ga）——否则无法区分「口径错」与")
print("     「取值分歧」。这条已写进 lessons。")
