# -*- coding: utf-8 -*-
"""第 203 轮 · F31 严格化学正解（只用引擎自己的数据，独立于引擎代码）。

投料：1 mol GaCl3 + 3 mol NaOH, V = 1 L, 298.15 K。
未知量：pH。Ga 总量守恒 1 M；Na 3 M；Cl 3 M。

平衡关系（全部取自 chemkit/data，数值已与文献核对通过）：
    logβ₁ = 11.4   [Ga(OH)]^{2+} = β₁·[Ga³⁺]·[OH⁻]
    logβ₄ = 37.6   [Ga(OH)₄]⁻   = β₄·[Ga³⁺]·[OH⁻]⁴
    logKsp = -35.1 Ga(OH)₃(s) 饱和时 [Ga³⁺] = 10^logKsp / [OH⁻]³
电荷：3[Ga³⁺] + 2[GaOH²⁺] + [H⁺] + [Na⁺]
      = [Ga(OH)₄⁻] + [OH⁻] + [Cl⁻]
Na⁺ = Cl⁻ = 3 ⟹ 3[Ga³⁺] + 2[GaOH²⁺] + [H⁺] = [Ga(OH)₄⁻] + [OH⁻]

判据：**先判固相是否饱和**（溶液侧 Q vs Ksp），再解电荷条件。
用法： python tools/f31_exact.py [pKw，默认 14.0]
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

PKW = float(sys.argv[1]) if len(sys.argv) > 1 else 14.0
LB1, LB4, LKSP = 11.4, 37.6, -35.1
GA_TOT, NA, CL = 1.0, 3.0, 3.0


def state(poh):
    """给定 pOH，返回溶液侧各量与饱和判据。"""
    oh = 10.0 ** (-poh)
    h = 10.0 ** (-(PKW - poh))
    ga3_sat = 10.0 ** LKSP / oh ** 3          # 饱和线要求的游离 Ga³⁺
    ga3_free = GA_TOT                          # 无固相时全部在溶液
    # 无固相假设下，用总量守恒反解 [Ga³⁺]：1 = G·(1 + β₁·oh + β₄·oh⁴)
    den = 1.0 + 10.0 ** LB1 * oh + 10.0 ** LB4 * oh ** 4
    ga3_nosol = GA_TOT / den
    ga1_nosol = 10.0 ** LB1 * ga3_nosol * oh
    ga4_nosol = 10.0 ** LB4 * ga3_nosol * oh ** 4
    q = ga3_nosol * oh ** 3                    # 溶液侧离子积
    ksp = 10.0 ** LKSP
    return {"oh": oh, "h": h, "ga3_sat": ga3_sat, "q": q, "ksp": ksp,
            "ga3": ga3_nosol, "ga1": ga1_nosol, "ga4": ga4_nosol,
            "sat": q > ksp}


print(f"=== F31 严格解（pKw = {PKW}）===")
print(f"β₁=10^{LB1}  β₄=10^{LB4}  Ksp=10^{LKSP}  "
      f"Ga(OH)₃(s)+OH⁻→[Ga(OH)₄]⁻ logK={LKSP + LB4:+.2f}")

print("\n--- ① 无固相假设：溶液能否容纳全部 1 M Ga？---")
print(f"  {'pH':>6} {'[Ga3+]':>11} {'[GaOH2+]':>11} {'[Ga(OH)4-]':>12} "
      f"{'Q/Ksp':>10}  饱和?")
for poh in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0):
    st = state(poh)
    print(f"  {PKW - poh:6.2f} {st['ga3']:11.3g} {st['ga1']:11.3g} "
          f"{st['ga4']:12.3g} {st['q'] / st['ksp']:10.3g}  "
          f"{'★过饱和' if st['sat'] else '未饱和'}")

print("\n--- ② 溶液恰好饱和的 pH（Q = Ksp 的边界）---")
lo, hi = 0.0, PKW
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if state(mid)["sat"]:
        lo = mid
    else:
        hi = mid
poh_sat = lo
st = state(poh_sat)
print(f"  饱和边界 pH = {PKW - poh_sat:.3f}")
print(f"  该点 [Ga³⁺] = {st['ga3']:.6g} M   [Ga(OH)₄⁻] = {st['ga4']:.6g} M")
print(f"  溶解总 Ga = {st['ga3'] + st['ga1'] + st['ga4']:.6g} M")

print("\n--- ③ 固相存在时的严格电荷解（[Ga³⁺] 由 Ksp 钉住）---")
best = None
poh = 0.0
while poh <= PKW:
    oh = 10.0 ** (-poh)
    h = 10.0 ** (-(PKW - poh))
    g3 = 10.0 ** LKSP / oh ** 3
    g1 = 10.0 ** LB1 * g3 * oh
    g4 = 10.0 ** LB4 * g3 * oh ** 4
    chg = 3 * g3 + 2 * g1 + h - g4 - oh      # Na⁺−Cl⁻ 抵消
    if best is None or abs(chg) < abs(best[1]):
        best = (poh, chg, g3, g1, g4, h, oh)
    poh += 0.0005
poh, chg, g3, g1, g4, h, oh = best
sol = g3 + g1 + g4
print(f"  电荷平衡点: pH = {PKW - poh:.3f}   残余电荷 = {chg:+.3e} M")
print(f"    [Ga³⁺]      = {g3:.6g} M")
print(f"    [Ga(OH)²⁺]  = {g1:.6g} M")
print(f"    [Ga(OH)₄⁻]  = {g4:.6g} M")
print(f"    [H⁺]={h:.3g}  [OH⁻]={oh:.3g}")
print(f"    溶液总 Ga  = {sol:.6g} M")
print(f"    Ga(OH)₃(s) 析出 = {GA_TOT - sol:.6f} mol   (标准要求 ≥ 0.9)")

print("\n=== ④ 结论 ===")
print(f"  严格化学解：pH ≈ {PKW - poh:.2f}，Ga(OH)₃(s) 析出 ≈ "
      f"{GA_TOT - sol:.3f} mol")
print(f"  引擎实测  ：pH = 1.62，Ga(OH)₃(s) 析出 = 0.000 mol，"
      f"余额为 [Ga(OH)₄]⁻ 0.743 + Ga³⁺ 0.253")
print(f"  ⟹ 偏差：pH 差 {abs((PKW - poh) - 1.62):.2f} 个 pH 单位；"
      f"固相差 {GA_TOT - sol:.3f} mol。")
print("  ⟹ 数据经核对与文献一致（β₄ 37.6 / pKsp 35.1 / 两性 logK +2.5）")
print("     ⟹ 矛盾**只能**归给引擎走步，不是数据。")
