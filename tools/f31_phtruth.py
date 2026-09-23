# -*- coding: utf-8 -*-
"""第 203 轮 · F31 终态 pH 的"电荷自洽真值"（判定 estimate_pH 的错幅）。

终态账本：{Ga³⁺ 0.253447, [Ga(OH)₄]⁻ 0.743147, [Ga(OH)]²⁺ 0.003406,
           Na⁺ 3, Cl⁻ 3}
Na⁺ 与 Cl⁻ 相消，故电荷条件给出游离 OH⁻/H⁺ 的唯一取值：
    Σz·n(除 H⁺/OH⁻) = +0.0172 ⟹ [OH⁻] ≈ 0.0172 M ⟹ pH ≈ 12.24
引擎自报 1.62（差 ~10.6 个 pH 单位）。

本脚本把该账本在 **He 轴上**逐点求 pH（沿用引擎自己的 estimate_pH），
并给出"电荷自洽 pH"（由账本直接决定，与 He 无关），两条线一起看，
即可看出引擎的 pH 是被 He 单变量决定的、而账本对该 He 并不自洽。

用法： python tools/f31_phtruth.py
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
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402
from chemkit.speciation import pKw_of                       # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith("F31 ")][0]
probe = {}
r = eng.judge([{"name": k, "mol": m} for k, m in c["subs"]],
              c.get("cond") or {"V_L": V}, T, _probe=probe)
led = dict(probe.get("ledger") or {})
led["H_2O"] = 55.6

print("=== 终态账本 ===")
net = 0.0
for sp, m in sorted(led.items(), key=lambda kv: -abs(kv[1])):
    if sp == "H_2O" or m <= 0:
        continue
    z = spec._charge_cached(sp)
    net += z * m
    print(f"  {sp:22s} {m:12.6g}  z={z:+d}  z·n={z * m:+12.6g}")
print(f"  {'Σz·n（除 H⁺/OH⁻）':22s} {'':12s}      {net:+12.6g} mol")

pKw = pKw_of(T_K)
print(f"\n  pKw(298.15) = {pKw}")
if net > 0:
    oh = net / V
    poh = -__import__("math").log10(oh)
    print(f"  净正电荷 ⟹ 必须由 OH⁻ 平衡： [OH⁻] = {oh:.6g} M")
    print(f"  ⟹ 电荷自洽 pH = {pKw - poh:.3f}  (pOH {poh:.3f})")
    truth = pKw - poh
else:
    h = -net / V
    truth = -__import__("math").log10(h)
    print(f"  净负电荷 ⟹ 必须由 H⁺ 平衡： [H⁺] = {h:.6g} M")
    print(f"  ⟹ 电荷自洽 pH = {truth:.3f}")

print("\n=== 引擎 estimate_pH 沿 He 轴（同一账本）===")
print("  注意：账本固定，只有 He 变——若 estimate_pH 同时依赖两者，")
print("  它会给出与账本电荷条件矛盾的答案。")
for he in (-3.0, -1.0, -0.5, -0.1, -0.05, -0.024007, -0.02, -0.01,
           -0.001, 0.0, 0.001, 0.1, 1.0):
    try:
        ph = eng.estimate_pH(led, he, V, T, T_K)
        print(f"  He={he:+11.6g}  estimate_pH = {ph:8.4f}")
    except Exception as exc:                                # noqa: BLE001
        print(f"  He={he:+11.6g}  异常 {type(exc).__name__}: {exc}")

print(f"\n=== 判定 ===")
print(f"  电荷自洽 pH = {truth:.3f}")
print(f"  引擎自报   = {probe.get('pH')}  (He={probe.get('H_excess')})")
print(f"  走步口径 He=-0.024007 对应 estimate_pH = "
      f"{eng.estimate_pH(led, -0.024007, V, T, T_K):.3f}")
print(f"\n  ⟹ 账本自身要求 pH≈{truth:.2f}（强碱侧）；引擎按 (ledger, He) 给出")
print("     He<0 时 ≈12.4（方向对但幅度错），He≥0 时 ≈1.6（方向都错）。")
print("     **同一个账本，pH 随 He 在 1.6 与 12.4 间跳变 10.8 个单位**")
print("     ⟹ `estimate_pH` 不是账本的自洽解，而是 (ledger, He) 的启发式。")

print("\n=== 若把 pH 钉在电荷自洽值，S 会变成多少？===")
from chemkit.templates import enumerate_candidates           # noqa: E402
cands = enumerate_candidates(led, 0.0, truth, V, T_K, T, True)
print(f"  {'S(真pH)':>10} {'kind':>10}  反应")
rows = []
for cd in cands:
    ps = cd.pres_specs
    if not (all(led.get(x, 0.0) > eng.X_MIN for x in ps[0])
            and all(led.get(x, 0.0) > eng.X_MIN for x in ps[1])):
        continue
    try:
        S = eng.S_of(cd, led, V, truth, T_K, T, frozenset(), eng.P_EXT_KPA,
                     True, {})
    except Exception:                                       # noqa: BLE001
        continue
    if abs(S) > 0.5:
        rows.append((S, cd.kind, f"{cd.r} -> {cd.pr}"))
rows.sort(key=lambda t: -abs(t[0]))
for S, kind, eq in rows[:10]:
    print(f"  {S:+10.3f} {kind:>10}  {eq[:76]}")
print("\n  （若 Ga(OH)₃ 沉淀候选在此 pH 上 S 为强正 ⟹ 真平衡就是'析出固相'，")
print("    引擎的终态在化学上是错的，且错因是 pH 而非沉淀机制。）")
