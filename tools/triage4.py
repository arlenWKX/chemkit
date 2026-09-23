# -*- coding: utf-8 -*-
"""第 205 轮 · 那 3 例是否与 F31 同类？（决定"修复面"到底几例）

`tools/narrowband.py` 用合取条件筛出 109 例，其中 gap>2 的只有 4 例：
`Se33` / `F31` / `D43` / `X12`。若那 3 例与 F31 不同类，则**修复面 = 1 例**
（可精确修）；若同类，则判据要能区分它们。

本脚本对每例打印：账本全量、He、净电荷、引擎 pH、电荷自洽 pH，
以及**账本里有没有强酸/强碱的痕迹**（分子态强酸、游离 H⁺/OH⁻ 量级）。

判据设想：F31 的账本是"**碱被记成配体**"（净电荷 +0.024 来自 Ga 物种的
零头，**没有强酸**）；而 `Se33`（H₂SeO₄ 强酸）、`X12`（醋酸弱酸）账本里
应当有**大量未解离的酸或游离质子**的痕迹。

用法： python tools/triage4.py
"""
import io
import json
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

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases, run_case           # noqa: E402

T = load_tables()
cases = {c["name"]: c for c in load_cases(None)}
PREFIXES = ["Se33", "F31", "D43", "X12"]

for pre in PREFIXES:
    hit = [v for n, v in cases.items() if n.startswith(pre + " ") or n == pre]
    if not hit:
        hit = [v for n, v in cases.items() if n.startswith(pre)]
    if not hit:
        print(f"### {pre}: 未找到")
        continue
    c = hit[0]
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    T_K = float(cond.get("T_K", 298.15))
    probe = {}
    ok = run_case(c, T, verbose=False)
    eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
              {"V_L": V, "T_K": T_K}, T, _probe=probe)
    led = dict(probe.get("ledger") or {})
    net = 0.0
    for sp, m in led.items():
        if sp == "H_2O" or m <= 0 or sp.startswith("__"):
            continue
        net += spec.charge_of(sp) * m
    pKw = spec.pKw_of(T_K)
    if net > 1e-6:
        pH_c = pKw + math.log10(net / V)
    elif net < -1e-6:
        pH_c = -math.log10(-net / V)
    else:
        pH_c = None
    print(f"\n{'=' * 74}")
    print(f"### {c['name']}   套件={'PASS' if ok else 'FAIL'}")
    print(f"  投料 = {c['subs']}")
    print(f"  引擎 pH = {probe.get('pH_solver')}   He(呈现) = "
          f"{probe.get('H_excess')}   max|S| = {probe.get('max_abs_S')}"
          f"   iters = {probe.get('iters')}  exit = {probe.get('exit')}")
    print(f"  净电荷 Σz·n = {net:+.6e}   ⟹ 电荷自洽 pH = "
          f"{'None(已平衡)' if pH_c is None else round(pH_c, 3)}")
    print(f"  账本（>1e-6）:")
    for sp, m in sorted(led.items(), key=lambda kv: -abs(kv[1])):
        if m > 1e-6 and not sp.startswith("__"):
            z = spec.charge_of(sp)
            print(f"    {sp:26s} {m:12.6g}  z={z:+d}  z·n={z * m:+12.5g}")
    # 未解离酸 / 游离质子痕迹
    free_h = 10 ** (-probe.get("pH_solver", 7.0)) * V
    print(f"  游离 [H⁺]·V = {free_h:.3g} mol（引擎 pH 下）")
    for sp in ("HNO_3", "H_2SO_4", "H_2SeO_4", "CH_3COOH", "HC_2H_3O_2"):
        if led.get(sp, 0.0) > 1e-9:
            print(f"    ★ 未解离酸 {sp} = {led[sp]:.6g}")

print("\n=== 判读 ===")
print("  F31 的账本应当**没有**未解离强酸、净电荷只是 Ga 物种的零头；")
print("  若 Se33/X12/D43 账本里有大量未解离酸或游离质子痕迹 ⟹ **不同类**，")
print("  判据可用「账本无未解离强酸 ∧ 净电荷 ≪ 1 mol」把它们分开。")
