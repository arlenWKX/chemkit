# -*- coding: utf-8 -*-
"""第 206 轮 · 量化"同一数值由多套算法算出且互相冲突"的实际暴露面。

使用者提问：引擎是否大量存在"同一数值、不同计算方法、互相冲突"？
本脚本只读地回答**"用户实际看到的数与引擎内部搜索用的数是否一致"**，
这是该问题最要紧的一面（内部口径分歧若只在诊断字段里，危害有限；
若 `final_pH` 与求解器自用 pH 不一致，则**引擎报告的答案与它搜索的答案不是同一个**）。

逐例采集（`judge(..., _probe={})`）：
  · `final_pH`      —— 用户看到的答案
  · `probe["pH"]`   —— 呈现口径（同上但经探针）
  · `probe["pH_solver"]` —— **走步实际用于判 S 的 pH**
  · `pH_charge`     —— 账本电荷自洽值（若账本完整则应等于真值）
并统计：
  ① `|final_pH − pH_solver|` 分布（答案与搜索是否脱钩）
  ② `|pH_charge − pH_solver|` 分布（搜索是否用了与账本矛盾的 pH）
  ③ 分歧 >1 / >3 的例数与清单

⚠️ 注意口径：`charge_pH` 类比较**必须排除账本欠定的用例**（引擎把未解离
强酸按中性分子记账 ⟹ Σz·n 看似不平衡，其实账本有缺口）。故本脚本把
"账本里有没有未解离酸/碱储备"一并打印，供判读时剔除假阳性。

用法： python tools/conflict_audit.py [步长，默认 1]
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
from chemkit.testsuit import load_cases                     # noqa: E402

STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 1
T = load_tables()
SOLIDS = frozenset(T.solids)
KSP_CATS = frozenset(e["pair"][0] for e in T.ksp if e["pair"][1] == "OH^-")
# 未解离酸/碱储备（引擎按中性分子记账 ⟹ 会使 Σz·n 看似不平衡）
UNDISSOC = frozenset(("HNO_3", "H_2SO_4", "H_2SeO_4", "HCl", "HBr", "HI",
                      "HClO_4", "CH_3COOH", "HC_2H_3O_2", "HF", "H_2CO_3",
                      "H_2S", "H_3PO_4", "H_2C_2O_4"))

cases = load_cases(None)
rows = []
scanned = 0
for i, c in enumerate(cases):
    if i % STEP:
        continue
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    T_K = float(cond.get("T_K", 298.15))
    probe = {}
    try:
        r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                      {"V_L": V, "T_K": T_K}, T, _probe=probe)
    except Exception:                                       # noqa: BLE001
        continue
    if not probe:
        continue                                    # OVERRIDE 直出，无探针
    scanned += 1
    led = dict(probe.get("ledger") or {})
    ph_final = r.get("final_pH")
    ph_pres = probe.get("pH")
    ph_solver = probe.get("pH_solver")
    if not isinstance(ph_solver, float):
        continue
    net = 0.0
    for sp, m in led.items():
        if sp == "H_2O" or m <= 0 or sp.startswith("__"):
            continue
        try:
            net += spec.charge_of(sp) * m
        except Exception:                                   # noqa: BLE001
            pass
    pKw = spec.pKw_of(T_K)
    if net > 1e-6:
        pH_ch = pKw + math.log10(net / V)
    elif net < -1e-6:
        pH_ch = -math.log10(-net / V)
    else:
        pH_ch = None
    has_solid = any(m > spec.X_MIN and sp in SOLIDS for sp, m in led.items())
    has_cat = any(m > spec.X_MIN and sp in KSP_CATS for sp, m in led.items())
    has_und = any(led.get(sp, 0.0) > 1e-6 for sp in UNDISSOC)
    d_ans = (abs(ph_final - ph_solver)
             if isinstance(ph_final, float) else None)
    d_chg = (abs(pH_ch - ph_solver) if pH_ch is not None else None)
    rows.append({"case": c["name"], "final": ph_final, "pres": ph_pres,
                 "solver": ph_solver, "charge": pH_ch, "d_ans": d_ans,
                 "d_chg": d_chg, "solid": has_solid, "cat": has_cat,
                 "undissoc": has_und, "resid": probe.get("max_abs_S")})

print(f"扫描 {scanned} 例（有探针者）\n")

print("=" * 76)
print("① **用户看到的答案 vs 走步搜索用的 pH**（|final_pH − pH_solver|）")
print("=" * 76)
ds = sorted((r["d_ans"] for r in rows if r["d_ans"] is not None))
if ds:
    n = len(ds)
    print(f"  n = {n}")
    for q, tag in ((0.5, "p50"), (0.9, "p90"), (0.99, "p99")):
        print(f"  |Δ| {tag} = {ds[min(int(n * q), n - 1)]:.4f}")
    print(f"  |Δ| max = {ds[-1]:.4f}")
    print(f"  |Δ| > 0.001 的: {sum(1 for d in ds if d > 0.001)}")
    print(f"  |Δ| > 1.0   的: {sum(1 for d in ds if d > 1.0)}")
    print(f"  |Δ| > 3.0   的: {sum(1 for d in ds if d > 3.0)}")
    bad = sorted((r for r in rows if (r["d_ans"] or 0) > 1.0),
                 key=lambda r: -(r["d_ans"] or 0))
    if bad:
        print(f"\n  脱钩最严重的用例（前 15）:")
        print(f"   {'用例':34s} {'final_pH':>9} {'solver':>8} {'Δ':>7} "
              f"{'max|S|':>9} 固相")
        for r in bad[:15]:
            print(f"   {r['case'][:34]:34s} {r['final']:9.3f} "
                  f"{r['solver']:8.3f} {r['d_ans']:7.3f} "
                  f"{str(r['resid'])[:9]:>9} {r['solid']}")

print()
print("=" * 76)
print("② **搜索用的 pH vs 账本电荷自洽值**（|pH_charge − pH_solver|）")
print("=" * 76)
ds2 = sorted((r["d_chg"] for r in rows if r["d_chg"] is not None))
n_all = sum(1 for r in rows if r["d_chg"] is not None)
# 剔除"账本欠定"（有未解离酸/碱 ⟹ 电荷不平衡是账本缺口，非 pH 错）
clean = [r for r in rows if r["d_chg"] is not None and not r["undissoc"]]
ds2c = sorted(r["d_chg"] for r in clean)
print(f"  有电荷自洽值可比: {n_all} 例")
if ds2:
    print(f"    |Δ| > 1.0: {sum(1 for d in ds2 if d > 1.0)}")
    print(f"    |Δ| > 3.0: {sum(1 for d in ds2 if d > 3.0)}")
print(f"  剔除「账本含未解离酸/碱」（欠定）后: {len(clean)} 例")
if ds2c:
    print(f"    |Δ| > 1.0: {sum(1 for d in ds2c if d > 1.0)}  ← 真冲突面")
    print(f"    |Δ| > 3.0: {sum(1 for d in ds2c if d > 3.0)}")
    print(f"    |Δ| p90 = {ds2c[int(len(ds2c) * 0.9)]:.4f}")

print()
print("=" * 76)
print("③ 判读")
print("=" * 76)
print("  ① 若 max 很大 ⟹ **引擎报告的答案与它自己搜索的答案不是同一个**")
print("     （这是最严重的一类'同一数值多套算法'）。")
print("  ② 才是「搜索用了与账本矛盾的 pH」，须先剔除欠定用例再计数。")

out = os.path.join(ROOT, "logs", "conflict_audit.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"scanned": scanned, "rows": rows}, f, ensure_ascii=False,
              indent=1, default=str)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
