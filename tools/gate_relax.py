# -*- coding: utf-8 -*-
"""第 213 轮 · 因果模拟：放宽退化区接管的两个子条件，F31 会怎样？

第 212 轮定影（`tools/gatecheck.py`）：`estimate_state` L1006 的入口
```
if not _res and _nf >= 2 and abs(_net + He_res) <= 1e-6:
```
在 F31 上 `_nf = 0`（镓梯在 beta 表）、`abs(_net+He_res) = 0.024` ⟹ 接管被挡。
其中 `abs(_net + He_res) <= 1e-6` 是**循环论证**（要求账本先自洽，
才允许调用"解出那个让账本自洽的 pH"的函数）。

本脚本用**不改源码**的方式测三个变体（monkey-patch `speciation` 内的判定
不可行——判定是 `estimate_state` 内联的；故改为**直接复算**：
在 F31 终态上按三种放宽各自调用 `charge_pH`，看会得到什么 pH），
并在**整套件**上量"放宽后有多少例的 `estimate_pH` 会改变"。

变体：
  V0 现状                       ：_nf>=2 ∧ |net+He|<=1e-6
  V1 去循环条款                 ：_nf>=2
  V2 去循环条款 + 放宽 nfam     ：(_nf>=2 **或** 账本含 beta 梯物种)

对每个变体统计：全库多少例的 estimate_pH 会与现状不同、差多大、
其中当前 PASS/FAIL 各多少（风险面）。

用法： python tools/gate_relax.py [步长，默认 1]
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
import chemkit.acidbase as ab                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 1
T = load_tables()
FAM = frozenset(ab.build_families(T))
SOLIDS = frozenset(T.solids)
# beta 梯物种（含 center 与 complex），用于 V2 的 "含 beta 梯物种" 判据
BETA_SP = set()
for b in T.beta:
    if b.get("ligand") == "OH^-":
        BETA_SP.add(b["complex"])
        BETA_SP.add(b["center"])
BETA_SP = frozenset(BETA_SP)

passmap = {}
try:
    d = json.load(io.open(os.path.join(ROOT, "logs",
                                       "suite-parallel-latest.json"),
                          encoding="utf-8"))
    passmap = {c["name"]: c["ok"] for c in d["cases"]}
except Exception:                                           # noqa: BLE001
    pass

cases = load_cases(None)
res = {v: [] for v in ("V0", "V1", "V2")}
scanned = 0
for i, c in enumerate(cases):
    if i % STEP:
        continue
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    T_K = float(cond.get("T_K", 298.15))
    probe = {}
    try:
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  {"V_L": V, "T_K": T_K}, T, _probe=probe)
    except Exception:                                       # noqa: BLE001
        continue
    led = dict(probe.get("ledger") or {})
    if not led:
        continue
    ph_est = probe.get("pH_solver")
    if not isinstance(ph_est, float):
        continue
    scanned += 1
    he = probe.get("H_excess")
    if not isinstance(he, float):
        continue
    _nf = sum(1 for sp, m in led.items() if m > 0.0 and sp in FAM)
    _res = any(m > spec.X_MIN and sp in SOLIDS for sp, m in led.items()
               if sp != spec.WATER)
    _has_beta = any(m > spec.X_MIN and sp in BETA_SP for sp, m in led.items())
    _net = 0.0
    try:
        for sp, m in led.items():
            if m > 0.0 and sp != spec.WATER and not sp.startswith("__"):
                _net += spec.charge_of(sp) * m
    except Exception:                                       # noqa: BLE001
        continue
    if _res:
        continue                                    # 有固相一律跳过（保持现状）
    cand = {"V0": (_nf >= 2 and abs(_net + he) <= 1e-6),
            "V1": (_nf >= 2),
            "V2": (_nf >= 2 or _has_beta)}
    for v, fire in cand.items():
        if not fire:
            continue
        try:
            ex = ab.charge_pH(led, V, T, T_K, fast=True)
        except Exception:                                   # noqa: BLE001
            continue
        if not isinstance(ex, float):
            continue
        ex = min(max(ex, -1.0), spec.pKw_of(T_K) + 1.0)
        if abs(ex - ph_est) > 1e-9:
            res[v].append((abs(ex - ph_est), c["name"], ph_est, ex,
                           passmap.get(c["name"]), _nf, _has_beta))

print(f"扫描 {scanned} 例（有固相的已排除）\n")
print(f"{'变体':6s} {'会改变':>7} {'其中PASS':>9} {'其中FAIL':>9} "
      f"{'|Δ|max':>9} {'|Δ|>1':>7}")
for v in ("V0", "V1", "V2"):
    rows = res[v]
    npass = sum(1 for r in rows if r[4] is True)
    nfail = sum(1 for r in rows if r[4] is False)
    mx = max((r[0] for r in rows), default=0.0)
    big = sum(1 for r in rows if r[0] > 1.0)
    print(f"{v:6s} {len(rows):>7} {npass:>9} {nfail:>9} {mx:>9.3f} {big:>7}")

print("\n=== V2（拟采用）中 |Δ|>1 的用例（前 25）===")
print(f"  {'用例':34s} {'est':>8} {'exact':>8} {'Δ':>7} {'套件':>5} "
      f"{'nfam':>5} beta")
for d_, nm, pe, ex, ok, nf, hb in sorted(res["V2"],
                                         key=lambda r: -r[0])[:25]:
    okk = {True: "PASS", False: "FAIL", None: "?"}[ok]
    print(f"  {nm[:34]:34s} {pe:8.3f} {ex:8.3f} {d_:7.3f} {okk:>5} "
          f"{nf:>5} {hb}")

print("\n=== 判读 ===")
print("  · V2 的『其中PASS』若很大 ⟹ 会打红现行用例，需再收紧。")
print("  · 若 V2 明显优于 V1 ⟹ beta 梯物种那一条是必要的。")
print("  · F31 是否在 V2 列表中 ⟹ 该变体至少**触到了**目标。")

out = os.path.join(ROOT, "logs", "gate_relax.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({v: [{"case": r[1], "est": r[2], "exact": r[3], "d": r[0],
                    "ok": r[4], "nfam": r[5], "beta": r[6]}
                   for r in res[v]] for v in res},
              f, ensure_ascii=False, indent=1)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
