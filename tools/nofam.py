# -*- coding: utf-8 -*-
"""第 210 轮 · 量化「账本无 pKa 质子族」条件的命中面（动手前必做）。

第 209 轮定死：F31 死锁源于**二分起点上的 pH 双模态**——同一账本在
`x=0` 取酸侧 1.620，挪 1e-6 就翻碱侧 5.764 ⟹ 二分判"无过零点"、返回 x*=0。

拟修法（窄）：**仅当账本里没有 pKa 质子族成员时**，让 `estimate_pH` 改用
`acidbase.charge_pH`。理由：这种账本上**没有 Henderson/族分布语义**，
诱导式分支（弱酸式/直读/比大小）本就是唯一依据，而 `charge_pH` 是精确的
（非族物种直接进电荷方程，族分布无歧义）。

本脚本回答三个数（决定该条件能不能用）：
  ① 全库有多少用例的**终态账本**落在这个条件下（无 pKa 族成员）？
  ② 这些例里 `charge_pH` 与 `estimate_pH` 差多少？差大的那些是不是
     就是"双模态"受害者（即改动的目标面）？
  ③ 这些例**当前套件通过性**如何（改动的风险面：PASS 的例不能被改坏）。

用法： python tools/nofam.py [步长，默认 1]
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
# pKa 质子族成员集（与 exact_proton_pH 的 fam_set 同源）
FAM = frozenset(__import__("chemkit.acidbase", fromlist=["x"]).build_families(T))

# 已存档的通过性（读档，不重跑）
passmap = {}
try:
    d = json.load(io.open(os.path.join(ROOT, "logs",
                                       "suite-parallel-latest.json"),
                          encoding="utf-8"))
    passmap = {c["name"]: c["ok"] for c in d["cases"]}
    print(f"载入通过性留档：{len(passmap)} 例")
except Exception as exc:                                    # noqa: BLE001
    print(f"（无通过性留档：{exc}）")

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
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  {"V_L": V, "T_K": T_K}, T, _probe=probe)
    except Exception:                                       # noqa: BLE001
        continue
    led = dict(probe.get("ledger") or {})
    if not led:
        continue
    scanned += 1
    nfam = sum(1 for sp, m in led.items()
               if m > spec.X_MIN and sp in FAM)
    if nfam:
        continue                                # 有 pKa 族 ⟹ 不在条件内
    ph_est = probe.get("pH_solver")
    if not isinstance(ph_est, float):
        continue
    try:
        ph_ch = ab.charge_pH(led, V, T, T_K)
    except Exception:                                       # noqa: BLE001
        ph_ch = None
    gap = (abs(ph_ch - ph_est) if isinstance(ph_ch, float) else None)
    rows.append({"case": c["name"], "ph_est": ph_est, "ph_ch": ph_ch,
                 "gap": gap, "ok": passmap.get(c["name"]),
                 "resid": probe.get("max_abs_S"),
                 "solid": any(m > spec.X_MIN and sp in T.solids
                              for sp, m in led.items())})

rows.sort(key=lambda r: -(r["gap"] if r["gap"] is not None else -1))
print(f"\n扫描 {scanned} 例；**终态账本无 pKa 质子族**的: {len(rows)} 例")
n_ok = sum(1 for r in rows if r["ok"] is True)
n_fail = sum(1 for r in rows if r["ok"] is False)
print(f"  其中当前 PASS {n_ok} 例 / FAIL {n_fail} 例 / 未知 "
      f"{sum(1 for r in rows if r['ok'] is None)} 例")
gaps = sorted(r["gap"] for r in rows if r["gap"] is not None)
if gaps:
    print(f"  |charge_pH − estimate_pH|: n={len(gaps)}  "
          f"p50={gaps[len(gaps)//2]:.3f}  max={gaps[-1]:.3f}  "
          f">1 的 {sum(1 for g in gaps if g > 1)} 例")

print(f"\n{'用例':34s} {'est':>8} {'charge':>8} {'Δ':>7} {'套件':>5} "
      f"{'max|S|':>9} 固相")
for r in rows[:30]:
    g = f"{r['gap']:7.3f}" if r["gap"] is not None else "      -"
    ch = f"{r['ph_ch']:8.3f}" if isinstance(r["ph_ch"], float) else "       -"
    ok = {True: "PASS", False: "FAIL", None: "?"}[r["ok"]]
    print(f"{r['case'][:34]:34s} {r['ph_est']:8.3f} {ch} {g} {ok:>5} "
          f"{str(r['resid'])[:9]:>9} {r['solid']}")
if len(rows) > 30:
    print(f"  … 另 {len(rows) - 30} 例")

print("\n=== 判读 ===")
print("  · 命中面小（≤30）+ 其中 PASS 的少 ⟹ 条件可用，改动窄。")
print("  · 若命中面里 **PASS 的很多**且 gap 大 ⟹ 改动会波及现行通过的例，")
print("    必须逐例核对（或再加窄条件）。")
print("  · 若 gap 大且 FAIL ⟹ 正是双模态受害者，是修复目标。")

out = os.path.join(ROOT, "logs", "nofam.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"n_nofam": len(rows), "n_pass": n_ok, "n_fail": n_fail,
               "rows": rows}, f, ensure_ascii=False, indent=1, default=str)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
