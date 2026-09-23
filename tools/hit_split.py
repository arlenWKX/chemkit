# -*- coding: utf-8 -*-
"""第 231 轮 · 32 例签名命中例的「固相在场」分布 —— 决定修法形态。

分支 4 的双模态要修，但修法形态取决于这些态的共同结构。已知：
  · `charge_pH`（无分支精确解）在**无固相**时可用（F31 实测对）；
  · 有固相时它必错（第 213 轮实测：202 例会变、171 例 PASS 被打红）。
  · 引擎已有 `pinned` 储库机制（把 Ksp 写进方程）为"有固相"这一半而生。
故先量：**32 例命中例里，有固相 / 无固相各多少？** —— 它决定：
  · 若多为无固相 ⟹ 可用 `charge_pH` 直接接管（但要能区分开"CuSO4 类"）；
  · 若多为有固相 ⟹ 必须走 `pinned` 路线。

用法： python tools/hit_split.py
"""
import io
import json
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
import chemkit.acidbase as ab                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
SOLIDS = frozenset(T.solids)
cen = json.load(io.open(os.path.join(ROOT, "logs",
                                     "signature_census.json"),
                        encoding="utf-8"))
hits = [h[0] for h in cen["hits"]]
cases = {c["name"]: c for c in load_cases(None)}
V, T_K = 1.0, 298.15

n_solid = n_free = 0
rows = []
for nm in hits:
    c = cases.get(nm)
    if c is None:
        continue
    probe = {}
    eng.judge([{"name": a, "mol": b} for a, b in c["subs"]],
              c.get("cond") or {"V_L": V}, T, _probe=probe)
    led = dict(probe.get("ledger") or {})
    solids = {sp: round(m, 6) for sp, m in led.items()
              if sp in SOLIDS and m > eng.X_MIN}
    ph_e = probe.get("pH_solver")
    try:
        ph_c = ab.charge_pH(led, V, T, T_K)
    except Exception:                                       # noqa: BLE001
        ph_c = None
    if solids:
        n_solid += 1
    else:
        n_free += 1
    rows.append((nm, probe.get("resid_live") or probe.get("max_abs_S"),
                 solids, ph_e, ph_c, led))

print(f"32 例命中例：**有固相 {n_solid} 例 / 无固相 {n_free} 例**\n")
print(f"{'用例':32s} {'max|S|':>8} {'pH_eng':>8} {'pH_chg':>8}  固相")
for nm, rv, solids, phe, phc, _led in rows:
    mark = "★" if solids else " "
    sv = ",".join(f"{k}:{v:.4g}" for k, v in list(solids.items())[:2])
    pc = f"{phc:8.3f}" if isinstance(phc, float) else "     None"
    print(f"{mark}{nm[:32]:32s} {str(round(rv,2) if rv else rv):>8} "
          f"{str(phe)[:8]:>8} {pc}  {sv[:34]}")

print("\n=== 判读 ===")
print("  · 有固相占比高 ⟹ 修法主体是 **pinned 路线**（把 Ksp 写进方程），")
print("    因为 charge_pH 在有固相时已知会错。")
print("  · 无固相那部分才可能用 charge_pH 直接接管。")
print("  · 两类都多 ⟹ 修法是**两个入口**：无固相走 charge_pH，")
print("    有固相走 pinned；共同点是**都不做 h_c/o_c 比大小**。")
