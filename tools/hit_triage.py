# -*- coding: utf-8 -*-
"""第 232 轮 · 32 例命中例的二次细分：**修分支 4 真正能治几例？**

第 231 轮发现 32 例里有一批 `pH_eng == pH_chg`（引擎 pH 与电荷自洽解一致）
⟹ 它们的残差**不是 pH 错**，修 pH 治不了。

本脚本（一次运行，写全部产物）把 32 例分成三类并给出**收益上界**：
  A **pH 错**（`|pH_eng − pH_chg| > 1` 且无固相）⟹ `charge_pH` 可治
  B **pH 错但有固相** ⟹ 需 `pinned` 路线
  C **pH 对**（`|Δ| ≤ 1`）⟹ 残差另有来源，**pH 修法治不了**
并给出"若 A+B 全修好"的残差减少量（= 修分支 4 的收益上界）。

用法： python tools/hit_triage.py
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

A, B, C = [], [], []
for nm in hits:
    c = cases.get(nm)
    if c is None:
        continue
    probe = {}
    eng.judge([{"name": a, "mol": b} for a, b in c["subs"]],
              c.get("cond") or {"V_L": V}, T, _probe=probe)
    led = dict(probe.get("ledger") or {})
    solid = any(sp in SOLIDS and m > eng.X_MIN for sp, m in led.items())
    phe = probe.get("pH_solver")
    try:
        phc = ab.charge_pH(led, V, T, T_K)
    except Exception:                                       # noqa: BLE001
        phc = None
    rv = probe.get("resid_live") or probe.get("max_abs_S") or 0.0
    d = (abs(phe - phc) if isinstance(phe, float)
         and isinstance(phc, float) else None)
    rec = (nm, rv, phe, phc, d, solid)
    if d is not None and d > 1.0:
        (B if solid else A).append(rec)
    else:
        C.append(rec)

print(f"32 例细分：")
print(f"  A 无固相且 pH 错（charge_pH 可治） : {len(A):2d} 例  "
      f"resid 合计 {sum(r[1] for r in A):.1f}")
print(f"  B 有固相且 pH 错（需 pinned）      : {len(B):2d} 例  "
      f"resid 合计 {sum(r[1] for r in B):.1f}")
print(f"  C pH 本就对（残差另有来源）        : {len(C):2d} 例  "
      f"resid 合计 {sum(r[1] for r in C):.1f}")

print(f"\n=== A + B（修分支 4 的收益上界）===")
print(f"  例数 {len(A) + len(B)}，残差合计 "
      f"{sum(r[1] for r in A) + sum(r[1] for r in B):.1f}")
for tag, grp in (("A", A), ("B", B), ("C", C)):
    if not grp:
        continue
    print(f"\n  --- {tag} ---")
    for nm, rv, phe, phc, d, solid in sorted(grp, key=lambda r: -r[1]):
        dd = f"{d:6.2f}" if d is not None else "  None"
        print(f"    {nm[:32]:32s} resid={rv:7.3f} pH={str(phe)[:7]:>7} "
              f"chg={str(phc)[:7]:>7} Δ={dd} 固相={solid}")

out = os.path.join(ROOT, "logs", "hit_triage.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"A": A, "B": B, "C": C}, f, ensure_ascii=False, indent=1,
              default=str)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
print("\n=== 判读 ===")
print("  `len(A)+len(B)` = 修分支 4 能触及的例数；其 resid 合计 = 收益上界。")
print("  `len(C)` 若不小 ⟹ 这 32 例里近半的残差**与 pH 无关**，")
print("  必须单独记账，否则会把它们的残差误算进修法收益。")
