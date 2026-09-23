# -*- coding: utf-8 -*-
"""第 203 轮 · 窄化 `exact_proton_pH` 闸的**触发面测量**（改 speciation.py L1191 前必做）。

现状闸门（speciation.py L1185-1201）：
    res_set = 全部固相 ∪ 全部 Ksp-OH 阳离子
    for sp, m in ledger: if m > X_MIN and sp in res_set: return None
    if not nfam or abs(net + H_excess) > 1e-6: return None
    return charge_pH(ledger, V, T, T_K, fast=True)

拟改动：`res_set` 的"Ksp-OH 阳离子"一半收窄为"**该阳离子的固相在场**"。
理由：闸的原始理由是"pH 由 Ksp 决定"——那只在**固相在场**时成立；
固相不在场时 pH 是自由变量，账本质子条件恰好是权威解（F31 实测：
闸后 charge_pH 给 12.3803 正确，闸前被拒，退回启发式给 1.62 错）。

本脚本（只读）逐例取**终态账本**，判定该例属于哪一类：
  A 有固相在场                      → 收窄后**仍被拒**（行为不变，安全）
  B 无固相但在场 Ksp-OH 阳离子      → 收窄后**放行**（行为会变 ⟹ 风险面）
  C 两者皆无                        → 收窄前后都放行（行为不变）
并对 B 类逐例给出：现在 exact_proton_pH / estimate_pH / charge_pH 的值，
以便看清"放行后 pH 会移动多少、方向对不对"。

用法： python tools/gateface.py [样本步长，默认 1]
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
import chemkit.speciation as spec                           # noqa: E402
import chemkit.acidbase as ab                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 1
T = load_tables()
cases = load_cases(None)
V_DEF = 1.0

# B3 静态集（与 exact_proton_pH 内部同源）
solids = frozenset(T.solids)
ksp_cats = frozenset(e["pair"][0] for e in T.ksp if e["pair"][1] == "OH^-")
fams = frozenset(__import__("chemkit.acidbase", fromlist=["x"]).build_families(T))
print(f"固相 {len(solids)} 种；Ksp-OH 阳离子 {len(ksp_cats)} 种；质子族成员 {len(fams)}")
print(f"Ksp-OH 阳离子: {sorted(ksp_cats)}")

A = B = C = 0
blist = []
print(f"\n扫描 {len(cases)} 例（步长 {STEP}）…")
for i, c in enumerate(cases):
    if i % STEP:
        continue
    probe = {}
    try:
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": V_DEF}, T, _probe=probe)
    except Exception:                                       # noqa: BLE001
        continue
    led = dict(probe.get("ledger") or {})
    if not led:
        continue
    V = float((c.get("cond") or {}).get("V_L", V_DEF))
    T_K = float((c.get("cond") or {}).get("T_K", 298.15))
    has_solid = any(m > spec.X_MIN and sp in solids for sp, m in led.items())
    cats = [sp for sp, m in led.items()
            if m > spec.X_MIN and sp in ksp_cats]
    if has_solid:
        A += 1
        continue
    if cats:
        B += 1
        try:
            est = eng.estimate_pH(led, probe.get("H_excess") or 0.0, V, T, T_K)
        except Exception as exc:                            # noqa: BLE001
            est = f"EXC {type(exc).__name__}"
        try:
            chg = ab.charge_pH(led, V, T, T_K)
        except Exception as exc:                            # noqa: BLE001
            chg = f"EXC {type(exc).__name__}"
        try:
            ex = spec.exact_proton_pH(led, probe.get("H_excess") or 0.0, V, T, T_K)
        except Exception as exc:                            # noqa: BLE001
            ex = f"EXC {type(exc).__name__}"
        blist.append({"case": c["name"], "cats": cats,
                      "probe_pH": probe.get("pH"), "est": est,
                      "charge": chg, "exact": ex,
                      "d": (abs(chg - est) if isinstance(chg, float)
                            and isinstance(est, float) else None)})
    else:
        C += 1

print(f"\n=== 分类 ===")
print(f"  A 有固相在场（收窄后仍拒，行为不变）: {A}")
print(f"  B **无固相但在场 Ksp-OH 阳离子**（收窄后放行 ⟹ 风险面）: {B}")
print(f"  C 两者皆无（行为不变）: {C}")

if blist:
    blist.sort(key=lambda r: -(r["d"] if r["d"] is not None else -1))
    print(f"\n=== B 类逐例（按 |charge_pH − estimate_pH| 降序，前 30）===")
    print(f"  {'用例':34s} {'阳离子':22s} {'现在pH':>8} {'est':>8} "
          f"{'charge_pH':>10} {'Δ':>8}")
    for r in blist[:30]:
        est = f"{r['est']:.4f}" if isinstance(r["est"], float) else str(r["est"])[:8]
        chg = f"{r['charge']:.4f}" if isinstance(r["charge"], float) else str(r["charge"])[:8]
        d = f"{r['d']:+.3f}" if r["d"] is not None else "-"
        print(f"  {r['case'][:34]:34s} {','.join(r['cats'])[:22]:22s} "
              f"{str(r['probe_pH'])[:8]:>8} {est:>8} {chg:>10} {d:>8}")
    big = [r for r in blist if r["d"] is not None and r["d"] > 1.0]
    print(f"\n  B 类中 |Δ|>1.0 的: {len(big)} / {B}"
          f"  ⟹ 这些例的 pH 会**移动**（可能是纠正，也可能是回归）")

print(f"\n=== 判读 ===")
print("  · A 类完全不受影响 ⟹ 闸的核心用途（防 L08 型储库误判）被保留。")
print("  · B 类是会变的用例：需逐例判'移动后的 pH 是更对还是更错'，")
print("    不能只看数量。判据：与账本电荷自洽值一致者为纠正。")
print("  · 若 B 类里多数 |Δ| 很小(<0.01) ⟹ 改动近乎零风险（两口径本就一致）。")

out = os.path.join(ROOT, "logs", "gateface.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"A_solid": A, "B_cation_no_solid": B, "C_neither": C,
               "B_detail": blist}, f, ensure_ascii=False, indent=1, default=str)
print(f"\n已写 {out}")
