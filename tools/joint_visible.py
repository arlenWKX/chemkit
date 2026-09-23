# -*- coding: utf-8 -*-
"""第 238 轮 · 验证新线索：`_joint_collect` 在 B 组终态到底"看得见"几个通道？

第 237 轮留下的怀疑：`_joint_collect` 只收**两侧物种都在场**的方程
（`all(ledger.get(s,0) > X_MIN for s in ps[0] + ps[1])`），
而 B 组的**主驱动通道是单侧的**（如 `Al(OH)₃+Cl⁻+3H⁺→[AlCl]²⁺`
的产物侧全缺席）⟹ 联立看不到它们 ⟹ `joint_solve` 无方程可解。

本脚本**在终态账本上重建 `_joint_collect` 的过滤**（不改引擎），
对每个候选按同一判据分类，给出：
  · 两侧在场（联立**看得见**）
  · 单侧在场（联立**看不见**）—— 并列出其 |S|
从而判定：B 组停滞时，联立"看得见"的方程里有没有强驱动的。

用法： python tools/joint_visible.py [用例前缀...]
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
import chemkit.joint as joint                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
PREFIXES = sys.argv[1:] or ["H45", "F31", "H43", "T52", "E41", "E55", "M01"]
cases = {c["name"]: c for c in load_cases(None)}

for pre in PREFIXES:
    hit = [v for n, v in cases.items() if n.startswith(pre)]
    if not hit:
        print(f"{pre}: 未找到")
        continue
    c = hit[0]
    probe = {}
    eng.judge([{"name": a, "mol": b} for a, b in c["subs"]],
              c.get("cond") or {"V_L": V}, T, _probe=probe)
    led = dict(probe.get("ledger") or {})
    pH = probe.get("pH_solver")
    He = probe.get("H_excess")
    print(f"\n{'=' * 74}")
    print(f"{c['name']}   pH={pH}  max|S|={probe.get('max_abs_S')}")
    cands = enumerate_candidates(led, He or 0.0, pH or 7.0, V, T_K, T, True)
    vis, inv = [], []
    for cd in cands:
        ps = cd.pres_specs
        r_ok = all(led.get(s, 0.0) > eng.X_MIN for s in ps[0])
        p_ok = all(led.get(s, 0.0) > eng.X_MIN for s in ps[1])
        if not (r_ok or p_ok):
            continue
        try:
            S = eng.S_of(cd, led, V, pH, T_K, T, frozenset(),
                         eng.P_EXT_KPA, True, {})
        except Exception:                                   # noqa: BLE001
            continue
        if r_ok and p_ok:
            vis.append((abs(S), cd, S))
        else:
            inv.append((abs(S), cd, S, "反应物侧" if r_ok else "产物侧"))
    vis.sort(key=lambda t: -t[0])
    inv.sort(key=lambda t: -t[0])
    print(f"  **两侧在场（联立看得见）**: {len(vis)} 条")
    for ab, cd, S in vis[:5]:
        print(f"    |S|={ab:7.3f} {cd.kind:9s} {cd.r} -> {cd.pr}")
    print(f"  **单侧在场（联立看不见）**: {len(inv)} 条")
    for ab, cd, S, side in inv[:5]:
        print(f"    |S|={ab:7.3f} {cd.kind:9s} [{side}在场] "
              f"{cd.r} -> {cd.pr}")
    vmax = vis[0][0] if vis else 0.0
    imax = inv[0][0] if inv else 0.0
    print(f"  ⟹ 看得见的最大 |S| = {vmax:.3f}；看不见的最大 |S| = {imax:.3f}")
    print(f"     JOINT_MIN_M = {joint.JOINT_MIN_M}（需 ≥ 此数才问联立）")

print("\n=== 判读 ===")
print("  · 若'看得见'的条数 < JOINT_MIN_M 或其中最大 |S| 很小 ⟹ **联立无方程可解**，")
print("    证实第 237 轮线索 ⟹ 修法是**放宽 `_joint_collect` 的在场判据**")
print("    （允许单侧通道进入，因为它同样携带耦合约束）。")
print("  · 若'看得见'里有强驱动者 ⟹ 线索不成立，联立是'有方程却解不动'。")
