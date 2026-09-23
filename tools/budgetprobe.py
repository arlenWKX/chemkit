# -*- coding: utf-8 -*-
"""第 204 轮 · 直接测 `solve_extent` 的酸碱预算是否对**非 redox** 生效。

第 203 轮的判断需要更正：F31 被闸挡住有**两个**理由——
`exact_solid_reason = "no_family"`（镓羟合梯在 beta 表里，账本里**没有
pKa 质子族**）＋ 曾记的"Ksp-OH 阳离子在场"。而 `charge_pH` 在全库里与
引擎现值差 >1 的有 **156 例**（E55 6.26 vs 10.39、B26 1.36 vs 12.64…），
**那些例里引擎是对的**（charge_pH 把总量当溶解量，忽略固相）⟹
"放行精确解"不是通用答案。

回到机制层：F31 首步把 0.744 mol 的 `Ga³⁺+4H₂O→[Ga(OH)₄]⁻+4H⁺` 走成，
需要 2.976 mol 碱（He=−3.0 有），但沉淀步 `Ga³⁺+3H₂O→Ga(OH)₃+3H⁺`
（需 3.0 mol 碱，恰好够）**从未执行**。本脚本直接问引擎：
在**起始态**上，这两条候选各自的 `solve_extent` 给出什么？
若沉淀步的 (ext, x_max) 被酸碱预算正确封顶而仍未执行，则病根在
**pick 的择优顺序**；若 x_max 异常（预算没施加），则病根在预算闸。

用法： python tools/budgetprobe.py
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
import chemkit.normalize as _norm                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith("F31 ")][0]
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
led, He, _st, _un = _norm.normalize(
    subs, {"V_L": V, "T_K": T_K, "c_H": None, "c_OH": None, "pH": None,
           "p_kpa": eng.P_EXT_KPA}, T, {})
pH = eng.estimate_pH(led, He, V, T, T_K)
print(f"起始态: He = {He}  pH = {pH:.4f}  账本 = "
      f"{ {k: v for k, v in led.items() if k != 'H_2O'} }")
print(f"碱储备 = -He = {-He:.4f} mol\n")

cands = enumerate_candidates(led, He, pH, V, T_K, T, True)
print(f"{'kind':10s} {'nu_H':>5} {'S':>9} {'ext':>12} {'x_max':>10} "
      f"{'预算(-He/nu_H)':>14}  反应")
rows = []
for cd in cands:
    ps = cd.pres_specs
    pres = all(led.get(x, 0.0) > eng.X_MIN for x in ps[0])
    if not pres:
        continue
    kinds = set(cd.r) | set(cd.pr)
    if not any("Ga" in k for k in kinds):
        continue
    nu_H = cd.pr.get(eng.H_ION, 0) - cd.r.get(eng.H_ION, 0)
    S = eng.S_of(cd, led, V, pH, T_K, T, frozenset(), eng.P_EXT_KPA, True, {})
    d = 1 if S > 0 else -1
    try:
        ext, x_max = eng.solve_extent(cd, d, led, He, V, T_K, T, frozenset())
    except Exception as exc:                                # noqa: BLE001
        ext, x_max = None, f"EXC {type(exc).__name__}"
    budget = (He / (-nu_H)) if nu_H < 0 else ((-He / nu_H) if nu_H > 0
                                              else float("inf"))
    rows.append((abs(S), cd.kind, nu_H, S, ext, x_max, budget,
                 f"{cd.r} -> {cd.pr}", d))
rows.sort(key=lambda t: -t[0])
for ab, kind, nu_H, S, ext, x_max, budget, eq, d in rows:
    _e = f"{ext:.6g}" if isinstance(ext, float) else str(ext)
    _x = f"{x_max:.6g}" if isinstance(x_max, float) else str(x_max)
    print(f"{kind:10s} {nu_H:+5d} {S:+9.3f} {_e:>12} {_x:>10} "
          f"{budget:>14.4g}  {eq[:56]}")

print("\n=== 判读 ===")
print("  · 若沉淀步（nu_H=+3）的 x_max ≫ 预算值 ⟹ 酸碱预算**没施加**在")
print("    非 redox 上 ⟹ 走步可以让产酸步超出碱储备，制造自相矛盾的态。")
print("  · 若沉淀步 ext 或 x_max 被正确限制 ⟹ 预算没问题，病根在 pick 择优：")
print("    引擎按 |S| 选步，配体腿 S 更大（+50.5 vs +35.0），先走配体腿；")
print("    这本身合理，**问题在于走完配体腿后沉淀腿再没被评估过**。")
