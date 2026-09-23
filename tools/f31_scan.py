# -*- coding: utf-8 -*-
"""第 203 轮 · F31 首步 f(x) 细扫（机制收口：走步方向为何是错的）。

已定影事实：
  · 起始态 ledger={Ga3+ 1, Na+ 3, Cl- 3}, He=-3.0, estimate_pH=14.477 ✓ 化学正确
  · 第一个执行步 = `Ga^{3+} → [Ga(OH)_4]^- + 4H^+`，记录 logK=-18.4, S=+50.46,
    extent=0.743966 —— 这一步把 3 mol 强碱用光并造出 2.98 mol 强酸
  · 终态 pH=1.62 / 固相 0；严格化学解 pH=5.764 / 固相 1.000 mol
  · 终态该反应 S=+13.305（在 pH 1.62 上评）但走步口径 He 把同一账本映到 pH 12.38

本脚本把**起始态**上每一步的 f(x)=S(x) 扫出来（直接用引擎的 S_of +
_exec 记账语义），回答：首步的二分根在哪、为什么根值会落在 0.744（化学
上不可能），以及 S(x) 曲线是否在 x 轴上有多个根（跨 pH 分支）。

用法： python tools/f31_scan.py [用例前缀，默认 F31]
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

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith(PRE + " ")][0]
V, T_K = 1.0, 298.15
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
led0, He0, _st, _un = _norm.normalize(
    subs, {"V_L": V, "T_K": T_K, "c_H": None, "c_OH": None, "pH": None,
           "p_kpa": eng.P_EXT_KPA}, T, {})
print(f"起始 He = {He0}  ledger = "
      f"{ {k: round(v, 6) for k, v in led0.items() if k != 'H_2O'} }")
pH0 = eng.estimate_pH(led0, He0, V, T, T_K)
print(f"起始 estimate_pH = {pH0:.4f}")

cands = enumerate_candidates(led0, He0, pH0, V, T_K, T, True)
print(f"起始候选 {len(cands)} 条")

# 目标反应族：Ga 羟合/沉淀（含固相生成）
TARGETS = ("[Ga(OH)_4]^-", "Ga(OH)_3", "[Ga(OH)]^{2+}")


def _nu_H(cd):
    return cd.pr.get(eng.H_ION, 0) - cd.r.get(eng.H_ION, 0)


def _exec_ledger(led, cd, d, x):
    """按 _exec 的记账语义推进一步（不改原账本）。"""
    out = dict(led)
    rr = cd.r if d > 0 else cd.pr
    pp = cd.pr if d > 0 else cd.r
    he = 0.0
    for s_, nu in rr.items():
        if s_ == eng.H_ION:
            he -= nu * x
        elif s_ != eng.WATER:
            out[s_] = out.get(s_, 0.0) - nu * x
    for s_, nu in pp.items():
        if s_ == eng.H_ION:
            he += nu * x
        elif s_ != eng.WATER:
            out[s_] = out.get(s_, 0.0) + nu * x
    return out, he


sel = []
for cd in cands:
    names = set(cd.r) | set(cd.pr)
    if any(t in names for t in TARGETS):
        sel.append(cd)
print(f"含 Ga 羟合/固相的候选 {len(sel)} 条\n")

for cd in sel:
    _d = 1
    ps = cd.pres_specs
    pres = all(led0.get(x, 0.0) > eng.X_MIN for x in ps[0])
    print(f"--- {cd.kind:9s} {cd.r} -> {cd.pr}")
    print(f"    反应物在场={pres}  nu_H={_nu_H(cd):+d}")
    if not pres:
        continue
    S0 = eng.S_of(cd, led0, V, pH0, T_K, T, frozenset(), eng.P_EXT_KPA,
                  True, {})
    print(f"    S(x=0) = {S0:+.3f}")
    print(f"    {'x':>10} {'pH_est':>9} {'He':>12} {'S(x)':>12}")
    for x in (0.0, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5,
              0.743, 0.75, 1.0):
        lx, he_x = _exec_ledger(led0, cd, _d, x)
        He_x = He0 + he_x
        try:
            ph = eng.estimate_pH(lx, He_x, V, T, T_K)
            Sx = eng.S_of(cd, lx, V, ph, T_K, T, frozenset(), eng.P_EXT_KPA,
                          True, {})
            print(f"    {x:10.5g} {ph:9.3f} {He_x:+12.5g} {Sx:+12.3f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"    {x:10.5g} 异常 {type(exc).__name__}: {exc}")
    try:
        ext, x_max = eng.solve_extent(cd, _d, led0, He0, V, T_K, T,
                                      frozenset())
        print(f"    ⟹ solve_extent 给出 ext={ext:.6g} x_max={x_max:.6g}")
    except Exception as exc:                                # noqa: BLE001
        print(f"    ⟹ solve_extent 异常 {type(exc).__name__}: {exc}")
    print()

print("=== 判读 ===")
print("  ① 若 S(x) 在 x∈(0, x_max) 内单调穿零一次，根位置即平衡程度；")
print("     根 ≈0.744 表示'引擎认为把 0.744 Ga 变成镓酸根后 S 才归零'——")
print("     但该态 pH≈1.6、游离 H⁺ 0.024 mol，与'释放 2.98 mol H⁺'矛盾；")
print("  ② 若 S(x) 因 pH 分支切换出现**多根/不连续**，则二分落在跨分支的")
print("     伪根上（D14 接缝类的 pH 分支悬崖），修法在分支一致性而非预算。")
