# -*- coding: utf-8 -*-
"""第 217 轮 · 通用版 `S_of(x)` 扫描：判定"零推进"家族病是否同机制。

第 216 轮发现 F31 与 H45 的最后一环同构（候选被正确挑中、pH 也对、
S 强正，但 `solve_extent` 返回 `x*=0`），共同特征 = **被选通道需多质子
（4 / 3 个）而账本处于强碱态（游离 H⁺ ≈1e-12）**。

本脚本把第 209 轮 `extzero.py` 的手法**通用化**：对任意用例
  ① 跑到终态，取账本与 He；
  ② 找出"残留最大"的那个两侧候选（引擎自报 `resid_src_eq`）；
  ③ 按 `_exec` 语义手工推账本，扫 `S_of(x)`，打印每点的 **pH 与 S**。
判据：
  · 若 `x=0` 处 pH 取酸侧（低）而挪动极小 x 就翻碱侧（高）⟹ **pH 双模态**，
    与 F31 同机制；
  · 若 pH 单调而 S 仍不自洽 ⟹ 另一机制。

用法： python tools/xscan.py <用例前缀> [点数，默认 12]
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
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "H45"
NPTS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
T = load_tables()
V, T_K = 1.0, 298.15
hit = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
       if n.startswith(PRE + " ") or n == PRE]
if not hit:
    print(f"未找到用例 {PRE}")
    sys.exit(1)
c = hit[0]
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
probe = {}
eng.judge(subs, c.get("cond") or {"V_L": V}, T, _probe=probe)
led = dict(probe.get("ledger") or {})
He = probe.get("H_excess")
pH = probe.get("pH_solver")
print(f"用例 {c['name']}")
print(f"  终态 pH_solver={pH}  He={He}  max|S|={probe.get('max_abs_S')}")
print(f"  账本={ {k: round(v, 6) for k, v in led.items() if k != 'H_2O'} }")
src = probe.get("resid_src") or {}
print(f"  引擎自报残差源: {src.get('eq')}  S={src.get('S')}")
print(f"\n=== 逐候选的 (S, ext, x_max) 在终态 ===")
cands = enumerate_candidates(led, He or 0.0, pH or 7.0, V, T_K, T, True)
rows = []
for cd in cands:
    ps = cd.pres_specs
    pres = (all(led.get(x, 0.0) > eng.X_MIN for x in ps[0])
            and all(led.get(x, 0.0) > eng.X_MIN for x in ps[1]))
    if not pres:
        continue
    try:
        S = eng.S_of(cd, led, V, pH or 7.0, T_K, T, frozenset(),
                     eng.P_EXT_KPA, True, {})
    except Exception:                                       # noqa: BLE001
        continue
    if abs(S) < 1.0:
        continue
    d = 1 if S > 0 else -1
    try:
        ext, xm = eng.solve_extent(cd, d, led, He or 0.0, V, T_K, T,
                                   frozenset())
    except Exception as exc:                                # noqa: BLE001
        ext, xm = None, f"EXC {type(exc).__name__}"
    rows.append((abs(S), cd, d, S, ext, xm))
rows.sort(key=lambda r: -r[0])
for ab, cd, d, S, ext, xm in rows[:5]:
    _e = f"{ext:.6g}" if isinstance(ext, float) else str(ext)
    _x = f"{xm:.6g}" if isinstance(xm, float) else str(xm)
    print(f"  |S|={ab:8.3f} d={d:+d} ext={_e:>12} x_max={_x:>12}  "
          f"{cd.r} -> {cd.pr}")

if not rows:
    print("  （没有 |S|>1 的两侧候选）")
    sys.exit(0)

TARGET = rows[0][1]
D = rows[0][2]
rr = TARGET.r if D > 0 else TARGET.pr
x_max = min(led.get(s, 0.0) / nu for s, nu in rr.items()
            if s not in (eng.WATER, eng.H_ION))
print(f"\n=== 扫描最大残差通道的 S(x)："
      f"{TARGET.r} -> {TARGET.pr}（d={D:+d}, x_max={x_max:.6g}）===")


def step(ledger, cd, d, x):
    out = dict(ledger)
    he = 0.0
    _rr = cd.r if d > 0 else cd.pr
    _pp = cd.pr if d > 0 else cd.r
    for s_, nu in _rr.items():
        if s_ == eng.H_ION:
            he -= nu * x
        elif s_ != eng.WATER:
            out[s_] = out.get(s_, 0.0) - nu * x
    for s_, nu in _pp.items():
        if s_ == eng.H_ION:
            he += nu * x
        elif s_ != eng.WATER:
            out[s_] = out.get(s_, 0.0) + nu * x
    return out, he


print(f"  {'x':>12} {'pH':>10} {'He':>13} {'S(x)':>12}")
xs = [0.0] + [x_max * (i / NPTS) for i in range(1, NPTS)] + [x_max]
prev_ph = None
jump = None
for x in xs:
    lx, dhe = step(led, TARGET, D, x)
    He_x = (He or 0.0) + dhe
    try:
        ph_x = eng.estimate_pH(lx, He_x, V, T, T_K)
        S_x = eng.S_of(TARGET, lx, V, ph_x, T_K, T, frozenset(),
                       eng.P_EXT_KPA, True, {})
    except Exception as exc:                                # noqa: BLE001
        print(f"  {x:12.6g} 异常 {type(exc).__name__}: {exc}")
        continue
    mark = ""
    if prev_ph is not None and abs(ph_x - prev_ph) > 1.0:
        mark = f"   <== pH 跳 {ph_x - prev_ph:+.2f}"
        if jump is None:
            jump = (x, prev_ph, ph_x)
    print(f"  {x:12.6g} {ph_x:10.4f} {He_x:+13.6f} {S_x:+12.3f}{mark}")
    prev_ph = ph_x

print("\n=== 判读 ===")
if jump:
    xj, p0, p1 = jump
    print(f"  **发现 pH 跳变**：x={xj:.6g} 处 {p0:.3f} -> {p1:.3f}"
          f"（跳 {p1 - p0:+.2f}）")
    print("  ⟹ 与 F31 同机制（estimate_pH 双模态）⟹ **修一处可兼治**。")
else:
    print("  **未发现 pH 跳变** ⟹ H45 与 F31 **不同机制**，是独立的第二缺陷。")
