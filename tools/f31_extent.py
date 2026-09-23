# -*- coding: utf-8 -*-
"""第 203 轮 · F31 终态逐候选 extent 定影（机制判定的最后一步）。

背景：F31 终态 pH=1.62, 固相 0，而严格化学解 pH=5.764, 固相 1.000 mol。
两个互斥解释必须先分开：
  (H1) `solve_extent` 对**非 redox** 候选不做酸碱预算 ⟹ 反向步被算出巨大
       x_max，实际 ext≈0 来自 f(x) 形状 / 在环 pH 错；
  (H2) 反向步确实被拦住（ext 大但被别处挡下）。

做法：把 F31 跑到终态，用**引擎自己的** `enumerate_candidates` +
`solve_extent`（micro_rel=None，即完整二分）在终态账本上逐候选求
(x*, x_max, S)，把 `logK`、`S`、H⁺ 计量、kind 并排打印。

用法： python tools/f31_extent.py [用例前缀，默认 F31]
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
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith(PRE + " ")][0]
V, T_K = 1.0, 298.15
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
import importlib                                            # noqa: E402
import chemkit.speciation as _spec                          # noqa: E402
import chemkit.normalize as _norm                           # noqa: E402
led0, He0, _st, _un = _norm.normalize(
    subs, {"V_L": V, "T_K": T_K, "c_H": None, "c_OH": None, "pH": None,
           "p_kpa": eng.P_EXT_KPA}, T, {})

# 终态：用真跑一遍拿到（probe 的 ledger 是受限视图，这里直接复跑）
probe = {}
r = eng.judge(subs, {"V_L": V}, T, _probe=probe)
led = dict(probe.get("ledger") or {})
for s_, m in led0.items():
    led.setdefault(s_, m)
led.setdefault("H_2O", 55.6)

_he_raw = probe.get("H_excess")
He = float(_he_raw if _he_raw is not None else 0.0)
# probe["H_excess"] 是**呈现口径**（幻影碱归零闸）。走步口径的真 He 需另取：
# 由实验记录 = -0.024007（`_presentation_He` 的 He_in）。此处两口径都算。
for tag, He_use in (("呈现 He=0.0", 0.0), ("走步 He=-0.024007", -0.024007)):
    print(f"\n{'=' * 74}\n### {tag}   ledger={ {k: round(v, 6) for k, v in led.items() if k != 'H_2O'} }")
    pH = eng.estimate_pH(led, He_use, V, T, T_K)
    print(f"  estimate_pH = {pH:.4f}")
    cands = enumerate_candidates(led, He_use, pH, V, T_K, T, True)
    print(f"  候选 {len(cands)} 条；以下为含 Ga 或含 H⁺ 且 |S|>0.1 者：")
    rows = []
    for cd in cands:
        try:
            S = eng.S_of(cd, led, V, pH, T_K, T, frozenset(), eng.P_EXT_KPA,
                         True, {})
        except Exception as exc:                            # noqa: BLE001
            rows.append((999.0, cd.kind, str(cd.r) + "->" + str(cd.pr),
                         f"S异常 {type(exc).__name__}: {exc}", None, None))
            continue
        ps = cd.pres_specs
        pres = (all(led.get(x, 0.0) > eng.X_MIN for x in ps[0])
                and all(led.get(x, 0.0) > eng.X_MIN for x in ps[1]))
        if not pres:
            continue
        if abs(S) < 0.1:
            continue
        nu_H = cd.pr.get(eng.H_ION, 0) - cd.r.get(eng.H_ION, 0)
        d = 1 if S > 0 else -1
        try:
            ext, x_max = eng.solve_extent(cd, d, led, He_use, V, T_K, T,
                                          frozenset())
        except Exception as exc:                            # noqa: BLE001
            ext, x_max = None, f"异常 {type(exc).__name__}: {exc}"
        rows.append((abs(S), cd.kind, f"{cd.r} -> {cd.pr}", S, ext, x_max))
        rows[-1] = rows[-1] + (nu_H, d)
    rows.sort(key=lambda t: -t[0])
    for ab, kind, eq, S, ext, x_max, nu_H, d in rows[:14]:
        _e = f"{ext:.6g}" if isinstance(ext, float) else str(ext)
        _x = f"{x_max:.6g}" if isinstance(x_max, float) else str(x_max)
        print(f"  |S|={ab:8.3f} d={d:+d} nu_H={nu_H:+d} {kind:9s} "
              f"ext={_e:>12s} x_max={_x:>12s}")
        print(f"      {eq[:96]}")

print("\n=== 判读 ===")
print("  H1 成立 ⟺ 反向步(镓酸根+4H⁺→Ga³⁺ / Ga³⁺→Ga(OH)₃+3H⁺) 的 x_max")
print("  远大于 'H_excess/|nu_H|'（酸碱预算）⟹ 预算根本没施加（非 redox 类）。")
print("  化学对比：反向步耗 4 mol H⁺/mol Ga，而账本游离 H⁺ 仅 0.024 mol；")
print("  镓酸根→固耗 1 mol H⁺/mol（Ga(OH)₄⁻+H⁺→Ga(OH)₃+H₂O），上界 0.743。")
