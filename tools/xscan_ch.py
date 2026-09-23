# -*- coding: utf-8 -*-
"""第 265 轮 · 扫**指定通道**的 `S(x)`（`xscan.py` 只扫 |S| 最大者）。

需求来源：第 264 轮把 `I31 InCl3+3NaOH` 推进了"零推进"家族
（残差 0 -> 2.816）。它的**计入残差的**通道是
`In(OH)₃ + H₂O -> [In(OH)₄]⁻ + H⁺`（S=+2.816、x_max=0.002067、
ext=8e-10），而 `xscan.py` 扫的是 |S| 最大的那条
（`In³⁺ + 3H₂O -> In(OH)₃ + 3H⁺`，x_max 只有 1.06e-5，被
`ANN_MIN_EXTENT` 排除、不计入残差）。两条不是同一条 ⟹ 必须能指名扫。

用法：
    python tools/xscan_ch.py I31 "In(OH)_3 + H_2O"        # 子串匹配
    python tools/xscan_ch.py I31 "In(OH)_3 + H_2O" 20
"""
from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.candidates import H_ION, WATER                 # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def _eq(c, d) -> str:
    rr = c.r if d > 0 else c.pr
    pp = c.pr if d > 0 else c.r
    f = lambda dd: " + ".join(f"{'' if n == 1 else n}{s}"            # noqa: E731
                              for s, n in dd.items() if s != WATER)
    return f"{f(rr)} -> {f(pp)}"


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    pre = argv[0]
    needle = argv[1] if len(argv) > 1 else ""
    npts = int(argv[2]) if len(argv) > 2 else 16
    T = load_tables()
    hit = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
           if n.startswith(pre)]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    subs = [{"name": n, "mol": m} for n, m in c["subs"]]
    pr = {}
    r = eng.judge(subs, c.get("cond") or {"V_L": 1.0}, T, _probe=pr)
    led = dict(pr.get("ledger") or {})
    He = pr.get("H_excess") or 0.0
    pH = pr.get("pH_solver") or 7.0
    V = float((c.get("cond") or {}).get("V_L", 1.0))
    T_K = float((c.get("cond") or {}).get("T_K", 298.15))
    print(f"用例 {c['name']}   pH={r.get('final_pH')} solver={pH} He={He:.6g}")
    print(f"  账本={ {k: round(v, 8) for k, v in led.items() if k != WATER} }")
    src = pr.get("resid_src") or {}
    print(f"  引擎自报残差源: {src.get('eq')}  S={src.get('S')}  "
          f"ext_max={src.get('ext_max')}")

    cands = enumerate_candidates(led, He, pH, V, T_K, T, True)
    pick = None
    print("\n=== 候选（|S|>0.5，两侧在场）===")
    for cd in cands:
        ps = cd.pres_specs
        if not (all(led.get(x, 0.0) > eng.X_MIN for x in ps[0])
                and all(led.get(x, 0.0) > eng.X_MIN for x in ps[1])):
            continue
        try:
            S = eng.S_of(cd, led, V, pH, T_K, T, frozenset(),
                         eng.P_EXT_KPA, True, {})
        except Exception:                                       # noqa: BLE001
            continue
        if abs(S) < 0.5:
            continue
        d = 1 if S > 0 else -1
        rr = cd.r if d > 0 else cd.pr
        xm = min((led.get(s, 0.0) / n for s, n in rr.items()
                  if s not in (WATER, H_ION)), default=0.0)
        txt = _eq(cd, d)
        mark = ""
        if (needle and needle in txt
                and "[In(OH)_4]^-" in txt):
            mark = "   <== 目标"
        print(f"  S={S:+8.3f} d={d:+d} x_max={xm:<12.6g} kind={cd.kind:<8} "
              f"{txt}{mark}")
        if needle and needle in txt and pick is None:
            pick = (cd, d, S, xm)
    if pick is None:
        print("\n[未匹配到指定通道]")
        return 1
    cd, d, S0, x_max = pick
    print(f"\n=== 扫描 S(x)：{_eq(cd, d)}（d={d:+d}, x_max={x_max:g}）===")
    nu_H = cd.pr.get(H_ION, 0) - cd.r.get(H_ION, 0)
    if d < 0:
        nu_H = -nu_H
    changing = [(s, (cd.pr.get(s, 0) - cd.r.get(s, 0)) * d)
                for s in set(cd.r) | set(cd.pr) if s not in (WATER, H_ION)]
    print(f"  每单位 x 净变（正=生成）: "
          + ", ".join(f"{s}{v:+g}" for s, v in changing)
          + f"   nu_H={nu_H:+g}")
    print(f"\n{'x':>14} {'pH_est':>9} {'pH_chg':>9} {'Δ':>8} {'He':>12} "
          f"{'S@est':>9} {'S@chg':>9}")
    prev = None
    for k in range(npts):
        x = x_max * k / (npts - 1)
        led2 = dict(led)
        for s, v in changing:
            led2[s] = led2.get(s, 0.0) + v * x
        hx = He + nu_H * x
        ph2 = spec.estimate_pH(led2, hx, V, T, T_K)
        try:
            phc = spec.charge_pH(led2, V, T, T_K, fast=True)
        except Exception:                                       # noqa: BLE001
            phc = None
        s2 = eng.S_of(cd, led2, V, ph2, T_K, T, frozenset(), eng.P_EXT_KPA,
                      True, None) * d
        s3 = (eng.S_of(cd, led2, V, phc, T_K, T, frozenset(), eng.P_EXT_KPA,
                       True, None) * d) if phc is not None else float("nan")
        jump = ""
        if prev is not None and abs(ph2 - prev) > 0.05:
            jump = f"   <== estimate pH 跳 {ph2 - prev:+.3f}"
        dc = f"{phc - ph2:+.4f}" if phc is not None else "—"
        print(f"{x:>14.6g} {ph2:>9.4f} "
              f"{(f'{phc:.4f}' if phc is not None else '—'):>9} {dc:>8} "
              f"{hx:>12.5g} {s2:>9.4f} {s3:>9.4f}{jump}")
        prev = ph2
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
