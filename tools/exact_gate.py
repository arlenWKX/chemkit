# -*- coding: utf-8 -*-
"""第 268 轮 · 逐条检查 `_exact_ok` 的三条闸（哪一条挡下了精确解）。

`estimate_state` 在 `_deg`（h_c ≈ o_c 退化区）成立时会尝试用 `charge_pH`
的精确解；若 `_exact_ok` 任一条件不满足就退回 `max(h_c, o_c)` 启发式。
`F31 GaCl3+3NaOH` 实测 `h_c == o_c == 0.024006826`（九位相同，D13）
但标签是 `酸侧max` ⟹ 精确解没被采用 ⟹ 必须查出是哪一条闸。

用法：python tools/exact_gate.py F31 H45 I31
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
import chemkit.speciation as _sp                            # noqa: E402
from chemkit.acidbase import build_families, charge_pH      # noqa: E402
from chemkit.candidates import H_ION, WATER, X_MIN          # noqa: E402
from chemkit.core import charge_of                          # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def main(argv: list[str]) -> int:
    pres = argv or ["F31"]
    T = load_tables()
    cases = {c["name"]: c for c in load_cases(None)}
    for pre in pres:
        hit = [v for k, v in cases.items() if k.startswith(pre)]
        if not hit:
            print(f"[跳过] {pre}")
            continue
        c = hit[0]
        cond = c.get("cond") or {}
        V = float(cond.get("V_L", 1.0))
        pr = {}
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T,
                  _probe=pr)
        led = dict(pr.get("ledger") or {})
        He = pr.get("H_excess") or 0.0
        print(f"\n=== {c['name']} ===  pH={pr.get('pH_solver')} He={He:.9g}")
        fams = build_families(T)
        fam_keys = frozenset(fams)
        solids = frozenset(T.solids)
        n_sol = [(s, m) for s, m in led.items()
                 if m > X_MIN and s != WATER and s in solids]
        n_fam = [(s, m) for s, m in led.items()
                 if m > 0.0 and s in fam_keys]
        net = sum(charge_of(s) * m for s, m in led.items()
                  if m > 0.0 and s != WATER and not s.startswith("__"))
        print(f"  ① 固相在场（>X_MIN={X_MIN:g}）：{n_sol or '无'}  "
              f"⟹ {'**挡下**' if n_sol else '通过'}")
        print(f"  ② 族成员数（≥2）：{len(n_fam)}  {n_fam}  "
              f"⟹ {'通过' if len(n_fam) >= 2 else '**挡下**'}")
        print(f"  ③ |Σz·n + He_res| = {abs(net + He):.3g}  "
              f"⟹ {'通过' if abs(net + He) <= 1e-6 else '**挡下**'}")
        print(f"  `_exact_ok` 综合 = {_sp._exact_ok(led, He, T, V)}")
        # 退化区闸
        _sp.PH_SRC = []
        _sp.PH_TAGS = []
        try:
            ph, _v, He_res = _sp.estimate_state(led, He, V, T, 298.15)
            src, tags = list(_sp.PH_SRC), list(_sp.PH_TAGS)
        finally:
            _sp.PH_SRC = None
            _sp.PH_TAGS = None
        d = {s: (k, v) for s, k, v in src}
        h_c = d.get("__branch4__", ("", 0))[1]
        o_c = abs(He_res) / V
        print(f"  h_c={h_c:.9g}  o_c={o_c:.9g}  "
              f"比值 h_c/o_c={h_c / o_c if o_c else float('inf'):.9f}")
        print(f"  `_deg`（0.1 ≤ h_c/o_c ≤ 10）= "
              f"{bool(h_c > 0 and o_c > 0 and 0.1 * o_c <= h_c <= 10.0 * o_c)}")
        print(f"  estimate_state -> pH {ph:.6f}  tags={tags}")
        print(f"  charge_pH(无 pinned) = {charge_pH(led, V, T, 298.15)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
