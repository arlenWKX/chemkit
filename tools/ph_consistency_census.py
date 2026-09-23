# -*- coding: utf-8 -*-
"""第 270 轮 · **呈现 pH 与「账本自身电荷平衡」的一致性普查**。

## 判据（只在**无固相**的终态上成立）

引擎不变量：`Σz·n(账本) + He = 0`
（`acidbase.ledger_charge` 的 docstring），而 `charge_pH` 按定义解
`V·h − V·oh + Σz·n(pH) = 0`。故对**完整账本**、**无 Ksp 储库在场**的状态，
`charge_pH` 给出的就是"这份账本在什么 pH 下电荷自洽"，可以直接与
**呈现层**给的 `final_pH` 对比。

⚠️ **两个必须排除的情形**（否则判据无效，第 263/265/268 轮的教训）：
① 终态账本里有固相 ⟹ pH 可能由 Ksp 储库决定，`charge_pH` 的 docstring
   明说此时它给的是"把账本拉回自洽"的 pH，**不是**真值 ⟹ 跳过；
② `charge_pH` 返回 `None`（无括号）⟹ 跳过。

**零侵入**（只读档 + 跑 `judge`），不改 `chemkit/` 一行。
输出 `logs/ph_consistency_census.json`。

用法：python tools/ph_consistency_census.py
      python tools/ph_consistency_census.py H45 ALU2
"""
from __future__ import annotations

import io
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.acidbase import charge_pH                      # noqa: E402
from chemkit.candidates import WATER, X_MIN                 # noqa: E402
from chemkit.data import load_tables                        # noqa: E402


def main(argv: list[str]) -> int:
    pres = tuple(a for a in argv if not a.startswith("-"))
    T = load_tables()
    solids = frozenset(T.solids)
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    from chemkit.engine import judge
    t0 = time.time()
    rows = []
    n_skip_solid = n_skip_none = 0
    for i, c in enumerate(cases, 1):
        nm = c["name"]
        if pres and not nm.startswith(pres):
            continue
        cond = c.get("cond") or {}
        V = float(cond.get("V_L", 1.0))
        T_K = float(cond.get("T_K", 298.15))
        pr = {}
        try:
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      cond, T, _probe=pr)
        except Exception:                                   # noqa: BLE001
            continue
        led = dict(pr.get("ledger") or {})
        if any(m > X_MIN and s != WATER and s in solids
               for s, m in led.items()):
            n_skip_solid += 1
            continue
        try:
            chg = charge_pH(led, V, T, T_K, fast=True)
        except Exception:                                   # noqa: BLE001
            chg = None
        if chg is None:
            n_skip_none += 1
            continue
        ph = pr.get("pH_solver")
        rows.append({
            "name": nm, "ph_solver": ph, "ph_chg": round(chg, 4),
            "d": (round(ph - chg, 4) if ph is not None else None),
            "ph_present": r.get("final_pH"),
            "resid": abs(pr.get("resid_live") or 0.0),
            "n_species": sum(1 for s, m in led.items()
                             if m > 1e-7 and s != WATER),
        })
        if i % 300 == 0:
            print(f"  {i}/{len(cases)} … {time.time() - t0:.0f}s", flush=True)
    bad = [r for r in rows if r["d"] is not None and abs(r["d"]) > 1.0]
    bad.sort(key=lambda r: -abs(r["d"]))
    print(f"\n可比对的例（无固相 + `charge_pH` 有解）：**{len(rows)}**"
          f"；跳过（固相在场）{n_skip_solid}、（无解）{n_skip_none}")
    print(f"其中 |pH_solver − charge_pH| > 1 的：**{len(bad)}** 例")
    print(f"\n{'用例':<40} {'solver':>8} {'charge':>8} {'Δ':>8} {'resid':>7}")
    for r in bad[:30]:
        print(f"   {r['name'][:40]:<40} {r['ph_solver']:>8} {r['ph_chg']:>8} "
              f"{r['d']:>+8.3f} {r['resid']:>7.3f}")
    out = {"n_compared": len(rows), "n_bad_gt1": len(bad),
           "skip_solid": n_skip_solid, "skip_none": n_skip_none,
           "rows": rows, "wall_s": round(time.time() - t0, 1)}
    with io.open(os.path.join(ROOT, "logs", "ph_consistency_census.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print("\n-> logs/ph_consistency_census.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
