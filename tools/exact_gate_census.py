# -*- coding: utf-8 -*-
"""第 268 轮 · 普查"`_deg` 成立但 `_exact_ok` 被第②条挡下"的规模。

## 背景

`F31 GaCl3+3NaOH`（当时 `resid_max` 榜首 13.305）实测：

    h_c = o_c = 0.024006826（九位相同，D13） ⟹ `_deg` **True**
    ① 固相在场      通过
    ② 族成员数 ≥2   **0**  ⟹ **挡下**
    ③ Σz·n + He ≈ 0  通过
    ⟹ 退回 `max(h_c, o_c)`，`h_c` 以**相对 1.9e-8** 胜出 ⟹ pH 1.62（而非 12.38）

**根因**：`build_families` 只从 `T.pka` 建族，而金属的羟合梯来自 `T.beta`（logβ）
⟹ 这些物种在 `charge_pH` 里被当作**固定电荷**、不参与再分配 ⟹
第②条把它们数成 0 ⟹ 精确解永不启用。

本工具量这个缺口的**规模**（多少例、多少次），据此决定是否值得把
`beta` 羟合梯并入 `build_families`。

## 口径

包 `speciation._exact_ok`：每次返回 False 时重算三条闸，归类是哪一条挡的。
零侵入（只包模块级函数），不改 `chemkit/` 一行。串行跑（~2 min），
运行期间**冻结 chemkit/**。输出 `logs/exact_gate_census.json`。

用法：python tools/exact_gate_census.py            # 全库
      python tools/exact_gate_census.py F31 H45    # 只看这些前缀
"""
from __future__ import annotations

import collections
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

import chemkit.speciation as _sp                            # noqa: E402
from chemkit.acidbase import build_families                 # noqa: E402
from chemkit.candidates import WATER, X_MIN                 # noqa: E402
from chemkit.core import charge_of                          # noqa: E402
from chemkit.data import load_tables                        # noqa: E402

_ORIG = _sp._exact_ok
REASON = collections.Counter()
PER = collections.Counter()
CALLS = [0]
CASE = [""]
T = None
FAMS: frozenset = frozenset()


def _hook(ledger, He_res, Tt, V):
    r = _ORIG(ledger, He_res, Tt, V)
    CALLS[0] += 1
    if not r:
        solids = Tt.solids
        n_sol = any(m > X_MIN and s != WATER and s in solids
                    for s, m in ledger.items())
        n_fam = sum(1 for s, m in ledger.items()
                    if m > 0.0 and s in FAMS)
        net = sum(charge_of(s) * m for s, m in ledger.items()
                  if m > 0.0 and s != WATER and not s.startswith("__"))
        bad = []
        if n_sol:
            bad.append("①固相在场")
        if n_fam < 2:
            bad.append(f"②族成员<2(={n_fam})")
        if abs(net + He_res) > 1e-6:
            bad.append("③电荷不变量")
        key = "+".join(bad) or "其它"
        REASON[key] += 1
        PER[(CASE[0], key)] += 1
    return r


def main(argv: list[str]) -> int:
    global FAMS
    pres = tuple(a for a in argv if not a.startswith("-"))
    T = load_tables()
    FAMS = frozenset(build_families(T))
    _sp._exact_ok = _hook
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    from chemkit.engine import judge
    t0 = time.time()
    rows = []
    for i, c in enumerate(cases, 1):
        nm = c["name"]
        if pres and not nm.startswith(pres):
            continue
        CASE[0] = nm
        before = sum(v for (kk, _r), v in PER.items() if kk == nm)
        try:
            judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                  c.get("cond") or {}, T)
        except Exception:                                   # noqa: BLE001
            pass
        after = sum(v for (kk, _r), v in PER.items() if kk == nm)
        rows.append({"name": nm, "n_blocked": after - before})
        if i % 150 == 0:
            print(f"  {i}/{len(cases)} … {time.time() - t0:.0f}s", flush=True)
    _sp._exact_ok = _ORIG
    hot = sorted((r for r in rows if r["n_blocked"]), key=lambda r: -r["n_blocked"])
    out = {"exact_ok_calls": CALLS[0],
           "blocked_total": sum(REASON.values()),
           "reasons": dict(REASON.most_common()),
           "n_cases_blocked": len(hot), "n_cases": len(rows),
           "hot": hot, "wall_s": round(time.time() - t0, 1)}
    with io.open(os.path.join(ROOT, "logs", "exact_gate_census.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(f"\n`_exact_ok` 调用 {CALLS[0]} 次；返回 False "
          f"**{out['blocked_total']}** 次")
    print("\n挡下原因分布：")
    for k, v in REASON.most_common(10):
        print(f"   {v:>7}  {k}")
    print(f"\n落在 {out['n_cases_blocked']} / {out['n_cases']} 例上")
    print("\n被挡最多的 25 例（**其中含②的**即为 build_families 缺口）：")
    for r in hot[:25]:
        rs = [k for (kk, k) in PER if kk == r["name"]]
        print(f"   {r['name'][:40]:<40} {r['n_blocked']:>6}  "
              f"{' / '.join(sorted(set(rs)))[:52]}")
    print("\n-> logs/exact_gate_census.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
