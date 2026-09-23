# -*- coding: utf-8 -*-
"""第 268 轮 · 量 `h_c/o_c` 比值的分布 —— 为"数值不可分"判据找**有证据的**分界。

## 为什么要量

`F31 GaCl3+3NaOH`（`resid_max` 榜首 13.305）卡在 `_exact_ok` 第②条
（族成员 <2；金属羟合梯来自 `T.beta` 不进 `build_families`）。
**无条件停用②** 可修好 `F31`（13.305 → 0、pH 1.62 → 5.76）但**10 例翻红**
（`La1`/`Ce1`/`Nd1`/`Gd1`/`Ho1` 稀土水解 pH、`MgN1`/`CdN1`/`MgI1`/`MgD1`/`NiE1`
阴离子不变性）。故须只在"两侧启发式**数值上不可分**"时绕过 —— 但分界不能凭感觉定。

## 判据（不依赖新的 chemkit 诊断开关）

分支 4 的返回式是 `pH = -log10(h_c) if h_c >= o_c else pKw + log10(o_c)`：

* 标签 `酸侧max` ⟹ `h_c = 10^(-pH)`，而 `o_c = max(1e-7, -He_res/V)` 可算；
* 标签 `碱侧max` ⟹ 只能算 `o_c`，`h_c` 未知 ⟹ 跳过（本工具只统计 `酸侧max`）。

`F31` 属 `酸侧max` 一档，正是要分开的那一档。
**读取方式**：包 `speciation.estimate_state` 拿 `(pH, He_res)`，读回 `PH_TAGS`
（引擎自己的标注，D2：不复制判据）。

用法：python tools/tie_census.py            # 全库
输出 `logs/tie_census.json`
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

import chemkit.speciation as _sp                            # noqa: E402
from chemkit.data import load_tables                        # noqa: E402

_ORIG = _sp.estimate_state
ROWS: list = []
CASE = [""]
V0 = [1.0]


def _wrap(ledger, H_excess, V, T, T_K, *a, **k):
    save = _sp.PH_TAGS
    _sp.PH_TAGS = []
    try:
        r = _ORIG(ledger, H_excess, V, T, T_K, *a, **k)
        tags = list(_sp.PH_TAGS)
    finally:
        _sp.PH_TAGS = save
    tag = tags[-1] if tags else ""
    if tag == "酸侧max":
        pH, _vled, He_res = r
        try:
            h_c = 10.0 ** (-pH)
            o_c = max(1e-7, -He_res / V)
            if h_c > 0 and o_c > 0 and 0.1 * o_c <= h_c <= 10.0 * o_c:
                ROWS.append((CASE[0], h_c, o_c, abs(h_c - o_c) / max(h_c, o_c)))
        except Exception:                                   # noqa: BLE001
            pass
    return r


def main(argv: list[str]) -> int:
    pres = tuple(a for a in argv if not a.startswith("-"))
    _sp.estimate_state = _wrap
    T = load_tables()
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    from chemkit.engine import judge
    t0 = time.time()
    done = set()
    for i, c in enumerate(cases, 1):
        nm = c["name"]
        if pres and not nm.startswith(pres):
            continue
        CASE[0] = nm
        try:
            judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                  c.get("cond") or {}, T)
        except Exception:                                   # noqa: BLE001
            pass
        done.add(nm)
        if i % 150 == 0:
            print(f"  {i}/{len(cases)} … {time.time() - t0:.0f}s", flush=True)
    _sp.estimate_state = _ORIG

    rel = sorted(r[3] for r in ROWS)
    print(f"\n`酸侧max` 且落在 `_deg` 带内的采样点：**{len(ROWS)}** 次，"
          f"覆盖 {len(set(r[0] for r in ROWS))} 例")
    if rel:
        import statistics
        print(f"  相对差 |h_c−o_c|/max 的分位："
              f"p10={rel[len(rel)//10]:.3g} p50={rel[len(rel)//2]:.3g} "
              f"p90={rel[9*len(rel)//10]:.3g} max={rel[-1]:.3g}")
        for band in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
            n = sum(1 for x in rel if x <= band)
            cs = len(set(r[0] for r in ROWS if r[3] <= band))
            print(f"    相对差 ≤ {band:<7g} ：{n:>7} 次 / {cs:>4} 例")
    tgt = ("F31", "La1", "Ce1", "Nd1", "Gd1", "Ho1",
           "MgN1", "CdN1", "MgI1", "MgD1", "NiE1")
    print("\n关注用例的比值（F31 应「小」、其余十个应「大」）：")
    per: dict = {}
    for nm, h, o, d in ROWS:
        per.setdefault(nm, []).append(d)
    for pre in tgt:
        hit = [k for k in per if k.startswith(pre)]
        for k in hit:
            v = sorted(per[k])
            print(f"   {k[:40]:<40} n={len(v):>5}  "
                  f"min={v[0]:.3g}  p50={v[len(v)//2]:.3g}  max={v[-1]:.3g}")
    out = {"n": len(ROWS), "n_cases": len(set(r[0] for r in ROWS)),
           "per_case": {k: sorted(v) for k, v in per.items()},
           "wall_s": round(time.time() - t0, 1)}
    with io.open(os.path.join(ROOT, "logs", "tie_census.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print("\n-> logs/tie_census.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
