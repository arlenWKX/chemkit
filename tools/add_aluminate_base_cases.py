# -*- coding: utf-8 -*-
"""第 277 轮 · **铝酸根碱侧锚**（本轮"角色注入被重复 `bases_map` 清空"修复的靶区）。

## 动因

第 276 轮把"含氧酸根 → 碱侧角色"写进了 `bases_map`，但**写了两遍
`bases_map: dict = {}`**，后一遍把前一遍整块清掉 ⟹ 角色**一条也没生效**。
`tools/branch4_src.py`（本轮改成直接读引擎真正的 `rolemap`）实测
`[Al(OH)_4]^-` 不在 182 项的 `rolemap` 里 ⟹ 分支 4 看不见 0.5 M 铝酸根
⟹ `H45` 那一路给 pH 7.0（纯水默认）。删掉重复定义后 `rolemap` 182 → 193 项。

## 本组锚锁的化学事实（与引擎无关）

* **铝酸根是弱碱**：`[Al(OH)₄]⁻ ⇌ Al(OH)₃(s) + OH⁻`，`Kb = 1/(β₄·Ksp)
  = 10^(pKsp − logβ₄) = 10^(33 − 34.5) = 0.0316`（`tools/ladder_base_census.py`）
  ⟹ `c = 0.1 M` 的铝酸根溶液 `[OH⁻] = (−Kb+√(Kb²+4Kb·c))/2 = 0.0426`
  ⟹ pOH 1.37 ⟹ **pH 12.63**；
* **两性氢氧化物溶于过量强碱**：`Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻`，`logK = +1.5`
  ⟹ 4 当量 `NaOH` 下 `Al(OH)₃` **全部**溶解。

## 用法

    python tools/add_aluminate_base_cases.py            # dry-run
    python tools/add_aluminate_base_cases.py --check
    python tools/add_aluminate_base_cases.py --write
"""
from __future__ import annotations

import io
import json
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

_TAG = "第 277 轮扩充 · 铝酸根碱侧锚。"


def build() -> list[dict]:
    return [
        {
            "name": "LT4 NaAlO2 0.1M 溶液 pH（铝酸根碱侧锚）",
            "subs": [["NaAlO_2", 0.1]],
            "ph": [12.2, 13.0],
            "has_not": {"Al(OH)_3": 0.01},
            "note": (
                _TAG + "**铝酸根是弱碱**：`[Al(OH)₄]⁻ ⇌ Al(OH)₃(s) + OH⁻`，"
                "`Kb = 1/(β₄·Ksp) = 10^(pKsp − logβ₄) = 10^(33 − 34.5) = 0.0316`"
                "（`tools/ladder_base_census.py`）⟹ 0.1 M 时"
                "`[OH⁻] = (−Kb+√(Kb²+4Kb·c))/2 = 0.0426` ⟹ pOH 1.37 ⟹ "
                "**pH 12.63**（区间取 ±0.4）。同时锁 `Al(OH)_3` 不得显著析出"
                "（0.1 M 铝酸根远未饱和）。"
                "**本条正是本轮修复的靶心**：角色注入被清空前，"
                "分支 4 看不见铝酸根 ⟹ 这种溶液会被算成 pH 7。"),
        },
        {
            "name": "LT5 Al(OH)3+4NaOH 过量复溶（铝酸根碱侧锚）",
            "subs": [["Al(OH)_3", 1], ["NaOH", 4]],
            "has": {"[Al(OH)_4]^-": 0.9},
            "note": (
                _TAG + "**两性氢氧化物在过量强碱中完全溶解**："
                "`Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻`（`logK = +1.5`）。"
                "4 当量 OH⁻ 远超沉淀所需的 3 当量 ⟹ 依 `x/(4−x) = 10^1.5` "
                "得 x = 3.88 > 1 ⟹ **1 mol Al(OH)₃ 全部溶解**"
                "（断言 `[Al(OH)_4]^- ≥ 0.9`）。"
                "与 `LT1`（1:1 恰好）配对，把\"恰好\"与\"过量\"两侧都锁住。"),
        },
    ]


def main(argv: list[str]) -> int:
    from chemkit.data import load_tables
    T = load_tables()
    new = build()
    path = os.path.join(ROOT, "chemkit", "data", "tests.json")
    with io.open(path, encoding="utf-8") as fh:
        db = json.load(fh)
    have = {c["name"] for c in db}
    for c in new:
        spec = c.get("has") or ({"ph": c["ph"]} if "ph" in c else {})
        ban = c.get("has_not") or {}
        print(f"{c['name'][:46]:<46} has={str(spec)[:24]:>24} "
              f"not={str(ban)[:20]:>20}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<46} {'引擎末态':>28} {'resid':>8} {'判读':>4}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            fin = {e["name"]: e["mol"] for e in (r.get("final") or [])}
            want = c.get("has") or {}
            ban = c.get("has_not") or {}
            got = " ".join(f"{k}={fin.get(k, 0.0):.4g}"
                           for k in list(want) + list(ban))
            ok = (all(fin.get(k, 0.0) >= v for k, v in want.items())
                  and all(fin.get(k, 0.0) < v for k, v in ban.items()))
            ph = r.get("final_pH")
            if "ph" in c:
                ok = ok and (c["ph"][0] <= (ph if ph is not None else -99)
                             <= c["ph"][1])
            print(f"{c['name'][:46]:<46} {got[:28]:>28} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>4}  pH={ph}")
        return 0

    if "--write" not in argv:
        print(f"\n[dry-run] 将新增 {sum(1 for c in new if c['name'] not in have)} 条"
              f"；加 --write 写入。")
        return 0
    add = [c for c in new if c["name"] not in have]
    db.extend(add)
    with io.open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(db, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"\n[写入] 新增 {len(add)} 条 ⟹ 共 {len(db)} 条")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
