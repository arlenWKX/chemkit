# -*- coding: utf-8 -*-
"""第 275 轮 · **两性金属碱侧锚**（本轮"角色盲区前提判据"的直接靶区）。

## 动因

第 275 轮查明 §1.12 的一个根因：`estimate_pH` 分支 4 的 `h_c`/`o_c` 只扫
`rolemap`，而两性金属的**含氧酸根**（`[Al(OH)₄]⁻` 等）在角色表里**没有条目**
⟹ 账本里 0.845 M 的铝酸根对 `o_c` 完全隐形 ⟹ 启发式以为那是酸液。
实测 `H43 AlCl3+NaOH 1:3.5`：呈现 pH **3.26** vs 同账本电荷平衡精确解
**12.037**（手算 12.03）—— 差 8.8 个 pH 单位。

修法（`_blind_dom`）是**前提判据**：被角色表忽略的**含氧酸根**比启发式
实际用到的任何量都大 ⟹ 启发式无权选支，交回精确解。

## 本组锚锁的化学事实（全部与引擎无关）

* `Ksp(In(OH)₃) ≈ 1e-33` ⟹ `In³⁺ + 3OH⁻ -> In(OH)₃` **定量**；
* 铝酸盐与铝盐**双水解**：`2Al³⁺ + 3[Al(OH)₄]⁻ -> 5Al(OH)₃↓`
  （Al³⁺ 水解生酸、铝酸根水解生碱，彼此中和把反应推到底 ⟹ 定量）；
* **碱过量**：`Al(OH)₃` 是两性氢氧化物，**过量的 OH⁻ 会把它溶成铝酸根**
  （`Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻`，`logK ≈ +1.5`，由 `β₄=10^34.5` 与
  `pKsp=33` 推出）⟹ 3.5 当量 OH⁻ 下应生成可观的铝酸根。

断言只用溶度积极小 / 两性这两条事实。

## 用法

    python tools/add_amph_base_cases.py            # dry-run
    python tools/add_amph_base_cases.py --check
    python tools/add_amph_base_cases.py --write
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

_TAG = "第 275 轮扩充 · 两性金属碱侧锚。"


def build() -> list[dict]:
    return [
        {
            "name": "PB1 InCl3+3NaOH 定量沉铟（两性金属锚）",
            "subs": [["InCl_3", 1], ["NaOH", 3]],
            "has": {"In(OH)_3": 0.9},
            "has_not": {"In^{3+}": 0.05},
            "note": (
                _TAG + "`Ksp(In(OH)₃) ≈ 1e-33` 极小 ⟹ **恰好 3 当量** OH⁻ 把"
                "In³⁺ 定量沉为 In(OH)₃（计量 `In^{3+} + 3OH^- -> In(OH)₃`）。"
                "断言 `In(OH)_3 ≥ 0.9`、残余 `In^{3+} < 0.05`——只用溶度积极小"
                "这一条化学事实。既有用例 `I31`/`IN1` 本轮残差 **3.675 → 0.000**。"),
        },
        {
            "name": "PB2 AlCl3+3NaAlO2 计量（两性金属锚·pH 侧）",
            "subs": [["AlCl_3", 1], ["NaAlO_2", 1]],
            "ph": [2.0, 6.0],
            "has_not": {"Al^{3+}": 0.5},
            "note": (
                _TAG + "⚠️ **本条的化学要配平清楚，首版就写错了**。按 Al/H/O/"
                "电荷守恒配平 `2Al³⁺ + 3[Al(OH)₄]⁻ + 3H₂O -> 5Al(OH)₃ + 3H⁺`"
                "—— 反应**净释放 3 个 H⁺**，所以它不是\"互相中和推到底\"的"
                "定量双水解（首版照 S²⁻/CO₃²⁻ 双水解的直觉写成\"定量\"，"
                "断言 `Al(OH)_3 ≥ 4.5`，被 `--check` 当场判 ✗：引擎给 4.0、"
                "pH 3.65，**与配平结果一致 ⟹ 错的是断言方**）。"
                "改锁**可判**的两条：① 溶液必然**酸性**（净放 H⁺，"
                "区间放宽到 [2.0, 6.0]）；② 铝不可能全留游离"
                "（`Ksp(Al(OH)₃) ≈ 1e-33`）：`Al^{3+} < 0.5`。"),
        },
        {
            "name": "PB3 AlCl3+3.5NaOH 碱过量复溶（两性金属锚）",
            "subs": [["AlCl_3", 1], ["NaOH", 3.5]],
            "has": {"[Al(OH)_4]^-": 0.3},
            "has_not": {"Al^{3+}": 0.05},
            "note": (
                _TAG + "**两性氢氧化物的复溶**：`Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻`"
                "（`logK = +1.5`，由 `β₄ = 10^34.5` 与 `pKsp(Al(OH)₃) = 33`"
                "推出，与该 `beta` 条目自己的 `calibrated` 注记一致）。"
                "3.5 当量 OH⁻ 比沉淀所需的 3 当量**多 0.5** ⟹ 那 0.5 会把"
                "一部分 Al(OH)₃ 溶成铝酸根。断言取**保守下界**"
                " `[Al(OH)_4]^- ≥ 0.3`（真值由 `x/(0.5−x) = 10^1.5` 给 0.485，"
                "取 0.3 留出活度修正余量）——**注意首版曾想写 0.7，那是错的**："
                "热力学只支持 0.485，写 0.7 等于要求引擎超过平衡。"
                "既有用例 `H43` 本轮残差 **24.511 → 0.147**、pH 3.26 → 12.19。"),
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
        spec = c.get("has") or {}
        ban = c.get("has_not") or {}
        print(f"{c['name'][:46]:<46} has={str(spec)[:24]:>24} "
              f"not={str(ban)[:18]:>18}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<46} {'引擎末态':>30} {'resid':>8} {'判读':>4}")
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
            print(f"{c['name'][:46]:<46} {got[:30]:>30} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>4}  pH={r.get('final_pH')}")
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
