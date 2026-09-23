# -*- coding: utf-8 -*-
"""第 276 轮 · **梯储备锚**（本轮"逐级终止于固相"修复的直接靶区）。

## 动因

`_buffer_titration` 的配离子碱储备原为**一步脱光到自由中心**
（`[M(OH)_k] + k·H⁺ -> M^{z+} + k·H₂O`），于是滴定器无法停在中途的
氢氧化物固相上。第 276 轮改成**逐级（终止于固相）**：

    [M(OH)_k] + (k−z)·H⁺ -> M(OH)_z(s) + (k−z)·H₂O
    logK = (k−z)·pKw + pKsp − logβ_k      （11 个配离子，全部库内可推）

## 本组锚锁的化学事实（与引擎无关）

* **两性氢氧化物溶于强碱**（铝酸根/锌酸根的生成）——`Al(OH)₃`、
  `Zn(OH)₂` 是教科书两性氢氧化物；
* **不溶于弱碱**：`Al(OH)₃` 不溶于氨水（`NH₃·H₂O` 碱强度不足以生成
  可观铝酸根）——这是定性分析的经典事实；
* **计量**：溶成 `[Al(OH)₄]⁻` 需要 **1 个** `OH⁻`（不是 3 个、更不是 4 个
  H⁺ 的逆过程）。

## 用法

    python tools/add_ladder_reserve_cases.py            # dry-run
    python tools/add_ladder_reserve_cases.py --check
    python tools/add_ladder_reserve_cases.py --write
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

_TAG = "第 276 轮扩充 · 梯储备锚。"


def build() -> list[dict]:
    return [
        {
            "name": "LT1 Al(OH)3+NaOH 恰好溶解（梯储备锚）",
            "subs": [["Al(OH)_3", 1], ["NaOH", 1]],
            "has": {"[Al(OH)_4]^-": 0.3},
            "note": (
                _TAG + "`Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻` 的 `logK = +1.5`"
                "（由 `β₄ = 10^34.5` 与 `pKsp(Al(OH)₃) = 33` 推出，"
                "`tools/ladder_reserve_census.py`），且**计量是 1 个 OH⁻**。"
                "1 mol Al(OH)₃ + 1 mol NaOH ⟹ 依 `x/(1−x) = 10^1.5` 得"
                "铝酸根 **0.969**（断言取保守下界 0.3）。"
                "**本条的靶心**：旧实现把 `[Al(OH)₄]⁻` 的质子储备当成"
                "`+4H⁺ -> Al³⁺`，滴定器因此看不见\"先停在固相\"这一步。"),
        },
        {
            "name": "LT2 Zn(OH)2+2NaOH 过量溶解（梯储备锚）",
            "subs": [["Zn(OH)_2", 1], ["NaOH", 2]],
            "has": {"[Zn(OH)_4]^{2-}": 0.3},
            "note": (
                _TAG + "`Zn(OH)₂ + 2OH⁻ -> [Zn(OH)₄]²⁻`（锌酸根）："
                "**两性氢氧化物在过量强碱中溶解**。"
                "`logK = 2·pKw + pKsp − logβ₄`（`tools/ladder_reserve_census.py`"
                " 给 27.200，即溶回 logK = 14×2 − 27.2 = **+0.8**）。"
                "2 mol OH⁻ 对 1 mol Zn(OH)₂ 恰是生成锌酸根所需 ⟹ 依"
                "`x/(2−x)²` 型平衡可得可观锌酸根（断言取保守下界 0.3）。"
                "与 `19 ZnCl2+2NaOH`（恰好 2 当量 ⟹ Zn(OH)₂ 定量）配对，"
                "把\"恰好\"与\"过量\"两侧都锁住。"),
        },
        {
            "name": "LT3 Al(OH)3+氨水 不溶（梯储备锚·负）",
            "subs": [["Al(OH)_3", 1], ["NH_3", 1]],
            "has_not": {"[Al(OH)_4]^-": 0.1},
            "note": (
                _TAG + "**负方向锚**：`Al(OH)₃` **不溶于氨水** —— "
                "`NH₃·H₂O` 的碱强度（`Kb = 1.8e-5`）远不足以把 `Al(OH)₃`"
                "溶成铝酸根（`logK(Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻) = +1.5`"
                " 需要 `[OH⁻]` 足够高）。这是定性分析的经典事实。"
                "断言 `[Al(OH)_4]^- < 0.1`。"
                "⚠️ 与既有用例 `U04 Al(OH)3 不溶氨水` 是**同一化学**；"
                "本条只做**上限**断言，不锁净方程形式（引擎给的是等价的"
                "2 倍配平式，锁形式会把表示差异当成化学差异）。"),
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
        print(f"{c['name'][:44]:<44} has={str(spec)[:26]:>26} "
              f"not={str(ban)[:20]:>20}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<44} {'引擎末态':>30} {'resid':>8} {'判读':>4}")
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
            print(f"{c['name'][:44]:<44} {got[:30]:>30} "
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
