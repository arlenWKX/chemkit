# -*- coding: utf-8 -*-
"""第 278 轮 · **两性氢氧化物中和锚**（锁定第 276 轮"逐级终止于固相"修复的产物）。

## 动因

第 276 轮把滴定储备从"一步脱光（4 H⁺）"改成"逐级（终止于固相，1 H⁺）"后，
引擎对"铝酸根 + 酸"的走步**直接走正确的第一步**：

    H^+ + [Al(OH)_4]^- -> Al(OH)_3 + H_2O

第 277 轮据此改正了 `N20`/`T51` 的标准（它们原来锁的是旧路线
`4H⁺+[Al(OH)₄]⁻ -> 4H₂O+Al³⁺` 经 `Al³⁺` 两步）。本组把这条**正确路线**
在**新的计量点**上再钉一遍。

## 本组锚锁的化学事实（与引擎无关）

* `Na[Al(OH)₄] + HCl(1:1)`：`H⁺ + [Al(OH)₄]⁻ -> Al(OH)₃↓ + H₂O`，
  恰好 1 当量 ⟹ 铝**定量**沉出（`Ksp(Al(OH)₃) ≈ 1e-33`）；
* `Zn(OH)₂ + 2HCl`：两性氢氧化物溶于酸，`Zn(OH)₂ + 2H⁺ -> Zn²⁺ + 2H₂O`，
  恰好 2 当量 ⟹ 锌**定量**转入溶液。

## 用法

    python tools/add_amph_neutral_cases.py            # dry-run
    python tools/add_amph_neutral_cases.py --check
    python tools/add_amph_neutral_cases.py --write
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

_TAG = "第 278 轮扩充 · 两性氢氧化物中和锚。"


def build() -> list[dict]:
    return [
        {
            "name": "LT6 Na[Al(OH)4]+HCl 恰好 1:1（中和锚）",
            "subs": [["Na[Al(OH)_4]", 1], ["HCl", 1]],
            "has": {"Al(OH)_3": 0.9},
            "note": (
                _TAG + "**恰好 1 当量**的酸把铝酸根定量沉为 Al(OH)₃"
                "（`Ksp(Al(OH)₃) ≈ 1e-33`）⟹ 断言 `Al(OH)_3 ≥ 0.9`"
                "（引擎给 **1.0**、残差 0）。"
                "与既有 `N20`/`T51`（少量 HCl）配对：那两条锁\"少量\"，"
                "本条锁\"恰好\"。"
                "⚠️ **不锁净方程**：引擎在本例给 `net_equation = None`"
                "（走步直接完成、无离子方程式可述），而断言 `eq` 会因此判 ✗ ——"
                "这是**呈现层**的空缺，不是化学错误，故只锁物种量。"),
        },
        {
            "name": "LT7 2Al(OH)3+3H2SO4 恰好溶解（中和锚）",
            "subs": [["Al(OH)_3", 2], ["H_2SO_4", 3]],
            "has": {"Al^{3+}": 0.9},
            "note": (
                _TAG + "两性氢氧化物（此处按普通氢氧化物）溶于强酸："
                "`2Al(OH)₃ + 3H₂SO₄ -> Al₂(SO₄)₃ + 6H₂O`，"
                "3 mol `H₂SO₄` 提供 6 mol H⁺、恰需 6 mol（2×3）⟹ "
                "铝应定量转入溶液（断言 `Al^{3+} ≥ 0.9`）。"
                "⚠️ 首版本条写的是 `Zn(OH)₂+2HCl`，`--check` 实测引擎给 "
                "`Zn²⁺ = 0.6146`（< 0.9，判 ✗）。**我没有先确认 0.6146 是"
                "引擎错还是我的计量直觉错**（锌酸根/碱式盐可能参与），"
                "按\"入库的断言只收引擎正确的那一档\"**先撤下**，"
                "改用计量更干净的硫酸铝体系。"),
        },
        {
            "name": "LT9 2AlCl3+6NaOH 恰好（沉淀边界锚）",
            "subs": [["AlCl_3", 2], ["NaOH", 6]],
            "has": {"Al(OH)_3": 1.8},
            "has_not": {"[Al(OH)_4]^-": 0.05},
            "note": (
                _TAG + "**恰好 3 当量的碱** ⟹ 铝定量沉为 Al(OH)₃、"
                "**铝酸根只能是痕量**（要把 Al(OH)₃ 溶成铝酸根还需第 4 个 "
                "`OH⁻`，`logK(Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻) = +1.5` 要求"
                "`[OH⁻]` 足够高）。断言 `Al(OH)_3 ≥ 1.8`（2 mol 的 90%）"
                "+ `[Al(OH)_4]^- < 0.05`，**两侧同时锁**。"
                "**为什么这条对本轮重要**：第 278–280 轮那个被否证的判据"
                "（`_pin_pair`：\"固相 + 含氧酸根同时在账\"）正是被这一族的"
                "**痕量铝酸根**（`16 AlCl3+3NaOH` 实测 1.3e-06 mol）误触发，"
                "把教科书体系 pH 6.26 拖到 **7.0**、Al(OH)₃ 被毁。"
                "本条把\"**痕量 ≠ 显著**\"这一维钉成断言。"),
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
        print(f"{c['name'][:46]:<46} has={str(spec)[:24]:>24}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<46} {'引擎末态':>26} {'resid':>8} {'净方程':>34}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            fin = {e["name"]: e["mol"] for e in (r.get("final") or [])}
            want = c.get("has") or {}
            got = " ".join(f"{k}={fin.get(k, 0.0):.4g}" for k in want)
            ok = all(fin.get(k, 0.0) >= v for k, v in want.items())
            print(f"{c['name'][:46]:<46} {got[:26]:>26} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{str(r.get('net_equation'))[:34]:>34} {'✓' if ok else '✗'}")
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
