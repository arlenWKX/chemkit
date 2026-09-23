# -*- coding: utf-8 -*-
"""第 268 轮 · **Ga 的「酸过量侧」锚** —— 把 Ga 数据缺口的边界围起来。

## 动因

第 268 轮把 `F31 GaCl3+3NaOH`（`resid_max` 13.305）查到根：**数据缺口**。
`Ga³⁺` 的羟合梯在库内只有 ν=1 与 ν=4（缺 ν=2,3）⟹ 梯子断成两段
⟹ `build_families` 连不成一条可再分配的族 ⟹ `_exact_ok` 第②条把精确解挡下
⟹ 两个启发式以**相对 1e-6** 分出胜负 ⟹ pH 落在酸侧 1.62（真值区 ~5.6）。

独立核查（`tools/ga_ladder_check.py`，只用库内 logβ 与 pKw）：
* 梯子**断开**时电荷平衡解出 **pH 1.88**（与引擎的 1.62 同侧）；
* 梯子**补全**后解出 **pH 5.61**（Ga(OH)₃ 两性最低溶解度区），
  且对中间级 ±0.5 个对数单位只敏感 ±0.18 pH ⟹ 结论稳健。

## 本脚本补什么

Ga 的三档里，**碱过量侧**已由第 266 轮的 `GA2` 覆盖、**恰好中和侧**是
已知缺陷（`F31`，按纪律「只收引擎正确的那一档」**不入库**）。
本轮补**酸过量侧**：强酸大大过量时 `Ga(OH)₃` 不可能存在（溶度积要求
`[OH⁻]³` 够大，而强酸把 pH 压到 ≪0）—— 断言用**上限**，纯化学事实。

## 用法

    python tools/add_ga_acid_cases.py            # dry-run
    python tools/add_ga_acid_cases.py --check    # 只测引擎
    python tools/add_ga_acid_cases.py --write
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


def _acid(tag: str, salt: str, solid: str, elem: str) -> dict:
    return {
        "name": f"{tag} {salt}+浓盐酸过量（酸过量侧）",
        "subs": [[salt, 1], ["HCl", 6]],
        "has_not": {solid: 0.01},
        "note": (
            f"第 268 轮扩充 · **{elem} 的「酸过量侧」锚**。`{solid}` 难溶，"
            "其存在要求 `[OH⁻]³` 足够大；强酸大大过量时 pH 被压到远低于 0，"
            f"`{solid}` **不可能析出** ⟹ 断言写成**上限**（纯化学事实，"
            "不含引擎输出）。本条与同族的「碱过量侧」"
            f"（`{'GA2' if elem == 'Ga' else 'AL4'}` 型）以及"
            "「恰好中和侧」（`ALU3`/`IN1`）一起，把两性氢氧化物"
            "**三档用量**都围起来；Ga 的恰好中和侧因**库内缺 ν=2,3 "
            "水解常数**而失败（`F31`，见 `tools/ga_ladder_check.py`），"
            "按纪律不写成锚。"),
    }


def build() -> list[dict]:
    return [
        _acid("GA3", "GaCl_3", "Ga(OH)_3", "Ga"),
        _acid("IN4", "InCl_3", "In(OH)_3", "In"),
        _acid("AL5", "AlCl_3", "Al(OH)_3", "Al"),
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
        spec = (c.get("has") or c.get("has_range") or c.get("has_not") or {})
        print(f"{c['name'][:46]:<46} {str(spec)[:22]:>22}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<46} {'引擎末态':>22} {'resid':>8} {'判读':>6}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            fin = {e["name"]: e["mol"] for e in (r.get("final") or [])}
            want = c.get("has") or {}
            rng = c.get("has_range") or {}
            ban = c.get("has_not") or {}
            keys = list(want) + list(rng) + list(ban)
            got = " ".join(f"{k}={fin.get(k, 0.0):.4g}" for k in keys)
            ok = (all(fin.get(k, 0.0) >= v for k, v in want.items())
                  and all(lo <= fin.get(k, 0.0) <= hi
                          for k, (lo, hi) in rng.items())
                  and all(fin.get(k, 0.0) < v for k, v in ban.items()))
            print(f"{c['name'][:46]:<46} {got[:22]:>22} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>6}  pH={r.get('final_pH')}")
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
