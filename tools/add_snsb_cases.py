# -*- coding: utf-8 -*-
"""第 270 轮 · **补齐 ν 空洞元素「恰好中和」锚**（Sn / Sb）。

## 动因

第 269 轮普查出 **10 个中心**的羟合梯有 ν 空洞，并为其中 5 个
（Zn/Cd/Cr/Be/Pb）补了「恰好中和」锚；`Al`/`Ga`/`In` 已由既有用例覆盖。
本轮补最后两个：**`Sn²⁺`（ν=1,3，缺 2）** 与 **`Sb³⁺`（ν=4，缺 1,2,3）**。

断言仍只用一条化学事实：**两性氢氧化物难溶 ⟹ 恰好按化学计量加碱时
金属定量成固**（`has` 取下限 0.9·n(M)，不含引擎输出）。

⚠️ `Sb³⁺` 的梯只有 ν=4 一项（连成不了任何相邻对）⟹ 它在
`build_families` 里**完全不在族中**；本条锚正是用来盯住"数据缺口最大的
那一个"在**产物**上仍须正确。

## 用法

    python tools/add_snsb_cases.py            # dry-run
    python tools/add_snsb_cases.py --check    # 只测引擎
    python tools/add_snsb_cases.py --write
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

# (tag, 盐, mol, NaOH mol, 固相, 中心, ν 空洞)
ROWS = [
    ("SN2", "SnCl_2", 1.0, 2.0, "Sn(OH)_2", "Sn^{2+}", "缺 2"),
    ("SB3", "SbCl_3", 1.0, 3.0, "Sb(OH)_3", "Sb^{3+}", "缺 1,2,3"),
]


def build() -> list[dict]:
    out = []
    for tag, salt, mol, naoh, solid, cat, hole in ROWS:
        out.append({
            "name": f"{tag} {salt}+{naoh:g}NaOH 恰好中和（ν 空洞元素）",
            "subs": [[salt, mol], ["NaOH", naoh]],
            "has": {solid: 0.9 * mol},
            "note": (
                f"第 270 轮扩充 · **羟合梯 ν 空洞最大的一档**"
                f"（`{cat}` {hole}）在恰好中和时的两性氢氧化物定量锚。"
                f"`{solid}` 难溶 ⟹ 恰好按化学计量加碱时金属应**定量**成固"
                f"（断言取 0.9·n(M)，纯化学事实，不含引擎输出）。"
                f"空洞的意义（`tools/ladder_gap_census.py`）："
                f"`build_families` 只能连**相邻** ν，缺口段的物种进不了质子化族"
                f"——`Sb³⁺` 只有 ν=4 一项，连一个相邻对都凑不出，"
                f"是 10 个中心里缺口最大的。本条锁住**产物**：即使中间级数据"
                f"缺失、族连接为零，'难溶 ⟹ 定量成固'仍必须成立。"),
        })
    return out


def main(argv: list[str]) -> int:
    from chemkit.data import load_tables
    T = load_tables()
    new = build()
    path = os.path.join(ROOT, "chemkit", "data", "tests.json")
    with io.open(path, encoding="utf-8") as fh:
        db = json.load(fh)
    have = {c["name"] for c in db}
    for c in new:
        print(f"{c['name'][:50]:<50} {str(c['has'])[:20]:>20}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<50} {'引擎末态':>20} {'resid':>8} {'判读':>6}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            fin = {e["name"]: e["mol"] for e in (r.get("final") or [])}
            want = c.get("has") or {}
            got = " ".join(f"{k}={fin.get(k, 0.0):.4g}" for k in want)
            ok = all(fin.get(k, 0.0) >= v for k, v in want.items())
            print(f"{c['name'][:50]:<50} {got[:20]:>20} "
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
