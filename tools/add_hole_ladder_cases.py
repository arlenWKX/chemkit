# -*- coding: utf-8 -*-
"""第 269 轮 · **「ν 空洞」两性金属的恰好中和锚**。

## 动因

第 269 轮把 `build_families` 扩到认识 `T.beta` 羟合梯（原先只读 `T.pka`），
`F31 GaCl3+3NaOH` **FAIL → PASS**（残差 13.305 → 0、pH 1.62 → 5.76），
`resid_max` 13.305 → **9.434（首次进 ≤10）**。

`tools/ladder_gap_census.py` 普查出 **10 个中心的羟合梯有 ν 空洞**：
`Al³⁺`(缺3) / `Be²⁺` / `Cd²⁺` / `Cr³⁺` / `Ga³⁺` / `In³⁺` / `Zn²⁺`(缺2,3) /
`Pb²⁺`(缺2) / `Sn²⁺`(缺2) / `Sb³⁺`(缺1,2,3)。
本轮为其中**恰好中和**这一档补锚（断言只用"两性氢氧化物难溶 ⟹ 金属定量成固"
这一条化学事实）。

## 已有覆盖（不重复）

`ALU3` = `AlCl₃+3NaOH`；`IN1` = `InCl₃+3NaOH`；
`F31` = `GaCl₃+3NaOH`（本轮起通过，作为既有用例已覆盖）。

## 用法

    python tools/add_hole_ladder_cases.py            # dry-run
    python tools/add_hole_ladder_cases.py --check    # 只测引擎
    python tools/add_hole_ladder_cases.py --write
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

# (tag, 盐投料, mol, NaOH mol, 固相, 金属离子, ν 空洞说明)
ROWS = [
    ("ZN2", "ZnCl_2", 1.0, 2.0, "Zn(OH)_2", "Zn^{2+}", "缺 2,3"),
    ("CD2", "CdCl_2", 1.0, 2.0, "Cd(OH)_2", "Cd^{2+}", "缺 2,3"),
    ("CR3", "CrCl_3", 1.0, 3.0, "Cr(OH)_3", "Cr^{3+}", "缺 2,3"),
    ("BE2", "BeCl_2", 1.0, 2.0, "Be(OH)_2", "Be^{2+}", "缺 2,3"),
    ("PB2", "Pb(NO_3)_2", 1.0, 2.0, "Pb(OH)_2", "Pb^{2+}", "缺 2"),
]


def build() -> list[dict]:
    out = []
    for tag, salt, mol, naoh, solid, cat, hole in ROWS:
        out.append({
            "name": f"{tag} {salt}+{naoh:g}NaOH 恰好中和（ν 空洞元素）",
            "subs": [[salt, mol], ["NaOH", naoh]],
            "has": {solid: 0.9 * mol},
            "note": (
                f"第 269 轮扩充 · **羟合梯有 ν 空洞的元素**（`{cat}` {hole}）"
                f"在**恰好中和**这一档的两性氢氧化物定量锚。"
                f"`{solid}` 难溶 ⟹ 恰好按化学计量加碱时金属应**定量**成固"
                f"（断言取 0.9·n(M)，纯化学事实，不含引擎输出）。"
                f"空洞的意义（第 269 轮 `tools/ladder_gap_census.py` 普查）："
                f"`build_families` 只能连长出来的相邻 ν，缺口段的物种**进不了"
                f"质子化族** ⟹ 该段的酸碱再分配缺失。本条锁住这族元素的"
                f"**产物**（即使中间级数据缺失，'难溶 ⟹ 定量成固'仍必须成立）。"
                f"碱过量侧由 `AL4`/`ZN4`/`GA2` 覆盖，酸过量侧由 "
                f"`GA3`/`IN4`/`AL5` 覆盖。"),
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
        print(f"{c['name'][:52]:<52} {str(c['has'])[:20]:>20}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<52} {'引擎末态':>20} {'resid':>8} {'判读':>6}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            fin = {e["name"]: e["mol"] for e in (r.get("final") or [])}
            want = c.get("has") or {}
            got = " ".join(f"{k}={fin.get(k, 0.0):.4g}" for k in want)
            ok = all(fin.get(k, 0.0) >= v for k, v in want.items())
            print(f"{c['name'][:52]:<52} {got[:20]:>20} "
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
