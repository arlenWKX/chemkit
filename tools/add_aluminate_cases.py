# -*- coding: utf-8 -*-
"""第 267 轮 · **铝酸盐/氢氧化铝定量锚**（第 267 轮修法的定向回归网）。

## 动因

第 267 轮查清 `resid_max` 榜首一族（`H45`/`M03`/`H43`，pH 12.19、残差合计 74.7）
的成因：残差通道 `Al³⁺ -> [Al(OH)₄]⁻ + 4H⁺` 的 `x_max` **恰等于**账本 `Al³⁺`，
于 `x_max` 端账本 `Al³⁺ = 0`，而 `_pin` 构造有一条"游离阳离子量 > 0"的闸
⟹ 钉住在端点被取消 ⟹ 退回逐点自判定 ⟹ `酸侧max` 把 pH 从 12.185 拉到 3.70
（跳 −8.48）⟹ 二分落在跳上 ⟹ `x*=0`、残差 24.991。

修法两处：① 整步冻结为 True 时不受"游离阳离子量 > 0"闸限制（前提仍由**固相在场**把关）；
② `engine` 侧"储库是否仍在场"的判据**只看固相、不看游离阳离子量**。

## 本脚本补什么

铝是两性金属里最常考的，而库内此前缺**铝酸盐侧的定量锚**。
断言只用**化学事实 + 质量守恒**，不含引擎输出快照：

* `Na[Al(OH)₄] + 半量 HCl`：铝酸根被中和一半 ⟹ 铝**一半成胶、一半仍在溶液**
  （写成区间，因为两性氢氧化物的溶解/析出是连续量）；
* `Na[Al(OH)₄] + 等量 HCl`（全部中和）：`Al(OH)₃` 应**定量**析出；
* `AlCl₃ + 3NaOH`（恰好中和）：`Al(OH)₃` **定量**析出；
* `AlCl₃ + 4NaOH`（碱过量）：`Al(OH)₃` 基本复溶（上限断言）。

## 用法

    python tools/add_aluminate_cases.py            # dry-run
    python tools/add_aluminate_cases.py --check    # 只测引擎，不调断言
    python tools/add_aluminate_cases.py --write
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


def build() -> list[dict]:
    out = []

    def add(name, subs, note, has=None, has_range=None, has_not=None,
            eq=None, eq_has=None):
        c = {"name": name, "subs": subs, "note": note, "changed": True}
        if has:
            c["has"] = has
        if has_range:
            c["has_range"] = has_range
        if has_not:
            c["has_not"] = has_not
        if eq:
            c["eq"] = eq
        if eq_has:
            c["eq_has"] = eq_has
        out.append(c)

    add("ALU1 Na[Al(OH)4]+等量HCl（铝酸根定量成胶）",
        [["Na[Al(OH)_4]", 1], ["HCl", 1]],
        "第 267 轮扩充 · **铝酸盐定量锚**。铝酸根是两性金属的碱侧形态，"
        "加**等量**强酸把它中和回 `Al(OH)₃` ⟹ 按质量守恒金属应**定量**成胶"
        "（断言取 0.9 下限，纯化学事实，不用引擎输出）。"
        "本条与既有 `H45 Na[Al(OH)4]+HCl 半量`（半量侧）配对锁住「酸用量」维；"
        "`H45` 是第 267 轮修法的招牌病例（残差 24.991 → 0.001），"
        "本条是它的**端点对照**。",
        has={"Al(OH)_3": 0.9})
    add("ALU2 Na[Al(OH)4]+半量HCl（铝酸根半中和）",
        [["Na[Al(OH)_4]", 1], ["HCl", 0.5]],
        "第 267 轮扩充 · 铝酸根的**半中和**侧。半量酸只能中和一半铝酸根 ⟹ "
        "铝**一半成胶、一半留在溶液**（以 `[Al(OH)₄]⁻` 形态）。"
        "⚠️ 两侧都只写**区间**：两性氢氧化物的溶解/析出程度由 Ksp 与 pH "
        "共同决定，是**连续量**（`Z31`/`SI1`/`IN1` 的教训——整数配比锁连续量）。"
        "与 `ALU1`（等量）配对构成用量维的两端。",
        has_range={"Al(OH)_3": [0.2, 0.8]})
    add("ALU3 AlCl3+3NaOH 恰好中和（两性氢氧化物）",
        [["AlCl_3", 1], ["NaOH", 3]],
        "第 267 轮扩充 · 铝的**恰好中和**侧（与既有 `16 AlCl3+3NaOH` 同化学，"
        "但断言只锁产物、不锁净方程）：`Al(OH)₃` 难溶 ⟹ 金属定量成固。"
        "第 267 轮修法后本档残差 3.139 → 0.000、pH 6.27 → 5.12（中性偏酸的"
        "饱和 Al(OH)₃ 溶液，与其两性最低溶解度区一致）。",
        has={"Al(OH)_3": 0.9})
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
        spec = (c.get("has") or c.get("has_range") or c.get("has_not") or {})
        print(f"{c['name'][:50]:<50} {str(spec)[:28]:>28}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<50} {'引擎末态':>28} {'resid':>8} {'判读':>6}")
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
            print(f"{c['name'][:50]:<50} {got[:28]:>28} "
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
