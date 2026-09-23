# -*- coding: utf-8 -*-
"""第 266 轮 · **两性金属「碱过量侧」锚**（补第 265 轮 IN2/GA2 的 Al/Cr 对照）。

## 动因

第 266 轮的锋利普查（`tools/cliff_census.py`）量出"口径不一致"家族
**145/1313 例**，且**高度集中在两性金属体系**
（`Amp14`/`M01`/`H42`/`H43`/`16`/`N20`/`T51`/`N34`/`U01`/`E55`/`I31`/`IN1`…）。

**库内现有用例几乎全是"碱不足/恰好"那一侧**（那一侧正是缺陷所在，故断言难写稳）；
**"碱过量 ⟹ 两性氢氧化物以配阴离子复溶"这一侧覆盖不足**。
本轮为 Al / Cr 补这一侧，与第 265 轮的 `IN2`/`GA2` 同形。

## 答案来源（只用化学事实，不含引擎输出）

* 两性氢氧化物在**强碱过量**时以 `[M(OH)₄]⁻` 复溶
  （`Al(OH)₃ + OH⁻ -> [Al(OH)₄]⁻`；`Cr(OH)₃ + OH⁻ -> [Cr(OH)₄]⁻`）
  ⟹ 1 mol 过量 OH⁻ 足以把 1 mol 氢氧化物基本溶完
  ⟹ 断言写成**上限**（`has_not`），不是点值；
* 复溶程度由 Ksp 与 pH 共同决定，是**连续量** ⟹ 不写点值、不写净方程
  （`Z31`/`SI1`/`IN1` 的教训）。

## 用法

    python tools/add_amph_alkali_cases.py            # dry-run
    python tools/add_amph_alkali_cases.py --check    # 只测引擎，不调断言
    python tools/add_amph_alkali_cases.py --write
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

    def add(name, subs, note, has=None, has_range=None, has_not=None):
        c = {"name": name, "subs": subs, "note": note, "changed": True}
        if has:
            c["has"] = has
        if has_range:
            c["has_range"] = has_range
        if has_not:
            c["has_not"] = has_not
        out.append(c)

    add("AL4 AlCl3+4NaOH 碱过量（两性复溶）",
        [["AlCl_3", 1], ["NaOH", 4]],
        "第 266 轮扩充 · **两性金属「碱过量侧」锚**（第 266 轮普查量出的家族，"
        "库内原本几乎只有「碱不足/恰好」那一侧）。"
        "`Al(OH)₃` 两性：1 mol 过量 OH⁻ 把它以 `[Al(OH)₄]⁻`（偏铝酸根）"
        "基本溶完 ⟹ 断言写成**上限**（`has_not`），不是点值——"
        "复溶程度由 Ksp 与 pH 共同决定，是**连续量**"
        "（`Z31`/`SI1`/`IN1` 的教训）。与既有用例 `16 AlCl3+3NaOH`"
        "（恰好中和侧）配对，锁住「碱用量」这一维。",
        has_not={"Al(OH)_3": 0.2})
    # ⚠️ **`CrCl3+4NaOH` 故意不入库**（第 266 轮实测）：引擎给 `Cr(OH)₃ = 0.5124`
    #    （约一半复溶，pH 13.71）。初版我按 `Al` 的形状写 `Cr(OH)₃ < 0.2`，
    #    被实测否掉 —— 复核化学：**亚铬酸根的生成常数远小于偏铝酸根**
    #    （`Al(OH)₃+OH⁻` 约 +1，`Cr(OH)₃+OH⁻` 约 −0.4 量级），故"1 当量过量碱
    #    基本溶完"对 Cr **不成立**，部分复溶才是化学上合理的。
    #    **但我不入库、也不写区间**：写区间需要一个"部分"的具体边界，
    #    而那要引用 Cr 的生成常数 —— 我手上只有**记忆值**，
    #    按纪律「记忆值不得当已核数据」**不据此写标准**。
    #    该观测记账（见 log 第 266 轮），等查到出处再补锚。
    add("ZN4 ZnCl2+4NaOH 碱过量（两性复溶）",
        [["ZnCl_2", 1], ["NaOH", 4]],
        "第 266 轮扩充 · 同族的 Zn 对照（`Zn(OH)₂` 两性，"
        "过量碱生成锌酸根 `[Zn(OH)₄]²⁻`）。2:1 化学计量下 4 mol NaOH "
        "= 恰好中和 2 mol + 2 mol 过量 ⟹ 氢氧化物基本溶完 ⟹ 上限断言。",
        has_not={"Zn(OH)_2": 0.2})
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
        spec = c.get("has") or c.get("has_range") or c.get("has_not") or {}
        print(f"{c['name'][:46]:<46} {str(spec)[:26]:>26}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<46} {'引擎末态':>26} {'resid':>8} {'判读':>6}")
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
            print(f"{c['name'][:46]:<46} {got[:26]:>26} "
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
