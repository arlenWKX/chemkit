# -*- coding: utf-8 -*-
"""第 265 轮 · **两性氢氧化物**族的定向锚（第 265 轮新诊断族）。

## 动因

第 265 轮把 `I31 InCl3+3NaOH` 的残差查到根：它是 §1.1 表示层族的成员
（同一迭代两个 pH 口径，见 log 第 265 轮）。该族的共同化学场景是
**两性氢氧化物在"恰好中和"附近**——产物是难溶氢氧化物，
溶液 pH 由「Ksp 钉住的游离金属」与「电荷平衡」共同要求，
而这两者在当前记账表示下不能同时满足。

本轮为该族补**答案可独立推导**的锚，作为将来表示层重构的回归网。
断言只用**化学事实 + 质量守恒**，不含引擎输出快照：

* 难溶氢氧化物的**定量性**：`M(OH)₃` 溶度积极小 ⟹ 恰好中和时金属
  基本全部成固（断言取下限 0.9×n(M)）；
* **两性**：强碱过量时氢氧化物以 `[M(OH)₄]⁻` 形式复溶 ⟹ 固相量下降
  （连续量，**只写区间、不写点值**——`Z31`/`SI1` 的教训）；
* **稀释不变性**：定量性由 Ksp 决定，与浓度无关（除非接近溶解度）。

## 用法

    python tools/add_amph_hydroxide_cases.py            # dry-run
    python tools/add_amph_hydroxide_cases.py --check    # 只测引擎，不调断言
    python tools/add_amph_hydroxide_cases.py --write
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

    def add(name, subs, note, has=None, has_range=None, eq=None,
            has_not=None):
        c = {"name": name, "subs": subs, "note": note, "changed": True}
        if has:
            c["has"] = has
        if has_range:
            c["has_range"] = has_range
        if has_not:
            c["has_not"] = has_not
        if eq:
            c["eq"] = eq
        out.append(c)

    # ---- In(OH)₃（本轮诊断的中心体系）----
    add("IN1 InCl3+3NaOH 恰好中和（两性氢氧化物）",
        [["InCl_3", 1], ["NaOH", 3]],
        "第 265 轮扩充 · **两性氢氧化物**族锚。恰好中和 ⟹ 产物 `In(OH)₃`；"
        "难溶 ⟹ 金属**定量**成固（断言取 0.9 下限，不用引擎输出）。"
        "本条是 `I31` 的化学基准：无论引擎内部走哪个 pH 口径，"
        "最终 `In(OH)₃` 都必须几乎全析出——这是**化学事实**，与算法无关。"
        "⚠️ **不写净方程断言**（第 265 轮自查修正）：恰好中和处固相已基本"
        "析完，净方程由痕量物种决定、且对 pKw 约定敏感 ⟹ 又一条"
        "「整数配比锁连续量」（`Z31`/`SI1` 同款）。",
        has={"In(OH)_3": 0.9})
    add("IN2 InCl3+4NaOH 碱过量（两性复溶）",
        [["InCl_3", 1], ["NaOH", 4]],
        "第 265 轮扩充 · 两性锚的**碱过量侧**：1 mol 过量 OH⁻ 把 `In(OH)₃` "
        "以 `[In(OH)₄]⁻`（铟酸根）**基本溶完** ⟹ 断言写成**上限**"
        "（`In(OH)₃` 必须很少），而不是下限。"
        "⚠️ 本条标准经**一次自查修正**：初版我写的是区间 [0.2, 0.95]"
        "（以为只溶一部分），实测引擎给 0.0040 —— 复核化学后确认"
        "**是我的标准错了**：`In(OH)₃ + OH⁻ -> [In(OH)₄]⁻` 在 1 当量过量下"
        "本就趋于完全（两性氢氧化物在强碱中溶解度大），故改为上限断言。"
        "与 IN1（恰好中和）配对锁住「碱用量」这一维。",
        has_not={"In(OH)_3": 0.2})

    # ⚠️ **0.1 M 稀档故意不入库**（第 265 轮实测）：该档锚定约定给 100% 成固、
    #    未锚定给 83.6%，**跨约定差 16 个百分点** ⟹ 落在某分支边界附近，
    #    不是稳健锚；而它想证的「难溶盐定量性与浓度无关」已由 IN1（1 mol 档）
    #    覆盖。**不把阈值从 0.09 调到 0.08**（那是为迁就引擎调标准）。

    # ---- Ga(OH)₃（同族第二个元素；库内 Ga 中间级有缺口，见 handoff §1.5）----
    # ⚠️ **GaCl3+3NaOH 这一档故意不入库**（第 265 轮实测）：引擎给
    #    `Ga(OH)₃ = 0`、pH **1.62**、残差 13.305 —— 与既有用例 `F31`
    #    同一条（`F31 GaCl3+3NaOH` 正是 §1.1 的招牌病例）。
    #    按既定纪律「入库的断言只收引擎正确的那一档，缺陷另行记账」，
    #    不把已知失败写成新锚（那只是重复 F31 的信息）。
    add("GA2 GaCl3+4NaOH 碱过量（两性复溶）",
        [["GaCl_3", 1], ["NaOH", 4]],
        "第 265 轮扩充 · Ga 的**碱过量侧**（同 IN2 的结构与标准形状）："
        "1 mol 过量 OH⁻ 把 `Ga(OH)₃` 以镓酸根基本溶完 ⟹ 上限断言。"
        "注意与 `F31`（GaCl₃+3NaOH **恰好**中和）的分工：那一档是已知的"
        "表示层缺陷（引擎给 pH 1.62），**不入库**；本条这一档引擎正确，"
        "故可入库，并同时锁住「碱过量侧」这一维。",
        has_not={"Ga(OH)_3": 0.2})

    # ---- Sc(OH)₃（第三元素；handoff §1.5 记 Sc³⁺ 库内无 OH⁻ 水解常数）----
    add("SC1 ScCl3+3NaOH 恰好中和（两性氢氧化物 · 数据缺口元素）",
        [["ScCl_3", 1], ["NaOH", 3]],
        "第 265 轮扩充 · `Sc(OH)₃` 难溶（Ksp 极小）⟹ 恰好中和时定量成固。"
        "⚠️ handoff §1.5 记 **`Sc³⁺` 库内无 OH⁻ 水解常数** ⟹ 这一条同时是"
        "**数据缺口的探针**：若断言通过，说明「难溶」这一条足以定产物，"
        "不必依赖缺失的一级水解常数（化学上合理：Ksp 已把游离金属钉死）；"
        "若不通过，就是缺数据的直接证据，按 §1.5 记账而非迁就引擎。",
        has={"Sc(OH)_3": 0.9},
        eq="Sc^{3+} + 3OH^- -> Sc(OH)_3")
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
        spec = c.get("has") or c.get("has_range") or {}
        print(f"{c['name'][:50]:<50} {str(spec)[:30]:>30}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<50} {'引擎末态':>30} {'resid':>8} {'判读':>6}")
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
            print(f"{c['name'][:50]:<50} {got[:30]:>30} "
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
