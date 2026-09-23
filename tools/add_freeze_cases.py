# -*- coding: utf-8 -*-
"""第 274 轮 · **冻结复核锚**（本轮"解冻不再有终身额度"的直接靶区）。

## 动因

第 274 轮查明：解冻（周期复核把仍是强驱动的冻结键作废）原实现有**终身额度**
（两个计数器各自全局封顶 12），`H35` 的 `CaCO_3 + 2H^+ -> Ca^{2+} + CO_2`
在强/弱之间来回、弱相位被冻结，额度**在状态长出 |S| = 11.094 之前就被别的键
花光** ⟹ 该通道永久冻结、走步不执行它、质量口径如实计入残差（X-38 病根一
"冻结在强驱动上属把问题藏起来"）。取消额度后全库质量不升反降、`n(|S|>0.1)`
196 → 182（详见 handoff §1.15）。

**本组锚锁的化学事实与引擎无关**，且都落在"本轮残差显著下降"的那批体系上
（把收益钉住，防将来悄悄回吐）：

* `Ksp(NiCO_3) ≈ 1.4e-7` ⟹ 碳酸氢盐能把 Ni²⁺ **定量**沉为 NiCO₃；
* `Ksp(Fe(OH)_3) ≈ 2.8e-39` ⟹ `Fe³⁺ + CO_3^{2-}` 双水解**定量**生成 Fe(OH)₃；
* `Ksp(Cr(OH)_3) ≈ 6.3e-31` ⟹ `Cr³⁺ + S^{2-}` 双水解**定量**生成 Cr(OH)₃
  （同时放出 H₂S，二者互相促进，方向无争议）。

断言一律取下限（计量的 90%）+ "游离离子基本耗尽"（上限 0.1），
**只用溶度积极小这一条化学事实**，不含引擎输出。

## 用法

    python tools/add_freeze_cases.py            # dry-run
    python tools/add_freeze_cases.py --check
    python tools/add_freeze_cases.py --write
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

_TAG = "第 274 轮扩充 · 冻结复核锚。"


def build() -> list[dict]:
    return [
        {
            "name": "FR1 NiCl2+2NaHCO3 定量沉镍（冻结复核锚）",
            "subs": [["NiCl_2", 1], ["NaHCO_3", 2]],
            "has": {"NiCO_3": 0.9},
            "has_not": {"Ni^{2+}": 0.1},
            "note": (
                _TAG + "`Ksp(NiCO_3) ≈ 1.4e-7` 极小 ⟹ 碳酸氢盐把 Ni²⁺"
                "**定量**沉为 NiCO₃（计量 `Ni^{2+} + 2HCO_3^- -> NiCO_3 + "
                "CO_2 + H_2O`）。断言 `NiCO_3 ≥ 0.9`、残余 `Ni^{2+} < 0.1`——"
                "只用\"溶度积极小\"这一化学事实。"
                "**为什么属于本轮**：既有用例 `Ni41`（同投料）本轮残差 "
                "4.528 → **0.000**，正是靠取消解冻终身额度把被永久冻结的"
                "通道放出来；本条把该结果钉成断言。"),
        },
        {
            "name": "FR2 2FeCl3+3Na2CO3 双水解定量（冻结复核锚）",
            "subs": [["FeCl_3", 2], ["Na_2CO_3", 3]],
            "has": {"Fe(OH)_3": 1.8},
            "has_not": {"Fe^{3+}": 0.1},
            "note": (
                _TAG + "`Ksp(Fe(OH)_3) ≈ 2.8e-39` ⟹ `2Fe^{3+} + 3CO_3^{2-} + "
                "3H_2O -> 2Fe(OH)_3 + 3CO_2` 双水解**定量**。"
                "断言 `Fe(OH)_3 ≥ 1.8`（2 mol 的 90%）、残余 `Fe^{3+} < 0.1`。"
                "与既有的 `E42 FeCl3+Na2CO3 双水解` 是**不同计量**的同族"
                "（那条本轮残差 1.402 → **0.000**）—— 本条锁计量侧。"),
        },
        {
            "name": "FR3 2CrCl3+3Na2S 双水解定量（冻结复核锚）",
            "subs": [["CrCl_3", 2], ["Na_2S", 3]],
            "has": {"Cr(OH)_3": 1.8},
            "has_not": {"Cr^{3+}": 0.1},
            "note": (
                _TAG + "`Ksp(Cr(OH)_3) ≈ 6.3e-31` ⟹ `2Cr^{3+} + 3S^{2-} + "
                "6H_2O -> 2Cr(OH)_3 + 3H_2S` 双水解**定量**（两个方向互相"
                "促进：Cr³⁺ 水解生酸、S²⁻ 水解生碱，彼此中和把反应推到底）。"
                "断言 `Cr(OH)_3 ≥ 1.8`、残余 `Cr^{3+} < 0.1`。"
                "既有的 `D41 CrCl3+Na2S 双水解` 本轮残差 2.551 → **0.503**；"
                "本条锁\"铬被定量沉出\"这一与走步细节无关的事实。"),
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
