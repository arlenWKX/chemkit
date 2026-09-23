# -*- coding: utf-8 -*-
"""第 272 轮 · **卤酸形态锚**（本轮落地的"共轭形态重建"机制的直接靶区）。

## 动因

第 272 轮查明 `S_of` 把 `ledger[s]/V` 当游离活度读，而账本是"引擎选定的
记账形态"：酸/碱对可以整池记在**另一侧**（`estimate_state` 的 `led_v` 把
`ClO_3^-` 整池搬成 `HClO_3`，`tools/x8gap.py` 实测），于是候选写作式里
的 `ClO_3^-` 槽位读成 0、落到 `ACT_FLOOR = 1e-12`，`logQ` 被凭空推低
22.5 个数量级 ⟹ 求解器 `f(0)≈0`、`x*=0`（零推进），而残差口径同态报
`|S| = 22.433`。修法是**按 Ka 在给定 pH 下从伙伴池重建游离分数**。

本组用例锁的化学事实**与引擎无关**，且正是该机制的靶区：**卤酸是强酸**
（`HClO_3`/`HBrO_3` pKa ≈ −1）⟹ 酸性液里卤酸根仍以接近分析浓度的
游离阴离子存在，**氧化性不因"记账成分子态"而消失**：

* `E°(ClO_3^-/Cl^-) = +1.45 V`、`E°(BrO_3^-/Br^-) = +1.44 V`（酸性）
  ≫ `E°(SO_4^{2-}/SO_3^{2-}) = +0.17 V`、`E°(I_2/I^-) = +0.54 V`
  ⟹ 卤酸根能把亚硫酸根、碘离子**定量**氧化；
* `E°(BrO_3^-/Br_2) = +1.52 V > E°(Br_2/Br^-) = +1.08 V`
  ⟹ 溴酸根 + 溴离子在酸中**归中**生成 Br₂（溴量法的发生反应）；
* 碱中 `3BrO^- -> 2Br^- + BrO_3^-`（次溴酸盐热力学不稳定，歧化完全）。

断言只用**氧化性强弱 + 计量（质量守恒）**两条，不含引擎输出；`has` 一律
取 90% 下限，`has_not` 锁"氧化剂基本耗尽"。

## 用法

    python tools/add_conjform_cases.py            # dry-run
    python tools/add_conjform_cases.py --check
    python tools/add_conjform_cases.py --write
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

_TAG = "第 272 轮扩充 · 卤酸形态锚。"


def build() -> list[dict]:
    return [
        {
            "name": "CH1 KClO3+3Na2SO3 酸性（卤酸形态锚）",
            "subs": [["KClO_3", 1], ["Na_2SO_3", 3], ["H_2SO_4", 1]],
            "has": {"SO_4^{2-}": 2.7},
            "has_not": {"ClO_3^-": 0.1},
            "note": (
                _TAG + "酸性液里氯酸根**仍以游离阴离子**参与氧化"
                "（HClO₃ pKa≈−1，是强酸）⟹ `ClO_3^- + 3SO_3^{2-} -> "
                "Cl^- + 3SO_4^{2-}` 定量（ΔE° = 1.45−0.17 = +1.28 V）。"
                "断言 `SO_4^{2-} ≥ 2.7`（3 mol 的 90%）+ 氯酸根耗尽 "
                "`< 0.1`。**这正是共轭形态重建的靶区**：账本在酸性下把"
                "氯酸记成 `HClO_3`，修前 `S_of` 读 `ClO_3^-` 槽位得 0、"
                "落到 ACT_FLOOR ⟹ 驱动力凭空消失。"),
        },
        {
            "name": "CH2 KBrO3+3Na2SO3 酸性（卤酸形态锚·强氧化剂）",
            "subs": [["KBrO_3", 1], ["Na_2SO_3", 3], ["H_2SO_4", 1]],
            "has": {"SO_4^{2-}": 2.7},
            "has_not": {"BrO_3^-": 0.1},
            "note": (
                _TAG + "溴酸根同理（HBrO₃ pKa≈−1）："
                "`BrO_3^- + 3SO_3^{2-} -> Br^- + 3SO_4^{2-}`，"
                "ΔE° = 1.44−0.17 = +1.27 V ⟹ 定量。与 `CH1` 配对，"
                "把「卤酸根在酸中不因记账形态失去氧化性」在 Cl/Br 两侧"
                "都锁住（两者走的是**同一条**共轭形态重建路径，"
                "但走步通道不同，互为对照）。"),
        },
        {
            "name": "CH3 NaClO3+3Na2SO3 中性（对照·无强酸）",
            "subs": [["NaClO_3", 1], ["Na_2SO_3", 3]],
            "has": {"SO_4^{2-}": 2.7},
            "has_not": {"ClO_3^-": 0.1},
            "note": (
                _TAG + "**对照组**：不加酸。氯酸根仍是强氧化剂"
                "（中性下 `ClO_3^- + 3SO_3^{2-} + 3H_2O -> Cl^- + "
                "3SO_4^{2-} + 6H^+` 的 ΔE° 仍 >1 V）⟹ 断言同 `CH1`。"
                "意义：把 `CH1`/`CH2` 的改善与\"pH 变了\"区分开——"
                "若只有酸性两例变好，说明机制只在 `led_v` 出现分子态时"
                "才起作用（符合预期）；本组三例都应定量。"),
        },
        {
            "name": "CH4 KBrO3+6KI 3H2SO4（溴量法滴定锚）",
            "subs": [["KBrO_3", 1], ["KI", 6], ["H_2SO_4", 3]],
            "has": {"I_2": 2.7},
            "has_not": {"BrO_3^-": 0.1},
            "note": (
                _TAG + "**溴量法的经典发生反应**（分析化学教科书事实）："
                "`BrO_3^- + 6I^- + 6H^+ -> Br^- + 3I_2 + 3H_2O`，"
                "ΔE° = 1.44−0.54 = +0.90 V ⟹ 6 mol I⁻ 被 1 mol 溴酸根"
                "**定量**氧化为 3 mol I₂。断言 `I_2 ≥ 2.7`（3 mol 的 90%）"
                "+ 溴酸根耗尽。本反应是容量分析里最标准的定量反应之一，"
                "化学事实无争议。"
                "**供酸量必须够**：反应要 6 H⁺，`H_2SO_4` 按 1:2 释出"
                "质子 ⟹ 取 3 mol（首版误取 1 mol，实测 `I_2 = 0.9999` —— "
                "恰好是 2/6 计量，**引擎是对的、断言方是错的**，"
                "这条错例本身反而验证了供质子与计量的耦合）。"),
        },
        {
            "name": "CH5 KBrO3+5KBr 3H2SO4（溴发生·归中锚）",
            "subs": [["KBrO_3", 1], ["KBr", 5], ["H_2SO_4", 3]],
            "has": {"Br_2": 2.7},
            "has_not": {"BrO_3^-": 0.1, "Br^-": 0.5},
            "note": (
                _TAG + "溴量法的**溴发生反应**："
                "`BrO_3^- + 5Br^- + 6H^+ -> 3Br_2 + 3H_2O`。"
                "方向判据：`E°(BrO_3^-/Br_2) = +1.52 V > "
                "E°(Br_2/Br^-) = +1.08 V` ⟹ 归中方向成立（反之不成立）。"
                "计量：1 mol 溴酸根 + 5 mol 溴离子 ⟹ **3 mol Br₂**，"
                "断言 `Br_2 ≥ 2.7`，且两个反应物都基本耗尽。"
                "本反应把 Br 的 0/+5 两端在**同一体系**里锁住。"
                "供酸同为 3 mol `H_2SO_4`（= 6 H⁺）；首版取 1 mol 实测"
                "`Br_2 = 0.9989、BrO_3^- = 0.6663`，正是 2/6 计量，"
                "再次印证**引擎按质子限量正确截断**。"),
        },
        {
            "name": "CH6 Br2+6NaOH 歧化（次溴酸盐不稳定）",
            "subs": [["Br_2", 3], ["NaOH", 6]],
            "has": {"BrO_3^-": 0.9, "Br^-": 4.5},
            "note": (
                _TAG + "碱中溴的**歧化**：先生成次溴酸盐"
                "`Br_2 + 2OH^- -> Br^- + BrO^- + H_2O`，"
                "而 `3BrO^- -> 2Br^- + BrO_3^-` 在热力学上完全"
                "（`E°(BrO^-/Br^-) = +0.76 V > E°(BrO_3^-/BrO^-) = "
                "+0.54 V`）⟹ 总反应 `3Br_2 + 6OH^- -> 5Br^- + "
                "BrO_3^- + 3H_2O`，断言 `Br^- ≥ 4.5`、`BrO_3^- ≥ 0.9`"
                "（各取计量的 90%）。"
                "**注意**：冷稀碱因**动力学**停在次溴酸盐，本引擎是"
                "热力学引擎、无动力学，故锁定平衡侧产物（BrO₃⁻）——"
                "这是化学事实的平衡表述，不是拟合引擎。"),
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
        print(f"{c['name'][:46]:<46} has={str(spec)[:26]:>26} "
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
