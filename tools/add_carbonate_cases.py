# -*- coding: utf-8 -*-
"""第 273 轮 · **碳酸族锚**（本轮"族级准入判据"的直接靶区）。

## 动因

第 273 轮把第 272 轮的"共轭形态重建"准入判据从**伙伴级**推广到**族级**：
`HCO_3^-` 与 `CO_3^{2-}` 在同一个酸碱族里，`CO_2 -> CO_3^{2-} + 2H^+`
的伙伴（`HCO_3^-`）不在反应里、但族内另一形态（`CO_3^{2-}`）在 ⟹
旧闸放行 ⟹ 两侧都从同一个 `HCO_3^-` 池重建 ⟹ 这条"碳酸族形态互变"的
`S` 变成恒等式，`Ni41` 因此走进 **333 步极限环**（残差 0.859 → 11.964）。

本组锚锁的是**碳酸族最独立可推导的化学事实** —— 全部来自
`pKa1 = 6.35`（`CO_2`/`HCO_3^-`）、`pKa2 = 10.33`（`HCO_3^-`/`CO_3^{2-}`），
**手算即可复现**，不含引擎输出：

* **两性盐**：`c = 0.1 M NaHCO_3` 的 pH = ½(pKa1 + pKa2) = **8.34**
  （两性物质近似式，且 `c` 远大于 `Ka1`、远小于 `Ka2` 的适用区间）；
* **弱酸强碱盐**：`c = 0.1 M Na_2CO_3` 的 `Kb = Kw/Ka2 = 2.0e-4`，
  `[OH^-] = sqrt(Kb·c) = 4.47e-3` ⟹ pOH 2.35 ⟹ pH = **11.65**；
* **计量**：`Ca^{2+} + 2HCO_3^- -> CaCO_3 + CO_2 + H_2O`（钙的碳酸盐
  定量沉淀，`Ksp(CaCO_3) = 3.3e-9` 极小）；
* **溶解方向**：`CaCO_3` 在 `CO_2` 溶液中溶解（石灰岩溶洞 / 硬水），
  `CaCO_3 + CO_2 + H_2O -> Ca^{2+} + 2HCO_3^-` —— 只断言**方向**
  （`Ca^{2+}` 增加、pH 落在酸性侧），不锁绝对量。

## 用法

    python tools/add_carbonate_cases.py            # dry-run
    python tools/add_carbonate_cases.py --check
    python tools/add_carbonate_cases.py --write
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

_TAG = "第 273 轮扩充 · 碳酸族锚。"


def build() -> list[dict]:
    return [
        {
            "name": "CA1 NaHCO3 0.1M 两性 pH（碳酸族锚）",
            "subs": [["NaHCO_3", 0.1]],
            "ph": [8.15, 8.55],
            "has_not": {"CO_3^{2-}": 0.02},
            "note": (
                _TAG + "**两性物质 pH 的独立可推导值**：0.1 M `NaHCO_3` 的"
                "pH = ½(pKa1 + pKa2) = ½(6.35 + 10.33) = **8.34**"
                "（适用条件 `c ≫ Ka1`、`c ≪ Ka2` 均满足）。区间取 ±0.2。"
                "同时锁死**一件不该发生的事**：`CO_3^{2-}` 不得 ≥ 0.02"
                "（真值 ≈ 0.0011，由 `[CO_3^{2-}] = c·Ka2/(Ka2+a_H)` 手算）。"
                "**本条正是本轮族级判据的靶心**：碳酸族在账本里以"
                "`HCO_3^-` 为主，若允许跨形态重建把 `CO_3^{2-}` 从"
                "`HCO_3^-` 池凭空放大，就会在这里造出不存在的碱度。"),
        },
        {
            "name": "CA2 Na2CO3 0.1M 水解 pH（碳酸族锚）",
            "subs": [["Na_2CO_3", 0.1]],
            "ph": [11.45, 11.85],
            "has": {"CO_3^{2-}": 0.05},
            "note": (
                _TAG + "**弱酸强碱盐的独立可推导值**：`Kb = Kw/Ka2 = "
                "1.0e-14/5.01e-11 = 2.0e-4`，`[OH^-] = sqrt(Kb·c) = "
                "sqrt(2.0e-4 × 0.1) = 4.47e-3` ⟹ pOH 2.35 ⟹ "
                "pH = **11.65**（区间 ±0.2）。"
                "计量侧锁**下界**：水解度只有 4.5% ⟹ `CO_3^{2-}` 至少还剩"
                "一半（断言 `≥ 0.05`，真值 0.0955）。"
                "⚠️ 这里**不能用上限断言**：真值 0.0955 大于任何"
                "形如 `< 0.05` 的上限，写上限等于要求引擎做错。"
                "首版正是这么写的，`--check` 立刻暴露 —— 记下来："
                "**上限断言要先手算真值再定阈，否则是在给引擎下错误的命令**。"),
        },
        {
            "name": "CA3 CaCl2+2NaHCO3 定量沉钙（碳酸盐计量锚）",
            "subs": [["CaCl_2", 1], ["NaHCO_3", 2]],
            "has": {"CaCO_3": 0.9},
            "has_not": {"Ca^{2+}": 0.1},
            "note": (
                _TAG + "**计量 + 溶解度**：`Ca^{2+} + 2HCO_3^- -> "
                "CaCO_3 + CO_2 + H_2O`。`Ksp(CaCO_3) = 3.3e-9` 极小 ⟹ "
                "钙基本定量沉出（断言 `CaCO_3 ≥ 0.9`、残余 `Ca^{2+} < 0.1`）。"
                "两条都只用**溶度积极小**这一化学事实，不含引擎输出。"
                "与 `CA4` 配对：本条锁\"碳酸氢盐能把钙沉下来\"，"
                "`CA4` 锁\"碳酸能把碳酸钙溶回去\"，方向相反、互为正负锚。"
                "⚠️ **已知问题区（如实记账，不作断言）**：引擎在本例给出 "
                "`pH = 12.17`、`resid = 6.853`。开放体系里 CO₂ 逸出后溶液"
                "应接近中性（~8–9），12.17 偏高 ⟹ 该态的**呈现 pH 不可信**。"
                "本条只锁\"钙被定量沉出\"这一条与 pH 无关的事实；"
                "pH 侧另立 P0 线索（与 §1.12 的\"呈现 pH 与账本矛盾\"同族）。"),
        },
        {
            "name": "CA4 CaCO3+CO2 溶解（碳酸族·方向锚）",
            "subs": [["CaCO_3", 1], ["CO_2", 1]],
            "has": {"Ca^{2+}": 0.005},
            "has_not": {"CO_2": 0.99},
            "note": (
                _TAG + "**方向锚**（石灰岩溶洞 / 硬水的成因）："
                "`CaCO_3 + CO_2 + H_2O -> Ca^{2+} + 2HCO_3^-`。"
                "断言取**保守下界** `Ca^{2+} ≥ 0.005`（即至少溶解 0.5%，"
                "真值远大于此），只锁\"碳酸能把碳酸钙溶解\"这一方向，"
                "**不锁绝对量**（绝对量依赖 `Ksp`/`Ka` 的具体取值，"
                "而那些是数据不是事实）。"
                "`CO_2` 侧只锁\"确有消耗\"（`< 0.99`，即至少 1% 被用掉）。"
                "⚠️ 首版写的是 `< 0.5`，`--check` 实测引擎给 `CO_2 = 0.9438` "
                "⟹ **是断言方凭空定了个没有依据的阈值**（0.0255 mol 钙溶出"
                "只需 0.0255 mol CO₂）。**上限断言必须先手算真值再定阈。**"),
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
        spec = c.get("has") or ({"ph": c["ph"]} if "ph" in c else {})
        ban = c.get("has_not") or {}
        print(f"{c['name'][:44]:<44} has={str(spec)[:26]:>26} "
              f"not={str(ban)[:18]:>18}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<44} {'引擎末态':>34} {'resid':>8} {'判读':>4}")
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
            ph = r.get("final_pH")
            if "ph" in c:
                ok = ok and (c["ph"][0] <= (ph if ph is not None else -99)
                             <= c["ph"][1])
            print(f"{c['name'][:44]:<44} {got[:34]:>34} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>4}  pH={ph}")
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
