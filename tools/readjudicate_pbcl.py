# -*- coding: utf-8 -*-
"""第 292 轮（续）· **PbCl₂ 溶解度族重裁 + ZN6 断言键纠正**。

## 一、ZN6：**断言键用错层**（用例 bug，非引擎缺陷）

`testsuit.py` 有**两套量值视图**（第 258–276 行）：

* `has` 读 **`amt`（池口径）**：池成员（`T.pool_ligands`）**不按物种名**记账，
  而是聚合到锚上；
* `has_not` / `has_range` 读 **`amt_sp`（物种级）**。

`T.pool_ligands = {Cu²⁺: [Cl⁻], Sn²⁺: [Cl⁻], Zn²⁺: [Cl⁻]}` —— **Cd²⁺ 不在池里**，
所以 `CD3`（Cd）用 `has` 断言 `[CdCl_4]^{2-}` 能过，而 `ZN6`（Zn）用同一个键断言
`[ZnCl_3]^-` 就恒读到 **0**。ZN6 本是本轮新加的锚，**引擎实测账本完全正确**
（`[ZnCl₃]⁻ = 0.072118`、`Zn²⁺ = 0.010509`，与库值独立不动点推导一致），
失败纯属"物种级下界写进了池级键" ⟹ 改用 `has_range`。

## 二、PbCl₂ 族 5 例：旧标是**简单模型快照**，与库内的氯络合数据**不能同时成立**

旧标全部来自

    s = (Ksp/4)^(1/3) ≈ 0.0158 M        （Ksp = 10^-4.8，库值）

该式的**前提是"溶解的铅全部是游离 Pb²⁺"**（电荷平衡 `2[Pb²⁺] = [Cl⁻]`）。
但库内已有 Pb²⁺–Cl⁻ 阶梯（`β₁ = 1.6` 一氯 IUPAC；`β₃ = 2.0` 三氯 Luo 2007），
而 PbCl₂ 自身溶解就提供 Cl⁻（饱和液 `[Cl⁻] ≈ 2s`）⟹ 络合必然参与：

* 络合把**游离 Pb²⁺ 压低**（铅被 `[PbCl]⁺` 分走，理想下占 64%）；
* 络合把**总溶解度抬高**（Cl⁻ 被络合消耗）。

`tools/pbcl_solubility.py`（**只读库内常数、独立写质量/电荷平衡**）实测：

| 体系 | 简单模型 s₀ | 含络合总溶 s | 游离 Pb²⁺ | `[PbCl]⁺` 占比 |
|---|---|---|---|---|
| 纯水 298 K | 0.0158 | **0.0223** | **0.0080** | 63.7% |
| 0.2 M NaCl | — | **0.0037** | **0.0004** | 81.4% |

引擎（SIT/DH 口径，I ≈ 0.05 / 0.21）实测总溶 **0.02746 / 0.00680**、
游离 **0.01086 / 0.00162** —— 比理想值高 23%/84%，方向与"Δz² = +6 的溶解反应
被活度层增强"一致（络合 Δz² = −4 被削弱，两者反向）。

⟹ **重裁**：断言改锁**游离 Pb²⁺ 与余固相**（物种级、区间），
撤掉单步总方程式断言（`PbCl_2 -> 2Cl^- + Pb^{2+}` 在络合参与下**本就不是**净方程；
按本轮 Cu 族既定裁定"产物复杂可不查总方程式"）。
区间一律取 **[理想下界, 引擎+余量]**，并**特意卡掉简单模型值**
（0.0158 游离 / 0.1842、0.9842 余固相）⟹ 旧标若回归立刻翻红。

## 三、诚实标注的数据待核项

**PbCl₂ 的实验总溶解度**文献分歧（教科书习题按 Ksp 给 ~0.016 M；手册给
~10.8 g/L ≈ 0.039 M）。本轮**不裁决**该分歧，只裁决"**旧标的游离=总溶**这一等式
与库内络合数据不能同时成立"；实验值待核项记入 handoff。

## 用法

    python tools/readjudicate_pbcl.py            # dry-run
    python tools/readjudicate_pbcl.py --check
    python tools/readjudicate_pbcl.py --write
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

_TAG = "第 292 轮重裁 · "
_DERIV = (
    "**判据（与旧标无关地重推）**：库内 `Ksp(PbCl₂) = 10^-4.8`、"
    "`β₁(PbCl⁺) = 1.6`（IUPAC）、`β₃(PbCl₃⁻) = 2.0`（Luo 2007，第 292 轮入库）。"
    "旧标用简单模型 `s = (Ksp/4)^(1/3) = 0.0158 M`，其前提是"
    "**游离 Pb²⁺ = 总溶解度**；而饱和液 `[Cl⁻] ≈ 2s ≈ 0.045 M` 下"
    "`[PbCl⁺]/[Pb²⁺] = β₁[Cl⁻] ≈ 1.8` ⟹ 该前提**不成立**。"
    "独立手算（`tools/pbcl_solubility.py`，只读库内常数、自带质量/电荷平衡）："
    "**总溶 0.0223 M、游离 Pb²⁺ 0.0080 M**。引擎（SIT 口径）总溶 0.0275、"
    "游离 0.0109 —— 高 23%，与 `Δz² = +6` 的溶解反应被活度层增强一致。"
    "区间取 `[理想下界, 引擎+余量]`，并**卡掉简单模型值**（游离 0.0158）⟹ 旧标回归即翻红。"
    "**撤掉单步 `eq`**：络合参与后 `PbCl_2 -> 2Cl^- + Pb^{2+}` 本就不是净方程"
    "（按本轮 Cu 族既定裁定「产物复杂可不查总方程式」）。"
    "**待核**：PbCl₂ 实验总溶解度文献分歧（教科书 ~0.016 vs 手册 ~0.039 M）——"
    "本轮不裁决，只裁决「游离=总溶」这一等式与库内络合数据不能同时成立。"
)

PATCHES: list[dict] = [
    {
        "name": "ZN6 ZnCl2 0.1+NaCl 6 三氯主导锚（与 CD3 对照）",
        "del": ["has"],
        "set": {"has_range": {"[ZnCl_3]^-": [0.05, 0.12],
                              "[ZnCl_2]": [0.004, 0.03],
                              "Zn^{2+}": [0.004, 0.02]}},
        "note_suffix": (
            " ⚠️ **第 292 轮修正断言键**：原用 `has: {'[ZnCl_3]^-': 0.05}` ⟹ 恒读到 0。"
            "根因在 `testsuit.py` 第 258–276 行：`has` 读**池口径** `amt`，而 "
            "`T.pool_ligands = {Cu²⁺, Sn²⁺, Zn²⁺} → Cl⁻` ⟹ **锌的氯合形态被聚合进 "
            "Zn²⁺ 池**，物种名下恒为 0；`has_not`/`has_range` 才读**物种级** `amt_sp`。"
            "（对照：`CD3` 用同一键能过，是因为 **Cd²⁺ 不在池里**。）"
            "引擎账本本身是对的：`[ZnCl₃]⁻ = 0.072118`、`[ZnCl₂] = 0.012104`、"
            "`[ZnCl]⁺ = 0.005268`、`Zn²⁺ = 0.010509` —— 与库值独立不动点推导的"
            "「三氯主导」逐项一致。**故改用 `has_range`（物种级区间）**。"),
    },
    {
        "name": "D45 PbCl2 冷水微溶",
        "del": ["has", "eq"],
        "set": {"has_range": {"Pb^{2+}": [0.006, 0.014],
                              "PbCl_2": [0.163, 0.181]}},
        "note": _TAG + _DERIV + " 本条（PbCl₂ 0.2 mol/L）：引擎游离 Pb²⁺ **0.01086**、"
                "余固相 **0.17255**（总溶 0.02746），与 J05 的 0.02744 **逐位一致**"
                "（饱和溶解度与投料量无关 ⟹ 自洽性内证）。",
    },
    {
        "name": "J05 PbCl2 溶解度@298K",
        "del": ["has_range", "eq"],
        "set": {"has_range": {"Pb^{2+}": [0.006, 0.014],
                              "PbCl_2": [0.963, 0.982]}},
        "note": _TAG + _DERIV + " 本条（1 mol/L 投料）：引擎游离 Pb²⁺ **0.01086**、"
                "余固相 **0.97256**（总溶 0.02744）。"
                "旧标 `Pb²⁺ ∈ [0.012, 0.02]` 是简单模型值，已由本裁取代。",
    },
    {
        "name": "J06 PbCl2 溶解度@363K",
        "del": ["has_range", "eq", "eq_has"],
        "set": {"has_range": {"Pb^{2+}": [0.03, 0.06],
                              "PbCl_2": [0.925, 0.955]}},
        "note": _TAG + "363.15 K。旧标 `Pb²⁺ ∈ [0.014, 0.045]` 同样是简单模型口径，"
                "而引擎值 **0.04493** 恰好压在旧上沿（差 0.00007）—— 这种「贴边过」"
                "说明旧窗不是按含络合的真值定的。同族的 298 K 情形已由 "
                "`tools/pbcl_solubility.py` 独立给出真值（总溶 0.0223、游离 0.0080），"
                "363 K 只把 `Ksp(T)` 抬高，络合比例同量级 ⟹ 区间改取 "
                "`[0.03, 0.06]`（卡掉简单模型的 ~0.027）与余固相 `[0.925, 0.955]`"
                "（引擎 0.93493）。**撤掉 `eq_has`**：络合参与后单步式不是净方程。"
                "**待核**：`Ksp(363.15)` 的独立推导（库内由 Hess 从 thermo.json 派生）"
                "本轮未做，故 363 K 两例的区间比 298 K 宽 —— 记入 handoff 待核项。",
    },
    {
        "name": "J14 PbCl2 溶解度@363K",
        "del": ["has_range", "eq", "eq_has"],
        "set": {"has_range": {"Pb^{2+}": [0.015, 0.05],
                              "PbCl_2": [0.0, 0.002]}},
        "note": _TAG + "363.15 K、投料仅 0.05 mol/L。引擎**全溶**（余固相 0）、"
                "游离 Pb²⁺ **0.02374**。判据：旧注自身给的趋势是"
                "「298 K 0.016 → 363 K 0.027」，再乘上络合因子（298 K 实测 F = 2.78）"
                "⟹ 363 K 饱和总溶 ≈ 0.06 M > 0.05 ⟹ **全溶是化学事实**，"
                "故余固相上界取 0.002（「基本全溶」）。"
                "**撤掉 `eq_has`**：络合参与后单步式不是净方程。",
    },
    {
        "name": "B24 PbCl2+0.2M NaCl 同离子",
        "del": ["has", "has_not", "eq"],
        "set": {"has_range": {"Pb^{2+}": [0.0002, 0.003],
                              "PbCl_2": [0.010, 0.018]}},
        "note": _TAG + "**同离子效应**：0.2 M Cl⁻ 下 `[Pb²⁺] = Ksp/[Cl⁻]² ≈ 4e-4 M`"
                "（旧注的 ~5e-4 是这个口径，对的）。**但旧标把「总溶解量」也当成了 "
                "4e-4**（故写 `PbCl₂ ≥ 0.018`、`eq: None`「不构成显著反应」）——"
                "而 0.2 M Cl⁻ 下 `[PbCl⁺]` 占总溶 **81%**（独立手算），"
                "总溶 ≈ 0.0037 M（理想）/ **0.00680 M**（引擎）。"
                "⟹ 游离 Pb²⁺ 仍然极低（旧标的**定性结论对**），"
                "但「总溶解 < 0.002」不成立 ⟹ 余固相改取 `[0.010, 0.018]`"
                "（引擎 0.01320、理想 0.0163）。**撤掉 `eq: None`**："
                "氯合形态生成是真反应，不再是「无净方程」。",
    },
]


def main(argv: list[str]) -> int:
    path = os.path.join(ROOT, "chemkit", "data", "tests.json")
    with io.open(path, encoding="utf-8") as fh:
        db = json.load(fh)
    index = {c["name"]: c for c in db}

    for p in PATCHES:
        c = index.get(p["name"])
        if c is None:
            print(f"✗ 未找到用例：{p['name']}")
            return 2
        print(f"{p['name'][:44]:<46} 删{len(p['del'])}键 "
              f"设{list(p['set'])} {'(改写note)' if 'note' in p else '(追加note)'}")

    if "--check" in argv:
        from chemkit.data import load_tables
        from chemkit.engine import judge
        T = load_tables()
        print(f"\n{'name':<46} {'判读':>4}   关键量")
        bad = 0
        for p in PATCHES:
            c = index[p["name"]]
            pr: dict = {}
            r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            amt_sp: dict = {}
            for e in r["production"] + r["final"]:
                amt_sp[e["name"]] = max(amt_sp.get(e["name"], 0.0), e["mol"])
            ok = True
            msgs = []
            for sp, (lo, hi) in p["set"].get("has_range", {}).items():
                m = amt_sp.get(sp, 0.0)
                good = lo <= m <= hi
                ok = ok and good
                msgs.append(f"{sp}={m:.5g}{'' if good else ' ✗'}")
            print(f"{p['name'][:44]:<46} {'✓' if ok else '✗':>4}   "
                  + "  ".join(msgs))
            bad += 0 if ok else 1
        print(f"\n不通过 {bad} 条")
        return 1 if bad else 0

    if "--write" not in argv:
        print("\n[dry-run] 加 --write 写入")
        return 0
    for p in PATCHES:
        c = index[p["name"]]
        for k in p["del"]:
            c.pop(k, None)
        c.update(p["set"])
        if "note" in p:
            c["note"] = p["note"]
        else:
            c["note"] = (c.get("note") or "") + p["note_suffix"]
    with io.open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(db, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"\n[写入] 重裁 {len(PATCHES)} 条 ⟹ 共 {len(db)} 条")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
