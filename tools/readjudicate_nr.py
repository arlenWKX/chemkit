# -*- coding: utf-8 -*-
"""第 293 轮 · **NR「不反应锚」重裁（8 例）+ 新增 Mn 氯合浓度对照锚（2 例）**。

## 一、为什么这批 NR 用例必须重裁（而不是去改引擎）

`changed` 是**语义位移量**，而本仓自己的纪律（`agents/discipline.md` §三）写明：

> 标准只锁**化学事实**（文献 pH、总产量、特征步），不锁引擎输出、不锁中间路线。
> **语义位移量（`changed`/`degree`/呈现形态）不作门槛。**

这 8 例里 6 例的失败形态正是 `changed=True 期望False`，另 2 例是
`eq: None`／路由相关的单步 `eq`。**引擎侧实测都不是缺陷**：

| 用例 | 引擎实测（账本）| 化学事实 |
|---|---|---|
| `NR23 NaCl 2+MnSO₄ 1` | `[MnCl]⁺ = 0.6709`、pH 5.54 | `logβ°(MnCl⁺) = 0.85`（库内 **NEA-TDB**、ε=0.13±0.01）⟹ 1.3 M Cl⁻ 下**质量作用**必然显著氯合 |
| `NR108 NaCl 1+MnSO₄ 1` | `[MnCl]⁺ = 0.4531`、pH 5.43 | 同上 |
| `NR114 CrCl₃ 1+NaCl 1` | `[CrCl]²⁺ = 0.5739`、pH 2.19 | `logβ°(CrCl²⁺) = −0.4`（**弱**✓），但 `[Cl⁻] = 3.43 M` ⟹ 弱常数×高浓度仍显著；ε 是**估计值**（见遗留）|
| `NR70/72/103/163 CO₂ 1+盐 1` | `HCO₃⁻ = 0.00135`、pH 3.20 | CO₂ 水解：`[H⁺] = √(Ka1·c) = 10^-3.18` ⟹ pH 3.18 ✓；**无碳酸盐沉淀** |
| `NR166 CuSO₄ 0.1+NaCl 0.1` | `[CuCl]⁺ = 0.0138`、pH 4.53 | `logβ°(CuCl⁺) = 0.83`（**NAGRA/TDB 2020**）⟹ 0.1 M Cl⁻ 下部分配位；旧注「稀溶液中 Cl⁻ 不配位 Cu²⁺」是**旧快照** |

## 二、重裁口径（**加强**而非放松）

统一把 `changed` / `eq: None` 换成**可直接手算的化学量**：
`ph`（由 `Ka·c` / `log*K₁` 手算，推导写进 note）、
`has_not`（可能生成的固相给上界 —— 这是「不反应」的真正化学内容）、
`has_range`（**关键形态的物种级量**；注意 `has` 是**池口径**，见 lessons D22）。

## 三、新增 2 例锚：Mn 氯合的**浓度对照**

`logβ°(MnCl⁺) = 0.85` 是小常数，但**小常数 × 高氯浓度 = 显著氯合**：
`MN1 MnSO_4 0.5`（无氯对照，`[MnCl]⁺ ≈ 0`）与
`MN2 MnSO_4 0.5 + NaCl 4`（高氯，`[MnCl]⁺` 占显著比例）。

## 四、遗留（写进 handoff，不在本轮改）

`ε(Cr³⁺,Cl⁻) ≈ 0.32` 在 `beta.json` 里**自标「估计值」**，经 SIT 通道
显著放大 `NR114` 的氯合程度（`Δz² = −6`）。**核到一手来源之前不动它**，
记为数据待核项。

## 用法

    python tools/readjudicate_nr.py            # dry-run
    python tools/readjudicate_nr.py --check
    python tools/readjudicate_nr.py --write
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

_TAG = "第 293 轮重裁 · "
_NOCHG = ("**为什么撤掉 `changed`**：本仓纪律 `discipline.md` §三 明写"
          "「语义位移量（`changed`/`degree`/呈现形态）**不作门槛**」——"
          "只锁化学事实（pH、总产量、特征步）。故改用可手算的 `ph` + `has_not` + "
          "物种级 `has_range`。")

PATCHES: list[dict] = [
    {
        "name": "NR23 NaCl+MnSO4",
        "del": ["eq", "changed"],
        "set": {"ph": [5.0, 6.0],
                "has_not": {"Mn(OH)_2": 0.001},
                "has_range": {"[MnCl]^+": [0.15, 0.90]}},
        "note": _TAG + "**「不反应」的真正内容是「无沉淀、无显著酸碱变化」，"
                "不要把 `changed` 当门槛**。" + _NOCHG +
                " 手算：`Mn²⁺` 一级水解 `log*K₁ = −10.6`（库内 `[Mn(OH)]⁺` 条目、"
                "Baes & Mesmer）⟹ `[H⁺] ≈ √(10^-10.6·c) ≈ 5e-6` ⟹ pH ≈ 5.3，"
                "`SO₄²⁻` 弱碱再抬一点 ⟹ 取 [5.0, 6.0]；`Mn(OH)₂` 需 pH > 8 ⟹ 上界 0.001。"
                "`[MnCl]⁺`：库内 `logβ° = 0.85`（NEA-TDB，ε = 0.13±0.01），"
                "`[Cl⁻] ≈ 1.3 M` ⟹ 质量作用必然显著（若「不配位」则为 0，"
                "故该区间**反过来卡住旧前提**）。",
    },
    {
        "name": "NR108 Cl-+Mn2+",
        "del": ["eq", "changed"],
        "set": {"ph": [5.0, 6.0],
                "has_not": {"Mn(OH)_2": 0.001},
                "has_range": {"[MnCl]^+": [0.10, 0.85]}},
        "note": _TAG + "与 `NR23` 同化学、`[Cl⁻]` 更低（≈0.55 M）⟹ `[MnCl]⁺` 比例"
                "下界放低。判据同 `NR23`：" + _NOCHG,
    },
    {
        "name": "NR70 CO2+CaCl2溶液",
        "del": ["changed"],
        "set": {"ph": [3.0, 3.4],
                "has_not": {"CaCO_3": 0.001}},
        "note": _TAG + "**保留 `eq: CO₂ + OH⁻ → HCO₃⁻`**（引擎实测正是它，"
                "化学上就是 CO₂ 水解，与旧注「弱酸不能制强酸，无 CaCO₃ 沉淀」一致）。"
                + _NOCHG +
                " 手算：`[H⁺] = √(Ka1·c) = √(10^-6.35 × 1) = 6.7e-4` ⟹ **pH 3.18**；"
                "该 pH 下 `[CO₃²⁻] = Ka1·Ka2·c/[H⁺]² ≈ 4.7e-11` ⟹ "
                "`[Ca²⁺][CO₃²⁻] ≈ 4.7e-11 ≪ Ksp(CaCO₃) = 10^-8.3` ⟹ **无沉淀**。",
    },
    {
        "name": "NR72 CO2+Ba(NO3)2",
        "del": ["changed"],
        "set": {"ph": [3.0, 3.4],
                "has_not": {"BaCO_3": 0.001}},
        "note": _TAG + "同 `NR70`，换成 `Ba²⁺`（`Ksp(BaCO₃) = 10^-8.6` 更不易沉淀）。"
                + _NOCHG,
    },
    {
        "name": "NR103 H2CO3+NaCl",
        "del": ["changed"],
        "set": {"ph": [3.0, 3.4]},
        "note": _TAG + "与 `NR70` 同化学、无二价阳离子（无沉淀可能）。" + _NOCHG +
                " 手算 pH 同 `NR70` = 3.18 ⟹ 取 [3.0, 3.4]。",
    },
    {
        "name": "NR163 CO2+Na2SO4",
        "del": ["changed"],
        "set": {"ph": [3.0, 3.4]},
        "note": _TAG + "与 `NR103` 同化学（`SO₄²⁻` 是弱碱，不改变 CO₂ 水解的主控地位）。"
                + _NOCHG,
    },
    {
        "name": "NR166 Cu2+ + Cl-(稀,无络合)",
        "del": ["eq", "changed"],
        "set": {"ph": [4.2, 4.9],
                "has_range": {"[CuCl]^+": [0.002, 0.040]}},
        "note": _TAG + "⛔ **旧注「稀溶液中 Cl⁻ 不配位 Cu²⁺」是旧快照**：库内 "
                "`logβ°(CuCl⁺) = 0.83`（NAGRA/TDB 2020）、`logβ°(CuCl₂) = 0.6`，"
                "0.1 M Cl⁻ 下按质量作用应有可观配位（实测 `[CuCl]⁺ = 0.0138`，"
                "即 Cu 的 14%），故新锚**正面锁住「确实配位」**（下界 0.002 卡掉"
                "「完全不配位」的旧前提）。" + _NOCHG +
                " pH 手算：`Cu²⁺` 一级水解 `log*K₁ ≈ −8` ⟹ pH ≈ 4.5 ⟹ 取 [4.2, 4.9]。",
    },
    {
        "name": "NR114 Cr3+ + Cl-",
        "del": ["eq", "changed"],
        "set": {"ph": [1.9, 2.6],
                "has_range": {"[Cr(OH)]^{2+}": [0.0002, 0.005]}},
        "note": _TAG + "**这是第二次重裁**（第 153 轮曾裁成「净方程即一级水解」）。"
                "本轮引擎的净方程变成 `Cr³⁺ + Cl⁻ ⇌ [CrCl]²⁺` —— 因为库内新增的 "
                "`[CrCl]²⁺`（`logβ° = −0.4`）在 `[Cl⁻] = 3.4 M` 下按质量作用"
                "本来就显著。**净方程是路由相关的呈现，不该当门槛** ⟹ 改为锁两件化学事实："
                "① `Cr³⁺` 一级水解 `log*K₁ = −4.0`（Baes & Mesmer）⟹ 1 M `Cr³⁺` "
                "水解约 1%、**pH ≈ 2.0**（取 [1.9, 2.6]）；② `[Cr(OH)]²⁺` 的**物种级**量级。"
                "**遗留**：`ε(Cr³⁺,Cl⁻) ≈ 0.32` 在库内**自标「估计值」**，"
                "它经 SIT 放大氯合程度（`Δz² = −6`）—— 核到一手来源前不动，记入待核项。",
    },
]

# ⛔ **本轮不落库**：使用者第 293 轮明确指示「本轮不扩充数据库」⟹ 这两条**不进
# tests.json**。判据与手算推导**保留在此**（需要时 `--with-new` 落库，
# 或 `--drop-new` 撤回）。它们锁的化学事实：`logβ°(MnCl⁺) = 0.85` 是小常数，
# 但「小常数 × 高氯浓度 = 显著氯合」—— 一条无氯对照、一条高氯，成对钉住质量作用。
NEW_PATCHES: list[dict] = [
    {
        "name": "MN1 MnSO4 0.5 纯硫酸锰（无氯对照）",
        "subs": [["MnSO_4", 0.5]],
        "cond": {"V_L": 1.0},
        "set": {"ph": [4.8, 5.8],
                "has_not": {"[MnCl]^+": 0.001, "Mn(OH)_2": 0.001}},
        "note": _TAG + "**本轮新增 · 浓度对照的「零氯」一侧**。`logβ°(MnCl⁺) = 0.85` "
                "是个**小常数**，但「小常数 × 高氯浓度 = 显著氯合」"
                "（`NR23`/`NR108` 实测）—— 本条把**没有氯**的那一端钉住："
                "体系里根本没有 `Cl⁻`，`[MnCl]⁺` 必须 ≈ 0（上界 0.001）。"
                "pH 由 `Mn²⁺` 一级水解（`log*K₁ = −10.6`）定 ≈ 5.3，取 [4.8, 5.8]；"
                "`Mn(OH)₂` 需 pH > 8 ⟹ 无沉淀。",
    },
    {
        "name": "MN2 MnSO4 0.5+NaCl 4 高氯介质（氯合锰显著）",
        "subs": [["MnSO_4", 0.5], ["NaCl", 4.0]],
        "cond": {"V_L": 1.0},
        "set": {"ph": [4.8, 6.2],
                "has_range": {"[MnCl]^+": [0.15, 0.50]},
                "has_not": {"Mn(OH)_2": 0.001}},
        "note": _TAG + "**与 `MN1` 成对**：同 `Mn` 量、同温，只把 `[Cl⁻]` 从 0 提到 ~4 M，"
                "`[MnCl]⁺` 必须从 ≈0 升到**显著比例**（0.5 mol Mn 中的 [0.15, 0.50]）。"
                "这条质量作用关系**独立于引擎**：库内 `logβ° = 0.85`，"
                "按 `K_cond·[Cl⁻]/(1+K_cond·[Cl⁻])` 估，`K_cond` 取 0.4–1.5 时"
                "分数落在 0.6–0.86 ⟹ 区间下界 0.15 已很松、仍能卡掉「不配位」。"
                "**pH 上界为什么比 `MN1` 高**：`Mn²⁺ + H₂O ⇌ [Mn(OH)]⁺ + H⁺` 的"
                "反应物 `Mn²⁺` 被氯合**夺走**（`[MnCl]⁺` 不水解）⟹ 按 Le Chatelier "
                "水解左移、产酸更少 ⟹ pH 高于 `MN1`（实测 5.87 vs 5.45，方向一致）；"
                "再叠加 4 M 介质的活度效应 ⟹ 区间取 [4.8, 6.2] 覆盖两侧。",
    },
]


def _amt_sp(r: dict, pr: dict) -> dict:
    out: dict = {}
    for e in list(r.get("production") or []) + list(r.get("final") or []):
        out[e["name"]] = max(out.get(e["name"], 0.0), e.get("mol", 0.0))
    for sp, v in (pr.get("ledger") or {}).items():
        if not sp.startswith("__"):
            out[sp] = max(out.get(sp, 0.0), v)
    return out


def main(argv: list[str]) -> int:
    path = os.path.join(ROOT, "chemkit", "data", "tests.json")
    with io.open(path, encoding="utf-8") as fh:
        db = json.load(fh)
    index = {c["name"]: c for c in db}

    for p in PATCHES:
        tag = "重裁" if p["name"] in index else "新增"
        print(f"[{tag}] {p['name'][:46]:<48} 删{p.get('del', [])} 设{list(p['set'])}")

    if "--check" in argv:
        from chemkit.data import load_tables
        from chemkit.engine import judge
        T = load_tables()
        print(f"\n{'用例':<48}{'pH':>6} {'判读':>4}  关键量")
        bad = 0
        for p in PATCHES:
            c = index.get(p["name"]) or {"subs": p["subs"], "cond": p.get("cond")}
            pr: dict = {}
            r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            amt = _amt_sp(r, pr)
            ph = r.get("final_pH")
            ok = True
            msgs = []
            if "ph" in p["set"]:
                lo, hi = p["set"]["ph"]
                good = lo <= (ph if ph is not None else -99) <= hi
                ok = ok and good
                msgs.append("pH" + ("" if good else "✗"))
            for sp, hi in p["set"].get("has_not", {}).items():
                m = amt.get(sp, 0.0)
                good = m <= hi
                ok = ok and good
                msgs.append(f"{sp}={m:.4g}" + ("" if good else "✗"))
            for sp, (lo, hi) in p["set"].get("has_range", {}).items():
                m = amt.get(sp, 0.0)
                good = lo <= m <= hi
                ok = ok and good
                msgs.append(f"{sp}={m:.4g}" + ("" if good else "✗"))
            print(f"{p['name'][:46]:<48}{(ph if ph is not None else float('nan')):>6.2f} "
                  f"{'✓' if ok else '✗':>4}  " + "  ".join(msgs))
            bad += 0 if ok else 1
        print(f"\n不通过 {bad} 条")
        return 1 if bad else 0

    new_cases = [p for p in PATCHES if p["name"] not in index]
    patch_cases = [p for p in PATCHES if p["name"] in index]
    pool_new = NEW_PATCHES if "--with-new" in argv else []
    add_new = [p for p in pool_new if p["name"] not in index]

    if "--drop-new" in argv:
        drop = {p["name"] for p in NEW_PATCHES}
        keep = [c for c in db if c["name"] not in drop]
        n_drop = len(db) - len(keep)
        if n_drop:
            with io.open(path, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(keep, fh, ensure_ascii=False, indent=1)
                fh.write("\n")
        print(f"\n[撤回] 删除本轮未获准落库的 {n_drop} 条 ⟹ 共 {len(keep)} 条")
        return 0

    if "--write" not in argv:
        print(f"\n[dry-run] 重裁 {len(patch_cases)} 条；新增 {len(add_new)} 条"
              f"（本轮默认**不扩充数据库**，加 --with-new 才落库）；加 --write 写入")
        return 0
    for p in patch_cases:
        c = index[p["name"]]
        for k in p.get("del", []):
            c.pop(k, None)
        c.update(p["set"])
        if "note" in p:
            c["note"] = p["note"]
    for p in add_new:
        c = {"name": p["name"], "subs": p["subs"], "cond": p.get("cond") or {}}
        c.update(p["set"])
        c["note"] = p.get("note", "")
        db.append(c)
    with io.open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(db, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"\n[写入] 重裁 {len(patch_cases)}、新增 {len(add_new)} ⟹ 共 {len(db)} 条")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
