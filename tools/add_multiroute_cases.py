# -*- coding: utf-8 -*-
"""第 264 轮 · 拾取修法的定向锚：**同产物多路线**体系的配比扫描。

## 动因

第 264 轮查明走步的拾取是**纯 S 贪心**，不看"这一步能走多少"。
`N15 Na₂SiO₃+足量CO₂` 里两条通往同一固相 `H₂SiO₃` 的路线：

    A `2H⁺ + SiO₃²⁻ -> H₂SiO₃`（precip）  S=4.162  可动量 1.99e-4
    B `HSiO₃⁻ + H⁺ -> H₂SiO₃`（proton）   S=4.077  可动量 0.5605

纯 S 贪心每轮选 A ⟹ 等量乒乓 3000 轮。修法 = 同产物通道按可动量回退。

## 答案来源（**不依赖引擎**，只用质量守恒 + 一条化学事实）

**化学事实**：碳酸强于硅酸（`CO₂` pKa₁ 6.4 < `HSiO₃⁻` pKa 12.0 / `H₂SiO₃` pKa 9.8）
⟹ `CO₂` 把硅酸根**定量**赶出溶液成为 `H₂SiO₃`。
**质量守恒**给死上下限：

    CO₂ : Si = 1 : 1   ⟹ `CO₂ + H₂O + SiO₃²⁻ -> CO₃²⁻ + H₂SiO₃`
                          （电荷守恒逼出 CO₃²⁻：LHS −2、RHS −2）
    CO₂ : Si = 2 : 1   ⟹ `2CO₂ + 2H₂O + SiO₃²⁻ -> 2HCO₃⁻ + H₂SiO₃`
    CO₂ : Si ≥ 2 : 1   ⟹ 同上（余量 CO₂ 留在溶液/逸出），Si 仍**全部**成固
    CO₂ : Si ≤ 1 : 1   ⟹ **CO₂ 是限量**，H₂SiO₃ 产量 ≤ CO₂ 量（本题取足量档）

故断言写成 `H₂SiO₃ ≥ 0.9 × n(Si)`（"定量析出"）与碳去向 `CO₃²⁻/HCO₃⁻ ≥ 0.9 × n(C)`；
**不含任何引擎输出的快照**。

**不选 1:0.5 这类"CO₂ 不足"档写硬断言**：此时产量由 `H₂SiO₃` 溶度积与 pH
共同决定（连续量），把它固化成具体数值正是 `Z31` 的教训（整数配比锁连续量）。

## 用法

    python tools/add_multiroute_cases.py            # dry-run
    python tools/add_multiroute_cases.py --check    # 只测引擎（不调断言）
    python tools/add_multiroute_cases.py --write
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

WHY = ("碳酸强于硅酸（CO₂ pKa₁ 6.4 < H₂SiO₃ 9.8 / HSiO₃⁻ 12.0）"
       "⟹ CO₂ 把硅酸根**定量**赶出溶液为 H₂SiO₃；"
       "碳的去向由 CO₂:Si 配比按**质量守恒 + 电荷守恒**定死")


def build() -> list[dict]:
    out = []

    def add(name, subs, has, note, eq=None, eq_has=None, ph=None,
            has_range=None):
        c = {"name": name, "subs": subs, "note": note, "changed": True}
        if has:
            c["has"] = has
        if has_range:
            c["has_range"] = has_range
        if eq:
            c["eq"] = eq
        if eq_has:
            c["eq_has"] = eq_has
        if ph:
            c["ph"] = ph
        out.append(c)

    # ---- 硅酸盐 / CO₂ 配比扫描（本轮修法的核心体系）----
    # ⚠️ 1:1 档**不写硬下限**：1 当量 CO₂ 把体系留在 pH ≈ 11（CO₃²⁻ 水解），
    #    而硅酸在碱性区会以 HSiO₃⁻/SiO₃²⁻ 重新溶解 ⟹ H₂SiO₃ 产量是
    #    「溶度积 × 该 pH」共同决定的**连续量**，把它固化成具体数值正是
    #    `Z31` 的教训（整数配比锁连续量）。故用 `has_range`。
    add("SI1 Na2SiO3+1当量CO2（同产物双路线·碱区复溶）",
        [["Na_2SiO_3", 1], ["CO_2", 1]],
        {},
        f"第 264 轮扩充 · 同产物双路线锚。{WHY}。"
        f"但 1:1 恰好把体系留在**碱区**（CO₃²⁻ 水解 ⟹ pH ≈ 11），"
        f"而硅酸在碱性区以 HSiO₃⁻/SiO₃²⁻ 复溶 ⟹ H₂SiO₃ 产量**不是**定量值，"
        f"而由「H₂SiO₃ 溶度积 × 该 pH」共同决定（连续量）。"
        f"故只写**区间**（大量析出但未到底）；碳侧同理不写硬断言。"
        f"这条故意留下「碱区复溶」这一维，与 SI2（2:1，pH ≈ 8，定量成固）配对。"
        f"⚠️ **不写净方程断言**（第 264 轮自查修正）：沉淀程度是连续量，"
        f"整数配比会把它锁死（`Z31` 的「整数配比锁连续量」教训）——"
        f"实测引擎净方程为 `12CO₂ + 11SiO₃²⁻ + 10H₂O + 2OH⁻ -> "
        f"12CO₃²⁻ + 11H₂SiO₃`（配比 12:11，元素与电荷逐项守恒），"
        f"而非理想化的 1:1。",
        has_range={"H_2SiO_3": [0.5, 0.95]})
    add("SI2 Na2SiO3+2当量CO2（同产物双路线）",
        [["Na_2SiO_3", 1], ["CO_2", 2]],
        {"H_2SiO_3": 0.9, "HCO_3^-": 1.8},
        f"第 264 轮扩充 · 同产物双路线锚。{WHY}。"
        f"2:1 恰好把碳全部送到 HCO₃⁻ ⟹ 净方程 "
        f"`2CO₂ + 2H₂O + SiO₃²⁻ -> 2HCO₃⁻ + H₂SiO₃`（与 N15 同式，"
        f"但 N15 是 2.5 当量即过量档，本条是**恰好**档）。"
        f"产物 pH ≈ 8（HCO₃⁻ 缓冲）⟹ 硅酸不再复溶，可写定量断言。",
        eq="2CO_2 + 2H_2O + SiO_3^{2-} -> 2HCO_3^- + H_2SiO_3")
    add("SI2b Na2SiO3+1.5当量CO2（同产物双路线·两可之间）",
        [["Na_2SiO_3", 1], ["CO_2", 1.5]],
        {"H_2SiO_3": 0.9},
        f"第 264 轮扩充 · 同产物双路线锚（1:1 与 2:1 之间的**中间档**）。"
        f"{WHY}。碳按 1.5 当量分摊在 CO₃²⁻/HCO₃⁻ 两形态上（比例随 pH 连续变化，"
        f"**不写碳的硬断言**），但硅的归宿仍是定量成固（pH 落在中性偏碱、"
        f"硅酸不复溶）⟹ 只锁 Si 侧。本条把「CO₂ 用量 → pH → 硅酸是否复溶」"
        f"这条连续链条取一个中间点。")
    add("SI3 Na2SiO3 0.1+2.5当量CO2（同产物双路线·稀档）",
        [["Na_2SiO_3", 0.1], ["CO_2", 0.25]],
        {"H_2SiO_3": 0.09},
        f"第 264 轮扩充 · 同产物双路线锚（**N15 的十分之一稀档**）。{WHY}。"
        f"稀档的意义：痕量族成员（SiO₃²⁻）的**绝对量**随浓度等比下降，"
        f"而 A 路线的可动量正是 [SiO₃²⁻] ⟹ 若拾取仍按纯 S 贪心，"
        f"爬行的**相对**代价与浓度无关（同样走不满），本条即该判据的稀档对照。")
    add("SI4 Na2SiO3 0.5+10当量CO2（同产物双路线·大过量）",
        [["Na_2SiO_3", 0.5], ["CO_2", 5]],
        {"H_2SiO_3": 0.45},
        f"第 264 轮扩充 · 同产物双路线锚（**大过量 CO₂**）。{WHY}。"
        f"CO₂ 大过量时 pH 被压到碳酸的酸性区，硅酸仍定量成固；"
        f"本条锁「过量程度不影响 Si 的归宿」这一条，"
        f"与 SI2（恰好档）、N15（2.5 当量）构成用量三维。")

    # ---- 同产物双路线的第二个体系：BaSO₃（Y05 的配比/浓度维）----
    add("BS1 BaCl2+Na2SO3 等摩尔（同产物双路线）",
        [["BaCl_2", 1], ["Na_2SO_3", 1]],
        {"BaSO_3": 0.9},
        "第 264 轮扩充 · 同产物双路线锚（第二个体系）。"
        "亚硫酸钡难溶，等摩尔即定量沉淀（Y05 的对照档："
        "Y05 在修法前后由 6.811 残差 → 0.0，本条锁其化学结论）。",
        eq="Ba^{2+} + SO_3^{2-} -> BaSO_3")
    add("BS2 BaCl2+Na2SO3 稀档 0.01M（同产物双路线）",
        [["BaCl_2", 0.01], ["Na_2SO_3", 0.01]],
        {"BaSO_3": 0.009},
        "第 264 轮扩充 · 同产物双路线锚（稀档）。同上，"
        "锁「沉淀完全程度与浓度无关」（难溶盐的定量性由 Ksp 决定，"
        "0.01 M 下 BaSO₃ 仍 ≥90% 成固）。",
        eq="Ba^{2+} + SO_3^{2-} -> BaSO_3")
    add("BS3 BaCl2过量+Na2SO3（同产物双路线·阳离子过量）",
        [["BaCl_2", 1], ["Na_2SO_3", 0.5]],
        {"BaSO_3": 0.45},
        "第 264 轮扩充 · 同产物双路线锚（**阴离子限量档**）。"
        "SO₃²⁻ 是限量者 ⟹ BaSO₃ 产量由它定死（0.5 mol 的 90%）；"
        "与 BS1（等摩尔）配对锁住「谁限量」这一维。",
        eq="Ba^{2+} + SO_3^{2-} -> BaSO_3")
    return out


def main(argv: list[str]) -> int:
    from chemkit.data import load_tables
    T = load_tables()
    new = build()
    path = os.path.join(ROOT, "chemkit", "data", "tests.json")
    with io.open(path, encoding="utf-8") as fh:
        db = json.load(fh)
    have = {c["name"] for c in db}
    print(f"{'name':<48} {'has':>34}  推导")
    for c in new:
        mark = "已存在" if c["name"] in have else "新增"
        spec = c.get("has") or c.get("has_range") or {}
        print(f"{c['name'][:48]:<48} {str(spec)[:34]:>34}  {mark}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<48} {'引擎末态产量':>28} {'resid':>7} {'判读':>6}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            fin = {e["name"]: e["mol"] for e in (r.get("final") or [])}
            want = c.get("has") or {}
            rng = c.get("has_range") or {}
            got = " ".join(f"{k}={fin.get(k, 0.0):.4g}"
                           for k in list(want) + list(rng))
            ok = (all(fin.get(k, 0.0) >= v for k, v in want.items())
                  and all(lo <= fin.get(k, 0.0) <= hi
                          for k, (lo, hi) in rng.items()))
            print(f"{c['name'][:48]:<48} {got[:28]:>28} "
                  f"{abs(_live(pr.get('active') or [])):>7.4f} "
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
