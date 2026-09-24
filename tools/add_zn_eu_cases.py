# -*- coding: utf-8 -*-
"""第 287 轮 · **Eu³⁺/Zn 还原度锚 + 锌的两性/氨合溶解锚**。

## 动因（第 287 轮查明的 `EU01` 空转机制）

第 286 轮拆掉气体活度地板后 `EU01`（`Zn 0.1 + Eu(NO3)3 0.1`）iters 76 → **2819**、
耗时 990 ms → 32.4 s。第 287 轮 `dev.py case EU01` + `CHEM_TRACE=1` 查明：

* 2815 步里绝大多数是 **`Zn²⁺ ⇌ [Zn(OH)₄]²⁻` 的结构性乒乓**，每步幅度仅
  **~3e-4 mol**（`0.000292→0.000290→0.000286→0.000280→…` 缓慢衰减），
  拣选时两个方向的 `S` 都很大（+15.9 / +10.8，**但这是先后两次拣选的值，
  不是同一态的矛盾** —— `tools/cand_audit.py` 全库 18014 条候选实测
  `logK(正)+logK(反)==0` **0 例违规**）；
* 四条检测器**全部**过不了线，因为阈值都是**绝对**的：近窗 `gross >= 0.05`、
  短周期 `Σext < 1e-3`（实测 ~2e-3 恰在带外）、物种周转 `>= 0.05`、
  以及"账本签名精确复现"（幅度缓慢衰减 ⟹ 签名**从不精确重现**）；
* 真正在跑的是一个 **`[freeze-limit-cycle]` → `强驱动跳过` → `[freeze-revive]`
  的活锁**：签名检测器**打中了**，但 `_freeze` 依 X-38 **拒绝冻结强驱动键**，
  而 `_try_revive`（第 274 轮改为**无终身额度**）每 16 步又把它们解冻 ⟹
  周而复始。

## 本组锚锁的化学事实（与引擎无关）

1. **Zn 能还原 Eu³⁺**：`E°(Zn²⁺/Zn) = −0.76 V` vs `E°(Eu³⁺/Eu²⁺) = −0.35 V`
   ⟹ `Zn + 2Eu³⁺ → Zn²⁺ + 2Eu²⁺` 自发（`E°cell = +0.41 V`）。投料 Zn 0.1 /
   Eu³⁺ 0.1 ⟹ Zn **过量**（只需 0.05）⟹ Eu²⁺ 应接近**全部** 0.1 mol。
2. **Zn(OH)₂ 是两性氢氧化物**：溶于过量强碱生成锌酸根
   `Zn(OH)₂ + 2OH⁻ → [Zn(OH)₄]²⁻`（教科书；这正是 `Zn²⁺` 与 `Fe³⁺`/`Al³⁺`
   分离的定性分析依据）。
3. **Zn(OH)₂ 也溶于氨水**生成氨合锌 `Zn(OH)₂ + 4NH₃ → [Zn(NH₃)₄]²⁺ + 2OH⁻`
   —— **这是 Zn²⁺ 与 Al³⁺ 的关键区别**（Al(OH)₃ 不溶于氨水）。

**为什么锁这些**：`EU01` 的终态正确性（Eu²⁺ 产量）**依赖走步能不能从锌的
形态乒乓里走出来**；而锌的形态化学正是第 287 轮定位的那一族。这三条锚
锁的是**产量/形态的化学事实**，不锁路线 —— 将来若改用别的机制（联立、
冻结策略调整）修好空转，它们仍然成立。

## 用法

    python tools/add_zn_eu_cases.py            # dry-run
    python tools/add_zn_eu_cases.py --check
    python tools/add_zn_eu_cases.py --write
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

_TAG = "第 287 轮扩充 · 锌/Eu 锚。"


def build() -> list[dict]:
    return [
        {
            "name": "EU6 Zn 0.1+Eu(NO3)3 0.1 还原度锚（Zn 过量）",
            "subs": [["Zn", 0.1], ["Eu(NO_3)_3", 0.1]],
            "has": {"Eu^{2+}": 0.08},
            "note": (
                _TAG + "**锁「还原到了多少」，不锁路线**。`E°(Zn²⁺/Zn) = −0.76 V`、"
                "`E°(Eu³⁺/Eu²⁺) = −0.35 V` ⟹ `Zn + 2Eu³⁺ → Zn²⁺ + 2Eu²⁺` "
                "自发（`E°cell = +0.41 V`）；Zn 0.1 对 Eu³⁺ 0.1 **过量**"
                "（只需 0.05）⟹ `Eu²⁺` 应接近全部 0.1 mol，下界取 0.08（80%）。"
                "**为什么本轮要加**：第 287 轮查明本条用例的走步在 "
                "`Zn²⁺ ⇌ [Zn(OH)₄]²⁻` 的**结构性乒乓**里空转 **2815 步**"
                "（幅度仅 ~3e-4 mol/步；四条检测器因阈值是**绝对**的而全部"
                "看不见；真正在跑的是 `[freeze-limit-cycle] → 强驱动跳过 → "
                "`[freeze-revive]` 的活锁）。**将来任何针对该族的修法都必须"
                "先过这一条** —— 它盯的是「Eu³⁺ 有没有被还原完」这个化学事实，"
                "而不是引擎的路线或步数。"),
        },
        {
            "name": "ZN1 Zn(OH)2 0.05+NaOH 0.3 两性溶解（锌酸根锚）",
            "subs": [["Zn(OH)_2", 0.05], ["NaOH", 0.3]],
            "ph": [12.8, 13.8],
            "has": {"[Zn(OH)_4]^{2-}": 0.035},
            "has_not": {"Zn(OH)_2": 0.01},
            "note": (
                _TAG + "**两性氢氧化物溶于过量强碱**："
                "`Zn(OH)₂ + 2OH⁻ → [Zn(OH)₄]²⁻`（教科书；`Zn²⁺` 与 `Fe³⁺` 的"
                "定性分离依据）。0.05 mol `Zn(OH)₂` 只需 0.1 mol OH⁻，投料 "
                "0.3 mol ⟹ OH⁻ **三倍过量** ⟹ 应几乎全溶。"
                "下界 `[Zn(OH)₄]²⁻ ≥ 0.035`（70%）取松；`Zn(OH)₂(s) ≤ 0.01` "
                "留 20% 余量。pH 由游离 OH⁻（0.3 − 0.1 = 0.2 mol/L 量级）定，"
                "取 `[12.8, 13.8]`（0.2 M 强碱 ≈ 13.3）。"
                "**为什么本轮要加**：锌酸根正是 `EU01` 乒乓的另一端 —— "
                "**它该不该存在、以多少量存在**必须有独立的化学锚。"),
        },
        {
            "name": "ZN2 Zn(OH)2 0.05+氨水 0.5 只微溶（**否掉的锚·留作数据待核**）",
            "subs": [["Zn(OH)_2", 0.05], ["NH_3", 0.5]],
            "has_not": {"[Zn(NH_3)_4]^{2+}": 0.01},
            "note": (
                _TAG + "⛔ **这一条原本想断言「Zn(OH)₂ 溶于氨水」，第 287 轮手算"
                "**否掉了那个断言**，故降级为「只微溶」的上界锚。**手算**："
                "`Ksp(Zn(OH)₂) ≈ 3e-17`、`β₄([Zn(NH₃)₄]²⁺) ≈ 1e9` ⟹ "
                "`Zn(OH)₂(s) + 4NH₃ ⇌ [Zn(NH₃)₄]²⁺ + 2OH⁻` 的 "
                "`K = Ksp·β₄ ≈ 3e-8`；设溶解 x，`4x³/(0.5)⁴ = 3e-8` ⟹ "
                "`x ≈ 7.8e-4 M` ⟹ **0.5 M 氨水下只有万分之八的 Zn(OH)₂ 溶解**。"
                "⟹ 教科书那条「Zn(OH)₂ 溶于氨水」是**浓氨水**事实（且通常从"
                "新制无定形沉淀出发），**不能按 0.5 M 断言** —— 这正是"
                "「先手算一个例子」纪律拦下的一次过强断言。"
                "**引擎实测**：`[Zn(NH₃)₄]²⁺ = 3e-05`、`Zn(OH)₂ = 0.04945`、"
                "pH 11.96 ⟹ 定性方向（几乎不溶）与手算一致，但比手算的 "
                "7.8e-4 还小 **26 倍**。**该 26 倍缺口记入数据待核**（"
                "`Ksp(Zn(OH)₂)` 的晶质/无定形取值 与 `β₄(氨合锌)` 的库值口径），"
                "**在核清之前不把氨合溶解写成正向断言**。"),
        },
        {
            "name": "ZN3 ZnCl2 0.05+NaOH 0.2 过量碱锌酸根（清洁体系）",
            "subs": [["ZnCl_2", 0.05], ["NaOH", 0.2]],
            "has": {"[Zn(OH)_4]^{2-}": 0.045},
            "has_not": {"Zn(OH)_2": 0.005},
            "note": (
                _TAG + "**与 `ZN1` 同化学、换成可溶锌盐投料**（`ZN1` 从固相出发）。"
                "`Zn²⁺ + 4OH⁻ → [Zn(OH)₄]²⁻`；投料 OH:Zn = **4:1** ⟹ "
                "应**全溶、无固相**。"
                "⛔ **不设 pH 区间（有意）**：4:1 时 0.2 mol OH⁻ **恰好被 0.05 mol "
                "Zn 全耗**（4 × 0.05 = 0.2）⟹ **没有游离碱**，pH 由锌酸根自身的"
                "水解定，而不是像 `ZN1` 那样由游离 OH⁻ 定。第 288 轮初版照抄了"
                "`ZN1` 的区间 [12.8, 13.8]，实测 **12.650** 翻红 —— **我的注文"
                "当时就写着'游离 OH⁻ ≈ 0'，区间却抄了另一条**，是自相矛盾。"
                "两例 pH 的**差**（`ZN1` 13.300 vs `ZN3` 12.650）恰好印证"
                "「游离碱有无」这件事，比给本条硬套一个区间更有信息量。"),
        },
        {
            "name": "ZN4 ZnCl2 0.05+NaOH 0.1 恰好沉淀（计量的另一侧）",
            "subs": [["ZnCl_2", 0.05], ["NaOH", 0.1]],
            "has": {"Zn(OH)_2": 0.045},
            "has_not": {"Zn^{2+}": 0.01},
            "note": (
                _TAG + "**与 `ZN3` 构成计量的两侧**：OH:Zn = **2:1** 恰好对应 "
                "`Zn²⁺ + 2OH⁻ → Zn(OH)₂(s)` ⟹ 应**定量沉淀**，游离 Zn²⁺ 极低"
                "（由 Ksp 与 pH 共同决定，取上界 0.01 = 投料的 20%，留两性复溶的"
                "余量）。**与 `ZN3` 一起**把「碱量决定产物是锌酸根还是氢氧化物」"
                "这条**两性化学的计量关系**钉住 —— 这正是 `EU01` 里锌物种"
                "反复进出账本的那一族化学。"),
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
        spec = ({"ph": c["ph"]} if "ph" in c else (c.get("has") or {}))
        print(f"{c['name'][:50]:<50} {str(spec)[:18]:>18}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<50} {'pH':>7} {'resid':>8} {'判读':>4}")
        for c in new:
            pr: dict = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            ph = r.get("final_pH")
            amt: dict = {}
            for e in (r.get("production") or []) + (r.get("final") or []):
                amt[e["name"]] = max(amt.get(e["name"], 0.0), e["mol"])
            ok = True
            if "ph" in c:
                ok = c["ph"][0] <= (ph if ph is not None else -99) <= c["ph"][1]
            for sp, lo in (c.get("has") or {}).items():
                ok = ok and amt.get(sp, 0.0) >= lo
            for sp, hi in (c.get("has_not") or {}).items():
                ok = ok and amt.get(sp, 0.0) <= hi
            print(f"{c['name'][:50]:<50} "
                  f"{(ph if ph is not None else float('nan')):>7.3f} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>4}")
            for sp in list((c.get("has") or {})) + list((c.get("has_not") or {})):
                print(f"      {sp:<22} {amt.get(sp, 0.0):.5g}")
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
