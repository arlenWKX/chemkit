# -*- coding: utf-8 -*-
"""第 291 轮 · **氯合一阶络合的 SIT-ε 锚（Cd/Zn）**。

## 动因（第 291 轮查明的 [MCl] 族摆锤机制）

第 291 轮查明：同一净反应 `[MCl]ⁿ⁺ ⇌ M + Cl⁻` 存在两个候选对象 ——
`templates.py` §4.4 的 `complex`/`decomplex`（此前 meta 无 eps，SIT 只加
Δz² 主项）与 `candidates.py` `build_derived` 的 `beta_direct`
（meta eps=−eps_rxn，加 Δz² **和** ε·I），两者零点相差整整 **ε·I**，
max(S) 拣选在两个候选各自的平衡态之间形成 0.6 mol 级大摆锤
（`E25` 实测 0.33×6.0=1.98 逐位吻合）。修法：**ε_rxn 属于平衡本身，
不属于候选风味** —— §4.4 的 direct 候选同步挂上 eps。

本组锚锁的是**带 ε 的条件稳定常数本身**（与引擎无关的库值推导）：

    logK_c(I) = logβ° + Δz²·D + ε·I，  D = 0.51·√I/(1+1.5√I)

其中 I 为引擎第 192 轮的**反应相关离子强度** `I_eff`（只计该反应的参与
物种，旁观强电解质不进项；纯盐体系中 `I_eff` 等于总 I）。

* `CD1`/`CD2` 用 **ε=+0.15** 的 `[CdCl]⁺`（IUPAC Part 4，logβ°=1.98），
  同一化学在 I≈0.14 与 I≈3.3 两个离子强度下份额从 ~80% 走到 ~99% ——
  **差分锚**，锁的是 D 与 ε·I 两项的**相对符号和量级**；
* `ZN5` 用 **ε=−0.14** 的 `[ZnCl]⁺`（IUPAC Part 5，logβ°=0.4），
  反号覆盖：ε<0 时条件常数被**压低**，2 M NaCl 中游离 Zn²⁺ 仍占多数。

**为什么锁这些**：修法之后 direct/derived 两候选的 S 在平衡态上逐位一致，
但「条件常数对不对」只能由库值独立判定。将来任何改动（SIT 口径、候选
构造、eps 传播）若让这三条的份额漂移，就是回归。

## 用法

    python tools/add_chloro_sit_cases.py            # dry-run
    python tools/add_chloro_sit_cases.py --check
    python tools/add_chloro_sit_cases.py --write
"""
from __future__ import annotations

import io
import json
import math
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

_TAG = "第 291 轮扩充 · 氯合 SIT-ε 锚。"


def D(I: float) -> float:
    """SIT 的 Debye-Hückel 项（298 K，与 tests.json 第 291 轮重裁 note 同口径）。"""
    s = math.sqrt(I)
    return 0.509 * s / (1 + 1.5 * s)


def build() -> list[dict]:
    return [
        {
            "name": "CD1 CdCl2 0.1 自络合锚（[CdCl]+ 占八成）",
            "subs": [["CdCl_2", 0.1]],
            "has": {"[CdCl]^+": 0.06},
            "has_range": {"Cd^{2+}": [0.005, 0.04]},
            "note": (
                _TAG + "**独立推导（只用库值）**：`Cd²⁺ + Cl⁻ ⇌ [CdCl]⁺`，"
                "logβ°=1.98（IUPAC Part 4）、ε=+0.15、Δz²=1−(4+1)=−4。"
                "设络合 x：`Cd²⁺=0.1−x`、`[CdCl]⁺=x`、`Cl⁻=0.2−x`。"
                "纯盐体系无旁观者，`I_eff`=总 I=½(4(0.1−x)+x+(0.2−x))=0.3−2x。"
                "迭代：x≈0.08 处 I≈0.14、D≈0.121 ⟹ logK_c=1.98−4×0.121"
                "+0.15×0.14≈1.52 "
                "⟹ K_c≈33 ⟹ x/(0.1−x)=33×(0.2−x)≈33×0.12≈4 ⟹ x≈0.08 "
                "（自洽）。⟹ **[CdCl]⁺≈0.08（八成）、游离 Cd²⁺≈0.02**。"
                "纯 CdCl₂ 溶液里一氯合镉就是主体形态 —— 这不是引擎结论，"
                "是 logβ°=1.98 这个库值的直接算术。"),
        },
        {
            "name": "CD2 CdCl2 0.1+NaCl 3 浓氯近全络合（ε>0 方向）",
            "subs": [["CdCl_2", 0.1], ["NaCl", 3.0]],
            "has": {"[CdCl]^+": 0.09},
            "has_range": {"Cd^{2+}": [0.0, 0.005]},
            "note": (
                _TAG + "**与 CD1 构成离子强度差分锚**。反应参与物种为 "
                "Cd²⁺/Cl⁻/[CdCl]⁺（Na⁺ 是旁观者，按 §192 `I_eff` 口径不进项）："
                "I_eff≈½(4×0.001+3.1+0.099)≈1.6、D≈0.223 ⟹ "
                "logK_c=1.98−4×0.223+0.15×1.6≈1.98−0.89+0.24=1.33 "
                "⟹ K_c≈21 ⟹ [CdCl]⁺/Cd²⁺=21×[Cl⁻≈3.1]≈66 "
                "⟹ **游离份额 ~1.5%**（Cd²⁺≈0.0015、[CdCl]⁺≈0.0985）。"
                "对比 CD1（I_eff≈0.14，游离 ~20%）：I 增大时 Δz²·D 项压低 K_c "
                "而 ε·I 项抬高 K_c，ε=+0.15 下后者胜 —— **两项的相对符号"
                "与量级被这对差分钉死**。"),
        },
        {
            "name": "ZN5 ZnCl2 0.1+NaCl 2 游离锌仍多数（ε<0 方向）",
            "subs": [["ZnCl_2", 0.1], ["NaCl", 2.0]],
            "has": {"Zn^{2+}": 0.025},
            "has_range": {"[ZnCl]^+": [0.005, 0.035]},
            "note": (
                _TAG + "**反号 ε 覆盖**：`Zn²⁺ + Cl⁻ ⇌ [ZnCl]⁺` logβ°=0.40、"
                "ε=−0.14（IUPAC Part 5）。Na⁺ 是旁观者（§192 `I_eff` 口径）："
                "I_eff≈½(4×0.035+2.07+0.018)≈1.1、D≈0.208 ⟹ "
                "logK_c=0.40−4×0.208−0.14×1.1≈0.40−0.83−0.15=−0.59 "
                "⟹ K_c≈0.26 ⟹ [ZnCl]⁺/Zn²⁺≈0.26×2.07≈0.54。"
                "二级/三级同口径可算：`ZnCl₂(aq)` β₂°=0.69、ε=−0.20、Δz²=−6 "
                "⟹ logK_c≈0.69−1.25−0.22=−0.78 ⟹ 比≈0.17×[Cl]²≈0.71；"
                "`[ZnCl₃]⁻` β₃°=0.18、ε=−0.05 ⟹ logK_c≈−1.12 ⟹ 比≈0.66。"
                "联立物料守恒迭代：四形态比 1 : 0.54 : 0.71 : 0.66 "
                "⟹ **游离 Zn²⁺≈0.034（~35%）、[ZnCl]⁺≈0.018**。"
                "**若无 ε 项**，一氯比会到 0.84（+55%）—— 本条盯的正是 "
                "ε<0 把条件常数**压低**这个方向。断言留量：Zn²⁺≥0.025、"
                "[ZnCl]⁺∈[0.005, 0.035]（推导值 ±ε 不确定度 ±0.02 对应的"
                "幅度）。"),
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
            for sp, lo in (c.get("has") or {}).items():
                ok = ok and amt.get(sp, 0.0) >= lo
            for sp, (lo, hi) in (c.get("has_range") or {}).items():
                ok = ok and lo <= amt.get(sp, 0.0) <= hi
            print(f"{c['name'][:50]:<50} "
                  f"{(ph if ph is not None else float('nan')):>7.3f} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>4}")
            for sp in list((c.get("has") or {})) + list((c.get("has_range") or {})):
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
