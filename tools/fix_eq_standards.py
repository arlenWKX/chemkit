"""修正**错误标准**的 eq（v0.5.0 数据修正）。

第一类（6 条）原文**元素与电荷不守恒**（tools/eqcheck.py 精确有理核验）：
  21 AgBr       240NH_3 + 91AgBr + …（银氨族缺 6H/2N，round-1 审计已列）
  N30/D32       99SCN^- + 52Fe^{3+} -> …（Fe 52≠51、电荷 +57≠+54）
  NR95          497CaCO_3 + …（C+6、O+18、q−12）
  Y14/U10       65H^+ + 32Cl^- + …（同族巨系数，不守恒）
新值 = 引擎当前输出的守恒形式（逐条过 eqcheck）；原值作为历史保留在 note。

第二类（**叙述化石**，见 REASON）：原文守恒、但化学情境被旧 pH 的
**错误值**固化成了另一种叙述。Q05（H₂S 半中和）：旧引擎 pH 10.50 ⟹ 走
"碱性 ⟹ OH⁻ 叙述"分支；精确质子条件给出 7.00 = pKa₁（1:1 H₂S/HS⁻ 缓冲对）
⟹ 按引擎既定政策（终态 pH > 7 用 OH⁻、否则用 H⁺）落到**解离叙述**——
与同情境的 B12/H78（醋酸半中和，pH 4.76）**一致**；同一化学的 AB03
（NaOH 更少，pH 9.84）仍走 OH⁻ 叙述。政策讨论见 architecture §7 O-5。

用法：python tools/fix_eq_standards.py [--write]
"""
from __future__ import annotations

import json
import sys

PATH = "chemkit/data/tests.json"

FIX = {
    "21 AgBr+浓氨水（合并D24）": (
        "2NH_3 + AgBr -> Br^- + [Ag(NH_3)_2]^+",
        "240NH_3 + 91AgBr + 56H_2O -> 91Br^- + 91[Ag(NH_3)_2]^+ + 56NH_4^+ + 56OH^-"),
    "N30 FeCl3+KSCN": (
        "4SCN^- + 2Fe^{3+} -> [Fe(SCN)]^{2+} + [Fe(SCN)_3]",
        "99SCN^- + 52Fe^{3+} -> 21[Fe(SCN)]^{2+} + 18[Fe(SCN)_3] + 12[Fe(SCN)_2]^+"),
    "D32 FeCl3+NH4SCN": (
        "4SCN^- + 2Fe^{3+} -> [Fe(SCN)]^{2+} + [Fe(SCN)_3]",
        "99SCN^- + 52Fe^{3+} -> 21[Fe(SCN)]^{2+} + 18[Fe(SCN)_3] + 12[Fe(SCN)_2]^+"),
    "NR95 AgCl+CaCO3": (
        "9CaCO_3 + 6H_2O + AgCl -> 9Ca^{2+} + 6HCO_3^- + 6OH^- + 3CO_3^{2-} + Ag^+ + Cl^-",
        "497CaCO_3 + 343H_2O + 56AgCl -> 497Ca^{2+} + 343HCO_3^- + 343OH^- + 160CO_3^{2-} + 56Ag^+ + 56Cl^-"),
    "Y14 CuO+盐酸 溶解": (
        "8H^+ + 4Cl^- + 4CuO -> 4H_2O + 3Cu^{2+} + [CuCl_4]^{2-}",
        "65H^+ + 32Cl^- + 32CuO -> 32H_2O + 24Cu^{2+} + 8[CuCl_4]^{2-}"),
    "U10 Cu+H2O2+盐酸 溶解": (
        "8H^+ + 4Cl^- + 4Cu + 4H_2O_2 -> 8H_2O + 3Cu^{2+} + [CuCl_4]^{2-}",
        "65H^+ + 32Cl^- + 32Cu + 32H_2O_2 -> 65H_2O + 24Cu^{2+} + 8[CuCl_4]^{2-}"),
    # 第二类：叙述化石（原文守恒，情境被旧 pH 的错误值固化）
    "Q05 H2S过量+NaOH": (
        "H_2S + OH^- -> HS^- + H_2O",
        "H_2S -> HS^- + H^+"),
}

# 逐条替换说明（缺省用第一类的守恒说明）
REASON = {
    "Q05 H2S过量+NaOH":
        "【v0.5.0 标准修正·叙述化石回正】本条的 eq 曾在 pH 修正的同一轮里"
        "被改成引擎当时的输出「H₂S → HS⁻ + H⁺」，随后**回正为化学上更贴"
        "试剂的原式**「H₂S + OH⁻ → HS⁻ + H₂O」（投料是 NaOH）。原因是同一轮"
        "给 `_eq_signature` 加了 **H⁺/OH⁻/H₂O 规范形**：两种写法只差水的"
        "自电离关系，属于**呈现规范**而非化学，断言不该 over-specify 它。"
        "回正后本式与引擎输出（pH 7.00 落在叙述政策的酸性侧，写解离式）"
        "在规范形下同一，断言通过。政策讨论见 architecture §7 O-5。",
}


def main(write: bool) -> None:
    raw = open(PATH, encoding="utf-8", newline="").read()
    ind = 1
    for line in raw.splitlines():
        if line.strip() == "{":
            ind = len(line) - len(line.lstrip(" "))
            break
    rows = json.loads(raw)
    n = 0
    for c in rows:
        f = FIX.get(c["name"])
        if f is None:
            continue
        new, old = f
        if c.get("eq") == new:
            continue
        c["eq"] = new
        add = REASON.get(c["name"]) or (
            f"【v0.5.0 标准修正】原 eq 为「{old}」——**元素与电荷不守恒**"
            f"（tools/eqcheck.py 精确有理核验；如 99:52 族 Fe 52≠21+18+12=51、"
            f"电荷 +57≠+54）。该错误式子是引擎旧美化路径的产物被固化下来的。"
            f"新值取引擎当前输出的守恒形式（已过 eqcheck）——守恒修复见 "
            f"architecture.md：净方程装配器改为单一守恒闸门裁决删项。")
        c["note"] = ((c.get("note") or "") + "；" + add).strip("；")
        n += 1
    print(f"修正 {n} 条标准")
    if not write:
        print("（--write 写回）")
        return
    with open(PATH, "w", encoding="utf-8", newline="") as fh:
        fh.write(json.dumps(rows, ensure_ascii=False, indent=ind) + "\n")
    print(f"已写回 {PATH}（indent={ind}）")


if __name__ == "__main__":
    main("--write" in sys.argv)
