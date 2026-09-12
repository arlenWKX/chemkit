"""把 6 条**错误标准**的 eq 修正为引擎的守恒输出（v0.5.0 数据修正）。

这 6 条的原文不守恒（tools/eqcheck.py 精确有理核验）：
  21 AgBr       240NH_3 + 91AgBr + …（银氨族缺 6H/2N，round-1 审计已列）
  N30/D32       99SCN^- + 52Fe^{3+} -> …（Fe 52≠51、电荷 +57≠+54）
  NR95          497CaCO_3 + …（C+6、O+18、q−12）
  Y14/U10       65H^+ + 32Cl^- + …（同族巨系数，不守恒）
新值 = 引擎当前输出的守恒形式（逐条过 eqcheck）；原值作为历史保留在 note。

用法：python tools/fix_eq_standards.py [--write]
"""
from __future__ import annotations

import json
import sys

PATH = "chemkit/data/tests.json"

FIX = {
    "21 AgBr+浓氨水（合并D24）": (
        "4.278NH_3 + 1.639AgBr + H_2O -> 1.639Br^- + 1.639[Ag(NH_3)_2]^+ + NH_4^+ + OH^-",
        "240NH_3 + 91AgBr + 56H_2O -> 91Br^- + 91[Ag(NH_3)_2]^+ + 56NH_4^+ + 56OH^-"),
    "N30 FeCl3+KSCN": (
        "84.602SCN^- + 44.72Fe^{3+} + 3OH^- -> 18.14[Fe(SCN)]^{2+} + 15.302[Fe(SCN)_3] + 10.278[Fe(SCN)_2]^+ + Fe(OH)_3",
        "99SCN^- + 52Fe^{3+} -> 21[Fe(SCN)]^{2+} + 18[Fe(SCN)_3] + 12[Fe(SCN)_2]^+"),
    "D32 FeCl3+NH4SCN": (
        "84.602SCN^- + 44.72Fe^{3+} + 3OH^- -> 18.14[Fe(SCN)]^{2+} + 15.302[Fe(SCN)_3] + 10.278[Fe(SCN)_2]^+ + Fe(OH)_3",
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
        add = (f"【v0.5.0 标准修正】原 eq 为「{old}」——**元素与电荷不守恒**"
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
