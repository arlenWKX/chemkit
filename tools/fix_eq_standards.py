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

# ⚠️ 第 243 轮修：**必须是有序 list of (name, new, old)，不能是 dict**。
# 原实现把四类修正合进一个 `FIX` dict，于是同一个用例若被两类修正覆盖
# （`NR95 AgCl+CaCO3` 就同时属于"整数倍放大"与"报告口径化石"两类），
# **后写的那条静默覆盖前一条**：
#   · 中间步骤丢失（只剩最终值）——对终态无害，但审计链断了；
#   · 更危险的是：若顺序反过来，"整数倍放大"那条就会把化石式当成最终值写回。
# 同时 `REASON` 只按 name 索引，导致 NR95 拿了通用文案、note 里写的是
# **99:52 Fe 族**的守恒问题（张冠李戴）。改成显式有序列表后两者都消失。
FIX: list[tuple[str, str, str]] = [
    # ── 第一类：整数倍放大（原文守恒，但配比无化学意义）────────────────
    ("21 AgBr+浓氨水（合并D24）",
     "2NH_3 + AgBr -> Br^- + [Ag(NH_3)_2]^+",
     "240NH_3 + 91AgBr + 56H_2O -> 91Br^- + 91[Ag(NH_3)_2]^+ + 56NH_4^+ + 56OH^-"),
    ("N30 FeCl3+KSCN",
     "4SCN^- + 2Fe^{3+} -> [Fe(SCN)]^{2+} + [Fe(SCN)_3]",
     "99SCN^- + 52Fe^{3+} -> 21[Fe(SCN)]^{2+} + 18[Fe(SCN)_3] + 12[Fe(SCN)_2]^+"),
    ("D32 FeCl3+NH4SCN",
     "4SCN^- + 2Fe^{3+} -> [Fe(SCN)]^{2+} + [Fe(SCN)_3]",
     "99SCN^- + 52Fe^{3+} -> 21[Fe(SCN)]^{2+} + 18[Fe(SCN)_3] + 12[Fe(SCN)_2]^+"),
    ("Y14 CuO+盐酸 溶解",
     "8H^+ + 4Cl^- + 4CuO -> 4H_2O + 3Cu^{2+} + [CuCl_4]^{2-}",
     "65H^+ + 32Cl^- + 32CuO -> 32H_2O + 24Cu^{2+} + 8[CuCl_4]^{2-}"),
    ("U10 Cu+H2O2+盐酸 溶解",
     "8H^+ + 4Cl^- + 4Cu + 4H_2O_2 -> 8H_2O + 3Cu^{2+} + [CuCl_4]^{2-}",
     "65H^+ + 32Cl^- + 32Cu + 32H_2O_2 -> 65H_2O + 24Cu^{2+} + 8[CuCl_4]^{2-}"),
    # ── 第二类：叙述化石（原文守恒，情境被旧 pH 的错误值固化）──────────
    ("Q05 H2S过量+NaOH",
     "H_2S + OH^- -> HS^- + H_2O",
     "H_2S -> HS^- + H^+"),
    # ── 第三类：报告口径化石（取自被 round(x,6) 量化过的净差）──────────
    ("P10 Fe(OH)3 不溶 NH4Cl",
     "3H^+ + Fe(OH)_3 -> 3H_2O + Fe^{3+}",
     "4H^+ + Fe(OH)_3 -> 3H_2O + Fe^{3+}"),
    # ── 第四类：标准自身内部不一致（两个比值互相矛盾）──────────────────
    # NR95 跨两类：先由"整数倍放大"落到 9:1 式，再由"报告口径化石"落到
    # 主通道守恒式。**两条都要保留**，顺序即审计链。
    ("NR95 AgCl+CaCO3",
     "9CaCO_3 + 6H_2O + AgCl -> 9Ca^{2+} + 6HCO_3^- + 6OH^- + 3CO_3^{2-} + Ag^+ + Cl^-",
     "497CaCO_3 + 343H_2O + 56AgCl -> 497Ca^{2+} + 343HCO_3^- + 343OH^- "
     "+ 160CO_3^{2-} + 56Ag^+ + 56Cl^-"),
    ("NR95 AgCl+CaCO3",
     "3CaCO_3 + 2H_2O -> 3Ca^{2+} + 2HCO_3^- + 2OH^- + CO_3^{2-}",
     "9CaCO_3 + 6H_2O + AgCl -> 9Ca^{2+} + 6HCO_3^- + 6OH^- + 3CO_3^{2-} + Ag^+ + Cl^-"),
    # ── 第四类：标准自身**内部不一致**（两个比值互相矛盾）──────────────
    ("F33 Na[Ga(OH)4]+CO2适量",
     "8[Ga(OH)_4]^- + 5CO_2 -> 8Ga(OH)_3 + 3CO_3^{2-} + 3H_2O + 2HCO_3^-",
     "11[Ga(OH)_4]^- + 7CO_2 -> 11Ga(OH)_3 + 4CO_3^{2-} + 4H_2O + 3HCO_3^-"),
]

# ⚠️ `FIX` 是**历史审计链**（含中间步），不是"当前应为何值"的目标表。
# 取每条的**最后一步**生成目标表：否则重复运行时会一直报"中间步需要修"
# （NR95 的 9:1 中间式就是这样被误报的）。
FINAL: dict[str, str] = {}
for _n, _new, _old in FIX:
    FINAL[_n] = _new

# 按 name 显式声明用的是**哪一类**修正的说明（缺省 = 第一类的守恒说明）。
# 用 class 名而不是裸文案：这样"某条用例跨两类"时不会张冠李戴。
CLASS2_NOTE = ("【v0.5.0 标准修正·叙述化石】原文守恒，但情境被旧 pH 的错误值固化，"
               "改回该情境下的真实反应式。")
CLASS3_NOTE = ("【v0.5.0 标准修正·报告口径化石】原式取自**被 round(x,6) 量化过**的"
               "净差——量化把迹量通道抬成了主通道；**其中有些原式还同时不守恒**"
               "（元素/电荷差见 tools/note_claim_audit.py 的逐条输出）。"
               "引擎改吃精确净差后按主通道给出守恒式（已过 eqcheck）。")
REASON = {
    "F33 Na[Ga(OH)4]+CO2适量":
        "【v0.5.0 标准修正·标准内部不一致】原 eq 的两个比值**互相矛盾**："
        "CO₃²⁻:HCO₃⁻ = 4:3 蕴含终态 pH = pKa₂ + log10(4/3) = **10.42**；"
        "而 CO₂:[Ga(OH)₄]⁻ = 7:11 蕴含溶解态镓 = 2 − 1.2×11/7 = **0.114 mol**，"
        "按本条自己的数据（pKsp 35.1、logβ₄ 37.6 ⟹ [Ga(OH)₄]⁻ = 10^2.5·[OH⁻]，"
        "与数据注记的 logK ≈ +2.5 一致）反推需要 pH ≈ **10.56**。两者不可同时成立。"
        "引擎末态（pH 10.46）四条独立关系同时成立：电荷 2 = x + [HCO₃⁻] + 2[CO₃²⁻] ✓、"
        "碳 1.2 = [HCO₃⁻] + [CO₃²⁻] ✓、碳酸平衡 [CO₃²⁻]/[HCO₃⁻] = 10^(pH−pKa₂) = 1.447 ✓、"
        "镓溶解度 [Ga(OH)₄]⁻ = 10^2.5·[OH⁻] = 0.0904 ✓。新值取引擎输出（已过 eqcheck）。",
    "Q05 H2S过量+NaOH": CLASS2_NOTE,
    "P10 Fe(OH)3 不溶 NH4Cl":
        "【v0.5.0 标准修正·报告口径化石】原 eq「4H⁺ + Fe(OH)₃ → 3H₂O + Fe³⁺」"
        "**本身电荷不平**（+4≠+3，tools/eqcheck.py 全库审计唯一的违规条），"
        "它是引擎旧路径吃 `round(x,6)` 后的净差印出来的：整条反应只有 "
        "1.17e-6 mol，H⁺ 净耗 3.515e-6 被量化成 4e-6。引擎现状：方程式装配改吃"
        "**精确净差**（`net_exact`），给出守恒的「3H⁺ + Fe(OH)₃ → 3H₂O + Fe³⁺」"
        "——教科书形式，也是本条应有的化学。",
    "NR95 AgCl+CaCO3":
        "【v0.5.0 标准修正·报告口径化石 + 整数倍放大，**跨两类**】"
        "原式 `497CaCO_3 + …` 是被 round(x,6) 量化后放大 56 倍印出来的"
        "（两条迹量通道都被量化）；第二步把 AgCl 通道（1.41e-5，占 CaCO₃ 通道 11%）"
        "按主通道呈现的接受线去掉，得 CaCO₃ 的溶解-水解守恒式（3:2:2:1）。"
        "两条通道都守恒，本条按引擎当前口径取主通道式。"
        "**注**：本条同时属于第一类（整数倍放大）与第三类（报告口径化石），"
        "审计链 = 497:343:56 → 9:6:1 → 3:2:2:1 式，中间步见 FIX 列表。",
}


def main(write: bool) -> None:
    raw = open(PATH, encoding="utf-8", newline="").read()
    ind = 1
    for line in raw.splitlines():
        if line.strip() == "{":
            ind = len(line) - len(line.lstrip(" "))
            break
    rows = json.loads(raw)
    by_name = {c["name"]: c for c in rows}
    n = 0
    skipped = []
    already = []
    # 按 FIX 的**顺序**逐条套用：同一条用例跨两类时，中间步也会被写进 note，
    # 审计链完整（旧 dict 实现会静默丢掉中间步）。
    for name, new, old in FIX:
        c = by_name.get(name)
        if c is None:
            print(f"!! 用例不存在: {name}")
            continue
        # ⚠️ 关键守卫（第 243 轮）：**只在用例确实有 `eq` 断言时才改**。
        # 有些用例的 `eq` 已被**有意移除**（改成只锁 `eq_has`/`has`，例如
        # `N30`/`D32`/`P10`）——此时再塞一个 `eq` 进去，等于**凭空新增一条
        # 断言**，而它并未经过"与引擎当前输出核对"这一步。
        if not c.get("eq"):
            skipped.append(name)
            continue
        if c.get("eq") == new:
            continue
        if FINAL.get(name) == c.get("eq"):
            # 已是该条的**最终**形态；本条只是审计链里的中间步，不重放
            already.append(name)
            continue
        c["eq"] = new
        add = REASON.get(name) or (
            f"【v0.5.0 标准修正】原 eq 为「{old}」——**元素与电荷不守恒**"
            f"（tools/eqcheck.py 精确有理核验；如 99:52 族 Fe 52≠21+18+12=51、"
            f"电荷 +57≠+54）。该错误式子是引擎旧美化路径的产物被固化下来的。"
            f"新值取引擎当前输出的守恒形式（已过 eqcheck）——守恒修复见 "
            f"architecture.md：净方程装配器改为单一守恒闸门裁决删项。")
        c["note"] = ((c.get("note") or "") + "；" + add).strip("；")
        n += 1
    print(f"修正 {n} 条标准")
    if already:
        print(f"已是最终形态、仅中间步未重放 {len(already)} 次: "
              + ", ".join(sorted(set(already))))
    if skipped:
        print(f"跳过 {len(skipped)} 条（`eq` 已被有意移除，仅锁 eq_has/has）: "
              + ", ".join(skipped))
    if not write:
        print("（--write 写回）")
        return
    with open(PATH, "w", encoding="utf-8", newline="") as fh:
        fh.write(json.dumps(rows, ensure_ascii=False, indent=ind) + "\n")
    print(f"已写回 {PATH}（indent={ind}）")


if __name__ == "__main__":
    main("--write" in sys.argv)
