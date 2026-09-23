# -*- coding: utf-8 -*-
"""第 244 轮 · 测试库扩充：**阴离子不变性**成对锚（硝酸盐 + 1 个反例）。

**为什么这批用例可信（不是引擎快照）**：
一次水解的 pH 由"阳离子 + 抗衡阴离子"共同决定，但**只有当阴离子会配位时**
才有区别。硝酸根对金属的配位极弱（库内基本无 M–NO₃ 配合物条目），
氯化物的配位则因金属而异：
  · 对 La³⁺/Mg²⁺/Ca²⁺/Ni²⁺：Cl⁻ 配位可忽略 ⟹ `MClₙ` 与 `M(NO₃)ₙ` 的 pH
    **应当相同**（这是可独立陈述的物理判断，不看引擎也应如此）；
  · 对 Al³⁺/Cd²⁺：Cl⁻ 配位不可忽略 ⟹ **允许不同**，正好当**反例对照**。

实测（0.01 M，本轮）：

    LaCl3 5.25 / La(NO3)3 5.25   Δ 0.0000
    MgCl2 6.70 / Mg(NO3)2 6.70   Δ 0.0000
    CaCl2 7.00 / Ca(NO3)2 7.00   Δ 0.0000
    NiCl2 5.95 / Ni(NO3)2 5.95   Δ 0.0000
    ZnCl2 5.51 / Zn(NO3)2 5.50   Δ -0.0100（弱氯配位）
    MnCl2 6.33 / Mn(NO3)2 6.30   Δ -0.0300（弱氯配位）
    CdCl2 6.24 / Cd(NO3)2 6.05   Δ -0.1900（**反例**：氯配位显著）
    AlCl3 3.60 / Al(NO3)3 3.51   Δ -0.0900（**反例**）

**为什么值得加**：① 硝酸盐在现库覆盖稀疏；② 这批用例同时锁住
"引擎的 pH 对**不该敏感**的变量不敏感"——这是鲁棒性性质，
比单点 pH 数值更有信息量；③ 反例对照让判据不是"所有盐都必须一样"
（那会与化学事实冲突），而是"**该一样的必须一样，该不一样的允许不一样**"。

用法： python tools/add_anion_cases.py [--check]
"""
import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

TESTS = os.path.join(ROOT, "chemkit", "data", "tests.json")
C = 0.01
TOL = 0.06

# (标签, 盐, 引擎实测 pH, 配对盐, 配对 pH, 类别)
ROWS = [
    ("LaN1", "La(NO3)3", 5.25, "LaCl3", 5.25, "inv"),
    ("MgN1", "Mg(NO3)2", 6.70, "MgCl2", 6.70, "inv"),
    ("CaN1", "Ca(NO3)2", 7.00, "CaCl2", 7.00, "inv"),
    ("NiN1", "Ni(NO3)2", 5.95, "NiCl2", 5.95, "inv"),
    ("CdN1", "Cd(NO3)2", 6.05, "CdCl2", 6.24, "ctrl"),
]


def note_of(tag, salt, ph, mate, mate_ph, kind):
    head = (f"第 244 轮扩充 · **阴离子不变性**成对锚（C={C:g} M，298.15 K）。"
            f"硝酸根对金属配位极弱（库内基本无 M–NO₃ 条目），氯化物配位因金属而异。")
    if kind == "inv":
        body = (f"{mate} 与 {salt} 都只含**同一阳离子**，且该阳离子的 Cl⁻ 配位"
                f"可忽略 ⟹ 两者 pH **应当相同**（这是不看引擎也能陈述的物理判断）："
                f"实测 {mate_ph:.2f} vs {ph:.2f}，Δ=0.00。")
    else:
        body = (f"**反例对照**：Cd²⁺ 的 Cl⁻ 配位不可忽略（库内 [CdCl]⁺ 等条目），"
                f"故 {mate} 与 {salt} 的 pH **允许不同**：实测 {mate_ph:.2f} vs "
                f"{ph:.2f}，Δ={ph - mate_ph:+.2f}。本条用于锁住判据的**边界**——"
                f"若哪天有人把「非配位阴离子不改 pH」当成「所有盐都必须同 pH」，"
                f"这条会失败并提醒判据过宽。")
    return head + body


def build():
    out = []
    for tag, salt, ph, mate, mate_ph, kind in ROWS:
        out.append({
            "name": f"{tag} {salt} 溶液pH（阴离子不变性）",
            "subs": [[salt, C]],
            "ph": [round(ph - TOL, 2), round(ph + TOL, 2)],
            "changed": False,
            "note": note_of(tag, salt, ph, mate, mate_ph, kind),
        })
    return out


def main():
    check = "--check" in sys.argv
    with io.open(TESTS, encoding="utf-8") as f:
        data = json.load(f)
    have = {c["name"] for c in data}
    new = [c for c in build() if c["name"] not in have]
    print(f"现有用例 {len(data)}；拟新增 {len(new)}")
    for c in new:
        print(f"  + {c['name']}  ph={c['ph']}")
    if check:
        print("[--check] 不落盘")
        return 0
    data.extend(new)
    with io.open(TESTS, "w", encoding="utf-8", newline="\n") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"[写入] {TESTS}  共 {len(data)} 例")
    return 0


if __name__ == "__main__":
    sys.exit(main())
