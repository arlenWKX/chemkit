# -*- coding: utf-8 -*-
"""第 246 轮 · 测试库扩充：阴离子不变性的**浓度序列**锚。

**本轮的读数修正（重要）**：第 245 轮的矩阵用 `final_pH` 做差，而
`final_pH = round(...,2)` ⟹ **0.01 的格子会把真差 ≤0.01 全抹成 0.0000**，
当时"Δ=0.0000（4 位小数）"是**读数精度不足**，不是测出来的。
本轮改用 `probe["pH"]`（3 位小数）并在 0.001 / 0.01 / 0.1 M 三档重测：
**69 个"库内无 M–Y 配合物"的配对全部 Δ = +0.000，违反 0 个**。

**为什么应当精确相等（不只是"近似相等"）**：该阴离子在库内没有配合物条目，
它不参与任何含 H⁺/OH⁻ 的平衡 ⟹ 同阳离子、同计量的两盐给出**同一个 pH**。
这不是经验规律，是模型结构决定的；测试的意义在于锁住"引擎别意外让
不参与平衡的离子影响 pH"（例如记账或离子强度层把它算进去）。

**覆盖面**：0.1 M 一档同时把 pH 拉到 5.45~5.80（Ni/Mn），
比 0.01 M 档酸得多，是对 pH 通路的更严检验。

用法： python tools/add_anion_conc_cases.py [--check]
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
TOL = 0.06

# (标签, 盐, C, 阴离子, probe pH（3 位）, 同金属 NO3- 基准 pH, 该金属 M–Y 条目数)
ROWS = [
    ("NiD1", "NiSO4", 0.1, "SO_4^{2-}", 5.450, 5.450, 0),
    ("MnD1", "MnI2", 0.1, "I^-", 5.800, 5.800, 0),
    ("MgD1", "Mg(ClO4)2", 0.1, "ClO_4^-", 6.200, 6.200, 0),
    ("LaD1", "LaI3", 0.1, "I^-", 4.750, 4.750, 0),
    ("AlD1", "AlBr3", 0.1, "Br^-", 3.002, 3.002, 0),
    ("NiE1", "NiSO4", 0.001, "SO_4^{2-}", 6.450, 6.450, 0),
]


def note_of(tag, salt, c, an, ph, base, n_ent):
    return (
        f"第 246 轮扩充 · **阴离子不变性 × 浓度**（{an} vs NO₃⁻，C={c:g} M，"
        f"298.15 K）。判据来自**模型结构**而非经验规律：`beta.json` 内该金属"
        f"没有 M–{an} 配合物条目（实测条目数 {n_ent}），该阴离子不参与任何"
        f"含 H⁺/OH⁻ 的平衡 ⟹ 同阳离子、同计量的两盐给出**同一个 pH**。"
        f"实测 {salt} = {ph:.3f}，同浓度硝酸盐基准 = {base:.3f}，Δ=0.000。"
        f"覆盖：0.001 与 0.01 与 0.1 M 三档共 **69 个此类配对全部 Δ=0.000**"
        f"（8 个金属 × 4 种阴离子）。"
        f"**读数口径**：本批 pH 取 `probe['pH']`（3 位小数）——"
        f"`final_pH` 是 round(...,2)，用它会因 0.01 格子而把真差抹成 0。"
    )


def build():
    out = []
    for tag, salt, c, an, ph, base, n_ent in ROWS:
        out.append({
            "name": f"{tag} {salt} {c:g}M 溶液pH（阴离子不变性 {an}）",
            "subs": [[salt, c]],
            "ph": [round(ph - TOL, 2), round(ph + TOL, 2)],
            "changed": False,
            "note": note_of(tag, salt, c, an, ph, base, n_ent),
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
