# -*- coding: utf-8 -*-
"""第 245 轮 · 测试库扩充（第二批）：**阴离子不变性**的多样本锚。

第 244 轮入了 5 条（Cl⁻ vs NO₃⁻）。本轮按 `tools/anion_matrix.py` 的
**配对矩阵**（基准 = NO₃⁻，因为库内无 M–NO₃ 配体）扩到更多阴离子：
Br⁻ / I⁻ / ClO₄⁻ / SO₄²⁻。

**判据来自数据库自身，不是来自引擎输出**：
`beta.json` 里若**没有**该 金属–阴离子 的配合物条目，该阴离子在模型里
完全不参与 ⟹ 与 NO₃⁻ 基准的 pH **应当逐位相同**。矩阵实测：
**30 个"库内无 M–Y 条目"的配对全部 ΔpH = +0.0000**（14 个金属中心）；
而有条目的配对（Zn–Br/Zn–I/Cd–Br/Cd–I/Cu–Br）给出 0 或小幅正值。

本轮取 6 条，覆盖 4 种新阴离子 + 3 个此前几乎没有该阴离子覆盖的体系：

    MgI1  MgI2       6.70   （I⁻）
    CaP1  Ca(ClO4)2  7.00   （ClO₄⁻）
    NiS1  NiSO4      5.95   （SO₄²⁻）
    LaI1  LaI3       5.25   （I⁻，稀土）
    ZnP1  Zn(ClO4)2  5.50   （ClO₄⁻）
    AlB1  AlBr3      3.51   （Br⁻，三价、强水解端）

用法： python tools/add_anion_cases2.py [--check]
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

# (标签, 对照盐, 阴离子, 实测 pH, 同金属 NO3- 基准 pH)
ROWS = [
    ("MgI1", "MgI2", "I^-", 6.70, 6.70),
    ("CaP1", "Ca(ClO4)2", "ClO_4^-", 7.00, 7.00),
    ("NiS1", "NiSO4", "SO_4^{2-}", 5.95, 5.95),
    ("LaI1", "LaI3", "I^-", 5.25, 5.25),
    ("ZnP1", "Zn(ClO4)2", "ClO_4^-", 5.50, 5.50),
    ("AlB1", "AlBr3", "Br^-", 3.51, 3.51),
]


def note_of(tag, salt, an, ph, base):
    return (
        f"第 245 轮扩充 · **阴离子不变性**（{an} vs NO₃⁻，C={C:g} M，298.15 K）。"
        f"判据来自**数据库自身**：`beta.json` 内没有该金属的 M–{an} 配合物条目，"
        f"故该阴离子在模型里完全不参与，与 NO₃⁻ 基准的 pH **应当逐位相同**。"
        f"实测 {salt} = {ph:.2f}，同金属硝酸盐基准 = {base:.2f}，"
        f"Δ=0.0000（配对矩阵覆盖 14 个金属 / 30 个此类配对，全部 Δ=0.0000）。"
        f"本条同时锁住「引擎 pH 对**不该敏感**的变量不敏感」这一鲁棒性性质。"
    )


def build():
    out = []
    for tag, salt, an, ph, base in ROWS:
        out.append({
            "name": f"{tag} {salt} 溶液pH（阴离子不变性 {an}）",
            "subs": [[salt, C]],
            "ph": [round(ph - TOL, 2), round(ph + TOL, 2)],
            "changed": False,
            "note": note_of(tag, salt, an, ph, base),
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
