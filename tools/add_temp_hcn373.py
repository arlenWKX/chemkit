# -*- coding: utf-8 -*-
"""第 251 轮 · 温度序列：**HCN 高温档**（唯一"高温仍精确"的弱酸）+ 记账。

**为什么只加 HCN**：同一批对拍里，第 247/249 轮收的体系在高温端表现不一：

    HBrO  pKa(373)=8.034  引擎 8.025  Δ-0.009   ✓（已收，249 轮）
    HCN   pKa(373)=7.668  引擎 7.664  Δ-0.004   ✓（本条补高温档）
    C6H5OH pKa(373)=9.950 引擎 9.576  Δ-0.374   ✗（**不收**）
    H2O2  pKa(373)=10.515 引擎 9.799  Δ-0.716   ✗（**不收**）

酚（dH 缺，pKa 恒 9.95）在 373 K 差 −0.374，而它的 `estimate_pH` 与
`exact_proton_pH` 有分歧 —— 说明**呈现层精确解 `exact_proton_pH` 在
部分弱酸 / 高温档不准**（见 log 第 251 轮法医核验）。该分歧**按 pKa 与
温度增长**，已作为缺陷候选记账，**不写断言、不凑容差**。

用法： python tools/add_temp_hcn373.py [--check]
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

ROWS = [
    ("TD5", 348.15, 8.106, 8.101),
    ("TD6", 373.15, 7.668, 7.664),
]


def note_of(tag, tk, expect, got):
    return (
        f"第 251 轮扩充 · **温度序列：HCN 高档**（补到 298/323/348/373 四档齐）。"
        f"判据同批（一元弱酸半中和 ⟹ pH = pKa(T) = pKa(298) − dH/(R·ln10)·"
        f"(1/298.15 − 1/T)，HCN 的 dH = +43.5 kJ/mol）。**断言以独立解为中心**："
        f"pKa({tk}) = {expect:.3f}（±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。"
        f"**选样理由**：同一批对拍里，酚（pKa298=9.95、无 dH）在 373 K 差 −0.374、"
        f"`H_2O_2` 差 −0.716（**均未收录**），而 HCN 在四个温度档全部 ≤0.005 —— "
        f"只收逐档都精确的体系，**不拿容差去盖住已知偏差**。")


def build():
    out = []
    for tag, tk, expect, got in ROWS:
        out.append({
            "name": f"{tag} HCN 半中和 {tk}K（温度序列高档）",
            "subs": [["HCN", 0.01], ["NaOH", 0.005]],
            "cond": {"V_L": 1.0, "T_K": tk},
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, tk, expect, got),
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
        print(f"  + {c['name']}  T={c['cond']['T_K']}  ph={c['ph']}")
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
