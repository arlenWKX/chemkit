# -*- coding: utf-8 -*-
"""第 248 轮 · 温度序列**中间档**（把线性检验做实）。

第 247 轮只做了 273.15 / 373.15 两端 —— 两端吻合**只能证明首尾对**，
中间若是曲线也能通过。本轮补 323.15 / 348.15 两档，**四点共线**才算证明
`pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)` 这条**直线**关系成立。

判据同上一批：一元弱酸半中和 ⟹ pH = pKa(T)（独立可算）。
实测（probe pH，3 位）：

    体系     pKa298  dH     273.15   298.15   323.15   348.15   373.15
    HClO     7.5    +13.8   7.721    7.500    7.313    7.152    7.013
    HCN      9.2    +43.5   9.896    9.197    8.607    8.101    7.664
    NH_4^+   9.25   +52.2  10.085    9.247    8.539    7.934    7.410

用法： python tools/add_temp_mid_cases.py [--check]
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

# (标签, 酸名, 投料, T_K, 描述, dH, pKa298, 独立解, 实测 probe pH)
ROWS = [
    ("TB3", "HClO", [["HClO", 0.01], ["NaOH", 0.005]], 323.15,
     "HClO 半中和 323.15K", 13.8, 7.5, 7.313, 7.313),
    ("TB4", "HClO", [["HClO", 0.01], ["NaOH", 0.005]], 348.15,
     "HClO 半中和 348.15K", 13.8, 7.5, 7.153, 7.152),
    ("TD3", "HCN", [["HCN", 0.01], ["NaOH", 0.005]], 323.15,
     "HCN 半中和 323.15K", 43.5, 9.2, 8.610, 8.607),
    ("TD4", "HCN", [["HCN", 0.01], ["NaOH", 0.005]], 348.15,
     "HCN 半中和 348.15K", 43.5, 9.2, 8.106, 8.101),
    ("TE3", "NH_4^+", [["NH_4Cl", 0.01], ["NaOH", 0.005]], 323.15,
     "NH4+ 半中和 323.15K", 52.2, 9.25, 8.543, 8.539),
    ("TE4", "NH_4^+", [["NH_4Cl", 0.01], ["NaOH", 0.005]], 348.15,
     "NH4+ 半中和 348.15K", 52.2, 9.25, 7.937, 7.934),
]
TOL = 0.06


def note_of(tag, desc, tk, dH, pk, expect, got):
    return (
        f"第 248 轮扩充 · **温度序列中间档**（第 247 轮只有 273/373 两端，"
        f"两端吻合只能证明首尾对，中间是曲线也会通过）。判据同批："
        f"一元弱酸半中和 ⟹ pH = pKa(T) = pKa(298) − dH/(R·ln10)·"
        f"(1/298.15 − 1/T)（van't Hoff 直线）。本条 {desc}："
        f"pKa(298)={pk}、dH={dH:+.1f} kJ/mol ⟹ 独立解 **{expect:.3f}**；"
        f"引擎实测 {got:.3f}，Δ={abs(got - expect):.4f}。"
        f"**四点共线**（273.15/298.15/323.15/348.15/373.15 K 五点全测）"
        f"才算证明这条关系在该条目上是**直线**、不是被两端凑出来的。")


def build():
    out = []
    for tag, acid, subs, tk, desc, dH, pk, expect, got in ROWS:
        out.append({
            "name": f"{tag} {desc}",
            "subs": [[n, m] for n, m in subs],
            "cond": {"V_L": 1.0, "T_K": tk},
            "ph": [round(got - TOL, 2), round(got + TOL, 2)],
            "note": note_of(tag, desc, tk, dH, pk, expect, got),
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
