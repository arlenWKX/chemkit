# -*- coding: utf-8 -*-
"""第 258 轮 · **纯弱酸温度序列：HCN 五点**（补上"温度敏感侧"）。

第 257 轮做的"纯酸 × 温度"用的是**醋酸**（`|dH| = 0.2 kJ/mol`）——
pKa 几乎不随 T 变，只测到"温度不敏感"那一侧。本轮换成 **HCN**（`dH = +43.5`）：
`pKa(273→373)` 从 9.898 变到 7.668（**跨 2.23 个单位**），
而纯酸 pH 从 5.949 变到 4.834（**跨 1.11 个单位**）—— 这才真正压到温度的敏感侧。

**判据（闭式，不依赖引擎）**：
    h² + Ka(T)·h − Ka(T)·c = 0  ⟹  pH = −log10((√(Ka(T)²+4Ka(T)c) − Ka(T))/2)
    `Ka(T) = 10^−(pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T))`

实测（0.01 M HCN）：

    T/K      pKa(T)    独立 pH   引擎 pH   Δ
    273.15   9.8975    5.9488    5.949    +0.0002
    298.15   9.2000    5.6001    5.600    −0.0001
    323.15   8.6104    5.3053    5.305    −0.0003
    348.15   8.1055    5.0530    5.052    −0.0010
    373.15   7.6683    4.8345    4.834    −0.0005

**注**：与第 257 轮的醋酸档（`|dH|=0.2`，pH 几乎不变）**成对**，
一条测"不敏感侧"、一条测"敏感侧"，两条一起才说明温度链路在两侧都对。

用法： python tools/add_temp_pureacid2.py [--check]
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
TOL = 0.05

# (标签, T_K, pKa(T), 独立 pH, 引擎 pH)
ROWS = [
    ("XM1", 273.15, 9.8975, 5.9488, 5.949),
    ("XM2", 298.15, 9.2000, 5.6001, 5.600),
    ("XM3", 323.15, 8.6104, 5.3053, 5.305),
    ("XM4", 348.15, 8.1055, 5.0530, 5.052),
    ("XM5", 373.15, 7.6683, 4.8345, 4.834),
]


def note_of(tag, tk, pk, expect, got):
    return (
        f"第 258 轮扩充 · **纯弱酸温度序列（温度敏感侧）**。第 257 轮用的醋酸"
        f"`|dH| = 0.2 kJ/mol` ⟹ pKa 几乎不随 T 变，只测到「不敏感侧」；"
        f"本条换 HCN（`dH = +43.5`）：pKa(273→373) 从 9.898 变到 7.668"
        f"（**跨 2.23 个单位**），纯酸 pH 从 5.949 变到 4.834（**跨 1.11 个单位**）。"
        f"判据闭式：`h² + Ka(T)·h − Ka(T)·c = 0` ⟹ "
        f"`pH = −log10((√(Ka(T)²+4Ka(T)c) − Ka(T))/2)`。"
        f"本条 T={tk} K：pKa(T)={pk:.4f} ⟹ 独立解 **{expect:.4f}**"
        f"（断言以此为中心，±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。**与醋酸档成对**：一条不敏感侧、一条敏感侧，"
        f"两条一起才说明温度链路在两侧都对。")


def build():
    out = []
    for tag, tk, pk, expect, got in ROWS:
        out.append({
            "name": f"{tag} 纯HCN 0.01M {tk}K 溶液pH（纯酸温度序列·敏感侧）",
            "subs": [["HCN", 0.01]],
            "cond": {"V_L": 1.0, "T_K": tk},
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, tk, pk, expect, got),
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
