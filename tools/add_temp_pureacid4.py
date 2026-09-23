# -*- coding: utf-8 -*-
"""第 260 轮 · 纯弱酸温度序列：**HClO 五点**（补齐敏感度梯度的低档）。

纯酸温度序列现有四个体系，**覆盖敏感度四档**：

| 体系 | dH (kJ/mol) | 273→373 的 pKa 跨度 | pH 跨度 |
|---|---|---|---|
| `CH₃COOH`（257 轮）| −0.2 | ~0.01 | ~0.005 |
| **`HClO`（本轮）** | **+13.8** | **0.71** | **0.35** |
| `H₂S`（259 轮）| +22.1 | 1.13 | 0.57 |
| `HCN`（258 轮）| +43.5 | 2.23 | 1.11 |

四档梯度比单点或双点更能说明 `Ka(T)` 链路**在整个敏感度范围内**都对。

**判据（闭式，不依赖引擎）**：
    h² + Ka(T)·h − Ka(T)·c = 0  ⟹  pH = −log10((√(Ka(T)²+4Ka(T)c) − Ka(T))/2)

实测（0.01 M HClO）：

    T/K      pKa(T)    独立 pH   引擎 pH   Δ
    273.15   7.7213    4.8609    4.861    +0.0001
    298.15   7.5000    4.7504    4.750    −0.0004
    323.15   7.3130    4.6570    4.657    +0.0000
    348.15   7.1528    4.5770    4.577    +0.0000
    373.15   7.0141    4.5077    4.508    +0.0003

用法： python tools/add_temp_pureacid4.py [--check]
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
    ("XO1", 273.15, 7.7213, 4.8609, 4.861),
    ("XO2", 298.15, 7.5000, 4.7504, 4.750),
    ("XO3", 323.15, 7.3130, 4.6570, 4.657),
    ("XO4", 348.15, 7.1528, 4.5770, 4.577),
    ("XO5", 373.15, 7.0141, 4.5077, 4.508),
]


def note_of(tag, tk, pk, expect, got):
    return (
        f"第 260 轮扩充 · **纯弱酸温度序列（低敏感档）**。纯酸温度序列现有"
        f"**四个体系覆盖敏感度四档**：醋酸 `dH=−0.2`（不敏感）、"
        f"**本条 HClO `dH=+13.8`（低敏感，273→373 跨 0.71）**、"
        f"H₂S `+22.1`（中等，跨 1.13）、HCN `+43.5`（高，跨 2.23）——"
        f"四档梯度比单点/双点更能说明 `Ka(T)` 链路在**整个敏感度范围**内都对。"
        f"判据闭式：`h² + Ka(T)·h − Ka(T)·c = 0` ⟹ "
        f"`pH = −log10((√(Ka(T)²+4Ka(T)c) − Ka(T))/2)`。"
        f"本条 T={tk} K：pKa(T)={pk:.4f} ⟹ 独立解 **{expect:.4f}**"
        f"（断言以此为中心，±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。")


def build():
    out = []
    for tag, tk, pk, expect, got in ROWS:
        out.append({
            "name": f"{tag} 纯HClO 0.01M {tk}K 溶液pH（纯酸温度序列·低敏感）",
            "subs": [["HClO", 0.01]],
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
