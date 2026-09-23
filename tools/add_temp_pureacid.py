# -*- coding: utf-8 -*-
"""第 257 轮 · **纯弱酸的温度序列**（把两条已有维度交叉起来）。

**新在哪**：现库有（a）弱酸的**温度序列**（但都是"半中和"，有 NaOH/Na⁺），
（b）**纯弱酸**（但只在 298 K）。**没有"纯弱酸 × 温度"这条交叉**
—— 它走的是**既无外加阳离子、又变温**的电荷平衡
（`[H⁺] = [A⁻] + [OH⁻]`，且 `Ka` 随 T 变）。

**判据（闭式，不依赖引擎）**：
    h² + Ka(T)·h − Ka(T)·c = 0
    ⟹ pH = −log10( (√(Ka(T)² + 4·Ka(T)·c) − Ka(T)) / 2 )
其中 `Ka(T) = 10^−(pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T))`。

实测（0.01 M 醋酸，dH = −0.2 kJ/mol ⟹ pKa 几乎不随 T 变，pH 变化也小）：

    T/K      pKa(T)    独立 pH   引擎 pH   Δ
    273.15   4.7568    3.3875    3.387    −0.0005
    298.15   4.7600    3.3891    3.389    −0.0001
    323.15   4.7627    3.3904    3.390    −0.0004
    348.15   4.7650    3.3915    3.392    +0.0005
    373.15   4.7670    3.3925    3.392    −0.0005

**注**：醋酸的 `|dH|` 很小，所以这条序列**测的是"温度不敏感"那一侧**；
它同时锁住"pKa 几乎不变"与"纯酸电荷平衡"两件事。要测温度敏感侧，
已有 HCN/NH₄⁺/H₂S 的半中和序列（第 247–252 轮）。

用法： python tools/add_temp_pureacid.py [--check]
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

# (标签, T_K, 独立 pH, 引擎 pH)
ROWS = [
    ("XL1", 273.15, 3.3875, 3.387),
    ("XL2", 323.15, 3.3904, 3.390),
    ("XL3", 373.15, 3.3925, 3.392),
]


def note_of(tag, tk, expect, got):
    return (
        f"第 257 轮扩充 · **纯弱酸的温度序列**（把两条已有维度交叉起来：现库有"
        f"「弱酸×温度」但都是半中和（有 Na⁺），有「纯弱酸」但只测 298 K）。"
        f"本条走**既无外加阳离子、又变温**的电荷平衡（`[H⁺] = [A⁻] + [OH⁻]`）。"
        f"判据闭式可解：`h² + Ka(T)·h − Ka(T)·c = 0` ⟹ "
        f"`pH = −log10((√(Ka(T)²+4Ka(T)c) − Ka(T))/2)`，"
        f"`Ka(T)` 由 van't Hoff（dH = −0.2 kJ/mol）给出。"
        f"独立解 **{expect:.4f}**（断言以此为中心，±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。**注**：醋酸 dH 很小 ⟹ 这条测的是"
        f"「温度不敏感」那一侧；温度敏感侧已由 HCN/NH₄⁺/H₂S 的半中和序列覆盖。")


def build():
    out = []
    for tag, tk, expect, got in ROWS:
        out.append({
            "name": f"{tag} 纯CH_3COOH 0.01M {tk}K 溶液pH（纯酸温度序列）",
            "subs": [["CH_3COOH", 0.01]],
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
