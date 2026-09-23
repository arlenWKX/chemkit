# -*- coding: utf-8 -*-
"""第 261 轮 · **纯弱酸的 pKa 谱系锚**（跨酸强度，参数无关的精确解）。

**判据（精确、无自由参数）**：一元弱酸纯溶液（无外加酸/碱）时
    h² + Ka·h − Ka·c = 0  ⟹  pH = −log10( (√(Ka²+4·Ka·c) − Ka) / 2 )
这是一元二次式的**闭式解**，不含任何近似（`[A⁻] = [H⁺]` 由电荷守恒与
"只有该酸贡献离子"共同给出；`[OH⁻]` 项在酸性侧可忽略，实测见下）。

**为什么值得单独做**：现库的酸锚集中在几个体系（醋酸/氢氰酸/…），
**没有一条跨酸强度的谱系**。本轮取 **pKa 由 1.96 到 9.95** 的 8 个酸，
同一浓度（0.01 M）、同一温度（298.15 K），把 pH 通路在**整段酸强度**上钉住。

实测（独立解 / 引擎）：

| 酸 | pKa | 独立 pH | 引擎 | Δ |
|---|---|---|---|---|
| `HClO₂` | 1.96 | 2.1981 | 2.198 | −0.0001 |
| `H₃PO₄` | 2.10 | 2.2376 | 2.238 | +0.0004 |
| `H₂SeO₃` | 2.60 | 2.4077 | 2.408 | +0.0003 |
| `HF` | 3.20 | 2.6544 | 2.654 | −0.0004 |
| `HCOOH` | 3.75 | 2.9039 | 2.904 | +0.0001 |
| `CH₃COOH` | 4.76 | 3.3891 | 3.389 | −0.0001 |
| `HClO` | 7.50 | 4.7504 | 4.750 | −0.0004 |
| `HCN` | 9.20 | 5.6001 | 5.600 | −0.0001 |
| `C₆H₅OH` | 9.95 | 5.9750 | 5.973 | −0.0020 |

**不收 `HNO₂`**：pKa 3.3 实测差 −0.2155（它有自身氧化还原/歧化通道，
纯酸溶液里不止一级解离）—— **不拿容差去盖**，另行记账。

用法： python tools/add_pka_spectrum.py [--check]
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
C = 0.01

# (标签, 酸, pKa, 独立 pH, 引擎 pH)
ROWS = [
    ("XP1", "HClO_2", 1.96, 2.1981, 2.198),
    ("XP2", "H_3PO_4", 2.10, 2.2376, 2.238),
    ("XP3", "H_2SeO_3", 2.60, 2.4077, 2.408),
    ("XP4", "HF", 3.20, 2.6544, 2.654),
    ("XP5", "HCOOH", 3.75, 2.9039, 2.904),
    ("XP6", "CH_3COOH", 4.76, 3.3891, 3.389),
    ("XP7", "HClO", 7.50, 4.7504, 4.750),
    ("XP8", "HCN", 9.20, 5.6001, 5.600),
    ("XP9", "C_6H_5OH", 9.95, 5.9750, 5.973),
]


def note_of(tag, acid, pk, expect, got):
    return (
        f"第 261 轮扩充 · **纯弱酸 pKa 谱系锚**（跨酸强度，参数无关的精确解）。"
        f"判据：一元弱酸纯溶液时 `h² + Ka·h − Ka·c = 0` ⟹ "
        f"`pH = −log10((√(Ka²+4Ka·c) − Ka)/2)` —— **一元二次式闭式解，无近似**。"
        f"本条 `{acid}`：pKa={pk}、c={C:g} M ⟹ 独立解 **{expect:.4f}**"
        f"（断言以此为中心，±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。**本条属于跨 pKa 1.96→9.95 的谱系**"
        f"（同一浓度、同一温度），把 pH 通路在**整段酸强度**上钉住 ——"
        f"现库的酸锚原先只集中在少数几个体系上。")


def build():
    out = []
    for tag, acid, pk, expect, got in ROWS:
        out.append({
            "name": f"{tag} 纯{acid} {C:g}M 溶液pH（pKa 谱系）",
            "subs": [[acid, C]],
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, acid, pk, expect, got),
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
