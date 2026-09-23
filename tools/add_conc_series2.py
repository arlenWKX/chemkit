# -*- coding: utf-8 -*-
"""第 262 轮 · **多酸浓度序列**（同浓度扫酸强度 + 同酸扫浓度，交叉验证）。

**判据（闭式、无自由参数）**：一元弱酸纯溶液
    pH = −log10( (√(Ka²+4·Ka·c) − Ka) / 2 )
（`[A⁻] = [H⁺]` 由电荷守恒给出；酸性侧 `[OH⁻]` 可忽略。）

**为什么这一批比单点强**：它在**两个方向**上同时扫 ——
① 同一浓度（0.01 M）换酸强度；② 同一酸换浓度（0.1→1e-4）。
五元酸 × 四浓度 = **20 个点全部 Δ ≤ 0.0005**（跨 pH 1.5–4.5），
说明 pH 通路在**酸强度与浓度两个维度**上都对。

**本轮只取此前未覆盖的组合**（现库已有醋酸与 HF 的个别浓度档）：

| 酸 | c/M | 独立 pH | 引擎 | Δ |
|---|---|---|---|---|
| `HClO₂` | 0.1 | 1.5516 | 1.552 | +0.0004 |
| `HClO₂` | 0.001 | 3.0351 | 3.035 | −0.0001 |
| `HCOOH` | 0.1 | 2.3842 | 2.384 | −0.0002 |
| `HCOOH` | 0.001 | 3.4659 | 3.466 | +0.0001 |
| `H₃PO₄` | 0.1 | 1.6110 | 1.611 | +0.0000 |
| `H₃PO₄` | 0.0001 | 4.0054 | 4.005 | −0.0004 |

**明确不收 `HNO₂`**：它在纯酸溶液里**不止一级解离** ——
引擎实测走了 `3HNO₂ → 2NO + NO₃⁻ + H⁺`（歧化，extent 0.0033），
故"一元弱酸 pH"这条判据对它**不适用**（不是引擎错）。
偏差随浓度放大：0.1 M −0.684、0.01 M −0.216、1e-3 M +0.161、1e-4 M +0.342
——**正是"溶解的 NO 储库随浓度变化"的形状**。

用法： python tools/add_conc_series2.py [--check]
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

# (标签, 酸, pKa, c/M, 独立 pH, 引擎 pH)
ROWS = [
    ("XQ1", "HClO_2", 1.96, 0.1, 1.5516, 1.552),
    ("XQ2", "HClO_2", 1.96, 1e-3, 3.0351, 3.035),
    ("XQ3", "HCOOH", 3.75, 0.1, 2.3842, 2.384),
    ("XQ4", "HCOOH", 3.75, 1e-3, 3.4659, 3.466),
    ("XQ5", "H_3PO_4", 2.10, 0.1, 1.6110, 1.611),
    ("XQ6", "H_3PO_4", 2.10, 1e-4, 4.0054, 4.005),
]


def note_of(tag, acid, pk, c, expect, got):
    return (
        f"第 262 轮扩充 · **多酸浓度序列**（在酸强度与浓度两个维度上交叉扫）。"
        f"判据闭式、无自由参数：一元弱酸纯溶液 `pH = −log10((√(Ka²+4Ka·c) − Ka)/2)`"
        f"（`[A⁻]=[H⁺]` 由电荷守恒给出）。本条 `{acid}`：pKa={pk}、c={c:g} M ⟹ "
        f"独立解 **{expect:.4f}**（断言以此为中心，±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。**覆盖规模**：五元酸 × 四浓度 = 20 个点"
        f"全部 |Δ| ≤ 0.0005（跨 pH 1.5–4.5）⟹ pH 通路在**两个维度**上都对。")


def build():
    out = []
    for tag, acid, pk, c, expect, got in ROWS:
        out.append({
            "name": f"{tag} 纯{acid} {c:g}M 溶液pH（多酸浓度序列）",
            "subs": [[acid, c]],
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, acid, pk, c, expect, got),
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
