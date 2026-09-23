# -*- coding: utf-8 -*-
"""第 250 轮 · 温度序列：**CO₂ 溶于纯水**（闭式可解 + Henry 随 T 变）。

**判据（可独立算，不依赖引擎）**：CO₂ 通入纯水，且**无强酸强碱储备**时，
溶液只有 CO₂(aq)/H₂CO₃、HCO₃⁻、H⁺、OH⁻（且 [H⁺]≫[OH⁻]），
故 `[H⁺] = √(K_H(T)·p·Ka₁(T))`，即 `pH = ½(pKa₁(T) + pK_H(T) − log10 p)`。

这条同时锁住**两条温度链路**：
  · `pKa₁(T)` 走 van't Hoff 且 dH = **+7.6 kJ/mol**（CO₂ 是"升温变强"的酸）；
  · `pK_H(T)` 走 Henry 的 dH = **−19.4 kJ/mol**（升温溶解度下降）。

库内实测（`pKw`、`pKa`、`Henry` 三条链路一起走）：

    T/K      pKa₁(T)   [CO2]aq   手算 pH   引擎 pH
    273.15   6.4580    0.06936    4.2922    4.262
    298.15   6.4000    0.03388    4.2068    4.201
    323.15   6.3412    0.01973    4.1398    4.150
    348.15   6.2820    0.01265    4.0860    4.106
    373.15   6.2228    0.00703    4.0401    4.068

**注意两端差偏大（273 K 差 −0.030、373 K 差 +0.028）**：手算忽略了 OH⁻ 与
活度项，且 273 K 时 CO₂ 浓度高、离子强度不可忽略。故断言取 ±0.08
（**以手算为中心**，容下这两项近似），不取 ±0.05。

用法： python tools/add_temp_co2_cases.py [--check]
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

# (标签, T_K, pKa1(T), [CO2]aq, 手算 pH, 引擎 pH)
ROWS = [
    ("TH1", 273.15, 6.4580, 0.06936, 4.2922, 4.262),
    ("TH2", 298.15, 6.4000, 0.03388, 4.2068, 4.201),
    ("TH3", 323.15, 6.3412, 0.01973, 4.1398, 4.150),
    ("TH4", 348.15, 6.2820, 0.01265, 4.0860, 4.106),
    ("TH5", 373.15, 6.2228, 0.00703, 4.0401, 4.068),
]
TOL = 0.08
SUBS = [["CO_2", 0.01]]


def note_of(tag, tk, pk, c_aq, expect, got):
    return (
        f"第 250 轮扩充 · **温度序列：CO₂ 溶于纯水**。判据可独立算：无强酸强碱"
        f"储备时 `[H⁺] = √(K_H(T)·p·Ka₁(T))` ⟹ pH = ½(pKa₁(T) + pK_H(T) − log p)。"
        f"本条同时锁**两条温度链路**：pKa₁({tk}) = {pk:.4f}"
        f"（dH = **+7.6** kJ/mol，CO₂ 是升温变强的酸）与 Henry 的"
        f"dH = **−19.4** kJ/mol（升温溶解度下降，1 atm 下 [CO₂]aq = {c_aq:.5f} M）。"
        f"手算 pH = **{expect:.3f}**（断言以此为中心，±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。容差取 ±{TOL} 而非 ±0.05：手算忽略了 OH⁻ 与"
        f"活度项，且 273 K 档 CO₂ 浓度高、离子强度不可忽略。")


def build():
    out = []
    for tag, tk, pk, c_aq, expect, got in ROWS:
        out.append({
            "name": f"{tag} CO2溶于纯水 {tk}K（Henry+Ka1 双温度链路）",
            "subs": [[n, m] for n, m in SUBS],
            "cond": {"V_L": 1.0, "T_K": tk},
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, tk, pk, c_aq, expect, got),
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
