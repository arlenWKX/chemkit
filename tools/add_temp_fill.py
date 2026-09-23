# -*- coding: utf-8 -*-
"""第 252 轮 · 温度序列**补齐空缺**（让每个体系都是完整五点）。

**为什么值得做**：第 247–251 轮陆续收了 7 个体系，但温度点参差
（H₂S 只做了 273/373 两端、HBrO 与 HClO 缺 298/348）。
**两点只能证明首尾对**（第 248 轮的教训），本轮统一补齐到
**273.15 / 298.15 / 323.15 / 348.15 / 373.15 五点**。

**判据**：一元弱酸半中和 ⟹ `pH = pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)`，
**断言区间以独立解为中心**（第 249 轮起口径）。

实测（独立解 / 引擎，全部 |Δ| ≤ 0.0095）：

    H2S   7.3544/7.354  7.0000/7.000  6.7005/6.700  6.4440/6.444  6.2218/6.222
    HBrO  9.0031/9.003  8.7000/8.699  8.4438/8.441  8.2245/8.219  8.0345/8.025
    HClO  7.7213/7.721  7.5000/7.500  7.3130/7.313  7.1528/7.152  7.0141/7.013

用法： python tools/add_temp_fill.py [--check]
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

# (标签, 酸, 盐名, T_K, dH, 独立 pKa(T), 引擎 pH)
ROWS = [
    ("TC3", "H_2S", "H_2S", 298.15, 22.1, 7.0000, 7.000),
    ("TC4", "H_2S", "H_2S", 323.15, 22.1, 6.7005, 6.700),
    ("TC5", "H_2S", "H_2S", 348.15, 22.1, 6.4440, 6.444),
    ("TG4", "HBrO", "HBrO", 298.15, 18.9, 8.7000, 8.699),
    ("TG5", "HBrO", "HBrO", 348.15, 18.9, 8.2245, 8.219),
    ("TB5", "HClO", "HClO", 298.15, 13.8, 7.5000, 7.500),
    ("TB6", "HClO", "HClO", 348.15, 13.8, 7.1528, 7.152),
]


def note_of(tag, acid, salt, tk, dH, expect, got):
    return (
        f"第 252 轮扩充 · **温度序列补齐**（把 {acid} 补成完整的五点："
        f"273.15/298.15/323.15/348.15/373.15 K）。判据同批：一元弱酸半中和 ⟹ "
        f"pH = pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)，"
        f"dH = {dH:+.1f} kJ/mol。**断言区间以独立解为中心**："
        f"pKa({tk}) = {expect:.4f}（±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。**补齐的理由**：两点只能证明首尾对，"
        f"中间若是曲线也会通过（第 248 轮教训）；五点共线才算证明这条直线。")


def build():
    out = []
    for tag, acid, salt, tk, dH, expect, got in ROWS:
        out.append({
            "name": f"{tag} {acid} 半中和 {tk}K（温度序列补齐）",
            "subs": [[salt, 0.01], ["NaOH", 0.005]],
            "cond": {"V_L": 1.0, "T_K": tk},
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, acid, salt, tk, dH, expect, got),
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
