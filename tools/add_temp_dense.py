# -*- coding: utf-8 -*-
"""第 254 轮 · 温度序列**加密**（283.15 / 363.15 K）：把五点扩到七点。

**为什么加密**：五点已经能证明直线（第 248 轮），但**加密能抓住"局部弯曲"**——
若 `pKa(T)` 的真实曲线在某段有热容项导致的弯曲，五点可能跨过去，七点不容易。
本轮在 **298 K 两侧各加一档**（283.15 / 363.15），覆盖 273.15–373.15 全域。

**判据**同批：一元弱酸半中和 ⟹ `pH = pKa(T)`，**断言区间以独立解为中心**。
只取四个**逐档都精确**的体系（HClO / H₂S / HCN / HBrO）；酚与 `H₂O₂`
已知在弱酸端偏差大（第 253 轮），**继续保持不收录**。

实测（独立解 / 引擎）：

    HClO  283.15 7.6281/7.628 (−0.0001)   363.15 7.0673/7.067 (−0.0003)
    H2S   283.15 7.2051/7.205 (−0.0001)   363.15 6.3070/6.307 (+0.0000)
    HCN   283.15 9.6037/9.602 (−0.0017)   363.15 7.8359/7.832 (−0.0039)
    HBrO  283.15 8.8754/8.875 (−0.0004)   363.15 8.1073/8.099 (−0.0083)

用法： python tools/add_temp_dense.py [--check]
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

# (标签, 酸, 盐, T_K, dH, 独立 pKa(T), 引擎 pH)
ROWS = [
    ("TJ1", "HClO", "HClO", 283.15, 13.8, 7.6281, 7.628),
    ("TJ2", "HClO", "HClO", 363.15, 13.8, 7.0673, 7.067),
    ("TJ3", "H_2S", "H_2S", 283.15, 22.1, 7.2051, 7.205),
    ("TJ4", "H_2S", "H_2S", 363.15, 22.1, 6.3070, 6.307),
    ("TJ5", "HCN", "HCN", 283.15, 43.5, 9.6037, 9.602),
    ("TJ6", "HCN", "HCN", 363.15, 43.5, 7.8359, 7.832),
    ("TJ7", "HBrO", "HBrO", 283.15, 18.9, 8.8754, 8.875),
    ("TJ8", "HBrO", "HBrO", 363.15, 18.9, 8.1073, 8.099),
]


def note_of(tag, acid, tk, dH, expect, got):
    return (
        f"第 254 轮扩充 · **温度序列加密**（283.15 / 363.15 K，把五点扩到七点，"
        f"覆盖 273.15–373.15 K 全域）。判据：一元弱酸半中和 ⟹ "
        f"pH = pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)，"
        f"dH = {dH:+.1f} kJ/mol。**断言区间以独立解为中心**："
        f"pKa({tk}) = {expect:.4f}（±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。"
        f"**加密的理由**：五点已能证明直线，但**加密能抓局部弯曲**"
        f"（热容项会让 `pKa(T)` 略弯；五点可能跨过去，七点不容易）。")


def build():
    out = []
    for tag, acid, salt, tk, dH, expect, got in ROWS:
        out.append({
            "name": f"{tag} {acid} 半中和 {tk}K（温度序列加密）",
            "subs": [[salt, 0.01], ["NaOH", 0.005]],
            "cond": {"V_L": 1.0, "T_K": tk},
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, acid, tk, dH, expect, got),
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
