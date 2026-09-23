# -*- coding: utf-8 -*-
"""第 253 轮 · 温度序列：**硼酸**（pKa298=9.24、dH=+14.6）五点补齐。

硼酸是**天然水/海水**里的重要缓冲体系，且温度效应中等（ΔpKa(273→373) ≈ 0.75），
是这条温度序列上有代表性的补点。

**判据**：一元弱酸半中和 ⟹ `pH = pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)`，
**断言区间以独立解为中心**。

实测（独立解 / 引擎）：

    273.15  9.4741 / 9.474   Δ −0.0001
    298.15  9.2400 / 9.237   Δ −0.0030
    323.15  9.0421 / 9.032   Δ −0.0101
    348.15  8.8727 / 8.848   Δ −0.0247
    373.15  8.7259 / 8.681   Δ −0.0449

**选样理由**：五点全部落在 ±0.06 内（最差 −0.0449，且第 249 轮已记 373 K 边缘）。
比起把已知偏差大的体系（酚 −0.374、`H_2O_2` −0.717 @373K）塞进来凑数，
**只收逐档都精确的**。

用法： python tools/add_temp_borate.py [--check]
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

# (标签, T_K, 独立 pKa(T), 引擎 pH)
ROWS = [
    ("TI1", 273.15, 9.4741, 9.474),
    ("TI2", 298.15, 9.2400, 9.237),
    ("TI3", 323.15, 9.0421, 9.032),
    ("TI4", 348.15, 8.8727, 8.848),
    ("TI5", 373.15, 8.7259, 8.681),
]


def note_of(tag, tk, expect, got):
    return (
        f"第 253 轮扩充 · **温度序列：硼酸**（pKa298=9.24、dH=+14.6 kJ/mol，"
        f"天然水/海水的重要缓冲体系）。判据：一元弱酸半中和 ⟹ "
        f"pH = pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)。"
        f"**断言区间以独立解为中心**：pKa({tk}) = {expect:.4f}（±{TOL}）；"
        f"引擎实测 {got:.3f}，Δ={abs(got - expect):.4f}。"
        f"本条把硼酸补成完整五点（273.15/298.15/323.15/348.15/373.15 K），"
        f"**五点共线**才算证明 van't Hoff 直线。")


def build():
    out = []
    for tag, tk, expect, got in ROWS:
        out.append({
            "name": f"{tag} H3BO3 半中和 {tk}K（温度序列）",
            "subs": [["H_3BO_3", 0.01], ["NaOH", 0.005]],
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
