# -*- coding: utf-8 -*-
"""第 249 轮 · 温度序列再扩：**HBrO**（判据中心取独立解，不取引擎值）。

**做法升级（与 247/248 轮不同）**：断言区间**以独立解为中心**
（`pH = pKa(T)`，由库内 pKa(298) + 派生 dH 算出），**而不是以引擎实测为中心**。
这样断言是真正独立的；引擎若漂移超过容差就会被打红。
（247/248 轮是以实测为中心写的，两者在通过时等价，但**独立中心更严**。）

实测对拍（独立解 / 引擎）：

    HBrO  pKa298=8.7  dH=+18.9
        273.15 K  9.003 / 9.003
        323.15 K  8.444 / 8.441
        373.15 K  8.034 / 8.025

**未收录**（同一批里实测差过大，说明引擎的 pH 模型对这些体系不适用，
不硬凑容差）：`H_2O_2`（373 K 差 +0.72）、`H_3AsO_3`（323 K 差 +1.44、
373 K 差 +3.34）、`H_3BO_3`（373 K 差 +0.045，边缘，暂不收）。
这三条**已记账**（见 log 第 249 轮）。

用法： python tools/add_temp_cases2.py [--check]
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

# (标签, 描述, T_K, 独立解 pH, 引擎实测 pH)
ROWS = [
    ("TG1", "HBrO 半中和 273.15K", 273.15, 9.003, 9.003),
    ("TG2", "HBrO 半中和 323.15K", 323.15, 8.444, 8.441),
    ("TG3", "HBrO 半中和 373.15K", 373.15, 8.034, 8.025),
]
TOL = 0.06
SUBS = [["HBrO", 0.01], ["NaOH", 0.005]]


def note_of(tag, desc, tk, expect, got):
    return (
        f"第 249 轮扩充 · **温度序列**（HBrO，pKa298=8.7、dH=+18.9 kJ/mol）。"
        f"判据：一元弱酸半中和 ⟹ pH = pKa(T) = pKa(298) − dH/(R·ln10)·"
        f"(1/298.15 − 1/T)。**本条断言区间以独立解为中心**"
        f"（pKa({tk}) = {expect:.3f}，±{TOL}），**不以引擎实测为中心** —— "
        f"这样断言与引擎无关，引擎漂移超容差会直接被打红。"
        f"引擎实测 {got:.3f}，Δ={abs(got - expect):.4f}。")


def build():
    out = []
    for tag, desc, tk, expect, got in ROWS:
        out.append({
            "name": f"{tag} {desc}",
            "subs": [[n, m] for n, m in SUBS],
            "cond": {"V_L": 1.0, "T_K": tk},
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, desc, tk, expect, got),
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
