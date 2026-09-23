# -*- coding: utf-8 -*-
"""第 259 轮 · 纯弱酸温度序列：**H₂S 五点**（第二个"敏感侧"体系）。

第 257 轮用醋酸（`|dH|=0.2`，不敏感侧）、第 258 轮用 HCN（`dH=+43.5`）。
本轮再补一个**中等敏感**的体系 **H₂S**（`dH=+22.1`）：
pKa(273→373) 从 7.3544 到 6.2218（跨 1.13），纯酸 pH 从 4.6776 到 4.1126（跨 0.57）。
**三个体系合起来覆盖"不敏感 / 中等 / 高敏感"三档**，比单点更能证明
`Ka(T)` 链路在全区间都对。

**判据（闭式，不依赖引擎）**：
    h² + Ka(T)·h − Ka(T)·c = 0  ⟹  pH = −log10((√(Ka(T)²+4Ka(T)c) − Ka(T))/2)

实测（0.01 M H₂S）：

    T/K      pKa(T)    独立 pH   引擎 pH   Δ
    273.15   7.3544    4.6776    4.678    +0.0004
    298.15   7.0000    4.5007    4.501    +0.0003
    323.15   6.7005    4.3512    4.351    −0.0002
    348.15   6.4440    4.2233    4.223    −0.0003
    373.15   6.2218    4.1126    4.113    +0.0004

用法： python tools/add_temp_pureacid3.py [--check]
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
    ("XN1", 273.15, 7.3544, 4.6776, 4.678),
    ("XN2", 298.15, 7.0000, 4.5007, 4.501),
    ("XN3", 323.15, 6.7005, 4.3512, 4.351),
    ("XN4", 348.15, 6.4440, 4.2233, 4.223),
    ("XN5", 373.15, 6.2218, 4.1126, 4.113),
]


def note_of(tag, tk, pk, expect, got):
    return (
        f"第 259 轮扩充 · **纯弱酸温度序列（中等敏感侧）**。纯酸温度序列已有"
        f"三个体系：醋酸 `dH=−0.2`（不敏感，第 257 轮）、HCN `dH=+43.5`"
        f"（高敏感，第 258 轮）、本条 **H₂S `dH=+22.1`**（中等敏感）——"
        f"三档合起来才说明 `Ka(T)` 链路在**整个敏感度范围**内都对。"
        f"判据闭式：`h² + Ka(T)·h − Ka(T)·c = 0` ⟹ "
        f"`pH = −log10((√(Ka(T)²+4Ka(T)c) − Ka(T))/2)`。"
        f"本条 T={tk} K：pKa(T)={pk:.4f}（273→373 跨 1.13）；"
        f"独立解 **{expect:.4f}**（断言以此为中心，±{TOL}）；"
        f"引擎实测 {got:.3f}，Δ={abs(got - expect):.4f}。")


def build():
    out = []
    for tag, tk, pk, expect, got in ROWS:
        out.append({
            "name": f"{tag} 纯H2S 0.01M {tk}K 溶液pH（纯酸温度序列·中等敏感）",
            "subs": [["H_2S", 0.01]],
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
