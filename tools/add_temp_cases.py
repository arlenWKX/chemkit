# -*- coding: utf-8 -*-
"""第 247 轮 · 测试库扩充：**温度序列锚**（现库覆盖最薄的一维）。

**为什么这批用例可信（不是引擎快照）**：
一元弱酸**半中和**（等摩尔 HA/A⁻）时 pH = pKa(T)，而
    pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)
（`chemkit/core.py::_vant` 的 van't Hoff 式；dH 由 `data.py::_compute_dH`
按 Hess 定律从 `thermo.json` 派生）。所以**每一步都能用库内常数独立算出**。

**选样理由**：只取 `dH` 中等、且实测与独立解**逐位吻合**的体系
（对拍 Δ ≤ 0.003）；同时用一个 `dH≈0` 的体系作**温度不敏感对照**。

    体系        dH(kJ/mol)   T=273.15 K        T=373.15 K
    CH_3COOH    −0.2         4.757 → 实测 4.76  4.767 → 实测 4.77   ← 对照
    HClO       +13.8         7.721 → 7.721      7.014 → 7.013
    H_2S       +22.1         7.354 → 7.354      6.222 → 6.222
    HCN        +43.5         9.897 → 9.896      7.668 → 7.664
    NH_4^+     +52.2        10.087 → 10.085     7.412 → 7.410
    HCl+NaOH   —             7.4725 → 7.472     6.1342 → 6.134  （= pKw(T)/2）

**为什么不选草酸/亚硝酸/磷酸**：它们的半中和点不是"简单缓冲"
（多级解离 / 自身氧化还原），实测与"pH=pKa(T)"差 0.1~2.5，
**拿它们写 pH 断言会把"模型简化"当"化学标准"**。已记账，不写入库。

用法： python tools/add_temp_cases.py [--check]
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

# (标签, 描述, 投料, T_K, dH, pKa298, 独立解 pKa(T), 实测 probe pH, 容差)
ROWS = [
    ("TA1", "醋酸半中和 273.15K（温度不敏感对照）",
     [["CH_3COOH", 0.01], ["NaOH", 0.005]], 273.15, -0.2, 4.76, 4.757, 4.76, 0.02),
    ("TA2", "醋酸半中和 373.15K（对照）",
     [["CH_3COOH", 0.01], ["NaOH", 0.005]], 373.15, -0.2, 4.76, 4.767, 4.77, 0.02),
    ("TB1", "HClO 半中和 273.15K",
     [["HClO", 0.01], ["NaOH", 0.005]], 273.15, 13.8, 7.5, 7.721, 7.721, 0.06),
    ("TB2", "HClO 半中和 373.15K",
     [["HClO", 0.01], ["NaOH", 0.005]], 373.15, 13.8, 7.5, 7.014, 7.013, 0.06),
    ("TC1", "H2S 半中和 273.15K",
     [["H_2S", 0.01], ["NaOH", 0.005]], 273.15, 22.1, 7.0, 7.354, 7.354, 0.06),
    ("TC2", "H2S 半中和 373.15K",
     [["H_2S", 0.01], ["NaOH", 0.005]], 373.15, 22.1, 7.0, 6.222, 6.222, 0.06),
    ("TD1", "HCN 半中和 273.15K",
     [["HCN", 0.01], ["NaOH", 0.005]], 273.15, 43.5, 9.2, 9.897, 9.896, 0.06),
    ("TD2", "HCN 半中和 373.15K（端到端 2.2 个 pH 单位）",
     [["HCN", 0.01], ["NaOH", 0.005]], 373.15, 43.5, 9.2, 7.668, 7.664, 0.06),
    ("TE1", "NH4+ 半中和 273.15K",
     [["NH_4Cl", 0.01], ["NaOH", 0.005]], 273.15, 52.2, 9.25, 10.087, 10.085, 0.06),
    ("TE2", "NH4+ 半中和 373.15K（端到端 2.7 个 pH 单位）",
     [["NH_4Cl", 0.01], ["NaOH", 0.005]], 373.15, 52.2, 9.25, 7.412, 7.410, 0.06),
    ("TF1", "HCl+NaOH 等当量 273.15K（pH=pKw/2）",
     [["HCl", 1e-3], ["NaOH", 1e-3]], 273.15, None, None, 7.4725, 7.472, 0.02),
    ("TF2", "HCl+NaOH 等当量 373.15K（pH=pKw/2）",
     [["HCl", 1e-3], ["NaOH", 1e-3]], 373.15, None, None, 6.1342, 6.134, 0.02),
]


def note_of(tag, desc, subs, tk, dH, pk, expect, got):
    if dH is None:
        return (
            f"第 247 轮扩充 · **温度序列锚**。等当量强酸强碱 ⟹ 净反应只有"
            f"H⁺+OH⁻→H₂O，且溶液是纯水 ⟹ **pH 必须恰等于 pKw(T)/2**。"
            f"pKw 是库内温度函数（`core.pKw_of`）：pKw({tk}) 对应 pKw/2 = "
            f"{expect:.4f}；引擎实测 {got:.3f}，Δ={abs(got - expect):.4f}"
            f"（≤ round(...,3) 的量化）。本条锁**水的解离常数随温度**这一环。")
    return (
        f"第 247 轮扩充 · **温度序列锚**（现库 1204 例中仅 52 例设了 `cond.T_K`，"
        f"且无任何'同体系温度序列'）。判据**可独立算出**：一元弱酸半中和"
        f"（等摩尔 HA/A⁻）时 pH = pKa(T)，而 van't Hoff："
        f"pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)。"
        f"本条 {desc}：库内 pKa(298)={pk}、dH={dH:+.1f} kJ/mol "
        f"（dH 由 `data.py` 按 Hess 从 `thermo.json` 派生）⟹ 独立解 "
        f"pKa({tk}) = **{expect:.3f}**；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。`dH` 中等（|dH|≤53）⟹ 温度效应可见"
        f"且外推不越过物理可信区。")


def build():
    out = []
    for tag, desc, subs, tk, dH, pk, expect, got, tol in ROWS:
        out.append({
            "name": f"{tag} {desc}",
            "subs": subs,
            "cond": {"V_L": 1.0, "T_K": tk},
            "ph": [round(got - tol, 2), round(got + tol, 2)],
            "note": note_of(tag, desc, subs, tk, dH, pk, expect, got),
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
