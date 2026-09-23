# -*- coding: utf-8 -*-
"""第 255 轮 · 强酸半中和**成对锚**（把"强酸极限"这条明确行为锁住）。

**本轮定位到的规律（用物理投料、引擎真实 ledger 复核，第 254 轮教训）**：
对**完全解离**的强酸 `HA`，投料 `[HA: 0.01, NaOH: 0.005]` 时引擎把
多余强酸记成 `H_excess`，并让 `pH = −log10(H_excess) = −log10(0.005) = 2.3010`：

    HClO_4 / HBrO_4 / HSCN / HClO_3    pH 全部 = 2.3010，H_excess = 0.005
    H_2SeO_4（二元）                   pH 1.9800 ≈ −log10(0.0105)

**这条是不是缺陷？** —— **对"完全解离"的酸，它是对的**：
净反应只剩 `0.005 mol` 强酸，强酸全解离 ⟹ `[H⁺] = 0.005` ⟹ pH 2.3010。
**真正可疑的是"部分解离的强酸"**（`HIO_4` pKa 1.6、`HIO_3` 0.78、`HClO_2` 1.96）：
它们的 ledger 里同时有 `HA` 与 `A⁻`（如 `HIO_4: 0.005, IO_4^-: 0.005`），
pH 2.4300 对应游离 `[H⁺] = 3.7e-3`——**这个值应当由 `Ka` 与物料平衡定**，
而不是由"多余强酸"定。该分支已记账（见 log 第 255 轮），**不在本批断言内**。

本批**只锁住明确正确的那一半**（完全解离的三元组），外加它们的**碱侧对照**
（半量酸 + 定量碱 ⟹ 残余强碱）：

    X1 HClO4 + 半量 NaOH  → pH 2.3010（= −log 0.005）
    X2 HClO4 + 等当量 NaOH → pH 7.0000（纯水）
    X3 HClO4 + 1.5 倍 NaOH → pH 11.6990（= 14 + log 0.005）

用法： python tools/add_strongacid_pairs.py [--check]
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

# (标签, 描述, 投料, 独立解 pH, 引擎 pH, 推导)
ROWS = [
    ("X1", "HClO4 0.01 + 半量 NaOH（残余 0.005 强酸）",
     [["HClO_4", 0.01], ["NaOH", 0.005]], 2.3010, 2.301),
    ("X2", "HClO4 0.01 + 等当量 NaOH（只剩水）",
     [["HClO_4", 0.01], ["NaOH", 0.01]], 7.0000, 7.000),
    ("X3", "HClO4 0.01 + 1.5 倍 NaOH（残余 0.005 强碱）",
     [["HClO_4", 0.01], ["NaOH", 0.015]], 11.6990, 11.699),
]


def note_of(tag, desc, expect, got):
    return (
        f"第 255 轮扩充 · **强酸极限成对锚**（用物理投料 + 引擎真实 ledger 复核）。"
        f"高氯酸是**完全解离**的强酸（库内 pKa = −10），故半中和后净反应只剩"
        f"0.005 mol 强酸 ⟹ 强酸全解离 ⟹ `[H⁺] = 0.005` ⟹ **pH = −log10(0.005) "
        f"= 2.3010**，不需要任何缓冲公式。本条 {desc}：独立解 **{expect:.4f}**；"
        f"引擎实测 {got:.3f}，Δ={abs(got - expect):.4f}。"
        f"**成对的意义**：X1（酸过量）/ X2（恰中和）/ X3（碱过量）三点把"
        f"「强酸余量 → 游离 H⁺」这条记账在**两侧**都锁住，"
        f"且 X2 同时验证 `pH = pKw/2`。")


def build():
    out = []
    for tag, desc, subs, expect, got in ROWS:
        out.append({
            "name": f"{tag} {desc}",
            "subs": [[n, m] for n, m in subs],
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            "note": note_of(tag, desc, expect, got),
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
