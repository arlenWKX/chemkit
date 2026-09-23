# -*- coding: utf-8 -*-
"""第 243 轮 · 测试库扩充：**反应计量/守恒**类锚（答案可精确推导，不依赖引擎）。

为什么选这一类（第 243 轮的方法论结论）：
    本轮本想继续用"独立解 pH"当判据，但在稀土/二价氯化物上**失败了**——
    那类体系的一级水解模型在"质量+电荷"联立下有**多根**，
    而引擎走的是"碱度平衡（虚拟强酸）"记账，两者可以相差几个 pH 单位
    （`LaCl₃` 0.01 M：引擎 5.25、弱根 6.88、实验观测 ≈5.2）。
    在这种"模型与数据都不确定"的地方写 pH 断言 = 把不确定写进标准。

    所以本轮改选**答案可精确推导**的断言：反应计量与守恒。
    这些量由投料直接决定，与活动度模型、与引擎实现都无关：
      · 水的自电离：pH 7 处 [H⁺]=[OH⁻]=1e-7，加酸必须让它**恰好**变多少；
      · 强酸强碱等当量：产物只有水，且 Na⁺/Cl⁻ 按投料守恒；
      · 弱酸盐完全中和：产物按计量守恒。

用法： python tools/add_stoich_cases.py [--check]
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

QW = 1e-14


def ph_of_strong_acid(c):
    """强酸（完全电离，无活动度修正）：[H+] = (c + sqrt(c²+4Kw))/2。"""
    return -__import__("math").log10((c + (c * c + 4 * QW) ** 0.5) / 2)


def ph_of_strong_base(c):
    """强碱对称式：[OH-] = (c + sqrt(c²+4Kw))/2 ⟹ pH = pKw − pOH。

    第 243 轮：只用来说明"1e-7 M NaOH 的精确解是 7.209"这一缺陷量级，
    **没有**对应的入库用例（引擎在该档报 7.00，属缺陷，不入库）。
    """
    oh = (c + (c * c + 4 * QW) ** 0.5) / 2
    return -__import__("math").log10(QW / oh)


CASES = [
    # ① 极稀强酸：**仅收录引擎正确的档位**（见文件头说明与 log 第 243 轮）
    #    1e-7 M 一档引擎报 pH 7.00（精确 6.79），已作为缺陷记账，**不入库**
    ("S2 1e-6 M HCl 极稀（水平衡起效）", [["HCl", 1e-6]],
     [round(ph_of_strong_acid(1e-6), 2) - 0.05,
      round(ph_of_strong_acid(1e-6), 2) + 0.05],
     {},
     "强酸极稀时 [H⁺] = (c+√(c²+4Kw))/2，不能直接取 c。"
     "1e-6 M 的精确解 pH = 5.9957（取 c 会得 6.00，差 0.004）；"
     "本档引擎正确。**更稀的 1e-7 M 是已知缺陷**（引擎报 7.00，精确 6.79），"
     "故此处只锁 1e-6 这一档——不为迁就引擎而放宽，也不把缺陷写成标准。"),
    # ② 强酸强碱等当量：产物只有水（计量守恒，与活动度无关）
    ("S4 HCl+NaOH 1e-3 等当量（只生成水）", [["HCl", 1e-3], ["NaOH", 1e-3]],
     [6.5, 7.5],
     {"Na^+": 1e-3, "Cl^-": 1e-3},
     "等当量强酸强碱 ⟹ 净反应只有 H⁺+OH⁻→H₂O；Na⁺/Cl⁻ 按投料**精确守恒**"
     "（各 1e-3 mol），pH 必须落回中性区。计量与守恒与活动度模型无关。"),
    # ③ 弱酸完全中和：产物按计量守恒
    ("S5 醋酸+NaOH 等当量（生成 NaAc）", [["CH_3COOH", 0.01], ["NaOH", 0.01]],
     [8.0, 9.2],
     {"Na^+": 0.01},
     "等当量中和 ⟹ 溶液是 0.01 M NaAc；pH 由 Ac⁻ 的 Kb = Kw/Ka 决定"
     "（Ka=1.8e-5 ⟹ Kb=5.6e-10 ⟹ 精确 pH≈8.37，区间取 8.0~9.2 容活动度）。"
     "Na⁺ 必须精确等于投料 0.01（守恒），引擎实测 8.38。"),
]


def build():
    out = []
    for name, subs, ph, has, note in CASES:
        c = {"name": name, "subs": subs, "ph": ph, "note": note}
        kept = {k: v for k, v in has.items() if v > 0}
        if kept:
            c["has"] = kept
        out.append(c)
    return out


def main():
    check = "--check" in sys.argv
    with io.open(TESTS, encoding="utf-8") as f:
        data = json.load(f)
    have = {c["name"] for c in data}
    new = [c for c in build() if c["name"] not in have]
    print(f"现有用例 {len(data)}；拟新增 {len(new)}")
    for c in new:
        print(f"  + {c['name']}  ph={c['ph']}  has={c.get('has')}")
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
