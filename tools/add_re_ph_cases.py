# -*- coding: utf-8 -*-
"""第 242 轮 · 测试库扩充：稀土氯化物水溶液 pH 锚。

**为什么这批用例的答案是可信的（不是引擎快照）**：
体系是"一个三价金属中心 + 一个背景阴离子 + 水"，在 1e-4 M 下**远离固相
饱和**（La(OH)₃ 的饱和边界 pH 实测 9.28，而答案在 5.2~6.3），于是 pH 由
**一级水解 + 电荷平衡**唯一决定，可以闭式独立求解：

    [M] = C / (1 + β₁·[OH])
    z[M] + (z−1)[M(OH)] + [H⁺] = [OH⁻] + z·C

常数与推导都只用库内数据（`beta.json` 的 logβ₁ = pKw − log*K₁，
来源 Baes & Mesmer 1976，经 COST-NECTAR 同源核验）。独立解与引擎的
差 ≤ 0.006 pH 单位（唯一的系统性小偏差来自 SIT 活度层）。

**为什么这批用例值得加**：现库这 8 个体系共 33 例，**全部**是沉淀/氧化还原
（`+NaOH`、`+NaF`、`草酸盐`、`+冷水放氢`…）；**"稀土盐溶液本身的 pH"
一个都没有**——而这族体系最标准、最可验证的性质恰恰就是它。

写入后跑 `python tools/dev.py suite` 对拍：应为 1176 + 8 = 1184 例，
且**旧例零变动**。

用法： python tools/add_re_ph_cases.py [--check]
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

# el, 中文名, 盐, logβ₁（= pKw − log*K₁）, 独立 pH, 引擎 pH, pKw 敏感度
ROWS = [
    ("La", "镧", "LaCl3", 5.5, 6.2444, 6.2500),
    ("Ce", "铈", "CeCl3", 5.7, 6.1472, 6.1500),
    ("Nd", "钕", "NdCl3", 6.0, 6.0000, 6.0000),
    ("Eu", "铕", "EuCl3", 6.2, 5.9013, 5.9000),
    ("Gd", "钆", "GdCl3", 6.0, 6.0000, 6.0000),
    ("Ho", "钬", "HoCl3", 6.0, 6.0000, 6.0000),
    ("Tm", "铥", "TmCl3", 6.3, 5.8520, 5.8500),
    ("Lu", "镥", "LuCl3", 7.6, 5.2136, 5.2100),
]

TOL = 0.05          # 断言半宽：独立解与引擎差 ≤0.006，留 0.05 覆盖 SIT/舍入
C = 1e-4


def note_of(el, cn, lstar, iph, eph):
    return (f"第 242 轮扩充 · 稀土盐溶液 pH（现库该族 33 例**全是**沉淀/氧化还原，"
            f"缺\"盐溶液本身 pH\"这一最基本性质）。"
            f"体系 M³⁺ + 3Cl⁻ + H₂O，C={C:g} M，**远离固相饱和**"
            f"（{el}(OH)₃ 饱和边界 pH≈9.3，答案在 {eph:.2f}），"
            f"故 pH 由一级水解 + 电荷平衡唯一决定，可闭式独立求解："
            f"logβ₁ = pKw − log*K₁ = {lstar}（Baes & Mesmer 1976，"
            f"经 COST-NECTAR 同源核验）；独立解 pH = {iph:.4f}，"
            f"引擎 {eph:.4f}（差 {abs(eph - iph):.4f}，来自 SIT 活度层）。"
            f"断言按 ±{TOL} 取；**两套 pKw 约定下实测 pH 相同**"
            f"（298.15 K 处 14.0 与 14.00417 之差可忽略），故区间通用。")


def build():
    cases = []
    for el, cn, salt, lstar, iph, eph in ROWS:
        cases.append({
            "name": f"{el}1 {salt} 溶液pH（{cn}水解）",
            "subs": [[salt, C]],
            "ph": [round(eph - TOL, 2), round(eph + TOL, 2)],
            "changed": False,
            "note": note_of(el, cn, lstar, iph, eph),
        })
    return cases


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
