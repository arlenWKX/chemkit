"""tools/gas_audit.py —— 气体活度约定对账（§7 X-21 的真因，v0.5.0）。

背景：`engine.S_of` 的气体分支是
    a = max(min(c_g, H(T)·p_ext), A_GAS)          # c_g = 账本摩尔 / V
其中两处**标准态混用**：
  · 上限 `min(c_g, H·p_ext)` 把"1 mol/L 溶解态"当标准态；物理上"1 atm 逸出"
    对应 `a = p/p° = 1`（H₂ 在 298 K 下差 ≈3 个 log 单位）；
  · 下限 `A_GAS`（≈1/101.3，**分压比**语义）被直接用作浓度活度的地板。
实测后果（NR57 Pb+ZnSO₄，锚定 pKw）：端点 `c_g = 7.9e-4 < A_GAS = 9.87e-3`
⟹ `a(H₂) = 9.87e-3` ⟹ 耦合平衡 `x_eq = 1 − sqrt(a/10^1.07) = 0.971`（实测端点
逐位吻合）；同一算式取 `a = p/p° = 1` 得 `x_eq = 0.708`——差 26 个百分点。

用法：
  python tools/gas_audit.py            # 常数对账（Henry/饱和浓度/A_GAS/标准态）
  python tools/gas_audit.py --scan     # 全库普查：产气用例的 a(引擎) vs a(p/p°)
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from chemkit.data import load_tables                        # noqa: E402
from chemkit.engine import A_GAS, P_EXT_KPA, henry_of, judge  # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402
from chemkit.candidates import X_MIN                        # noqa: E402

P0_KPA = 101.325          # 标准压力


def consts() -> None:
    T = load_tables()
    print(f"p_ext = {P_EXT_KPA} kPa   A_GAS = {A_GAS}（≈1/101.3 ⟹ 分压比语义）")
    print(f"{'气体':<8}{'H(T) mol/(L·kPa)':>20}{'c_sat=H·p_ext':>16}"
          f"{'a(p/p°) 折算':>16}")
    for g in sorted(T.gases):
        H = henry_of(T, g, 298.15)
        if H is None:
            print(f"{g:<8}{'（无 Henry 数据 ⟹ 恒定 A_GAS）':>38}")
            continue
        print(f"{g:<8}{H:>20.3e}{H * P_EXT_KPA:>16.3e}"
              f"{(P_EXT_KPA / P0_KPA):>16.4f}")


def scan(limit: int = 40) -> None:
    """全库普查：产气（账本残留 > X_MIN）用例的活度对照，按 |Δlog a| 排序。"""
    T = load_tables()
    rows = []
    for c in load_cases(None):
        cond = c.get("cond") or {"V_L": 1.0}
        V = float(cond.get("V_L", 1.0) or 1.0)
        try:
            r = judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T)
        except Exception:                                    # noqa: BLE001
            continue
        for e in r.get("final", []):
            s, m = e["name"], e["mol"]
            if s not in T.gases or m <= X_MIN:
                continue
            H = henry_of(T, s, 298.15)
            cg = m / V
            csat = None if H is None else H * P_EXT_KPA
            a_eng = A_GAS if H is None else max(min(cg, csat), A_GAS)
            a_pp = (A_GAS * (P_EXT_KPA / P0_KPA) if H is None
                    else max(min(cg / csat, 1.0), A_GAS) * (P_EXT_KPA / P0_KPA))
            rows.append((abs(__import__("math").log10(max(a_pp, 1e-30)
                                                     / max(a_eng, 1e-30))),
                         c["name"], s, cg, csat, a_eng, a_pp))
    rows.sort(reverse=True)
    print(f"产气用例 {len(rows)} 条（按 |Δlog a| 降序，前 {limit}）：")
    print(f"{'Δlog a':>8}  {'用例':<34}{'气':<7}{'c_g':>11}{'c_sat':>11}"
          f"{'a引擎':>11}{'a(p/p°)':>11}")
    for d, name, s, cg, csat, ae, ap in rows[:limit]:
        cs = "—" if csat is None else f"{csat:.3e}"
        print(f"{d:>8.3f}  {name[:34]:<34}{s:<7}{cg:>11.3e}{cs:>11}"
              f"{ae:>11.3e}{ap:>11.3e}")
    print("\n提示：Δlog a ≠ 0 的用例在'气体活度标准态统一'落地时会位移，")
    print("      须逐条按化学事实裁决（§7 X-21 下一轮主攻 ①）。")


if __name__ == "__main__":
    if "--scan" in sys.argv:
        scan()
    else:
        consts()
