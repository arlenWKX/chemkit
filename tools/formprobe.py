# -*- coding: utf-8 -*-
"""第 272 轮 · **质子化形态约定探针**：`estimate_state` 返回的虚拟账本
`led_v` 把物种记成什么形态，与真值（Ka 分裂）差多少。

## 为什么需要

`tools/x8gap.py` 实测：`BR2` 退出态上 `estimate_state` 把
`ClO_3^- = 0.1866103` **整池**搬成 `HClO_3 = 0.1866103`，同时
`He 1.1196825 -> 0.9330721`（pH 两口径同为 0.030085，自洽）。

于是 `S_of(c, led_v, …)` 按候选写作式去查 `ClO_3^-` 槽位，得 **0**
⟹ `_logc_of` 落到 `ACT_FLOOR = 1e-12` ⟹ `logQ` 的 `+2·log10(a)` 项被
凭空推低 `≈ 2×(12−0.73) = 22.5`，与实测口径差 `−22.58` 吻合。
后果：求解器 `f(0)=+7.6e-05 ≈ 0`、`x*=0`（零推进），而残差口径同一态
报 `|S|=22.433`。

本工具问的是**上游**问题：`led_v` 记成什么形态，是不是化学事实？
用单酸体系（无其它弱组分干扰）把引擎的形态约定逼出来，与 Ka 真值对照。

用法：python tools/formprobe.py
"""

from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import math                                                 # noqa: E402

import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402

V = 1.0
T_K = 298.15
T = load_tables()


def exact(ka: float, c_tot: float) -> tuple[float, float, float]:
    """强酸/弱酸单质子酸 HA（总浓度 c_tot，全电离记账 ⟹ 总 H = c_tot）
    的精确解：[H+]² + Ka·[H+] − Ka·c_tot = 0。"""
    h = (-ka + math.sqrt(ka * ka + 4.0 * ka * c_tot)) / 2.0
    return (-math.log10(h), c_tot - h, h)      # (pH, [HA], [A-])


def probe(tag: str, anion: str, cation: str, c: float, ka: float) -> None:
    led = {anion: c, cation: c}
    p, lv, _ = spec.estimate_state(dict(led), c, V, T, T_K)
    got = {k: round(v, 6) for k, v in sorted(lv.items())
           if v > 0 and k != "H_2O"}
    e_pH, e_ha, e_a = exact(ka, c)
    print(f"--- {tag}  c={c:g} M  Ka={ka:g} ---")
    print(f"  引擎 pH = {p:8.4f}   真值 pH = {e_pH:8.4f}   "
          f"差 {p - e_pH:+.4f}")
    print(f"  引擎 led_v = {got}")
    print(f"  真值      {{'HA': {e_ha:.6g}, 'A-': {e_a:.6g}}}   "
          f"（阴离子占比真值 {e_a / c:.4f}）")
    a_eng = lv.get(anion, 0.0)
    print(f"  阴离子槽位：引擎 {a_eng:.6g}（占 {a_eng / c:.4f}）  "
          f"真值 {e_a:.6g}（占 {e_a / c:.4f}）  "
          f"⟹ logQ 项偏移 {2 * (math.log10(max(a_eng, 1e-12)) - math.log10(e_a)):+.4f}"
          "（ν=2 时）\n")


def _consistency() -> None:
    """规范形是否 **Ka 自洽**：把 `led_v` 里 100% 分子态的那一池按 Ka 分裂
    回 (酸, 阴离子) 混合态（He 相应回补），再问引擎 pH。

    · 若返回**同一个 pH** ⟹ 规范形只是记账重排，HH 分裂就是引擎的模型
      ⟹ `S_of` 侧补 HH 是**重建**引擎自己的态；
    · 若 pH 变了 ⟹ 规范形是**有损归一化**，得改上游。
    """
    print("=== 规范形 Ka 自洽性（BR2 退出态）===")
    base = {"Cl^-": 1.4330516, "K^+": 0.5, "Br_2": 0.2499799,
            "Cl_2": 0.1901685, "HBrO": 1.9627e-05}
    He_v = 0.9330721
    p_can, _, _ = spec.estimate_state(dict(base, **{"HClO_3": 0.1866103}),
                                      He_v, V, T, T_K)
    print(f"  规范形  {{HClO_3: 0.1866103}}  He={He_v:.7f}  ->  pH={p_can:.6f}")
    ka = 10.0 ** -1.0
    for frac in (0.5, 0.9, 0.915, 1.0):
        a = 0.1866103 * frac
        m = 0.1866103 - a
        p, _, _ = spec.estimate_state(dict(base, **{"HClO_3": m, "ClO_3^-": a}),
                                      He_v + a, V, T, T_K)
        print(f"  分裂回 {frac:.3f} 阴离子  {{HClO_3:{m:.6f}, ClO_3^-:{a:.6f}}}"
              f"  He={He_v + a:.6f}  ->  pH={p:.6f}")
    # HH/精确解在 pH=0.0301 处应有的阴离子占比
    h = 10.0 ** -0.0301
    print(f"  Ka={ka:g} 在 a_H={h:.4f} 处精确分裂：阴离子占比 = "
          f"{ka / (ka + h):.4f}")
    print()


def main() -> int:
    print("=== 单质子酸形态约定（账本按全电离记账：Anion=c, Cation=c, He=c）===\n")
    probe("氯酸 HClO3", "ClO_3^-", "K^+", 0.18661, 0.1)
    probe("氯酸 HClO3 稀", "ClO_3^-", "K^+", 1e-3, 0.1)
    probe("溴酸 HBrO3", "BrO_3^-", "K^+", 0.00169584, 0.1)
    probe("亚硝酸 HNO2 对照", "NO_2^-", "Na^+", 0.1, 10.0 ** -3.3)
    probe("乙酸 HAc 对照", "CH_3COO^-", "Na^+", 0.1, 10.0 ** -4.76)
    probe("硝酸 HNO3（有 spec 曲线）", "NO_3^-", "K^+", 1.0, 10.0 ** 1.0)
    _consistency()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
