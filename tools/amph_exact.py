# -*- coding: utf-8 -*-
"""第 263 轮 · **独立于引擎**的二元酸碱精确解（判据），以及分支 4 合法性检验。

## 为什么需要它

第 263 轮定位到 `estimate_pH` 分支 4 的「两性中点」捷径
（`pH = (pKa1 + pKa2)/2`）在两类体系上给出错值：

* **P 类（纯两性盐）**：中点式由 `[H₂A] = [A²⁻]` 推出，**丢掉了 h 与 oh**；
  当 `oh` 与物种量可比时系统性偏碱（NaHS 0.01 M：中点 10.5、真值 9.4972）。
* **B 类（共轭缓冲对）**：账本里两性物种**与其共轭酸同时在账**时，
  `[H₂A]` 由投料定、不由歧化定 ⟹ 中点式的**推导前提整个不成立**
  （H₂S 半中和：中点 10.5、真值 7.0；CO₂/HCO₃⁻ 1:1：中点 8.35、真值 6.40）。

## 判据（本文件的核心，**不用引擎任何函数**）

二元酸 H₂A，总量 `C`、外加钠 `n_Na`、外加游离强酸 `n_H`，体积 `V`：

    D   = h² + Ka1·h + Ka1·Ka2
    [H₂A] = C·h²/D ;  [HA⁻] = C·Ka1·h/D ;  [A²⁻] = C·Ka1·Ka2/D
    电荷：n_Na + n_H + V·h  =  V·([HA⁻] + 2[A²⁻] + oh)

`f(pH)` 在 pH 上严格单调递减 ⟹ 唯一根，二分必收敛。
**只用库内 pKa 与 Kw**，不含引擎的任何启发式/分支。

## 用法

    python tools/amph_exact.py --calib     # 正/负对照（判据必须先能分开已知良例与已知病例）
    python tools/amph_exact.py --probe     # 与引擎三路 pH 对拍
"""
from __future__ import annotations

import math
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

from chemkit import speciation as _sp            # noqa: E402
from chemkit.acidbase import charge_pH           # noqa: E402
from chemkit.data import load_tables             # noqa: E402
from chemkit.speciation import estimate_state    # noqa: E402


def exact_pH(pKa1: float, pKa2: float, C: float, n_Na: float = 0.0,
             n_H: float = 0.0, V: float = 1.0, pKw: float = 14.0,
             lo: float = -2.0, hi: float = 18.0) -> float:
    """**引擎无关**的二元酸碱精确解（见模块 docstring 的方程）。"""
    Ka1, Ka2 = 10.0 ** (-pKa1), 10.0 ** (-pKa2)
    ct = C / V

    def f(pH: float) -> float:
        h = 10.0 ** (-pH)
        oh = 10.0 ** (pH - pKw)
        D = h * h + Ka1 * h + Ka1 * Ka2
        if D <= 0.0:
            return 1.0
        ha = ct * Ka1 * h / D
        a2 = ct * Ka1 * Ka2 / D
        return (n_Na + n_H) / V + h - ha - 2.0 * a2 - oh

    flo, fhi = f(lo), f(hi)
    if flo < 0.0:
        return lo
    if fhi > 0.0:
        return hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _pka(T, sp: str, T_K: float, which: int = 0) -> float:
    """库内第 `which` 级 pKa 在 T_K 的有效值（van't Hoff，与引擎同口径）。"""
    from chemkit.acidbase import eff_pka
    e = sorted(T.pka_acid[sp], key=lambda e: e["pka"])[which]
    pk = eff_pka((e["pka"],), (e.get("dH"),), T_K)[0]
    return pk / e["n"]


def calib(T) -> int:
    """正/负对照：判据必须判 NaHCO₃/NaH₂PO₄ 为"中点式够用"，
    判 NaHS / H₂S 半中和 / CO₂-HCO₃⁻ 缓冲为"中点式失效"。"""
    pKw = _sp.pKw_of(298.15)
    print(f"pKw(298.15) = {pKw}")
    rows = [
        # 标签, pKa1, pKa2, C, n_Na, 中点式预期, 说明
        ("NaHCO₃ 0.01 M（纯两性，教科书 pH≈8.35）",
         _pka(T, "CO_2", 298.15, 0), _pka(T, "HCO_3^-", 298.15, 0),
         0.01, 0.01, "正对照：oh 可忽略 ⟹ 中点式应准"),
        ("NaH₂PO₄ 0.01 M（纯两性）",
         _pka(T, "H_3PO_4", 298.15, 0), _pka(T, "H_2PO_4^-", 298.15, 0),
         0.01, 0.01, "正对照：Ka2/Ka1 ≈ 1e5 ⟹ 中点式应大致准"),
        ("NaHS 0.01 M（纯两性，**oh 不可忽略**）",
         _pka(T, "H_2S", 298.15, 0), _pka(T, "HS^-", 298.15, 0),
         0.01, 0.01, "负对照：oh=3e-5 与歧化量可比 ⟹ 中点式应偏碱"),
        ("H₂S 半中和 0.01/0.005（**1:1 缓冲**）",
         _pka(T, "H_2S", 298.15, 0), _pka(T, "HS^-", 298.15, 0),
         0.01, 0.005, "负对照：共轭对共存在 ⟹ 中点式前提不成立"),
        ("CO₂ + NaOH 1:1（**HCO₃⁻/CO₂ 缓冲**）",
         _pka(T, "CO_2", 298.15, 0), _pka(T, "HCO_3^-", 298.15, 0),
         0.01, 0.005, "负对照：同上（碳系）"),
        ("AB03 H₂S + 少量 NaOH（9:1 缓冲）",
         _pka(T, "H_2S", 298.15, 0), _pka(T, "HS^-", 298.15, 0),
         0.01, 0.001, "负对照：非 1:1 缓冲"),
    ]
    print(f"\n{'体系':<40} {'pKa1':>6} {'pKa2':>6} {'中点':>7} "
          f"{'精确':>7} {'中点−精确':>9}")
    for label, p1, p2, C, nNa, why in rows:
        mid = 0.5 * (p1 + p2)
        ex = exact_pH(p1, p2, C, n_Na=nNa, pKw=pKw)
        print(f"{label[:40]:<40} {p1:>6.3f} {p2:>6.3f} {mid:>7.3f} "
              f"{ex:>7.4f} {mid - ex:>+9.4f}   {why}")
    return 0


def probe(T) -> int:
    """与引擎三路 pH 对拍（引擎 ledger → 三条通路）。"""
    pKw = _sp.pKw_of(298.15)
    print(f"pKw(298.15) = {pKw}")
    cases = [
        # 标签, acid, amph, C, n_Na, engine ledger
        ("TC3 H₂S 半中和 298.15K", "H_2S", "HS^-", 0.01, 0.005,
         {"H_2S": 0.005, "HS^-": 0.005, "Na^+": 0.005}),
        ("纯 NaHS 0.01 M", "H_2S", "HS^-", 0.01, 0.01,
         {"HS^-": 0.01, "Na^+": 0.01}),
        ("AB03 H₂S+少量NaOH", "H_2S", "HS^-", 0.01, 0.001,
         {"H_2S": 0.009, "HS^-": 0.001, "Na^+": 0.001}),
        ("CO₂+NaOH 1:1", "CO_2", "HCO_3^-", 0.01, 0.005,
         {"CO_2": 0.005, "HCO_3^-": 0.005, "Na^+": 0.005}),
        ("纯 NaHCO₃ 0.01 M", "CO_2", "HCO_3^-", 0.01, 0.01,
         {"HCO_3^-": 0.01, "Na^+": 0.01}),
    ]
    print(f"\n{'体系':<26} {'精确(引擎无关)':>14} {'charge_pH':>10} "
          f"{'estimate_pH':>12} {'est−精确':>9} {'chg−精确':>9}")
    for label, acid, amph, C, nNa, led in cases:
        p1 = _pka(T, acid, 298.15, 0)
        p2 = _pka(T, amph, 298.15, 0)
        ex = exact_pH(p1, p2, C, n_Na=nNa, pKw=pKw)
        chg = charge_pH(led, 1.0, T, 298.15)
        _sp.PH_TAGS = []
        est = estimate_state(led, 0.0, 1.0, T, 298.15)[0]
        tags = list(_sp.PH_TAGS)
        _sp.PH_TAGS = None
        print(f"{label[:26]:<26} {ex:>14.4f} {chg:>10.4f} {est:>12.4f} "
              f"{est - ex:>+9.4f} {chg - ex:>+9.4f}   {'/'.join(tags)}")
    return 0


def main(argv: list[str]) -> int:
    T = load_tables()
    if "--probe" in argv:
        return probe(T)
    return calib(T)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
