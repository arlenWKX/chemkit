# -*- coding: utf-8 -*-
"""第 286 轮 · **溶解态气体锚（半中和 · c 落在旧 `A_GAS` 地板以下）**。

## 动因（第 286 轮查明的根因）

`S_of` 的气体分支（`typ=4`，只发给**气体产物**）原式是

    a = max(min(c_g, H(T)·p_ext), A_GAS)          # A_GAS = P_RES/P_STD = 1/101.325

`A_GAS = 0.0098692`（`log10 = −2.00579`）本意是"**自产气体逸出后**惰性环境的
残余分压约定"，却被写成了**溶质活度的下限**：只要 `c_g < A_GAS` 就把活度**抬高**
到 `A_GAS`。于是同一个物种、同一个状态，只因为写在箭头哪一边就拿到两个标准态：

* **产物**侧命中 `typ=4` 气体分支 ⟹ `a = max(·, A_GAS)`（被抬高）
* **反应物**侧命中 `typ=2` 走 `_logc_of` 溶质口径 ⟹ `a = c_g`（真实值）

⟹ `S_fwd ≠ −S_rev`（热力学硬性违反）。本组锚把**溶解态气体**的活度口径钉住。

## 本组锚锁的化学事实（与引擎无关）

半中和的弱酸缓冲液 **pH = pKa**（`n(HA) = n(A⁻)` 时
`pH = pKa − log([HA]/[A⁻]) = pKa`），**与总浓度无关**。取三档浓度
（`CO₂` 账上浓度 0.002 / 0.005 / 0.008 M）必须给**同一个 pH** —— 这一条
正好把"`c_g < A_GAS` 就把活度抬到 `A_GAS`"的旧口径顶出来：三档全在
`A_GAS = 0.00987` 以下，旧口径下**每一档的抬高量都相同**（+0.2952），
但真实活度不同，所以三档的"引擎答案"会被同一把错误的尺子量。

`CO₂` 的亨利常数 `H = 3.4e-4 mol/(L·kPa)`（`data.py::HENRY`）
⟹ 饱和浓度 `H·p_ext = 0.0344 M`，本组最高只到 0.008 M，**不触发扫气**
（扫气判据 `cap = H·p_ext·V`）⟹ 体系等价闭口，半中和恒等式成立。

## 用法

    python tools/add_gas_solute_cases.py            # dry-run
    python tools/add_gas_solute_cases.py --check
    python tools/add_gas_solute_cases.py --write
"""
from __future__ import annotations

import io
import json
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

_TAG = "第 286 轮扩充 · 溶解态气体锚。"
_PKA_CO2 = 6.35      # 教科书 pKa1(CO2 + H2O ⇌ HCO3^- + H^+), 298.15 K
_PKA_H2S = 7.02      # 教科书 pKa1(H2S ⇌ HS^- + H^+), 298.15 K
_TOL = 0.15

_CO2_NOTE = (
    "**半中和 ⟹ pH = pKa**：`HCO₃⁻ + H⁺ ⇌ CO₂ + H₂O`，"
    "`n(CO₂) = n(HCO₃⁻)` 时 `pH = pKa1(CO₂)`。教科书值 "
    f"`pKa1(CO₂, 298.15 K) = {_PKA_CO2}`、`pKa1(H₂S) = {_PKA_H2S}`（298.15 K）"
    "⟹ 区间取 ±0.15（与 `HN1`–`HN7` 同口径）。"
    "**与总浓度无关**是本组的关键：三档 `CO₂` 浓度（0.002/0.005/0.008 M）"
    "必须给**同一个 pH**。"
)


def build() -> list[dict]:
    hi = _PKA_CO2 + _TOL
    lo = _PKA_CO2 - _TOL
    return [
        {
            "name": "GS1 NaHCO3 0.01+HCl 0.005 半中和（CO2 0.005 M < A_GAS）",
            "subs": [["NaHCO_3", 0.01], ["HCl", 0.005]],
            "ph": [lo, hi],
            "note": (
                _TAG + _CO2_NOTE +
                "**本条取 `c(CO₂) = 0.005 M`，落在旧地板 "
                "`A_GAS = P_RES/P_STD = 1/101.325 = 0.0098692` 以下**："
                "旧口径 `a = max(min(c_g, H·p_ext), A_GAS)` 会把活度从真实的 "
                "0.005（`log10 = −2.30103`）**抬到 0.0098692**"
                "（`log10 = −2.00579`），Δlog = **+0.29524** —— 与第 286 轮"
                "`tools/sprobe.py` 实测的 `S_of` 不对称量 **±0.2953** 逐位吻合。"),
        },
        {
            "name": "GS2 NaHCO3 0.016+HCl 0.008 半中和（CO2 0.008 M 贴上沿）",
            "subs": [["NaHCO_3", 0.016], ["HCl", 0.008]],
            "ph": [lo, hi],
            "note": (
                _TAG + _CO2_NOTE +
                "**本条取 `c(CO₂) = 0.008 M`，紧贴旧地板 `A_GAS = 0.0098692` 的下沿**"
                "（再高一点就够到地板、旧口径不再抬高）。"
                "与 `GS1`/`GS3` 一起把 `A_GAS` 这条**人为阈值**两侧夹住："
                "三档真实活度不同，若引擎答案随浓度分档跳变，就是阈值在作祟。"),
        },
        {
            "name": "GS3 NaHCO3 0.004+HCl 0.002 半中和（CO2 0.002 M 深处）",
            "subs": [["NaHCO_3", 0.004], ["HCl", 0.002]],
            "ph": [lo, hi],
            "note": (
                _TAG + _CO2_NOTE +
                "**本条取 `c(CO₂) = 0.002 M`，深在地板以下**"
                "（旧口径抬高量同样是 +0.29524，但真实活度只有地板值的 1/5）。"
                "本条与 `GS1`/`GS2` 的 pH 必须相同 —— 半中和缓冲的 pH 只由"
                "`[CO₂]/[HCO₃⁻]` 之比定，与绝对浓度无关。"),
        },
        {
            "name": "GS4 Na2S 0.008+HCl 0.012 半中和（H2S 0.004 M < A_GAS）",
            "subs": [["Na_2S", 0.008], ["HCl", 0.012]],
            "ph": [_PKA_H2S - _TOL, _PKA_H2S + _TOL],
            "note": (
                _TAG + "**换成另一种溶解态气体（H₂S）**，把'地板'问题从 CO₂ "
                "一个体系扩到第二个体系。`Na₂S 0.008 + HCl 0.012`："
                "0.008 mol H⁺ 先把 S²⁻ 全变 HS⁻，余下 0.004 mol H⁺ 再把一半 "
                "HS⁻ 变 H₂S ⟹ `n(H₂S) = n(HS⁻) = 0.004` ⟹ `pH = pKa1(H₂S) = "
                f"{_PKA_H2S}`。`H(H₂S) = 1.0e-3 mol/(L·kPa)` ⟹ 饱和 "
                "`H·p_ext = 0.1013 M ≫ 0.004 M`，不触发扫气。"
                "**为什么与 CO₂ 不同**：`H₂S` 的 `pKa1` 落在中性附近，"
                "而 `CO₂` 落在酸性侧 —— 两个不同的 pH 区间同时钉住口径。"),
        },
        {
            "name": "GS5 NaHCO3 0.01+HCl 0.005 半中和（账本组成锚·CO2）",
            "subs": [["NaHCO_3", 0.01], ["HCl", 0.005]],
            "has_range": {"CO_2": [0.004, 0.006], "HCO_3^-": [0.004, 0.006]},
            "note": (
                _TAG + "**锁的不是 pH，而是账本的组成** —— 半中和的定义就是"
                "「`n(CO₂) = n(HCO₃⁻) = c/2`」，c = 0.01 ⟹ 两者都必须落在 "
                "**[0.004, 0.006]**。"
                "**为什么这条对本轮重要**：第 278–285 轮反复撞到"
                "「**pH 对、账本没落实**」（`TE1` 实测 `NH₃` 只有 0.000131、"
                "`He` 停在 −0.004869 未中和；`HN4`/`HN5` 实测 0 步、`He` 未动）。"
                "气体分支的活度错误正是通过**走步决策**改账本的，"
                "所以必须有一条断言直接盯住气体物种在账上的量。"),
        },
    ]


def main(argv: list[str]) -> int:
    from chemkit.data import load_tables
    T = load_tables()
    new = build()
    path = os.path.join(ROOT, "chemkit", "data", "tests.json")
    with io.open(path, encoding="utf-8") as fh:
        db = json.load(fh)
    have = {c["name"] for c in db}
    for c in new:
        spec = ({"ph": c["ph"]} if "ph" in c
                else (c.get("has_range") or c.get("has") or {}))
        print(f"{c['name'][:52]:<52} {str(spec)[:16]:>16}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<52} {'pH':>9} {'resid':>8} {'判读':>4}")
        for c in new:
            pr: dict = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            ph = r.get("final_pH")
            ok = True
            if "ph" in c:
                ok = c["ph"][0] <= (ph if ph is not None else -99) <= c["ph"][1]
            if c.get("has_range"):
                fin0 = {e["name"]: e["mol"] for e in (r.get("final") or [])}
                for sp, (lo, hi) in c["has_range"].items():
                    ok = ok and lo <= fin0.get(sp, 0.0) <= hi
            print(f"{c['name'][:52]:<52} {(ph if ph is not None else float('nan')):>9.3f} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>4}")
        return 0

    if "--write" not in argv:
        print(f"\n[dry-run] 将新增 {sum(1 for c in new if c['name'] not in have)} 条"
              f"；加 --write 写入。")
        return 0
    add = [c for c in new if c["name"] not in have]
    db.extend(add)
    with io.open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(db, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"\n[写入] 新增 {len(add)} 条 ⟹ 共 {len(db)} 条")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
