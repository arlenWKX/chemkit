# -*- coding: utf-8 -*-
"""**PbCl₂ 溶解度族的独立手算器**（只读库内数据，不复用引擎的任何求解逻辑）。

## 为什么需要它

`D45/J05/J06/J14/B24` 这五例的旧标全部来自**简单模型**

    s = (Ksp/4)^(1/3) ≈ 0.016 M

而该式的**前提是"溶解的铅全部以游离 Pb²⁺ 存在"**（电荷平衡 2[Pb²⁺] = [Cl⁻]）。
库内已有 Pb²⁺–Cl⁻ 络合阶梯（`beta.json`：β₁ = 1.6 一氯、β₃ = 2.0 三氯），
而 PbCl₂ 自身溶解就提供 Cl⁻（饱和液 [Cl⁻] ≈ 2s）⟹ 络合**必然**参与：

* 络合把**游离 Pb²⁺ 压低**（铅被 [PbCl]⁺/[PbCl₃]⁻ 分走）；
* 络合把**总溶解度抬高**（Cl⁻ 被络合消耗 ⟹ 需溶解更多固相才能维持 Ksp）。

⟹ 旧标的"**游离 Pb²⁺ ∈ [0.012, 0.02]**"与"**含络合的真值**"**不可能同时成立**。
本工具用**库内常数 + 独立写出的质量/电荷平衡**（二分自洽）算真值，
再把引擎实测并列 —— 供重裁判断"是引擎错还是标准是简单模型快照"。

## 口径声明（诚实标注）

* **理想（I→0）口径**：只施 Ksp 与累积 β，不含活度修正 ⟹ 与引擎的 SIT/DH
  口径**会有差异**；差异归因于活度层，不归因于"络合有没有参与"。
* 本工具**不读 `thermo.json` 的 dH** ⟹ 温度 ≠ 298.15 K 时它只把 `pKsp`
  当常数用，输出会标注 `(T≠298K: 未做 van't Hoff)`；T 例以引擎实测为准。

用法：

    python tools/pbcl_solubility.py              # 纯水 + 0.2M NaCl 同离子
    python tools/pbcl_solubility.py --cases      # 并列打印族内 5 例的引擎实测
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

DATA = os.path.join(ROOT, "chemkit", "data")


def load() -> tuple[float, list[tuple[int, float]]]:
    """库内 PbCl₂ 的 pKsp 与 Pb²⁺–Cl⁻ 累积 β 阶梯（只读 JSON）。"""
    with io.open(os.path.join(DATA, "ksp.json"), encoding="utf-8") as fh:
        ksp = [e for e in json.load(fh) if e["solid"] == "PbCl_2"]
    with io.open(os.path.join(DATA, "beta.json"), encoding="utf-8") as fh:
        beta = [e for e in json.load(fh)
                if e.get("center") == "Pb^{2+}" and e.get("ligand") == "Cl^-"]
    ladder = sorted((int(e["nu"]), float(e["logb"])) for e in beta)
    return 10.0 ** (-float(ksp[0]["pKsp"])), ladder


def solve(ksp: float, ladder: list[tuple[int, float]],
          cl_added: float = 0.0) -> dict:
    """独立质量/电荷平衡：cl = cl_added + 2s，[Pb²⁺] = Ksp/cl²，s = [Pb²⁺]·F。"""
    def f_of(s: float) -> tuple[float, float, float]:
        cl = cl_added + 2.0 * s
        F = 1.0 + sum(10.0 ** lb * cl ** nu for nu, lb in ladder)
        pb = ksp / cl ** 2
        return pb * F, pb, F
    lo, hi = 1e-9, 5.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f_of(mid)[0] > mid:
            lo = mid
        else:
            hi = mid
    s = 0.5 * (lo + hi)
    tot, pb_free, F = f_of(s)
    cl = cl_added + 2.0 * s
    dist = {f"[PbCl{nu if nu > 1 else ''}]": pb_free * 10.0 ** lb * cl ** nu
            for nu, lb in ladder}
    return {"s_total": s, "pb_free": pb_free, "cl": cl, "F": F,
            "dist": dist, "simple": (ksp / 4.0) ** (1.0 / 3.0),
            "cl_added": cl_added}


def show(tag: str, r: dict, ksp: float, T_K: float) -> None:
    print(f"\n== {tag}"
          + ("   (T≠298K: 未做 van't Hoff)" if abs(T_K - 298.15) > 0.01 else ""))
    print(f"   pKsp = {-__import__('math').log10(ksp):.2f}")
    print(f"   简单模型 s0 = (Ksp/4)^(1/3) = {r['simple']:.4f} M"
          f"   ← 旧标用的就是这个")
    print(f"   含络合 总溶 s = {r['s_total']:.4f} M   游离 Pb²⁺ = {r['pb_free']:.4f} M"
          f"   [Cl⁻] = {r['cl']:.4f} M   F = {r['F']:.3f}")
    print(f"   ⟹ 总溶/游离 = {r['s_total'] / r['pb_free']:.2f} 倍"
          f"（络合把游离压低、把总溶抬高）")
    for k, v in r["dist"].items():
        print(f"      {k:<12} {v:.4f} M   占总溶 {v / r['s_total'] * 100:5.1f}%")


def main(argv: list[str]) -> int:
    ksp, ladder = load()
    print(f"库内 Pb²⁺–Cl⁻ 阶梯（nu, logβ°）：{ladder}")
    show("纯水 · 298.15 K（D45/J05 的体系）", solve(ksp, ladder), ksp, 298.15)
    show("0.2 M NaCl 同离子（B24 的体系）", solve(ksp, ladder, 0.2), ksp, 298.15)

    if "--cases" in argv:
        from chemkit.data import load_tables
        from chemkit.engine import judge
        from chemkit.testsuit import load_cases
        T = load_tables()
        cases = load_cases(None)
        want = ("D45", "J05", "J06", "J14", "B24")
        print("\n" + "=" * 70)
        print("引擎实测（并列）：")
        for c in cases:
            if c["name"].split()[0] not in want:
                continue
            pr: dict = {}
            r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            led = dict(pr.get("ledger") or {})
            pb_spec = {k: v for k, v in led.items() if "Pb" in k and v > 1e-9}
            tot = sum(pb_spec.values())
            T_K = float((c.get("cond") or {}).get("T_K", 298.15))
            print(f"\n  {c['name']}   T={T_K}")
            print(f"    账本含铅物种："
                  + ", ".join(f"{k}={v:.5g}" for k, v in
                              sorted(pb_spec.items(), key=lambda x: -x[1])))
            print(f"    总溶铅 = {tot:.5f} M   游离 Pb²⁺ = {led.get('Pb^{2+}', 0.0):.5f}"
                  f"    Cl⁻ = {led.get('Cl^-', 0.0):.5f}"
                  f"    余固相 = {led.get('PbCl_2', 0.0):.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
