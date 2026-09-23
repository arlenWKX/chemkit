# -*- coding: utf-8 -*-
"""第 251 轮 · 弱酸半中和的**法医核验**：引擎的终态自洽吗？

**背景**：第 249 轮用"半中和 ⟹ pH = pKa(T)"对拍时，发现三个体系差得远
（`H_2O_2` 373 K 差 +0.72、`H_3AsO_3` 323 K 差 +1.44 / 373 K 差 +3.34、
`H_3BO_3` 373 K 差 +0.045）。当时只说"模型不适用"，没说清**是数据缺口
还是引擎缺陷**。

**本脚本判法（第 243 轮 LaCl₃ 用过的那套）**：不看引擎报的 pH，而是看它的
**终态物种分布**，逐条核三件事：
  ① **解离常数关系**：`[A⁻][H⁺]/[HA]` 是否等于 10^(−pKa)？不等则引擎用的 pKa
     与库内不一致（**口径问题**）；
  ② **质量守恒**：Σ(含该酸根骨架的物种) = 投料？（**记账问题**）
  ③ **电荷守恒**：Σz·n + [H⁺] − [OH⁻] = 0？（**电荷问题**）

三条都过 ⟹ 引擎内部自洽，差异来自"半中和公式在此不适用"（多级解离/其它形态）。
有一条不过 ⟹ 是**可定位的缺陷**。

用法： python tools/weakacid_forensic.py [--write]
"""
import io
import json
import math
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

from chemkit.core import charge_of, elements_of                # noqa: E402

R_LN10 = 8.31446261815324e-3 * math.log(10.0)

# (标签, 酸, 骨架元素, 投料, T_K)
CASES = [
    ("H2O2-298", "H_2O_2", "H_2O_2", 0.01, 298.15),
    ("H2O2-373", "H_2O_2", "H_2O_2", 0.01, 373.15),
    ("H3AsO3-298", "H_3AsO_3", "H_3AsO_3", 0.01, 298.15),
    ("H3AsO3-323", "H_3AsO_3", "H_3AsO_3", 0.01, 323.15),
    ("H3AsO3-373", "H_3AsO_3", "H_3AsO_3", 0.01, 373.15),
    ("H3BO3-373", "H_3BO_3", "H_3BO_3", 0.01, 373.15),
    ("CH3COOH-298", "CH_3COOH", "CH_3COOH", 0.01, 298.15),
]


def main():
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()

    print("=== 弱酸半中和终态法医核验 ===")
    rows = []
    for tag, acid, skel, c, tk in CASES:
        e = [x for x in T.pka if x["acid"] == acid]
        if not e:
            print(f"  {tag}: 库内无 {acid} 的 pKa 条目")
            continue
        e = e[0]
        dH = e.get("dH")
        pk = e["pka"] - (dH / R_LN10 * (1 / 298.15 - 1 / tk)
                        if dH is not None else 0.0)
        pr = {}
        r = judge([{"name": acid, "mol": c},
                   {"name": "NaOH", "mol": c / 2}],
                  {"V_L": 1.0, "T_K": tk}, T, _probe=pr)
        fin = {x["name"]: x["mol"] for x in (r.get("final") or [])}
        pH = pr.get("pH")
        h = 10.0 ** (-pH)
        # ① 解离常数关系：找共轭碱
        base = e["base"]
        ha = fin.get(acid, 0.0)
        a = fin.get(base, 0.0)
        k_ratio = (a * h / ha) if ha > 0 else float("nan")
        pk_eff = -math.log10(k_ratio) if k_ratio > 0 else float("nan")
        # ② 质量守恒：骨架元素总量
        skel_el = elements_of(skel)
        tot = 0.0
        for nm, mol in fin.items():
            try:
                el = elements_of(nm)
            except Exception:                               # noqa: BLE001
                continue
            if all(el.get(k, 0) >= v for k, v in skel_el.items()):
                tot += mol
        # ③ 电荷守恒
        q = sum(charge_of(nm) * mol for nm, mol in fin.items()) + h
        print(f"\n  {tag}  pKa({tk})={pk:.3f}  引擎 pH={pH}")
        print("     物种: " + ", ".join(
            f"{k}={v:.6g}" for k, v in sorted(fin.items(),
                                              key=lambda kv: -kv[1])[:6]))
        print(f"     ① 解离关系: [A-][H+]/[HA] = {k_ratio:.4g} "
              f"⟹ pKa_eff = {pk_eff:.3f}  (库内 {pk:.3f}, "
              f"差 {pk_eff - pk:+.3f})")
        print(f"     ② 质量守恒: 骨架总量 = {tot:.6g} (投料 {c:g}) "
              f"→ {'✓' if abs(tot - c) < 1e-6 else '✗ 差 %.2e' % (tot - c)}")
        print(f"     ③ 电荷守恒: Σz·n + [H+] = {q:+.3e} "
              f"→ {'✓' if abs(q) < 1e-5 else '✗'}")
        rows.append((tag, pk, pH, pk_eff, tot, c, q))
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "weakacid_forensic.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"tag": t, "pKa_T": pk, "engine_pH": ph,
                        "pKa_eff": pe, "mass": m, "feed": c, "charge": q}
                       for t, pk, ph, pe, m, c, q in rows], f,
                      ensure_ascii=False, indent=1)
        print(f"\n[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
