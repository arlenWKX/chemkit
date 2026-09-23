# -*- coding: utf-8 -*-
"""第 247 轮 · 温度维审计：**van't Hoff 路径是否真的在工作**。

**为什么查这里**：现库 1204 例里只有 **52 例**设了 `cond.T_K`（多为 350 K 的
氯气/二氧化锰族），**没有任何"同一体系的温度序列"**。而温度依赖在引擎里是
一条独立链路：`pka.dH`（只有 2/89 条直接给）或 **Hess 定律从 `thermo.json`
派生 dH** ⟹ `core._vant()`：`logK(T) = logK298 + dH/(R·ln10)·(1/298.15 − 1/T)`。

**判据（可独立算）**：
  · 一元弱酸**半中和缓冲**（等摩尔 HA/A⁻）⟹ pH = pKa(T)，
    而 pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)；
  · 等当量强酸强碱 ⟹ pH = pKw(T)/2（pKw 走 `core.pKw_of`）。

所以每条温度曲线都能**从库内常数直接算出来**，不依赖引擎输出。

用法： python tools/temp_series.py [--write]
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

R_LN10 = 8.31446261815324e-3 * math.log(10.0)          # kJ/mol/K × ln10
TS = (273.15, 283.15, 298.15, 313.15, 323.15, 333.15,
      343.15, 353.15, 363.15, 373.15)


def vant(dH, tk):
    return dH / R_LN10 * (1.0 / 298.15 - 1.0 / tk)


def main():
    from chemkit.data import load_tables
    from chemkit.engine import judge
    from chemkit import core
    T = load_tables()

    # 按 |dH| 排序，挑温度效应最大的弱酸（dH 大 ⟹ pKa(T) 变化大）
    acids = [e for e in T.pka
             if e.get("n") == 1 and e.get("dH") is not None]
    acids.sort(key=lambda e: -abs(e["dH"]))
    print(f"一元弱酸中带 dH 的：{len(acids)} 条；"
          f"|dH| 最大的 6 条：")
    for e in acids[:6]:
        print(f"   {e['acid']:>12s} / {e['base']:<12s} pKa298={e['pka']:<7} "
              f"dH={e['dH']:+.2f} kJ/mol")
    print()

    # 用"游离酸 + 半当量 NaOH"造缓冲（等摩尔 HA/A⁻ ⟹ pH = pKa(T)）
    CASES = [(e["acid"], e["base"], e["pka"], e["dH"]) for e in acids[:3]]
    CASES.append(("CH_3COOH", "CH_3COO^-", 4.76, -0.2))

    print("=== ① 半中和缓冲：pH 应 = pKa(T)（由库内 pKa298 + dH 独立算出）===")
    worst = []
    for acid, base, pk298, dH in CASES:
        print(f"\n  --- {acid}  pKa298={pk298}  dH={dH:+.2f} ---")
        print("     %-9s %-11s %-11s %-9s" % ("T/K", "应=pKa(T)",
                                              "引擎 pH", "Δ"))
        c = 0.01
        for tk in TS:
            expect = pk298 - vant(dH, tk)
            try:
                r = judge([{"name": acid, "mol": c},
                           {"name": "NaOH", "mol": c / 2}],
                          {"V_L": 1.0, "T_K": tk}, T)
                got = r.get("final_pH")
            except Exception as exc:                        # noqa: BLE001
                got = None
                print(f"     {tk:<9.2f} ERR {type(exc).__name__}")
                continue
            d = None if got is None else got - expect
            print("     %-9.2f %-11.4f %-11.4f %s"
                  % (tk, expect, got,
                     ("%+.4f" % d) if d is not None else "-"))
            if d is not None:
                worst.append((abs(d), acid, tk, expect, got))

    print("\n=== ② 等当量强酸强碱：pH 应 = pKw(T)/2 ===")
    print("     %-9s %-11s %-11s %-9s" % ("T/K", "pKw/2", "引擎 pH", "Δ"))
    for tk in TS:
        pkw = core.pKw_of(tk)
        r = judge([{"name": "HCl", "mol": 1e-3},
                   {"name": "NaOH", "mol": 1e-3}],
                  {"V_L": 1.0, "T_K": tk}, T)
        got = r.get("final_pH")
        d = got - pkw / 2
        print("     %-9.2f %-11.4f %-11.4f %+.4f" % (tk, pkw / 2, got, d))
        worst.append((abs(d), "HCl+NaOH", tk, pkw / 2, got))

    worst.sort(reverse=True)
    print("\n=== 汇总：|Δ| 最大的 6 条（阈值 0.05）===")
    bad = [w for w in worst if w[0] > 0.05]
    for a, name, tk, exp, got in worst[:6]:
        print(f"   |Δ|={a:.4f}  {name:<14s} T={tk:<8.2f} 应={exp:.4f} "
              f"引擎={got:.4f}")
    print(f"\n超过 0.05 的：{len(bad)} 条 / 共 {len(worst)} 条")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "temp_series.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"abs_d": a, "case": n, "T": tk, "expect": e,
                        "engine": g} for a, n, tk, e, g in worst],
                      f, ensure_ascii=False, indent=1)
        print(f"[写入] {os.path.relpath(p, ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
