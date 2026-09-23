# -*- coding: utf-8 -*-
"""第 261 轮 · 半中和的**精确解**审计（不依赖任何近似式）。

**判据（第 260 轮推得，闭式可解）**：一元酸 `HA` 投料 `c`、加 `c/2` 的 NaOH，
纯水、无其它酸碱储备时，由**电荷守恒 + 质量守恒**：

    [A⁻] = c/2 + [H⁺] − [OH⁻]
    [HA] = c/2 − [H⁺] + [OH⁻]
    Ka   = [A⁻]·[H⁺] / [HA]

（**注意**：`[A⁻]` 不是 Henderson 近似里的 `c/2` —— 当 `[H⁺]` 与 `c/2` 同量级时
两者差别巨大，这正是第 253–259 轮误判的根源。）

把 `h = 10⁻ᵖᴴ`、`oh = 10^(pH−pKw)` 代进去，对 h 是一元二次式，可直接求根，
**不需要任何近似**，也不需要迭代。

**本脚本做两件事**：
  1. 对库内**每个**一元酸做"精确解 vs 引擎"对拍（用物理投料，不手造账本）；
  2. 报告｜Δ｜分布，指出是否还有真偏差。

用法： python tools/halfneut_exact_audit.py [--write]
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

R_LN10 = 8.31446261815324e-3 * math.log(10.0)


def exact_pH(Ka: float, c: float, pKw: float = 14.0):
    """半中和精确解：解 (c/2 + h − oh)·h = Ka·(c/2 − h + oh)。"""
    Kw = 10.0 ** (-pKw)

    def f(pH):
        h = 10.0 ** (-pH)
        oh = Kw / h
        return (c / 2.0 + h - oh) * h - Ka * (c / 2.0 - h + oh)

    # f 随 pH 单调下降（h 降、oh 升）；找唯一过零
    lo, hi = 0.0, pKw
    flo = f(lo)
    if flo < 0.0:
        return lo
    if f(hi) > 0.0:
        return None
    a, b, fa = lo, hi, flo
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm < 0:
            b = m
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def main():
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()

    rows = []
    for e in T.pka:
        if e.get("n") != 1:
            continue
        acid, pk = e["acid"], e["pka"]
        dH = e.get("dH")
        for c in (0.01,):
            try:
                pr = {}
                judge([{"name": acid, "mol": c},
                       {"name": "NaOH", "mol": c / 2}],
                      {"V_L": 1.0, "T_K": 298.15}, T, _probe=pr)
            except Exception:                               # noqa: BLE001
                continue
            g = pr.get("pH")
            if g is None:
                continue
            ex = exact_pH(10.0 ** -pk, c)
            if ex is None:
                continue
            rows.append((abs(g - ex), acid, pk, c, ex, g, dH))

    rows.sort(reverse=True)
    print(f"可用'酸 + 半量 NaOH'跑通、且精确解存在的一元酸：{len(rows)} 条"
          f"（c = 0.01 M, 298.15 K）\n")
    print("  精确解：由电荷+质量守恒直接解，**不用 Henderson 近似**\n")
    print("  %-14s %-8s %-10s %-10s %s"
          % ("acid", "pKa298", "精确解", "引擎", "Δ"))
    for ad, acid, pk, c, ex, g, dH in rows[:20]:
        print("  %-14s %-8.3f %-10.4f %-10.4f %+.4f"
              % (acid, pk, ex, g, g - ex))
    for thr in (0.01, 0.05, 0.1, 0.3):
        k = sum(1 for ad, *_r in rows if ad > thr)
        print(f"  |Δ| > {thr}: {k} / {len(rows)}")
    bad = [r for r in rows if r[0] > 0.05]
    print(f"\n超过 0.05 的：{len(bad)} 条")
    for ad, acid, pk, c, ex, g, dH in bad[:15]:
        print(f"   {acid:<14} pKa={pk:<7.3f} 精确={ex:.4f} 引擎={g:.4f} "
              f"Δ={g - ex:+.4f} dH={dH}")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "halfneut_exact_audit.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"acid": r[1], "pKa": r[2], "c": r[3], "exact": r[4],
                        "engine": r[5], "d": r[5] - r[4], "dH": r[6]}
                       for r in rows], f, ensure_ascii=False, indent=1)
        print(f"\n[写入] {os.path.relpath(p, ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
