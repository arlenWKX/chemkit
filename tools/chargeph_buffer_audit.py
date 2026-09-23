# -*- coding: utf-8 -*-
"""第 253 轮 · `charge_pH` 在 1:1 缓冲上的偏差：**全孔径普查**。

**背景**：第 252 轮确认 `exact_proton_pH`（内部走 `charge_pH`）与
`estimate_pH` 在 253 个过闸状态里有 44 例 |Δ|>0.05。本轮把它收到一个
**可复现的最小窗口**：对库内**每个一元酸**构造"等摩尔 HA/A⁻ + Na⁺"的
理想缓冲账本，直接问 `charge_pH`。

**判据（数学上必然）**：1:1 缓冲时 `pH = pKa + log([A⁻]/[HA]) = pKa`。
所以偏离量就是**该函数的纯误差**，不需要引擎其它部分参与。

**为何不用引擎跑用例**：直接构造账本可以扫**全部**一元酸（含没有对应用例的），
不受"哪些体系有 test case"限制。

用法： python tools/chargeph_buffer_audit.py [--write]
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


def main():
    from chemkit.data import load_tables
    from chemkit.acidbase import charge_pH
    from chemkit.core import charge_of
    T = load_tables()

    rows = []
    for e in T.pka:
        if e.get("n") != 1:
            continue
        a, b, pk = e["acid"], e["base"], e["pka"]
        try:
            if charge_of(a) != 0 or charge_of(b) != -1:
                continue
        except Exception:                                   # noqa: BLE001
            continue
        led = {a: 0.005, b: 0.005, "Na^+": 0.005}
        try:
            cp = charge_pH(led, 1.0, T, 298.15, fast=True)
        except Exception:                                   # noqa: BLE001
            continue
        if cp is None:
            continue
        dH = e.get("dH")
        rows.append((abs(cp - pk), a, pk, cp, cp - pk, dH))

    rows.sort(reverse=True)
    print(f"可构造 1:1 缓冲的一元酸：{len(rows)} 条（298.15 K）")
    print("判据：pH 必然等于 pKa；Δ 即 charge_pH 的纯误差\n")
    print("  %-14s %-9s %-11s %-9s %s"
          % ("acid", "pKa", "charge_pH", "Δ", "dH"))
    for ad, a, pk, cp, d, dH in rows[:16]:
        print("  %-14s %-9.3f %-11.4f %+9.4f %s"
              % (a, pk, cp, d, ("%+.1f" % dH) if dH is not None else "-"))
    for thr in (0.05, 0.1, 0.3):
        k = sum(1 for ad, *_r in rows if ad > thr)
        print(f"  |Δ| > {thr}: {k} / {len(rows)}")
    # 温度依赖（以 H2O2 为例）
    e = [x for x in T.pka if x["acid"] == "H_2O_2"]
    if e and e[0].get("dH") is not None:
        dH, pk0 = e[0]["dH"], e[0]["pka"]
        print("\n  === 同一缓冲（H2O2）扫温度：误差随 T 变宽 ===")
        for tk in (273.15, 298.15, 323.15, 348.15, 373.15):
            pk = pk0 - dH / R_LN10 * (1 / 298.15 - 1 / tk)
            cp = charge_pH({"H_2O_2": 0.005, "HO_2^-": 0.005, "Na^+": 0.005},
                           1.0, T, tk, fast=True)
            print("     T=%-8.2f pKa(T)=%-8.4f charge_pH=%-9.4f Δ=%+.4f"
                  % (tk, pk, cp, cp - pk))
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "chargeph_buffer_audit.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"acid": a, "pKa": pk, "charge_pH": cp, "d": d,
                        "dH": dH} for _ad, a, pk, cp, d, dH in rows], f,
                      ensure_ascii=False, indent=1)
        print(f"\n[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
