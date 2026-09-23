# -*- coding: utf-8 -*-
"""第 255 轮 · 用**物理投料**重做弱酸的 pH 通路审计（修正第 253 轮的输入缺陷）。

**第 253 轮错在哪**（第 254 轮查明）：我手工构造
`{HA:0.005, A⁻:0.005, Na⁺:0.005}` 喂给 `charge_pH`，但它是**元素总量账**——
同族两个物种会被**相加**成 `M=0.01`。于是"双重记账"，强酸端结论整片作废。

**本轮做法（物理投料，不再手造账本）**：
对库内每个一元弱酸 `HA`，投料 `[HA: 0.01 mol, NaOH: 0.005 mol]`（半中和），
让**引擎自己**走出 ledger，然后：
  ① 从探针拿 `pH_engine`（`final_pH` 用的那条通路）与 ledger；
  ② 直接调 `estimate_pH` / `exact_proton_pH` 得两条通路的值；
  ③ 与独立解 `pKa(T)` 并排比。

**判据**：半中和时 `pH = pKa(T)`（前提：酸为**一元**、**只发生一级解离**、
且**无其它酸碱储备**）。第 254 轮学到：该前提**不是我能指定的**，
所以本轮**只把三条数并排打出来**，判定留给"哪条通路与 pKa 更近"这个经验事实，
不宣称"必然"。

用法： python tools/weakacid_paths.py [--write]
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
    from chemkit.engine import judge, estimate_pH, exact_proton_pH
    T = load_tables()

    rows = []
    for e in T.pka:
        if e.get("n") != 1:
            continue
        a, pk = e["acid"], e["pka"]
        # 只做"酸+半量 NaOH"这条路：要求酸本身能作为投料被引擎接受
        try:
            pr = {}
            r = judge([{"name": a, "mol": 0.01},
                       {"name": "NaOH", "mol": 0.005}],
                      {"V_L": 1.0, "T_K": 298.15}, T, _probe=pr)
        except Exception:                                   # noqa: BLE001
            continue
        led = pr.get("ledger")
        if not led:
            continue
        ep = pr.get("pH")
        est = estimate_pH(led, pr.get("H_excess", 0.0), 1.0, T, 298.15)
        ex = exact_proton_pH(led, pr.get("H_excess", 0.0), 1.0, T, 298.15)
        rows.append((a, pk, ep, est, ex))

    print(f"可用'酸 + 半量 NaOH'跑通的 一元弱酸：{len(rows)} 条\n")
    print("  %-14s %-8s %-9s %-9s %-9s %-8s %s"
          % ("acid", "pKa298", "引擎pH", "estimate", "exact", "∩", "谁更近 pKa"))
    n_est, n_ex, n_eng = 0, 0, 0
    tab = []
    for a, pk, ep, est, ex in rows:
        # 每条通路与 pKa 的距离（exact 可能是 None）
        de = abs(est - pk) if est is not None else None
        dx = abs(ex - pk) if ex is not None else None
        dg = abs(ep - pk) if ep is not None else None
        who = "-"
        if de is not None and dx is not None:
            if de < dx:
                who = "estimate"
                n_est += 1
            elif dx < de:
                who = "exact"
                n_ex += 1
            else:
                who = "tie"
        elif dx is None and de is not None:
            who = "exact=None"
        if dg is not None and de is not None and dx is not None:
            near = min(de, dx)
            # 引擎实际跟谁走
            if dg <= near + 0.02:
                n_eng += 1
        tab.append((a, pk, ep, est, ex, who, dg, de, dx))
    for a, pk, ep, est, ex, who, dg, de, dx in sorted(tab, key=lambda r: r[1]):
        print("  %-14s %-8.3f %-9s %-9s %-9s %-8s %s"
              % (a, pk,
                 ("%.3f" % ep) if ep is not None else "-",
                 ("%.3f" % est) if est is not None else "-",
                 ("%.3f" % ex) if ex is not None else "None",
                 who, ""))
    print(f"\n  estimate 更近 pKa 的：{n_est}；exact 更近的：{n_ex}")
    got = [r for r in tab if r[5] in ("estimate", "exact")]
    if got:
        big = [r for r in got if min(r[7] or 9, r[8] or 9) > 0.05]
        print(f"  两条通路都不近（min|Δ|>0.05）的：{len(big)}")
        for r in big:
            print(f"     {r[0]:<14} pKa={r[1]:.3f} est={r[3]} "
                  f"exact={r[4]}")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "weakacid_paths.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"acid": r[0], "pKa": r[1], "engine": r[2],
                        "estimate": r[3], "exact": r[4], "closer": r[5]}
                       for r in tab], f, ensure_ascii=False, indent=1)
        print(f"\n[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
