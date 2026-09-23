# -*- coding: utf-8 -*-
"""第 265 轮 · 取某用例**末态**上 `estimate_pH` 的分支归属 + 各条 pH 通路并排。

为什么需要它：`tools/xscan_ch.py` 里我用 `charge_pH(led, ...)`（**不带 `pinned`**）
当"精确解"，得到 pH 8.984 而 `estimate_pH` 给 6.168。但账本里
`In(OH)₃(s) = 0.997 mol` **确实在场** ⟹ 该阳离子应由 **Ksp 钉住**
（`estimate_state` 的 `pinned` 路径），不带 pinned 的解**不是**同一口径。
**判据必须同口径**（第 253/263 轮教训）。

本工具直接问引擎：① 赢下 pH 的是哪一支（`PH_TAGS`）；
② 参与竞争的来源（`PH_SRC`）；③ 退化区闸 `_deg` 是否成立。

用法：python tools/ph_branch.py I31
"""
from __future__ import annotations

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

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.acidbase import charge_pH                      # noqa: E402
from chemkit.candidates import H_ION, WATER                 # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def main(argv: list[str]) -> int:
    pre = argv[0] if argv else "I31"
    T = load_tables()
    hit = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
           if n.startswith(pre)]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    T_K = float(cond.get("T_K", 298.15))
    pr = {}
    r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T,
                  _probe=pr)
    led = dict(pr.get("ledger") or {})
    He = pr.get("H_excess") or 0.0
    print(f"用例 {c['name']}  pH={r.get('final_pH')} "
          f"solver={pr.get('pH_solver')} He={He:.8g}")
    led_s = {k: v for k, v in sorted(led.items()) if v and k != WATER}
    print("账本: " + "  ".join(f"{k}={v:.8g}" for k, v in led_s.items()))

    spec.PH_TAGS = []
    spec.PH_SRC = []
    try:
        ph, vled, He_res = spec.estimate_state(led, He, V, T, T_K)
        tags, src = list(spec.PH_TAGS), list(spec.PH_SRC)
    finally:
        spec.PH_TAGS = None
        spec.PH_SRC = None
    print(f"\nestimate_state -> pH = {ph:.6f}   tags = {tags}")
    print(f"  He_res = {He_res:.8g}")
    print("  来源（前 12，按贡献降序）:")
    for s, k, v in sorted(src, key=lambda t: -abs(t[2]))[:12]:
        print(f"     {s:<24} [{k}] = {v:.8g}")

    # 退化区闸（与 estimate_state 同式）：h_c/o_c 是否能从来源里复原
    print("\n-- 各 pH 通路并排（**同口径提醒**：不带 pinned 的那条不是同一模型）--")
    print(f"  estimate_pH                 = {ph:.6f}")
    for kw in ({}, {"fast": True}):
        try:
            x = charge_pH(led, V, T, T_K, **kw)
        except Exception as exc:                            # noqa: BLE001
            x = f"EXC {exc!r}"
        print(f"  charge_pH(无 pinned{', fast' if kw else ''})   = {x}")
    # 纯水/自由碱口径参照：He 直接给的下限
    pKw = spec.pKw_of(T_K)
    if He < 0:
        print(f"  「自由碱」纯水口径 pKw+log10(-He) = "
              f"{pKw + __import__('math').log10(-He):.6f}")
    elif He > 0:
        print(f"  「自由酸」纯水口径 -log10(He)     = "
              f"{-__import__('math').log10(He):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
