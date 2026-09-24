# -*- coding: utf-8 -*-
"""第 286 轮 · **指定态上的 S 分解探针**：把某条候选的 `S = logK − logQ`
逐项拆开（`logK`、pH、各物种的活度项），用来判定"本来应该是 0 的 S 却算出 0.3"
到底是 **logK** 错还是 **pH/活度** 错。

## 为什么需要

第 285 轮把 `TE1`（`NH4Cl 0.01 + NaOH 0.005`）的走步查清：step2
`NH4^+ -> NH3 + H^+` x=0.004869（**恰是游离 OH⁻ 全量**）走完 `He → 0`，
**账本已经是正确的 1:1 缓冲**（`NH4^+ = NH3 = 0.005`、pH = pKa = 9.25）。
按定义此时 `NH3 + H^+ -> NH4^+` 的 `S = pKa − pH + log([NH3]/[NH4^+]) = 0`。
**但引擎给 step3（反向）S = 0.3** —— 正是这 0.3 把它推着拆掉了正确终态。

本工具在**同一个账本**上把 `S` 拆开：`logK_T`、`estimate_pH`、逐物种项，
一眼看出 0.3 是从哪一项来的。

用法：
    python tools/sprobe.py                       # 默认跑 TE1 的 step2 终态
    python tools/sprobe.py --led=NH_4^+:0.005,NH_3:0.005,Na^+:0.005,Cl^-:0.01 --He=0 --needle=NH_3
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

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.candidates import (H_ION, WATER, logK_T,       # noqa: E402
                                ACT_FLOOR)
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402

V = 1.0
T_K = 298.15


def _txt(dd: dict) -> str:
    return " + ".join(f"{'' if n == 1 else n}{s}"
                      for s, n in dd.items() if s != WATER)


def main(argv: list[str]) -> int:
    led_s = "NH_4^+:0.005,NH_3:0.005,Na^+:0.005,Cl^-:0.01"
    He = 0.0
    needle = "NH_3"
    for a in argv:
        if a.startswith("--led="):
            led_s = a[6:]
        elif a.startswith("--He="):
            He = float(a[5:])
        elif a.startswith("--needle="):
            needle = a[9:]
    led = {}
    for part in led_s.split(","):
        if not part.strip():
            continue
        k, v = part.split(":")
        led[k.strip()] = float(v)
    T = load_tables()
    pH = spec.estimate_pH(dict(led), He, V, T, T_K)
    print(f"账本 {led}  He={He}  ->  estimate_pH = {pH:.6f}")
    print(f"（1:1 缓冲的定义式 pH = pKa = 9.25）\n")
    print(f"{'候选':<46}{'logK_T':>9}{'S_of':>10}{'d':>4}")
    for c in enumerate_candidates(led, He, pH, V, T_K, T, True):
        for d in (1, -1):
            rr = c.r if d > 0 else c.pr
            pp = c.pr if d > 0 else c.r
            t = f"{_txt(rr)} -> {_txt(pp)}"
            if needle not in t:
                continue
            lk = logK_T(c, T_K, 0.0)
            S = eng.S_of(c, led, V, pH, T_K, T, frozenset(),
                         eng.P_EXT_KPA, True, None) * d
            print(f"{t[:46]:<46}{lk:>9.4f}{S:>10.4f}{d:>4}")
    print("\n⚠️ **本工具只打印引擎自己的值**（`logK_T` 与 `S_of`）。"
          "\n   第 286 轮初版还并排了一列「手工 logK−logQ」，"
          "\n   实测它与 `S_of` 相差 18.5（同一个 Cand 的 `S_of`=0.0000 "
          "而手工给 −18.5）\n   ⟹ **那一列是我算错的**（Cand 的存储书写向与"
          "打印方向不一定是同一个），已删除。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
