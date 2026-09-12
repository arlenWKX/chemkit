"""`charge_pH` 的 Newton 快路径 vs 二分：全量账本上的逐例等价性与加速比。

用途（§7 O）：精确质子条件要进引擎热路径，前提是解算器够快且与二分
**同根**。本工具在真实退出账本上逐例对拍（最大 |ΔpH|）+ 计时。

用法：python tools/phcmp.py [-n 10]
"""
from __future__ import annotations

import sys
import time

sys.path.insert(0, ".")

from chemkit.acidbase import charge_pH                     # noqa: E402
from chemkit.data import load_tables                       # noqa: E402
from chemkit.engine import judge                           # noqa: E402
from chemkit.testsuit import load_cases                    # noqa: E402


def main() -> int:
    n = 10
    if "-n" in sys.argv:
        n = int(sys.argv[sys.argv.index("-n") + 1])
    T = load_tables()
    rows = []
    t_bi = t_nw = 0.0
    both = 0
    for c in load_cases(None):
        p: dict = {}
        cond = c.get("cond") or {"V_L": 1.0}
        judge([{"name": x, "mol": m} for x, m in c["subs"]], cond, T, _probe=p)
        if "ledger" not in p:
            continue
        led = p["ledger"]
        V = float(cond.get("V_L", 1.0) or 1.0)
        c_H = float(cond.get("c_H", 0.0) or 0.0)
        t0 = time.perf_counter()
        a = charge_pH(led, V, T, 298.15, c_H=c_H)
        t1 = time.perf_counter()
        b = charge_pH(led, V, T, 298.15, c_H=c_H, fast=True)
        t2 = time.perf_counter()
        t_bi += t1 - t0
        t_nw += t2 - t1
        if a is not None and b is not None:
            both += 1
            rows.append((abs(a - b), c["name"], a, b))
        elif (a is None) != (b is None):
            rows.append((99.0, c["name"] + " [None 不一致]", a, b))
    rows.sort(reverse=True)
    print(f"可对照 {both} 例；二分 {t_bi * 1000:.0f} ms，"
          f"Newton {t_nw * 1000:.0f} ms（×{t_bi / max(t_nw, 1e-9):.1f}）")
    bad = [r for r in rows if r[0] > 1e-9]
    print(f"|ΔpH| > 1e-9 的例：{len(bad)}")
    for d, nm, a, b in rows[:n]:
        print(f"   Δ{d:.3e}  二分 {a}  Newton {b}  {nm[:46]}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
