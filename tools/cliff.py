"""pH 机器的 He 连续性探针：扫 He，看 estimate_state 的 pH 有无跳变。

J06「酸侧悬崖」的可复现证据 + 修复后的回归闸门。

用法：python tools/cliff.py [用例名前缀 ...]
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from chemkit.data import load_tables                    # noqa: E402
from chemkit.engine import judge                        # noqa: E402
from chemkit.speciation import estimate_pH              # noqa: E402
from chemkit.testsuit import load_cases                 # noqa: E402

GRID = [-1e-2, -3e-3, -1.5e-3, -1.1e-3, -1e-3, -0.999e-3, -0.99e-3, -7e-4,
        -3e-4, -1e-4, -1e-5, -1e-6, 0.0,
        1e-6, 1e-5, 1e-4, 3e-4, 7e-4, 0.99e-3, 0.999e-3, 1e-3, 1.1e-3,
        1.5e-3, 3e-3, 1e-2]


def scan(prefixes: list[str]) -> None:
    T = load_tables()
    worst = []
    for c in load_cases(None):
        if not c["name"].startswith(tuple(prefixes)):
            continue
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        probe: dict = {}
        judge(subs, cond, T, _probe=probe)
        if "ledger" not in probe:
            continue
        V = float(cond.get("V_L", 1.0) or 1.0)
        T_K = float(cond.get("T_K", 298.15) or 298.15)
        led = probe["ledger"]
        print(f"\n== {c['name']}  V={V} T_K={T_K}  "
              f"（退出账本 He={probe['H_excess']}）")
        prev = None
        mx = 0.0
        for He in GRID:
            ph = estimate_pH(led, He * V, V, T, T_K)
            d = "" if prev is None else f"  Δ={ph - prev:+.4f}"
            if prev is not None:
                mx = max(mx, abs(ph - prev))
            print(f"   He={He:+10.3e}  pH={ph:8.4f}{d}")
            prev = ph
        worst.append((mx, c["name"]))
    print("\n== 相邻栅格最大跳变排行 ==")
    for mx, nm in sorted(worst, reverse=True)[:15]:
        print(f"   Δmax={mx:8.3f}  {nm}")


if __name__ == "__main__":
    a = sys.argv[1:] or ["J06"]
    scan(a)
