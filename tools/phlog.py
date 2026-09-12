"""pH 估计轨迹：逐次记录 solve_extent 二分内 `estimate_pH` 的输入/输出。

用途：定位"幻影碱/幻影酸"类走步偏移——某个 x 上 pH 估计偏离真实值，
S = logK − logQ 的符号随之翻转，走步就把反应推到不存在的产物上
（F35 NaH2PO4：H2PO4- + H+ → H3PO4 跑掉 0.5 mol，He 变 −0.5）。

用法：
    python tools/phlog.py F35           # 首尾各 12 次调用
    python tools/phlog.py F35 -n 40     # 首尾各 40 次
    python tools/phlog.py F35 --all     # 全部
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from chemkit.data import load_tables                 # noqa: E402
from chemkit import engine                            # noqa: E402
from chemkit.testsuit import load_cases               # noqa: E402

_ORIG = engine.estimate_pH
_CALLS: list = []


def _wrap(ledger, H_excess, V, T, T_K, cache=None, touch=None):
    out = _ORIG(ledger, H_excess, V, T, T_K, cache, touch)
    _CALLS.append((dict(ledger), H_excess, out, V, T_K, cache, touch))
    return out


engine.estimate_pH = _wrap


def _fmt(led, He, pH):
    parts = [f"{s}={m:.6g}" for s, m in sorted(led.items(), key=lambda kv: -kv[1])
             if s != "H_2O" and m > 1e-9]
    return f"He={He:+10.6f} → pH {pH:7.3f} | " + "  ".join(parts[:9])


def show(prefixes: list[str], n: int, all_: bool) -> None:
    T = load_tables()
    for c in load_cases(None):
        if not c["name"].startswith(tuple(prefixes)):
            continue
        _CALLS.clear()
        engine.judge([{"name": x, "mol": m} for x, m in c["subs"]],
                     c.get("cond") or {"V_L": 1.0}, T)
        print(f"\n{'='*78}\n{c['name']}  调用 {len(_CALLS)} 次")
        print(f"  ph 期望 {c.get('ph')}   注: {(c.get('note') or '')[:70]}")
        idx = range(len(_CALLS)) if all_ else list(range(min(n, len(_CALLS)))) \
            + list(range(max(n, len(_CALLS) - n), len(_CALLS)))
        last = -1
        for i in idx:
            if i == last:
                continue
            if i == n and not all_:
                print("   …")
            last = i
            led, He, pH = _CALLS[i][0], _CALLS[i][1], _CALLS[i][2]
            V, T_K = _CALLS[i][3], _CALLS[i][4]
            print(f"  #{i:<4d} {_fmt(led, He, pH)}   [V={V} T={T_K} "
                  f"keys={len(led)}]")


if __name__ == "__main__":
    argv = sys.argv[1:]
    nn = 12
    if "-n" in argv:
        i = argv.index("-n")
        nn = int(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    args = [x for x in argv if not x.startswith("-")]
    show(args or ["F35"], nn, "--all" in sys.argv)
