"""口径边界普查：changed / reacted / 净方程 三者的分布与"净方程类型"分类。

用途（用户提问：水解算不算 reacted？有意义的平衡过程该不该给方程？）：
把全库 1173 例按 (changed, reacted, 净方程是否存在) 分格，再把"有方程但
reacted=False"的 120 例按**方程里有没有固相、固相是不是投料**细分——
这个细分正是"痕量新相"这类边界判据的依据。

用法：python tools/eq_semantics.py [-n 15]
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from chemkit.data import load_tables                     # noqa: E402
from chemkit.engine import judge                          # noqa: E402
from chemkit.system import Reaction                       # noqa: E402
from chemkit.testsuit import load_cases                   # noqa: E402


def main() -> None:
    n = 15
    if "-n" in sys.argv:
        n = int(sys.argv[sys.argv.index("-n") + 1])
    T = load_tables()
    solids = frozenset(T.solids)
    grid: dict = {}
    cats: dict = {}
    for c in load_cases(None):
        r = judge([{"name": nm, "mol": m} for nm, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
        rxn = Reaction(r)
        ch, rx, eq = bool(r["changed"]), bool(r.get("reacted")), rxn.net_equation
        grid[(ch, rx, bool(eq))] = grid.get((ch, rx, bool(eq)), 0) + 1
        if eq and not rx:
            feed = {e["name"] for e in r.get("initial", [])}
            prod = {e["name"] for e in r.get("production", [])}
            cons = {e["name"] for e in r.get("consumption", [])}
            new_solid_p = sorted(s for s in prod if s in solids and s not in feed)
            feed_solid = sorted(s for s in cons | prod if s in solids and s in feed)
            if new_solid_p:
                k = "生成固相（投料里没有）"
            elif feed_solid:
                k = "投料自带固相的溶解/析出（溶解度表达）"
            else:
                k = "纯离子/分子平衡（解离、水解、配位、氧化还原）"
            cats.setdefault(k, []).append((c["name"], eq, ch))
    print("== (changed, reacted, 有净方程) 分布 ==")
    for k in sorted(grid, key=lambda k: -grid[k]):
        print(f"   changed={int(k[0])} reacted={int(k[1])} 有方程={int(k[2])}"
              f"  → {grid[k]:4d} 例")
    print("\n== 有方程但 reacted=False 的分类 ==")
    for k, v in sorted(cats.items(), key=lambda t: -len(t[1])):
        print(f"\n-- {k}：{len(v)} 例（其中 changed=True 的 "
              f"{sum(1 for _, _, ch in v if ch)} 例）")
        for nm, eq, ch in v[:n]:
            print(f"     changed={int(ch)} {nm[:30]:32s} {eq[:58]}")


if __name__ == "__main__":
    main()
