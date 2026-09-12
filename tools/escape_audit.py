"""逸出气体质量守恒审计：初态 − 终态 − 逸出 = 0？

gas escape 路径（`_sweep_gases` + `escaped` 账）如果丢量，会在**所有**
涉及气体逸出的用例上留下按元素计的缺口，而净方程装配器拿到的是一个
本身不平的净差——表现为"非整系数 + 大系数"（E42 实测 C 缺 0.0242 mol，
于是美化器永远找不到干净整数式）。

本工具对全库逐例核验：每个**逸出气体所含元素**的
`初始 − 终态 − 逸出`，报出缺口排行。

用法：python tools/escape_audit.py [-n 20]
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from chemkit.core import elements_of          # noqa: E402
from chemkit.data import load_tables          # noqa: E402
from chemkit.engine import judge              # noqa: E402
from chemkit.testsuit import load_cases       # noqa: E402


def main(n: int) -> None:
    T = load_tables()
    rows = []
    n_esc = 0
    for c in load_cases(None):
        r = judge([{"name": a, "mol": b} for a, b in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
        esc = r.get("escaped") or []
        if not esc:
            continue
        n_esc += 1
        els = set()
        for e in esc:
            els |= set(elements_of(e["name"]))
        # H/O 不走元素账：它们由水（溶剂，不在 initial/final 列表里）与
        # 质子账本 H_excess 承载，纳入只会得到"缺口 −12 H"这类假象。
        els -= {"H", "O"}
        if not els:
            continue
        ini = {el: 0.0 for el in els}
        fin = {el: 0.0 for el in els}
        escv = {el: 0.0 for el in els}
        for e in r["initial"]:
            for el, k in elements_of(e["name"]).items():
                if el in ini:
                    ini[el] += k * e["mol"]
        for e in r["final"]:
            for el, k in elements_of(e["name"]).items():
                if el in fin:
                    fin[el] += k * e["mol"]
        for e in esc:
            for el, k in elements_of(e["name"]).items():
                if el in escv:
                    escv[el] += k * e["mol"]
        for el in els:
            gap = ini[el] - fin[el] - escv[el]
            if abs(gap) > 1e-6:
                rows.append((abs(gap), gap, el, c["name"],
                             {e["name"]: round(e["mol"], 6) for e in esc},
                             round(escv[el], 6)))
    rows.sort(reverse=True)
    print(f"涉及气体逸出的用例 {n_esc} 例；**元素缺口 {len(rows)} 处**")
    for a, gap, el, nm, escd, tot in rows:
        print(f"  缺 {gap:+10.6f} {el}  (逸出 {el} 共 {tot:g})  {nm[:40]}")
        print(f"       escaped={escd}")
    if not rows:
        print("  逸出路径元素守恒 ✓")


if __name__ == "__main__":
    nn = 20
    a = sys.argv[1:]
    if "-n" in a:
        nn = int(a[a.index("-n") + 1])
    main(nn)
