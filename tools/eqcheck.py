"""期望方程式守恒核验（元素 + 电荷），逐条给出偏差。

用法：python tools/eqcheck.py "式子1" "式子2" ...
      python tools/eqcheck.py --cases N30 NR95        # 从 tests.json 取 eq
"""
from __future__ import annotations

import re
import sys
from fractions import Fraction as F

sys.path.insert(0, ".")

from chemkit.core import charge_of, elements_of   # noqa: E402


def _split_terms(side: str) -> list[str]:
    """按**带空格的**顶层 ` + ` 切项。

    不能按裸 `+` 切：`Fe^{3+}`、`[Fe(SCN)_2]^+` 的电荷号也是 `+`，且
    `]^+` 的 `+` 还在括号**外面**（深度判定救不了）。引擎的方程格式化
    一律用 `" + ".join(...)`，故以 ` + `（两侧空白）为分隔符是稳的。
    """
    return [t.strip() for t in re.split(r"\s+\+\s+", side.strip()) if t.strip()]


def parse(side: str) -> dict:
    """'99SCN^- + 52Fe^{3+}' -> {物种: 系数(有理)}"""
    out: dict = {}
    for term in _split_terms(side):
        m = re.match(r"^(\d+)\s*(.*)$", term)
        if m:
            k, sp = F(int(m.group(1))), m.group(2).strip()
        else:
            k, sp = F(1), term
        out[sp] = out.get(sp, F(0)) + k
    return out


def check(eq: str) -> None:
    if "->" in eq:
        lhs, rhs = eq.split("->")
    else:
        lhs, rhs = eq.split("=")
    L, R = parse(lhs), parse(rhs)
    bad = []
    els = set()
    for s in list(L) + list(R):
        els |= set(elements_of(s))
    for el in sorted(els):
        d = (sum(k * elements_of(s).get(el, 0) for s, k in R.items())
             - sum(k * elements_of(s).get(el, 0) for s, k in L.items()))
        if d:
            bad.append(f"{el}{float(d):+g}")
    dq = (sum(k * charge_of(s) for s, k in R.items())
          - sum(k * charge_of(s) for s, k in L.items()))
    if dq:
        bad.append(f"q{float(dq):+g}")
    print(f"  {'守恒 ✓' if not bad else '不守恒 ✗ ' + ', '.join(bad)}   {eq[:88]}")


def main(argv: list[str]) -> None:
    if argv and argv[0] == "--cases":
        from chemkit.testsuit import load_cases
        want = tuple(argv[1:])
        for c in load_cases(None):
            if c["name"].startswith(want) and c.get("eq"):
                print(f"## {c['name']}")
                check(c["eq"])
        return
    for e in argv:
        check(e)


if __name__ == "__main__":
    main(sys.argv[1:])
