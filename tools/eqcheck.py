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
        m = re.match(r"^([0-9]*\.?[0-9]+)\s*(.*)$", term)
        if m and m.group(1):
            k, sp = F(m.group(1)), m.group(2).strip()
        else:
            k, sp = F(1), term
        if not sp:
            continue
        out[sp] = out.get(sp, F(0)) + k
    return out


def check(eq: str) -> list[str]:
    """返回违规列表（空 = 守恒）。"""
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
    return bad


def main(argv: list[str]) -> int:
    """无参数 = 审计**整个用例库**的 eq/eq_has（非零退出 = 有违规）。

    这两种用法是审计的基本要求：默认模式必须覆盖全库（首版默认什么都不查，
    P10 的 `4H^+ + Fe(OH)_3 -> 3H_2O + Fe^{3+}`（电荷 +4≠+3）因此长期躺在
    标准里没被抓到），且**违规必须反映到退出码**（否则 CI/批量审计拿不到
    信号）。"""
    n_bad = 0
    if argv and argv[0] == "--cases":
        from chemkit.testsuit import load_cases
        want = tuple(argv[1:])
        for c in load_cases(None):
            if want and not c["name"].startswith(want):
                continue
            if c.get("eq"):
                print(f"## {c['name']}")
                n_bad += bool(check(c["eq"]))
            for e in (c.get("eq_has") or []):
                if check(e):
                    n_bad += 1
        print(f"\n合计违规 {n_bad} 条")
        return 1 if n_bad else 0
    if not argv:
        return main(["--cases"])
    for e in argv:
        n_bad += bool(check(e))
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
