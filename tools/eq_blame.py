"""对 9 条不守恒期望式：打印实际引擎输出与标准，判定责任方。"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chemkit.data import load_tables                    # noqa: E402
from chemkit.core import elements_of, charge_of         # noqa: E402
from chemkit.engine import judge                        # noqa: E402
from chemkit.system import Reaction, _parse_equation    # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FLAGGED = ("P10", "NR95", "Y14", "U10", "N30", "D32", "W11", "Fe31",
           "21 AgBr")


def bal(eq):
    if not isinstance(eq, str):
        return "n/a"
    try:
        r, p = _parse_equation(eq)
    except Exception:                                   # noqa: BLE001
        return "parse-fail"
    el = {}
    q = 0
    for d, sg in ((r, -1), (p, 1)):
        for sp, nu in d.items():
            q += sg * nu * charge_of(sp)
            for e, c in elements_of(sp).items():
                el[e] = el.get(e, 0) + sg * nu * c
    bad = {k: v for k, v in el.items() if v}
    return "OK" if not bad and not q else f"不守恒 {bad} q={q:+g}"


def main():
    T = load_tables()
    judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, T)
    cases = json.load(open(os.path.join(HERE, "chemkit/data/tests.json"),
                           encoding="utf-8"))
    for c in cases:
        if not any(c["name"].startswith(p) for p in FLAGGED):
            continue
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        r = judge(subs, c.get("cond") or {"V_L": 1.0}, T)
        got = Reaction(r).net_equation
        exp = c.get("eq")
        print(f"== {c['name'][:40]}")
        print(f"   cond      {c.get('cond')}")
        print(f"   标准 eq   {exp}")
        print(f"   标准平衡  {bal(exp) if isinstance(exp, str) else '-'}")
        print(f"   引擎实际  {got}")
        print(f"   引擎平衡  {bal(got)}")
        print(f"   一致?     {got == exp}")


if __name__ == "__main__":
    main()
