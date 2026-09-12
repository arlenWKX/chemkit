"""一次性：N15 净方程呈现链（_format_equation/_rationalize/_beautify 输入输出）。"""
import sys

sys.path.insert(0, ".")

import chemkit.equations as E
from chemkit.data import load_tables
from chemkit.engine import judge
from chemkit.testsuit import load_cases

_fmt = E._format_equation
_rat = E._rationalize
_bea = E._beautify_big_coeff
_drop = E._drop_budget_ok


def fmt(c, p):
    out = _fmt(c, p)
    print(f"  _format_equation → {out!r}\n     c={ {k: round(v, 8) for k, v in c.items()} }"
          f"\n     p={ {k: round(v, 8) for k, v in p.items()} }")
    return out


def rat(vals, validate=None):
    out = _rat(vals, validate)
    print(f"  _rationalize({[round(v, 8) for v in vals]}) → {out}")
    return out


def bea(c, p):
    out = _bea(c, p)
    print(f"  _beautify_big_coeff → {None if out is None else out}")
    return out


def drop(bc, bp, c, p):
    out = _drop(bc, bp, c, p)
    print(f"  _drop_budget_ok → {out}  候选 c={ {k: round(v, 8) for k, v in c.items()} }"
          f" p={ {k: round(v, 8) for k, v in p.items()} }")
    return out


E._format_equation = fmt
E._rationalize = rat
E._beautify_big_coeff = bea
E._drop_budget_ok = drop
T = load_tables()
c = [x for x in load_cases(None) if x["name"].startswith("N15")][0]
r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
          c.get("cond") or {"V_L": 1.0}, T)
print("consumption:", r["consumption"])
print("production:", r["production"])
from chemkit.system import Reaction  # noqa: E402
rxn = Reaction(r)
print("net_equation:", rxn.net_equation)
