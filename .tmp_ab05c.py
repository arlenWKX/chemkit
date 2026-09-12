"""一次性：AB05 中 pH ≤ −1.0 的调用栈（定位 −1.0 的来源）。"""
import sys
import traceback

sys.path.insert(0, ".")

from chemkit.data import load_tables
from chemkit import engine
from chemkit.testsuit import load_cases

_ORIG = engine.estimate_pH
N = [0]


def wrap(ledger, He, V, T, T_K, cache=None, touch=None):
    out = _ORIG(ledger, He, V, T, T_K, cache, touch)
    if out <= -1.0 and N[0] < 2:
        N[0] += 1
        print(f"\n!! estimate_pH → {out}  He={He!r}")
        print("   ledger:", {k: repr(v) for k, v in ledger.items()
                             if k != "H_2O" and abs(v) > 1e-12})
        print("   cache a:", None if not cache or "a" not in cache
              else cache["a"][3])
        traceback.print_stack()
    return out


engine.estimate_pH = wrap
T = load_tables()
c = [x for x in load_cases(None) if x["name"].startswith("AB05")][0]
engine.judge([{"name": n, "mol": m} for n, m in c["subs"]],
             c.get("cond") or {"V_L": 1.0}, T)
