"""一次性：AB05 中 pH ≤ −1.0 时 `estimate_state` 内部量（tit/He_res/分支）。"""
import sys

sys.path.insert(0, ".")

import chemkit.speciation as S
from chemkit.core import pKw_of
from chemkit.data import load_tables
from chemkit import engine
from chemkit.testsuit import load_cases

_ORIG = S.estimate_state
N = [0]


def wrap(ledger, H_excess, V, T, T_K, cache=None, touch=None):
    out = _ORIG(ledger, H_excess, V, T, T_K, cache, touch)
    if out[0] <= -1.0 and N[0] < 3:
        N[0] += 1
        pKw = pKw_of(T_K)
        t = S._buffer_titration(dict(ledger), H_excess, V, T, pKw, T_K=T_K)
        tml = S._buffer_titration(dict(ledger), H_excess, V, T, pKw,
                                  multilevel=True, T_K=T_K)
        print(f"\n!! pH={out[0]} He_in={H_excess!r} He_res={out[2]!r}")
        print(f"   ml=False: pH={t[0]} He_res={t[1]!r}")
        print(f"   ml=True : pH={tml[0]} He_res={tml[1]!r}")
        print("   ledger:", {k: round(v, 8) for k, v in ledger.items()
                             if k != "H_2O" and abs(v) > 1e-12})
        print("   vled(ml=F):", {k: round(v, 8) for k, v in t[2].items()
                                 if k != "H_2O" and abs(v) > 1e-12})
    return out


S.estimate_state = wrap
# estimate_pH 内部按模块全局名调用 ✓；engine 的 estimate_pH 也指向它
T = load_tables()
c = [x for x in load_cases(None) if x["name"].startswith("AB05")][0]
engine.judge([{"name": n, "mol": m} for n, m in c["subs"]],
             c.get("cond") or {"V_L": 1.0}, T)
