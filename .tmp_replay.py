"""一次性诊断：F35 每次 estimate_pH 调用的 cache 命中/未命中复算对照。"""
import copy
import sys

sys.path.insert(0, ".")

from chemkit.data import load_tables
from chemkit import engine
from chemkit.testsuit import load_cases
from chemkit.speciation import estimate_pH as raw

ORIG = engine.estimate_pH
CALLS = []


def wrap(ledger, He, V, T, T_K, cache=None, touch=None):
    out = ORIG(ledger, He, V, T, T_K, cache, touch)
    CALLS.append([dict(ledger), He, V, T, T_K, copy.deepcopy(cache), touch, out])
    return out


engine.estimate_pH = wrap
T = load_tables()
c = [x for x in load_cases(None) if x["name"].startswith("F35")][0]
engine.judge([{"name": n, "mol": m} for n, m in c["subs"]],
             c.get("cond") or {"V_L": 1.0}, T)
print("calls", len(CALLS))
for i, rec in enumerate(CALLS[:14]):
    led, He, V, _T, T_K, ca, touch, out = rec
    a = raw(dict(led), He, V, _T, T_K, None, None)
    b = (raw(dict(led), He, V, _T, T_K, copy.deepcopy(ca), touch)
         if ca is not None else None)
    print(f"#{i} He={He:+.6f} engine={out:.4f} nocache={a:.4f} "
          f"cached={None if b is None else round(b, 4)} nkeys={len(led)} "
          f"touch={None if touch is None else sorted(touch)} "
          f"cache={None if ca is None else list(ca)}")
    if ca:
        for k, v in ca.items():
            print(f"     cache[{k}] keytuple={v[0]} presence={v[1]} "
                  f"touchsp={v[2]} rows={v[3]} cnt={v[4]}")
