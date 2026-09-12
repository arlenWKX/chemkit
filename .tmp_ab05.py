"""一次性：AB05 末段 `_titration_heap` 命中/未命中与堆序对照。"""
import sys

sys.path.insert(0, ".")

import chemkit.speciation as S
from chemkit.data import load_tables
from chemkit import engine
from chemkit.testsuit import load_cases

ORIG = S._titration_heap
LOG = []
CALLS = []


def wrap(ledger, entries, cache, ckey, touch, nominal=None):
    ent = cache.get(ckey) if cache is not None else None
    hit = bool(ent is not None and ent[0] == tuple(ledger) and all(
        (ledger.get(sp, 0.0) > S.X_MIN) == ent[1][sp] for sp in ent[2]))
    why = ""
    if ent is not None and not hit:
        why = ("keydiff" if ent[0] != tuple(ledger) else
               "presence" + str([(sp, ledger.get(sp, 0.0) > S.X_MIN, ent[1][sp])
                                 for sp in ent[2]]))
    heap, cnt = ORIG(ledger, entries, cache, ckey, touch, nominal)
    LOG.append((ckey, hit, why, [(e[2], round(e[0], 2)) for e in heap],
                {k: round(v, 8) for k, v in ledger.items() if k != "H_2O"}))
    return heap, cnt


S._titration_heap = wrap
_ORIG_BT = S._buffer_titration


def bt(ledger, He, V, T, pKw, multilevel=False, T_K=298.15, cache=None,
       touch=None):
    out = _ORIG_BT(ledger, He, V, T, pKw, multilevel, T_K, cache, touch)
    CALLS.append((He, out[0], out[1], len(LOG)))
    return out


S._buffer_titration = bt
# engine 已 import 名字，需同步
engine._buffer_titration = bt

T = load_tables()
c = [x for x in load_cases(None) if x["name"].startswith("AB05")][0]
engine.judge([{"name": n, "mol": m} for n, m in c["subs"]],
             c.get("cond") or {"V_L": 1.0}, T)
print(f"titration 调用 {len(CALLS)}；heap 调用 {len(LOG)}")
print("== 末 10 次 heap 构建 ==")
for i, (ckey, hit, why, order, led) in enumerate(LOG[-10:]):
    print(f"  ck={ckey} hit={int(hit)} {why} 堆序={order}\n      {led}")
print("== 末 10 次 ti tration ==")
for He, pH, Her, n in CALLS[-10:]:
    print(f"  He={He:+.6f} → pH={pH} He_res={Her:+.6g}")
