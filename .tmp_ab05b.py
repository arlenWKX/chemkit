"""一次性：AB05 −1.0 pH 复现（同一账本、有无缓存、逐条 pop 追踪）。"""
import sys

sys.path.insert(0, ".")

import chemkit.speciation as S
from chemkit.candidates import X_MIN
from chemkit.data import load_tables
from chemkit import engine
from chemkit.testsuit import load_cases

_ORIG_BT = S._buffer_titration
HITS = []


def bt(ledger, He, V, T, pKw, multilevel=False, T_K=298.15, cache=None,
       touch=None):
    out = _ORIG_BT(ledger, He, V, T, pKw, multilevel, T_K, cache, touch)
    if out[0] is not None and out[0] <= -1.0 and len(HITS) < 3:
        alt = _ORIG_BT(dict(ledger), He, V, T, pKw, multilevel, T_K, None, None)
        print(f"\n!! pH={out[0]} He_res={out[1]:+.6g} ml={multilevel} "
              f"He_in={He:+.9g} freshcache→{alt[0]}")
        print("   ledger:", {k: repr(v) for k, v in ledger.items() if k != "H_2O"
                             and v > 1e-12})
        ent = cache.get("a") if cache else None
        print("   cache rows:", None if ent is None else ent[3])
        # 逐条 pop 追踪（multilevel 同参）
        he = -He if He < 0 else He
        led2 = dict(ledger)
        rows = [] if ent is None else list(ent[3])
        import heapq
        heapq.heapify(rows)
        print(f"   walk he0={he!r}")
        while rows and he > 0.0:
            pka, cnt, acid, base = heapq.heappop(rows)
            avail = led2.get(acid, 0.0)
            if avail <= 0.0:
                print(f"     pop {acid} pKa={pka} avail=0 → skip")
                continue
            take = min(he, avail)
            he -= take
            led2[acid] = avail - take
            led2[base] = led2.get(base, 0.0) + take
            mark = "  <== PLATEAU" if take < avail else ""
            print(f"     pop {acid} pKa={pka} avail={avail!r} take={take!r} "
                  f"→ {acid}={led2[acid]!r} {base}={led2[base]!r}{mark}")
        print(f"   walk 结束 he={he!r} 余堆={[r[2] for r in rows]}")
    HITS.append(1)
    return out


S._buffer_titration = bt
engine._buffer_titration = bt
T = load_tables()
c = [x for x in load_cases(None) if x["name"].startswith("AB05")][0]
engine.judge([{"name": n, "mol": m} for n, m in c["subs"]],
             c.get("cond") or {"V_L": 1.0}, T)
