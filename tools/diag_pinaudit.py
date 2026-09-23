"""第 198 轮诊断：pinned 接管的触发面测量 + A/B 对照（一次性探针，不入库）。

同一进程内跑两遍全套件（开/关 pinned），按用例记录：
  · pinned 分支命中次数（经 PH_TAGS 审计钩子，生产路径零成本）；
  · 用例通过/失败，求翻转集（pass->fail / fail->pass）。
输出 ASCII（控制台 GBK 纪律）。
"""
import json
import sys

import chemkit.speciation as sp
from chemkit.data import load_tables
from chemkit.engine import judge
from chemkit.testsuit import load_cases, run_case, RESULTS, TIMES


def sweep(T, cases, pinned_on):
    sp.PINNED_TAKEOVER = pinned_on
    out = {}
    TIMES.clear()
    RESULTS.clear()
    for c in cases:
        sp.PH_TAGS = []
        ok = run_case(c, T, verbose=False)
        tags = sp.PH_TAGS
        sp.PH_TAGS = None
        n_pin = sum(1 for t in tags if "pinned" in t)
        n_ex = sum(1 for t in tags if t == "电荷平衡精确解")
        out[c["name"]] = {"ok": ok, "pin": n_pin, "ex": n_ex}
    return out


def main():
    T = load_tables()
    judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, T)   # 预热
    cases = load_cases(None)

    on = sweep(T, cases, True)
    off = sweep(T, cases, False)
    sp.PINNED_TAKEOVER = True   # 还原（A/B 纪律：末尾断言已还原）
    assert sp.PINNED_TAKEOVER is True

    n_on = sum(1 for v in on.values() if v["ok"])
    n_off = sum(1 for v in off.values() if v["ok"])
    flipped_bad = sorted(n for n in on if on[n]["ok"] and not off[n]["ok"])
    flipped_good = sorted(n for n in on if off[n]["ok"] and not on[n]["ok"])
    hits = {n: v["pin"] for n, v in on.items() if v["pin"]}
    hits_pass = sorted(n for n in hits if on[n]["ok"])
    hits_fail = sorted(n for n in hits if not on[n]["ok"])

    print(f"suite pinned=ON : {n_on}/{len(on)}")
    print(f"suite pinned=OFF: {n_off}/{len(off)}")
    print(f"flipped pass->fail ({len(flipped_bad)}): {flipped_bad}")
    print(f"flipped fail->pass ({len(flipped_good)}): {flipped_good}")
    print(f"pinned-hit cases ({len(hits)}), total hits {sum(hits.values())}")
    print(f"  hit & passing ({len(hits_pass)}): {hits_pass}")
    print(f"  hit & failing ({len(hits_fail)}): {hits_fail}")
    with open(".tmp_pin_audit.json", "w", encoding="utf-8") as f:
        json.dump({"on": on, "off": off, "hits": hits,
                   "flipped_bad": flipped_bad, "flipped_good": flipped_good},
                  f, ensure_ascii=False, indent=1)
    print("written .tmp_pin_audit.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
