"""pinned 捕获探针（第 199 轮通用版）：python tools/dev.py run .tmp_pin_cap.py <用例前缀>

捕获 pinned 接管的每次调用：账本、pinned 项、返回 pH；
并统计 PH_TAGS 里 pinned 命中/耗尽回绝次数。输出 ASCII。
"""
import json
import sys

import chemkit.speciation as sp
from chemkit.acidbase import charge_pH as orig_charge_pH
from chemkit.data import load_tables
from chemkit.testsuit import load_cases, run_case

CAP = []


def spy(ledger, V, T, T_K, **kw):
    r = orig_charge_pH(ledger, V, T, T_K, **kw)
    if kw.get("pinned"):
        CAP.append({"ledger": {k: v for k, v in ledger.items() if v > 1e-9},
                    "V": V, "T_K": T_K,
                    "pinned": [list(p) for p in kw["pinned"]], "r": r})
    return r


def main():
    prefix = sys.argv[1] if len(sys.argv) > 1 else "E55"
    T = load_tables()
    sp.charge_pH = spy
    sp.PH_TAGS = []
    cases = [c for c in load_cases(None) if c["name"].startswith(prefix)]
    assert cases, f"{prefix} not found"
    ok = run_case(cases[0], T, verbose=False)
    tags = list(sp.PH_TAGS)
    sp.PH_TAGS = None
    sp.charge_pH = orig_charge_pH
    assert sp.charge_pH is orig_charge_pH

    n_hit = sum(1 for t in tags if t == "电荷平衡精确解(pinned)")
    n_rej = sum(1 for t in tags if t == "pinned耗尽回绝")
    print(f"case ok={ok}; pinned hits={n_hit}, exhaustion-rejects={n_rej}, "
          f"captured={len(CAP)}")
    # 只打印前 3 个与最后 2 个不同的 (账本签名, r)
    seen = []
    for cap in CAP:
        sig = (tuple(sorted((k, round(v, 4)) for k, v in cap["ledger"].items()
                            if k != "H_2O")), cap["r"])
        if not seen or seen[-1] != sig:
            seen.append(sig)
    print(f"distinct states: {len(seen)}")
    show = seen[:4] + ([("...", None)] if len(seen) > 6 else []) + seen[-2:]
    for sig, r in show:
        if sig == "...":
            print("  ...")
            continue
        print(f"  r={r} ledger={dict(sig)}")
    fn = f".tmp_pin_cap_{prefix}.json"
    with open(fn, "w", encoding="utf-8") as f:
        json.dump(CAP, f, ensure_ascii=False, indent=1)
    print(f"written {fn}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
