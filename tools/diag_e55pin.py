"""E55 深探（第 198 轮）：捕获 pinned 接管瞬间的真实账本并复算。

做法：包一层 speciation.charge_pH（speciation 内部按模块全局名调用，
patch 模块属性对这条调用链有效），只记录带 pinned 的调用；
跑完 E55 后对每个捕获点复算：
  A) 当前口径：pop 游离阳离子后 + pinned（即线上行为）；
  B) 对照：同一账本不带 pinned 的 charge_pH 精确解（若有根）；
  C) 钉住浓度 vs 该阳离子总量（固相是否本不该在场的回检）。
输出 ASCII。
"""
import json
import sys

import chemkit.speciation as sp
from chemkit.acidbase import charge_pH as orig_charge_pH
from chemkit.core import charge_of
from chemkit.data import load_tables
from chemkit.testsuit import load_cases, run_case

CAP = []


def spy(ledger, V, T, T_K, **kw):
    r = orig_charge_pH(ledger, V, T, T_K, **kw)
    if kw.get("pinned"):
        CAP.append({"ledger": dict(ledger), "V": V, "T_K": T_K,
                    "pinned": [list(p) for p in kw["pinned"]], "r": r})
    return r


def main():
    T = load_tables()
    sp.charge_pH = spy
    cases = [c for c in load_cases(None) if c["name"].startswith("E55")]
    assert cases, "E55 not found"
    run_case(cases[0], T, verbose=False)
    sp.charge_pH = orig_charge_pH
    assert sp.charge_pH is orig_charge_pH

    print(f"captured pinned calls: {len(CAP)}")
    for i, cap in enumerate(CAP[:8]):
        led = cap["ledger"]
        V, T_K = cap["V"], cap["T_K"]
        print(f"\n== hit {i}: pinned r={cap['r']}")
        print("  ledger:", {k: round(v, 6) for k, v in led.items()
                             if v > 1e-9})
        print("  pinned:", cap["pinned"])
        # B) 同账本不带 pinned
        r_np = orig_charge_pH(led, V, T, T_K, fast=True)
        print(f"  no-pinned exact on popped ledger: {r_np}")
        # C) 钉住浓度 vs 总量回检
        from chemkit.core import pKw_of
        pKw = pKw_of(T_K)
        for pin in cap["pinned"]:
            z, pKsp, x, y = pin[:4]
            lad = pin[4] if len(pin) > 4 else []
            if cap["r"] is not None:
                cm = 10.0 ** ((y * (pKw - cap["r"]) - pKsp) / x)
                print(f"  pinned conc at r: [M]={cm:.3e} M, ladder={[(c[3], c[1], c[2]) for c in lad]}")
    with open(".tmp_e55_pin.json", "w", encoding="utf-8") as f:
        json.dump(CAP, f, ensure_ascii=False, indent=1)
    print("\nwritten .tmp_e55_pin.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
