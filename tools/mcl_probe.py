# -*- coding: utf-8 -*-
"""第 291 轮 · [MCl] 族：重构 E25 走步中间态 B，直接对拍
"拾取时看到的 S" 与 "S_of 平推扫描" 的不一致。

背景：E25 末态 A（Al3+=0.0098, [AlCl]2+=1.99, Cl-=4.01）上
S(解离)=+1.98，xscan 平推 0.612 处过零；但走步在 B=A+0.61·解离 上
拾取**络合**方向 S=+2.55 并执行 0.586 ⟹ 账本大摆锤。
本探针重构 B 态账本，枚举候选，打印该通道两个方向的 S。

用法： python tools/mcl_probe.py [E25]
"""
from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

AL, CL, MCL = "Al^{3+}", "Cl^-", "[AlCl]^{2+}"


def main(argv: list[str]) -> int:
    pre = argv[0] if argv else "E25"
    T = load_tables()
    c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
         if n.startswith(pre)][0]
    subs = [{"name": n, "mol": m} for n, m in c["subs"]]
    pr = {}
    r = eng.judge(subs, c.get("cond") or {"V_L": 1.0}, T, _probe=pr)
    led = dict(pr.get("ledger") or {})
    He = pr.get("H_excess") or 0.0
    V = float((c.get("cond") or {}).get("V_L", 1.0))
    T_K = float((c.get("cond") or {}).get("T_K", 298.15))
    pH = pr.get("pH_solver") or 7.0
    print(f"用例 {c['name']}  末态 pH={r.get('final_pH')} He={He:.6g}")
    print(f"  账本A={ {k: round(v, 8) for k, v in led.items()} }")

    # 末态 A 上的枚举（对照 crawl2）
    for tag, ledx, pHx in (("A(末态)", led, pH),):
        _report(T, ledx, He, pHx, V, T_K, tag)

    # 重构 B：A 反向施加"最后一次络合步"（x=0.610095，即解离方向 +0.61）
    x_last = 0.610095
    ledB = dict(led)
    ledB[AL] = ledB.get(AL, 0.0) + x_last
    ledB[CL] = ledB.get(CL, 0.0) + x_last
    ledB[MCL] = ledB.get(MCL, 0.0) - x_last
    print(f"\n重构 B 账本（A + {x_last}·解离）: "
          f"{ {k: round(v, 8) for k, v in ledB.items()} }")
    import chemkit.speciation as spec
    pHB = spec.estimate_pH(ledB, He, V, T, T_K)
    print(f"  B 态 estimate_pH = {pHB:.4f}")
    _report(T, ledB, He, pHB, V, T_K, "B(重构)")

    # B 态上沿络合方向再平推扫描 S(x)
    cands = enumerate_candidates(ledB, He, pHB, V, T_K, T, True)
    cd = next(cd for cd in cands
              if AL in cd.r and CL in cd.r and MCL in cd.pr
              and cd.kind == "complex")
    print(f"\n=== B 态络合方向 S(x) 平推（x_max={ledB[AL]:.6g}）===")
    xm = min(ledB[AL], ledB[CL])
    for k in range(9):
        x = xm * k / 8
        led2 = dict(ledB)
        led2[AL] -= x
        led2[CL] -= x
        led2[MCL] = led2.get(MCL, 0.0) + x
        S = eng.S_of(cd, led2, V, pHB, T_K, T, frozenset(),
                     eng.P_EXT_KPA, True, None)
        print(f"  x={x:10.6f}  S_complex={S:+9.4f}")
    return 0


def _report(T, led, He, pH, V, T_K, tag):
    cands = enumerate_candidates(led, He, pH, V, T_K, T, True)
    print(f"\n--- {tag} 枚举候选中 AlCl 通道 ---")
    for cd in cands:
        if MCL not in cd.r and MCL not in cd.pr:
            continue
        S = eng.S_of(cd, led, V, pH, T_K, T, frozenset(),
                     eng.P_EXT_KPA, True, {})
        print(f"  kind={cd.kind:<9} S={S:+9.4f}  "
              f"r={cd.r} pr={cd.pr}")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
