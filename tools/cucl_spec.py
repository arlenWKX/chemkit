# -*- coding: utf-8 -*-
"""第 292 轮 · 高阶氯络合物族（Cu 6 例）重裁探针。

**用法**：`python tools/cucl_spec.py` —— 对 6 个 Cu 案例：
  1. 用**库值独立推导** Cu(II)-Cl⁻ 形态分布（不走引擎走步）：
     累积 β_k 的 SIT 条件常数 logβ_k(I)=logβ°+Δz²·D+ε·I（I=I_eff，
     反应参与物种口径，§192），与 Cl 物料守恒联立不动点；
  2. 并排打印引擎终态账本，核对两套独立计算是否互洽；
  3. 给出建议的 has_range 锚（推导值 ± 余量）。

依据（第 292 轮查明）：这 6 例的失败**不是走步缺陷**——账本里
[CuCl]⁺/[CuCl₂]/[CuCl₃]⁻/[CuCl₄]²⁻ 全在且满足平衡关系；失败在呈现层
**弱池折叠口径**（v0.4.0 事件口径：池成员净差并入 anchor Cu²⁺，净方程
只写 `2Fe³⁺+Cu→2Fe²⁺+Cu²⁺`）。使用者裁定："产物过于复杂可以不检查
总方程式" ⟹ 重裁 = 撤 eq、锁形态分布锚。
"""
from __future__ import annotations

import json
import math
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

# 库值（chemkit/data/beta.json，NAGRA/TDB 2020 校准）
#   (complex, logβ°, Δz²累积, ε, z²(络合物), 结合 Cl 数)
CU_LADDER = [
    ("[CuCl]^+",     0.83, -4, -0.05, 1, 1),
    ("[CuCl_2]",     0.60, -6, -0.10, 0, 2),
    ("[CuCl_3]^-",  -0.15, -6, None,  1, 3),   # 条目无 eps_rxn：只施 DH 主项
    ("[CuCl_4]^{2-}", -1.2, -4, -0.06, 4, 4),
]


def D(I: float) -> float:
    s = math.sqrt(I)
    return 0.51 * s / (1 + 1.5 * s)


def solve(cu_tot: float, cl_tot: float):
    """整体阻尼不动点（β 依赖 I_eff，I_eff 依赖形态，形态依赖 β）。"""
    c = cl_tot * 0.7
    fr = [0.0] * 4           # x_k / z
    for _ in range(500):
        z = cu_tot / (1 + sum(fr))
        new = []
        for k, (name, lb, dz2, eps, z2, ncl) in enumerate(CU_LADDER):
            xk = z * fr[k]
            # I_eff 参与者：Cu²⁺、Cl⁻、该络合物（§192 口径）
            I = 0.5 * (4 * z + c + z2 * xk)
            lg = lb + dz2 * D(I) + (eps * I if eps is not None else 0.0)
            new.append((10 ** lg) * c ** ncl)
        fr = [0.5 * a + 0.5 * b for a, b in zip(fr, new)]
        z = cu_tot / (1 + sum(fr))
        xs = [z * f for f in fr]
        c_new = cl_tot - sum(x * ncl for x, (_, _, _, _, _, ncl) in zip(xs, CU_LADDER))
        if abs(c_new - c) < 1e-12:
            c = c_new
            break
        c = 0.5 * c + 0.5 * c_new
    z = cu_tot / (1 + sum(fr))
    xs = [z * f for f in fr]
    return z, xs, c


CASES = {
    "39 FeCl3+Cu（合并E22/T29/H65）": (1.0, 6.0),
    "J12 Cu+FeCl3@350K 仍完全": (1.0, 6.0),   # sit_logK 的 A 不随 T 变、条目无 dH
    "O03 Cu+FeCl3 腐蚀": (1.0, 6.0),
    "RX1 FeCl3+Cu(过量)": (1.0, 6.0),
    "Y14 CuO+盐酸 溶解": (1.0, 2.0),
    "U10 Cu+H2O2+盐酸 溶解": (1.0, 2.0),
}


def engine_ledger(name):
    import json as j
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()
    cases = j.load(open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                       encoding="utf-8"))
    c = [x for x in cases if x["name"] == name][0]
    pr = {}
    r = judge([{"name": s, "mol": float(m)} for s, m in c["subs"]],
              c.get("cond") or {}, T, _probe=pr)
    amt = {}
    for e in (r.get("production") or []) + (r.get("final") or []):
        amt[e["name"]] = max(amt.get(e["name"], 0.0), e["mol"])
    return amt


def main():
    print(f"{'case':<34} {'species':<14} {'推导':>8} {'引擎':>8} {'比值':>6}")
    for name, (cu, cl) in CASES.items():
        z, xs, c = solve(cu, cl)
        amt = engine_ledger(name)
        rows = [("Cu^{2+}", z)] + [(CU_LADDER[i][0], xs[i]) for i in range(4)]
        rows.append(("Cl^-", c))
        print(f"-- {name}  (Cu={cu}, Cl={cl})")
        for sp, v in rows:
            e = amt.get(sp, 0.0)
            ratio = (e / v) if v > 1e-12 else float("nan")
            print(f"{'':<34} {sp:<14} {v:>8.4f} {e:>8.4f} {ratio:>6.2f}")


if __name__ == "__main__":
    main()
