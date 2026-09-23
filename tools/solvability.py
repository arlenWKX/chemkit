# -*- coding: utf-8 -*-
"""第 240 轮 · **判定二分**：B 组是"解不动"还是"无解"？

第 239 轮排除了"欠迭代"。本轮用**独立于引擎的鲁棒优化器**（Nelder-Mead，
纯标准库实现，不走 Jacobian）对同一残差向量做最小化：

  · 变量：`x_1..x_m`（各平衡的净程度）**加上 pH**（与 `_solve_ph` 同构）
  · 目标：`F(x, pH)` 的最大绝对值（与 `_resid` 同口径）
  · 若 Nelder-Mead 能把 `max|F|` 压到 `JOINT_TOL=0.05` 以下
    ⟹ 不动点**存在** ⟹ 属"**解不动**"（Newton 能力不足），值得换求解器；
  · 若压不下去（停在某个正值）⟹ 不动点**不存在** ⟹ 属"**无解**"，
    更强求解器也无用，应改"让走步不进入这种态"。

同时打印**逐方程残差**，指出是哪几条在互相冲突。

用法： python tools/solvability.py [用例前缀...]
"""
import io
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
import chemkit.joint as joint                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
PREFIXES = sys.argv[1:] or ["H45", "H43"]
cases = {c["name"]: c for c in load_cases(None)}


def build(pre):
    hit = [v for n, v in cases.items() if n.startswith(pre)]
    if not hit:
        return None
    c = hit[0]
    probe = {}
    eng.judge([{"name": a, "mol": b} for a, b in c["subs"]],
              c.get("cond") or {"V_L": V}, T, _probe=probe)
    led = dict(probe.get("ledger") or {})
    pH0 = probe.get("pH_solver") or 7.0
    He = probe.get("H_excess") or 0.0
    cands = enumerate_candidates(led, He, pH0, V, T_K, T, True)
    seen = {}
    for cd in cands:
        ps = cd.pres_specs
        if not all(led.get(s, 0.0) > eng.X_MIN for s in ps[0] + ps[1]):
            continue
        if cd.meta.get("slow") or cd.meta.get("deferred"):
            continue
        nk = (cd.netkey_fwd if cd.netkey_fwd <= cd.netkey_rev
              else cd.netkey_rev)
        try:
            S = eng.S_of(cd, led, V, pH0, T_K, T, frozenset(),
                         eng.P_EXT_KPA, True, {})
        except Exception:                                   # noqa: BLE001
            continue
        if nk not in seen or abs(S) > abs(seen[nk][1]):
            seen[nk] = (cd, S)
    act = sorted(((cd, 1 if S > 0 else -1, S) for cd, S in seen.values()),
                 key=lambda e: -abs(e[2]))[:joint.JOINT_MAX_M]
    return c, led, He, pH0, act


def make_F(led, He, act):
    m = len(act)

    def F(xs):
        led_w = dict(led)
        He_w = He
        for (cd, d, _S0), xj in zip(act, xs[:m]):
            if xj == 0.0:
                continue
            rr = cd.r if d > 0 else cd.pr
            pp = cd.pr if d > 0 else cd.r
            for s_, nu in rr.items():
                if s_ == eng.H_ION:
                    He_w -= nu * xj
                elif s_ != eng.WATER:
                    led_w[s_] = led_w.get(s_, 0.0) - nu * xj
            for s_, nu in pp.items():
                if s_ == eng.H_ION:
                    He_w += nu * xj
                elif s_ != eng.WATER:
                    led_w[s_] = led_w.get(s_, 0.0) + nu * xj
        ph = xs[m]
        out = []
        for cd, d, _S0 in act:
            S = eng.S_of(cd, led_w, V, ph, T_K, T, frozenset(),
                         eng.P_EXT_KPA, True, {})
            out.append(S if d > 0 else -S)
        # 闭合行与 _solve_ph 口径一致（启发式 pH 机器）
        out.append(joint.estimate_pH(led_w, He_w, V, T, T_K) - ph)
        return out
    return F, m


def nelder_mead(f, x0, step=1e-2, maxit=3000, tol=1e-10):
    """纯标准库 Nelder-Mead（不依赖 Jacobian）——用于判定"不动点是否存在"。"""
    n = len(x0)
    sim = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += step
        sim.append(p)
    vals = [f(p) for p in sim]
    for _ in range(maxit):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        sim = [sim[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol:
            break
        cen = [sum(sim[i][j] for i in range(n)) / n for j in range(n)]
        xr = [cen[j] + (cen[j] - sim[-1][j]) for j in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [cen[j] + 2.0 * (cen[j] - sim[-1][j]) for j in range(n)]
            fe = f(xe)
            sim[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            sim[-1], vals[-1] = xr, fr
        else:
            xc = [cen[j] + 0.5 * (sim[-1][j] - cen[j]) for j in range(n)]
            fc = f(xc)
            if fc < vals[-1]:
                sim[-1], vals[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    sim[i] = [sim[0][j] + 0.5 * (sim[i][j] - sim[0][j])
                              for j in range(n)]
                    vals[i] = f(sim[i])
    best = min(range(n + 1), key=lambda i: vals[i])
    return sim[best], vals[best]


for pre in PREFIXES:
    got = build(pre)
    if got is None:
        print(f"{pre}: 未找到")
        continue
    c, led, He, pH0, act = got
    F, m = make_F(led, He, act)
    print(f"\n{'=' * 74}")
    print(f"{c['name']}   m={m}  pH0={pH0}  actives 最大 |S|="
          f"{max(abs(a[2]) for a in act):.2f}")

    def obj(xs):
        try:
            return max(abs(v) for v in F(xs))
        except Exception:                                   # noqa: BLE001
            return 1e9

    x0 = [0.0] * m + [pH0]
    xb, fb = nelder_mead(obj, x0, step=0.05)
    print(f"  Nelder-Mead 最优 max|F| = **{fb:.4f}**"
          f"  （JOINT_TOL={joint.JOINT_TOL}）")
    print(f"  x(程度) = {[round(v, 5) for v in xb[:m]]}   pH* = {xb[m]:.4f}")
    Fr = F(xb)
    print(f"  逐方程残差（fail 时谁最大）:")
    for i, (cd, d, _S0) in enumerate(act):
        print(f"    F[{i}] = {Fr[i]:+9.3f}  d={d:+d} {cd.kind:9s} "
              f"{cd.r} -> {cd.pr}")
    print(f"    F[{m}] = {Fr[m]:+9.3f}  ← pH 闭合行")
    verdict = ("**解不动**（不动点存在，Newton 到不了）"
               if fb < joint.JOINT_TOL
               else "**无解**（不动点不存在或极远，需改走步策略）")
    print(f"  ⟹ 判定：{verdict}")
