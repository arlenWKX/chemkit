# -*- coding: utf-8 -*-
"""第 228 轮 · **签名普查**：`resid_live>1` 的例子里有多少命中"pH 跳变"签名。

第 227 轮抽检发现该签名与化学体系无关（Ga/Al/Ca/Ba 全命中），
故 53 例的旧估计**偏小**。本脚本对全库逐例测签名，得到真实命中数，
用它重新评估"修分支 4"的收益。

签名判据（与 `tools/xscan.py` 同源）：
  对某例终态上 |S| 最大的**两侧候选**，按 `_exec` 语义推 `x ∈ (0, x_max]`，
  测 `estimate_pH` 是否出现 >1 单位的跳变（且随后 `S` 变号）。

判据写成"**扫描中出现 pH 跳变 > 1 单位**"这一条即可——
第 209/217/227 轮三次实测都伴随 `S` 变号，故不必再加第二条件。

用法： python tools/signature_census.py [阈值 resid_live，默认 1.0] [步长，默认 1]
"""
import io
import json
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
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

THR = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
STEP = int(sys.argv[2]) if len(sys.argv) > 2 else 1
NPTS = 10
T = load_tables()

# 从留档取候选名单（resid_live > THR），避免重跑
d = json.load(io.open(os.path.join(ROOT, "logs", "suite-latest.json"),
                      encoding="utf-8"))
cand_names = {c["name"] for c in d
              if (c.get("resid_live") or 0) > THR}
resid_of = {c["name"]: c.get("resid_live") for c in d}
print(f"留档中 resid_live > {THR} 的用例：{len(cand_names)} 例")

cases = {c["name"]: c for c in load_cases(None)}
V, T_K = 1.0, 298.15


def step(ledger, cd, d_, x):
    out = dict(ledger)
    he = 0.0
    rr = cd.r if d_ > 0 else cd.pr
    pp = cd.pr if d_ > 0 else cd.r
    for s_, nu in rr.items():
        if s_ == eng.H_ION:
            he -= nu * x
        elif s_ != eng.WATER:
            out[s_] = out.get(s_, 0.0) - nu * x
    for s_, nu in pp.items():
        if s_ == eng.H_ION:
            he += nu * x
        elif s_ != eng.WATER:
            out[s_] = out.get(s_, 0.0) + nu * x
    return out, he


hits = []
miss = []
err = []
for i, nm in enumerate(sorted(cand_names)):
    if i % STEP:
        continue
    c = cases.get(nm)
    if c is None:
        miss.append((nm, "用例未找到"))
        continue
    try:
        probe = {}
        eng.judge([{"name": a, "mol": b} for a, b in c["subs"]],
                  c.get("cond") or {"V_L": V}, T, _probe=probe)
        led = dict(probe.get("ledger") or {})
        pH = probe.get("pH_solver")
        He = probe.get("H_excess")
        if not led or not isinstance(pH, float) or He is None:
            miss.append((nm, "无探针画像"))
            continue
        cands = enumerate_candidates(led, He, pH, V, T_K, T, True)
        best = None
        for cd in cands:
            ps = cd.pres_specs
            if not (all(led.get(x, 0.0) > eng.X_MIN for x in ps[0])
                    and all(led.get(x, 0.0) > eng.X_MIN for x in ps[1])):
                continue
            try:
                S = eng.S_of(cd, led, V, pH, T_K, T, frozenset(),
                             eng.P_EXT_KPA, True, {})
            except Exception:                               # noqa: BLE001
                continue
            if best is None or abs(S) > abs(best[2]):
                best = (cd, 1 if S > 0 else -1, S)
        if best is None:
            miss.append((nm, "无两侧候选"))
            continue
        cd, ddir, S0 = best
        rr = cd.r if ddir > 0 else cd.pr
        xm = min(led.get(s_, 0.0) / nu for s_, nu in rr.items()
                 if s_ not in (eng.WATER, eng.H_ION))
        jump = 0.0
        prev = None
        for k in range(NPTS + 1):
            x = xm * k / NPTS
            lx, dhe = step(led, cd, ddir, x)
            try:
                ph_x = eng.estimate_pH(lx, He + dhe, V, T, T_K)
            except Exception:                               # noqa: BLE001
                continue
            if prev is not None and abs(ph_x - prev) > 1.0:
                jump = max(jump, abs(ph_x - prev))
            prev = ph_x
        if jump > 1.0:
            hits.append((nm, resid_of.get(nm), round(jump, 2),
                         str(cd.r) + " -> " + str(cd.pr)))
        else:
            miss.append((nm, f"无跳变(max|ΔpH|<1)"))
    except Exception as exc:                                # noqa: BLE001
        err.append((nm, f"{type(exc).__name__}: {exc}"))

print(f"\n=== 结果 ===")
print(f"**命中签名: {len(hits)} 例**   未命中: {len(miss)}   异常: {len(err)}")
print(f"\n命中清单（按 resid 降序）:")
print(f"  {'用例':36s} {'resid':>8} {'jump':>6}  通道")
for nm, rv, j, ch in sorted(hits, key=lambda t: -(t[1] or 0)):
    print(f"  {nm[:36]:36s} {str(round(rv, 3) if rv else rv):>8} {j:6.2f}  "
          f"{ch[:50]}")

if miss:
    print(f"\n未命中清单（前 20）:")
    for nm, why in miss[:20]:
        print(f"  {nm[:40]:40s} {why}")

out = os.path.join(ROOT, "logs", "signature_census.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"thr": THR, "hits": hits,
               "miss": [[n, w] for n, w in miss],
               "err": err}, f, ensure_ascii=False, indent=1)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
print("\n=== 判读 ===")
print("  命中数 = '修分支 4 能直接改善的例数'。")
print("  若接近总数 ⟹ 该缺陷是**全库性**的，架构改动收益极高。")
