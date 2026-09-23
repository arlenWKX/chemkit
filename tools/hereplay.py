# -*- coding: utf-8 -*-
"""第 204 轮 · F31 走步**逐步记账回放**：He 与账本到底在哪儿对不上。

第 204 轮实测矛盾（必须查清，不能猜）：
  · 引擎终态账本 = {[Ga(OH)₄]⁻ 0.7431, Ga³⁺ 0.2534, [Ga(OH)]²⁺ 0.0034,
    Na⁺ 3, Cl⁻ 3}
  · 由该账本**独立**算净电荷（除 H⁺/OH⁻）= +0.024007
  · 引擎自报 `H_excess = 0.0`
  · 但 `_exec` 的记账语义下，仅第 1 步 `Ga³⁺ + 4H₂O → [Ga(OH)₄]⁻ + 4H⁺`
    （ext 0.743966）就应使 He 从 −3.0 增到 −0.0241（产 H⁺）
    ⟹ 与"自报 0.0"差 0.024，而"独立算净电荷 +0.024"正好是 −(−0.024)。
  到底哪个口径对、He 最终是多少，只能回放。

本脚本用 `judge` 拿到完整 `steps`，按 `_exec` 语义逐步重算 (ledger, He)，
并与引擎终态自报值对照，定位不一致发生的位置。

用法： python tools/hereplay.py [用例前缀，默认 F31]
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
import chemkit.normalize as _norm                           # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
V, T_K = 1.0, 298.15
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith(PRE + " ")][0]
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
led, He, _st, _un = _norm.normalize(
    subs, {"V_L": V, "T_K": T_K, "c_H": None, "c_OH": None, "pH": None,
           "p_kpa": eng.P_EXT_KPA}, T, {})
led0, He0 = dict(led), He
probe = {}
r = eng.judge(subs, c.get("cond") or {"V_L": V}, T, _probe=probe)

print(f"用例 {c['name']}")
print(f"起始: He = {He0:+.6f}  账本 = "
      f"{ {k: round(v, 6) for k, v in led0.items() if k != 'H_2O'} }")

steps = r.get("steps") or []
print(f"\n共 {len(steps)} 步（judge 回报）\n")
print(f"{'#':>3} {'extent':>12} {'ΔHe':>13} {'He 累计':>12}  反应")
cur, Hec = dict(led0), He0
for i, st in enumerate(steps, 1):
    if not isinstance(st, dict):
        print(f"{i:3d}  (非 dict: {str(st)[:60]})")
        continue
    eq = st.get("equation") or ""
    ex = st.get("extent")
    # 从方程字符串解析反应式不可靠 ⟹ 用 extent 与候选表对照
    # 简化：直接看 extent 与 He 的变化（引擎未直接给每步 ΔHe）
    print(f"{i:3d} {str(ex)[:12]:>12} {'?':>13} {'?':>12}  {eq[:64]}")

print("\n=== 引擎终态 vs 独立复算 ===")
ledf = dict(probe.get("ledger") or {})
print(f"  引擎 probe He      = {probe.get('H_excess')}")
print(f"  引擎 probe pH      = {probe.get('pH')}  (pH_solver="
      f"{probe.get('pH_solver')})")
net = 0.0
for sp, m in ledf.items():
    if sp == "H_2O" or m <= 0:
        continue
    z = spec.charge_of(sp)
    net += z * m
print(f"  终态账本 Σz·n     = {net:+.6f}   （除 H⁺/OH⁻）")
print(f"  ⟹ 若 He 与账本自洽，应有 He = {-net:+.6f}（He ≡ −Σz·n）")
print(f"  ⟹ 引擎自报 He = {probe.get('H_excess')}，差 "
      f"{(probe.get('H_excess') or 0.0) - (-net):+.6f}")

print("\n=== 逐步 extent 与理论 ΔHe（按反应式 H⁺ 计量）===")
# 由 extent 与反应的 H+ 计量算 ΔHe：需要反应的 nu_H。
# 用候选枚举在**起始态**上找同方程的候选取 nu_H 不可靠（多形态），
# 故这里直接从方程字符串数 H^+ 的系数——本用例方程简单，够用。
import re                                                   # noqa: E402


def nu_H_of(eq: str):
    """'A -> B' 形式下，右侧 H^+ 系数 − 左侧 H^+ 系数。"""
    if "->" not in eq:
        return None
    lhs, rhs = eq.split("->", 1)

    def coef(side, tok):
        # 匹配 '4H^+' 或 'H^+'（可带系数）
        m = re.findall(r"(\d*\.?\d*)\s*" + re.escape(tok), side)
        tot = 0.0
        for x in m:
            tot += float(x) if x else 1.0
        return tot
    return coef(rhs, "H^+") - coef(lhs, "H^+")


tot_dHe = 0.0
for i, st in enumerate(steps, 1):
    if not isinstance(st, dict):
        continue
    eq, ex = st.get("equation") or "", st.get("extent")
    if ex is None:
        continue
    n = nu_H_of(eq)
    d = (n * ex) if n is not None else None
    if d is not None:
        tot_dHe += d
    print(f"  {i:2d}. nu_H={str(n):>5} ext={ex:+.6f} ΔHe={d if d is None else round(d, 6)}  {eq[:52]}")
print(f"\n  各步 ΔHe 合计 = {tot_dHe:+.6f}")
print(f"  起始 He {He0:+.6f} + 合计 = {He0 + tot_dHe:+.6f}")
print(f"  引擎终态 He = {probe.get('H_excess')}")
print(f"  账本自洽 He = {-net:+.6f}")

out = os.path.join(ROOT, "logs", f"hereplay-{PRE}.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"case": c["name"], "He0": He0, "steps": steps,
               "probe_He": probe.get("H_excess"),
               "ledger_net": net, "sum_dHe": tot_dHe,
               "ledger": ledf}, f, ensure_ascii=False, indent=1, default=str)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
