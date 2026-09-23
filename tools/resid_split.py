# -*- coding: utf-8 -*-
"""第 215 轮 · resid 头部用例的"固相在场"二分 + **同族异果对比**。

第 214 轮归因：resid 第一名族 = 两性金属羟合形态（Al 占前 14 中 9 例），
与 F31 同源。本轮回答两个问题：

  ① 这些例**有固相在场**吗？（决定修法：有固相 -> pinned 储库精确解；
     无固相 -> charge_pH 适用面。第 213 轮已证 charge_pH 在有固相时必错。）
  ② **同为 Al + 强碱，为什么 `16 AlCl3+3NaOH` 做对了（Al(OH)₃ 0.998）
     而 `H43/M01/Amp14/H42` 卡住了？** 找出差异变量——这比逐个调参有价值。

一次运行打印全部所需字段（不再为看另一面重跑）。

用法： python tools/resid_split.py
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
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases, run_case           # noqa: E402

T = load_tables()
SOLIDS = frozenset(T.solids)
cases = {c["name"]: c for c in load_cases(None)}

# 第 214 轮 resid 头部 + 已知做对的同族对照
TOP = ["H45 Na[Al(OH)4]+HCl 半量", "M03 Na[Al(OH)4]+少量HCl",
       "H43 AlCl3+NaOH 1:3.5", "F31 GaCl3+3NaOH",
       "M01 AlCl3+少量NaOH", "Amp14 AlCl3+少量NaOH",
       "H42 AlCl3+NaOH 1:2.5", "E55 明矾+适量NaOH", "N34 明矾+适量Ba(OH)2",
       "36 Cl2+NaBr（合并H62）", "Y05 BaCl2+Na2SO3 白沉"]
GOOD = ["16 AlCl3+3NaOH（合并T46）", "N10 CrCl3+3NaOH",
        "H48 FeCl3+3NaOH", "I31 InCl3+3NaOH", "S32 SbCl3+NaOH"]


def probe(c):
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    p = {}
    ok = run_case(c, T, verbose=False)
    eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
              {"V_L": V}, T, _probe=p)
    led = dict(p.get("ledger") or {})
    solids = {sp: round(m, 6) for sp, m in led.items()
              if sp in SOLIDS and m > 1e-9}
    hydroxo = {sp: round(m, 5) for sp, m in led.items()
               if "(OH)" in sp and m > 1e-6}
    free = {sp: round(m, 5) for sp, m in led.items()
            if m > 1e-6 and "(OH)" not in sp and sp not in ("H_2O",)}
    return {"ok": ok, "pH": p.get("pH_solver"), "resid": p.get("max_abs_S"),
            "iters": p.get("iters"), "steps": p.get("steps_n"),
            "exit": p.get("exit"), "solids": solids, "hydroxo": hydroxo,
            "free": free}


print("=" * 78)
print("① resid 头部用例：**固相是否在场**")
print("=" * 78)
for nm in TOP:
    c = cases.get(nm)
    if c is None:
        print(f"  {nm}: 未找到")
        continue
    d = probe(c)
    print(f"\n  {nm}")
    print(f"    套件={'PASS' if d['ok'] else 'FAIL'}  pH={d['pH']}  "
          f"max|S|={d['resid']}  iters={d['iters']} steps={d['steps']}  "
          f"exit={d['exit']}")
    print(f"    **固相在场**: {d['solids'] if d['solids'] else '（无）'}")
    print(f"    羟合物: {d['hydroxo']}")
    print(f"    其它溶解: {d['free']}")

print()
print("=" * 78)
print("② 同族异果：Al/Cr/Fe/In/Sb + 强碱 —— 做对的 vs 卡住的")
print("=" * 78)
print(f"  {'用例':32s} {'套件':>5} {'pH':>8} {'max|S|':>8} {'it':>5} "
      f"{'固相':>10}  羟合物/游离")
for nm in GOOD + ["16 AlCl3+3NaOH（合并T46）"]:
    c = cases.get(nm)
    if c is None:
        continue
    d = probe(c)
    sm = sum(d["solids"].values()) if d["solids"] else 0.0
    print(f"  {nm[:32]:32s} {'PASS' if d['ok'] else 'FAIL':>5} "
          f"{str(d['pH'])[:8]:>8} {str(d['resid'])[:8]:>8} "
          f"{str(d['iters']):>5} {sm:>10.4f}  "
          f"{ {k: v for k, v in list(d['hydroxo'].items())[:2]} }")
for nm in ["H43 AlCl3+NaOH 1:3.5", "M01 AlCl3+少量NaOH",
           "Amp14 AlCl3+少量NaOH", "H42 AlCl3+NaOH 1:2.5",
           "E55 明矾+适量NaOH", "N34 明矾+适量Ba(OH)2",
           "H45 Na[Al(OH)4]+HCl 半量"]:
    c = cases.get(nm)
    if c is None:
        continue
    d = probe(c)
    sm = sum(d["solids"].values()) if d["solids"] else 0.0
    print(f"  {nm[:32]:32s} {'PASS' if d['ok'] else 'FAIL':>5} "
          f"{str(d['pH'])[:8]:>8} {str(d['resid'])[:8]:>8} "
          f"{str(d['iters']):>5} {sm:>10.6f}  "
          f"{ {k: v for k, v in list(d['hydroxo'].items())[:2]} }")

print("\n=== 判读 ===")
print("  · 若「卡住的」多为**无固相**或固相极微量，而「做对的」固相近 1 mol")
print("    ⟹ 差异变量 = **固相是否真的析出**；修法在'让析出发生'（形核/驱动），")
print("    而不是 pH 求解器。")
print("  · 若两者固相都近 1 mol 而残差差很多 ⟹ 差异在别处（碱量比/配体）。")
