# -*- coding: utf-8 -*-
"""第 207 轮 · fixsim 说"修 pH 修不好 F31" ⟹ 查清走步为何仍不动。

`tools/fixsim.py` 实测（改 estimate_pH 返回账本电荷自洽值）：
  · 旁路触发 3 次，值 1.6197 -> 12.3803（正确）
  · 但 `pH_solver` 仍 1.62、`iters` 仍 5、`exit=no-cands`、账本不变、固相 0
⟹ 第 204/205 轮"改 pH 就能修 F31"的因果链**不成立**。

本脚本在**不改源码**的前提下逼近真相：给 `estimate_pH` 装一个**记录型**
shim，打印每次调用时
  · 传入的 ledger 是否含镓羟合物（即是否已到终态）
  · 原值 vs 电荷自洽值
  · He
并在**终态账本 + 正确 pH** 上手算沉淀候选的驱动与可达程度，看
"若走步用 12.38，沉淀步到底会不会被选中"。

用法： python tools/whystuck.py
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
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith("F31 ")][0]
subs = [{"name": n, "mol": m} for n, m in c["subs"]]

# ---------- ① 记录每次 estimate_pH 调用（仅当账本含镓羟合物时） ----------
_orig = eng.estimate_pH
HITS = []


def _shim(ledger, H_excess, V_, T_, T_K_, *a, **k):
    val = _orig(ledger, H_excess, V_, T_, T_K_, *a, **k)
    if any("Ga(OH)" in sp and m > 1e-6 for sp, m in ledger.items()):
        HITS.append((round(val, 4), round(H_excess, 6),
                     {sp: round(m, 5) for sp, m in ledger.items()
                      if m > 1e-6 and sp != "H_2O"}))
    return val


eng.estimate_pH = _shim
try:
    probe = {}
    r = eng.judge(subs, c.get("cond") or {"V_L": V}, T, _probe=probe)
finally:
    eng.estimate_pH = _orig
    assert eng.estimate_pH is _orig

print(f"=== 走步中「账本含镓羟合物」时的 estimate_pH 调用：{len(HITS)} 次 ===")
for val, he, led in HITS[:12]:
    print(f"  pH={val:8.4f}  He={he:+10.6f}  {led}")

print(f"\n=== 终态 ===")
print(f"  iters={probe.get('iters')} steps={probe.get('steps_n')} "
      f"exit={probe.get('exit')} pH_solver={probe.get('pH_solver')}")
led_f = dict(probe.get("ledger") or {})
print(f"  账本={ {k: round(v,6) for k,v in led_f.items() if k != 'H_2O'} }")

# ---------- ② 在终态账本上，用"正确的" pH 评沉淀候选 ----------
print(f"\n=== 在终态账本上换 pH 重评候选（关键实验）===")
print(f"  {'pH':>8} {'Ga(OH)3 沉淀候选 S':>20} {'ext_max':>10} "
      f"{'在场?':>6}  反应物")
for pH in (1.62, 5.0, 7.0, 12.38):
    cands = enumerate_candidates(led_f, -0.024007, pH, V, T_K, T, True)
    for cd in cands:
        names = set(cd.r) | set(cd.pr)
        if "Ga(OH)_3" not in names or cd.kind != "precip":
            continue
        ps = cd.pres_specs
        pres = all(led_f.get(x, 0.0) > eng.X_MIN for x in ps[0])
        try:
            S = eng.S_of(cd, led_f, V, pH, T_K, T, frozenset(),
                         eng.P_EXT_KPA, True, {})
        except Exception as exc:                            # noqa: BLE001
            S = float("nan")
            print(f"  S 异常 {type(exc).__name__}: {exc}")
        lim = [led_f.get(x, 0.0) / nu for x, nu in cd.r.items()
               if x not in (eng.WATER, eng.H_ION) and nu > 0]
        ex = min(lim) if lim else float("inf")
        print(f"  {pH:8.2f} {S:20.3f} {ex:10.4g} {str(pres):>6}  {cd.r}")

print("\n=== 判读 ===")
print("  ① 若走步中「含镓羟合物」的 estimate_pH 调用**极少/不含终态账本**")
print("     ⟹ 走步在到达该账本后**几乎没再评估过候选**（no-cands 提前退出）")
print("     ⟹ 病根是**退出/枚举**，不是 pH。")
print("  ② 若换 pH 后沉淀候选 S 由负转正 ⟹ pH 确实能改变判据，")
print("     但**必须在退出之前**生效才有用。")
