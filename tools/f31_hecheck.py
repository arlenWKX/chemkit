# -*- coding: utf-8 -*-
"""第 203 轮 · F31 质子/碱量账目一致性核查（D15 形核帧缺陷的核心）。

化学事实：1 mol GaCl3 + 3 mol NaOH。NaOH 是强碱，其 OH^- 无论走哪条路
（沉淀成 Ga(OH)_3 固 / 配位成 [Ga(OH)_4]^-）都**必须**中和等当量的酸。
引擎终态 pH=1.62 ⟹ 加了 3 当量强碱却得酸性溶液，这是硬矛盾。

本脚本做两件事（只读）：
  ① 打印 probe 的 H_excess / pH / pH_solver；
  ② **独立复算**终态账本的质子过剩 He(ledger)（不依赖引擎），与引擎值比对。
     He(ledger) = Σ z_i·n_i  (z_i = 电荷，含 H^+ 自身 z=+1)
     再与"初始 He = -3.0（3 mol OH^- 相对中性水）"比对。
  ③ 逐物种列出对 He 的贡献，定位差额来自哪个物种。

用法： python tools/f31_hecheck.py [用例前缀，默认 F31]
"""
import io
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
from chemkit.testsuit import load_cases                     # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith(PRE + " ")][0]
print(f"用例: {c['name']}  subs={c['subs']}")

probe = {}
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
r = eng.judge(subs, c.get("cond") or {"V_L": 1.0}, T, _probe=probe)

V = float((c.get("cond") or {}).get("V_L", 1.0))
print(f"\n=== ① 引擎自报 ===")
for k in ("pH", "pH_solver", "exit", "iters", "steps_n", "H_excess",
          "max_abs_S", "frozen_n", "disabled_n"):
    print(f"  {k:12s} = {probe.get(k)}")

print(f"\n=== ② 独立复算终态账本质子过剩 He(ledger) ===")
led = probe.get("ledger") or {}
charged = []
for sp, m in led.items():
    if m <= 0:
        continue
    z = spec._charge_cached(sp)
    if z:
        charged.append((sp, m, z, z * m))
tot = sum(x[3] for x in charged)
print(f"  {'物种':24s} {'mol':>12s} {'z':>4s} {'z*mol':>12s}")
for sp, m, z, zm in sorted(charged, key=lambda x: -abs(x[3])):
    print(f"  {sp:24s} {m:12.6g} {z:+4d} {zm:+12.6g}")
print(f"  {'—':24s} {'':>12s} {'':>4s} {'—':>12s}")
print(f"  {'Σ z·n (净电荷，须=0)':24s} {'':>12s} {'':>4s} {tot:+12.4e}")

# H_excess 的物理意义：账本中"净质子过剩"= Σ_{i≠H2O,OH-,H+} z_i·n_i + n(H+) − n(OH-)
# 等价于 Σ z_i·n_i（含 H+/OH-；净电荷口径下二者抵消为 0 ⟹ 需单独列出）
nH = led.get("H^+", 0.0)
nOH = led.get("OH^-", 0.0)
print(f"\n  （账本未显式记录 H^+/OH^- 时，两者由 pH 隐含）")
print(f"  n(H^+)  = {nH:.6g} mol   n(OH^-) = {nOH:.6g} mol")

print(f"\n=== ③ 初始态 He（化学约定：3 mol NaOH ⟹ He = -3.0 mol）===")
init = {}
for n, m in c["subs"]:
    init[n] = init.get(n, 0.0) + m
he_init = 0.0
for n, m in sorted(init.items()):
    z = spec._charge_cached(n)
    print(f"  {n:24s} {m:12.6g} z={z:+3d}  z*m={z * m:+10.4g}")
    he_init += z * m
print(f"  Σ z·m(投料) = {he_init:+.6f} mol  （Cl^- 3×(-1) + Na^+ 3×(+1) "
      f"= 0；Ga^3+ 1×(+3) = +3 ⟹ +3.0）")
print("  说明：投料的 Σz·m 含 Ga^3+ 的 +3，而 'He' 约定只记**相对中性水的"
      "酸碱当量**⟹ 初始 He = -(强碱当量) = -3.0 mol。")

print(f"\n=== ④ 判读 ===")
print(f"  引擎 He = {probe.get('H_excess')}  vs  化学应得 He = -3.0")
_d = (probe.get("H_excess") or 0.0) - (-3.0)
print(f"  差 = {_d:+.6f} mol")
print("  残差 |Σz·n| 应≈0（电中性）；若非零则账本自身破坏电中性。")
print("  终态游离 H^+ ≈ V·10^(-pH)：pH=1.62 ⟹ ~0.024 mol **游离强酸**，"
      "与 'He=-3（碱过量）' 不可能共存 ⟹ 二者其一必错。")

print("\n=== ⑤ 走步事件（engine 自报 steps）===")
for i, st in enumerate((r.get("steps") or [])[:12], 1):
    print(f"  {i:2d}. {str(st)[:110]}")

out = os.path.join(ROOT, "logs", f"f31he-{PRE}.json")
import json                                                 # noqa: E402
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"case": c["name"], "probe": probe,
               "he_from_ledger": tot, "he_init_subs": he_init,
               "steps": r.get("steps")}, f, ensure_ascii=False,
              indent=1, default=str)
print(f"\n已写 {out}")
