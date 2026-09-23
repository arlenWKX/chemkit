# -*- coding: utf-8 -*-
"""第 203 轮 · F31 起始态账目定影：`normalize` 到底给了什么 (ledger, He, pH)。

已知硬矛盾：引擎 `H_excess = 0.0`（净零酸碱当量）却说 pH = 1.62（≈0.024 M
游离强酸）。化学上 0 净酸 ⟹ pH≈7（中性）。两者必有一错，且必须知道错在
normalize（投料记账）还是 pH 估计器。

用法： python tools/f31_init.py [用例前缀，默认 F31]
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
from chemkit.testsuit import load_cases                     # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith(PRE + " ")][0]
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
cond = {"V_L": 1.0, "T_K": 298.15, "c_H": None, "c_OH": None, "pH": None,
        "p_kpa": eng.P_EXT_KPA}
print(f"用例: {c['name']}  subs={subs}")

# ---- 勾 normalize（模块级函数，可勾）----
_orig = eng.normalize
REC = {}


def _wrapped(substances, conditions, TT, origins):
    led, he, steps, unk = _orig(substances, conditions, TT, origins)
    REC["ledger"] = dict(led)
    REC["He"] = he
    REC["unknown"] = unk
    REC["steps"] = [dict(x) if isinstance(x, dict) else str(x)
                    for x in (steps or [])]
    return led, he, steps, unk


eng.normalize = _wrapped
probe = {}
r = eng.judge(subs, {"V_L": 1.0}, T, _probe=probe)
eng.normalize = _orig

V, T_K = 1.0, 298.15
print("\n=== normalize 输出（起始账本）===")
led0 = REC.get("ledger") or {}
for sp, m in sorted(led0.items(), key=lambda kv: -abs(kv[1])):
    if abs(m) > 1e-9:
        print(f"  {sp:24s} {m:12.6g}  z={spec._charge_cached(sp):+d}")
print(f"  H_excess = {REC.get('He')!r}")
print(f"  unknown  = {REC.get('unknown')}")
print("  steps:")
for x in (REC.get("steps") or [])[:10]:
    print(f"    {str(x)[:120]}")

print("\n=== 由 normalize 起始账本独立求 pH（多个估计器）===")
led = dict(led0)
He = REC.get("He") or 0.0
for nm, fn in (("closed_pH", lambda: eng.closed_pH(led, He, V, T, T_K)),
               ("estimate_pH", lambda: eng.estimate_pH(led, He, V, T, T_K)),
               ("presentation_pH",
                lambda: eng.presentation_pH(led, He, V, T, T_K))):
    f = getattr(eng, nm, None)
    if f is None:
        print(f"  {nm:16s} [缺]")
        continue
    try:
        print(f"  {nm:16s} -> {fn()}")
    except Exception as exc:                                # noqa: BLE001
        print(f"  {nm:16s} 异常 {type(exc).__name__}: {exc}")

print("\n=== 对照：judge 终态 pH / He ===")
print(f"  probe pH={probe.get('pH')} pH_solver={probe.get('pH_solver')} "
      f"He={probe.get('H_excess')}")

print("\n=== 手工化学核算（以投料为准，不依赖引擎）===")
print("  GaCl_3 1 mol + NaOH 3 mol, V=1 L")
print("  强碱当量 = 3.0 mol ⟹ 起始 He 应 = -3.0 mol（碱过量）")
print("  终态若 Ga(OH)_3(s) ≈ 1 mol（Ksp 极小 → 溶解 Ga 由 pH 决定）")
print("  ⟹ pH 应 ≈ 11~12（饱和 Ga(OH)_3 的碱侧），绝不可能 1.62")

out = os.path.join(ROOT, "logs", f"f31init-{PRE}.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"case": c["name"], "init_ledger": led0, "init_He": He,
               "init_steps": REC.get("steps"), "probe": probe},
              f, ensure_ascii=False, indent=1, default=str)
print(f"\n已写 {out}")
