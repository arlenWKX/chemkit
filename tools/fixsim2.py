# -*- coding: utf-8 -*-
"""第 211 轮 · 因果模拟（第 210 轮三条件版）：改 pH 能否让 F31 动起来？

在**不改 chemkit/ 源码**的前提下，给 `engine.estimate_pH` 装 shim：当三条件
同时满足时返回 `acidbase.charge_pH` 的值，其余一律原样。

三条件（第 210 轮量化，爆破半径 = 1 例 = 只有 F31）：
  (a) 账本里**没有 pKa 质子族成员**
  (b) **无固相在场**
  (c) 引擎自认未达平衡：`max_abs_S > 3`
      —— (c) 在 `estimate_pH` 内部拿不到（那是探针量），故本模拟用
      **最接近的可判定替代**：要求"有 Ksp-OH 阳离子在场"（F31 的签名），
      并额外要求 He 与账本电荷不自洽。**实测后据结果再定最终条件。**

⚠️ 必须 patch `eng.estimate_pH`（engine 按名导入）。

用法： python tools/fixsim2.py [用例前缀，默认 F31]
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
import chemkit.acidbase as ab                               # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases, run_case           # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
SOLIDS = frozenset(T.solids)
KSP_CATS = frozenset(e["pair"][0] for e in T.ksp if e["pair"][1] == "OH^-")
FAM = frozenset(ab.build_families(T))

STATS = {"calls": 0, "bypass": 0, "no_fam": 0, "no_solid": 0, "cat": 0,
         "incons": 0}
SAMPLE = []


def _net(led):
    n = 0.0
    for sp, m in led.items():
        if sp == "H_2O" or m <= 0 or sp.startswith("__"):
            continue
        try:
            n += spec.charge_of(sp) * m
        except Exception:                                   # noqa: BLE001
            return None
    return n


def _make(orig):
    def _shim(ledger, H_excess, V, T_, T_K, *a, **k):
        STATS["calls"] += 1
        val = orig(ledger, H_excess, V, T_, T_K, *a, **k)
        # (a) 无 pKa 质子族成员
        if any(m > spec.X_MIN and sp in FAM for sp, m in ledger.items()):
            return val
        STATS["no_fam"] += 1
        # (b) 无固相
        if any(m > spec.X_MIN and sp in SOLIDS for sp, m in ledger.items()):
            return val
        STATS["no_solid"] += 1
        # (c) 有 Ksp-OH 阳离子在场（F31 签名；代替探针的 max|S|>3）
        if not any(m > spec.X_MIN and sp in KSP_CATS
                   for sp, m in ledger.items()):
            return val
        STATS["cat"] += 1
        # (d) 账本电荷与 He 不自洽
        net = _net(ledger)
        if net is None or abs(net + H_excess) <= 1e-6:
            return val
        STATS["incons"] += 1
        try:
            ph_ch = ab.charge_pH(ledger, V, T_, T_K)
        except Exception:                                   # noqa: BLE001
            return val
        if not isinstance(ph_ch, float):
            return val
        STATS["bypass"] += 1
        if len(SAMPLE) < 10:
            SAMPLE.append((round(val, 4), round(ph_ch, 4),
                           round(H_excess, 6)))
        return ph_ch
    return _shim


_orig = eng.estimate_pH
eng.estimate_pH = _make(_orig)
print("已装 shim（patch engine.estimate_pH）")
try:
    cases = {c["name"]: c for c in load_cases(None)}
    hit = [v for n, v in cases.items() if n.startswith(PRE + " ") or n == PRE]
    c = hit[0]
    probe = {}
    ok = run_case(c, T, verbose=False)
    eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
              c.get("cond") or {"V_L": 1.0}, T, _probe=probe)
    print(f"\n用例 {c['name']}")
    print(f"  断言: {'PASS' if ok else 'FAIL'}")
    print(f"  调用 {STATS['calls']} 次；逐级通过数: 无pKa族 "
          f"{STATS['no_fam']} -> 无固相 {STATS['no_solid']} -> 有Ksp阳离子 "
          f"{STATS['cat']} -> 电荷不自洽 {STATS['incons']} -> **旁路触发 "
          f"{STATS['bypass']}**")
    print(f"  旁路取值 (原 -> 新, He): {SAMPLE[:6]}")
    print(f"  终态 pH={probe.get('pH')}  pH_solver={probe.get('pH_solver')}")
    print(f"  iters={probe.get('iters')} steps={probe.get('steps_n')} "
          f"exit={probe.get('exit')} max|S|={probe.get('max_abs_S')}")
    led = dict(probe.get("ledger") or {})
    print("  终态账本:")
    for sp, m in sorted(led.items(), key=lambda kv: -abs(kv[1])):
        if m > 1e-6 and not sp.startswith("__"):
            print(f"    {sp:24s} {m:12.6g}")
    r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
    amt = {}
    for e in r["production"] + r["final"]:
        amt[e["name"]] = max(amt.get(e["name"], 0.0), e["mol"])
    print(f"  Ga(OH)_3 产量 = {amt.get('Ga(OH)_3', 0.0):.6g}（要 ≥ 0.9）")
finally:
    eng.estimate_pH = _orig
    assert eng.estimate_pH is _orig
    print("\n已还原 engine.estimate_pH（断言通过）")

print("\n=== 判读 ===")
print("  旁路触发 >0 且断言 PASS ⟹ 三条件有效，可改源码。")
print("  旁路触发 0 ⟹ 条件在**走步调用路径上**从不满足，需重新设计。")
print("  触发 >0 但仍 FAIL ⟹ 条件对但开关位置不对（要查 solve_extent 的调用点）。")
