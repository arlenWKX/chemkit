# -*- coding: utf-8 -*-
"""第 206 轮 · **因果验证**：把 F31 的 pH 改成账本电荷自洽值，能修好吗？

在动引擎源码之前先做这个实验——若"修 pH"不能修好 F31，则第 204/205 轮
锁定的因果链是错的，改 estimate_pH 就是白改（避免又一轮无效改动）。

做法（**不改 chemkit/ 源文件**）：
  1. 记录 `engine.estimate_pH` 原引用；
  2. 换成 shim：先算原值；若当次调用的 (ledger) 满足五条合取条件，
     则改返回 `acidbase.charge_pH(ledger, ...)`；
  3. 跑 F31，看断言与固相产量；
  4. 同进程内还原，并断言已还原（lessons 的 A/B 纪律）。

⚠️ 必须 patch `engine.estimate_pH`（engine.py 是**按名导入**，patch
`speciation.estimate_pH` 对热路径无效——lessons 已记两次前科）。

用法： python tools/fixsim.py [用例前缀，默认 F31]
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
import chemkit.acidbase as ab                               # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases, run_case           # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
SOLIDS = frozenset(T.solids)
KSP_CATS = frozenset(e["pair"][0] for e in T.ksp if e["pair"][1] == "OH^-")
# (c) 未解离酸/碱储备（引擎按中性分子记账 ⟹ 会使 Σz·n 看似不平衡）
UNDISSOC = frozenset(("HNO_3", "H_2SO_4", "H_2SeO_4", "HCl", "HBr", "HI",
                      "HClO_4", "CH_3COOH", "HC_2H_3O_2", "HF", "H_2CO_3",
                      "H_2S", "H_3PO_4", "H_2C_2O_4"))

STATS = {"calls": 0, "bypass": 0}
_detail = []


def _charge_consistent(led, V, T_K):
    """由账本净电荷求 pH（唯一确定）。"""
    net = 0.0
    for sp, m in led.items():
        if sp == "H_2O" or m <= 0 or sp.startswith("__"):
            continue
        try:
            net += spec.charge_of(sp) * m
        except Exception:                                   # noqa: BLE001
            return None
    pKw = spec.pKw_of(T_K)
    if net > 1e-6:
        return pKw + __import__("math").log10(net / V)
    if net < -1e-6:
        return -__import__("math").log10(-net / V)
    return None


def _make_shim(orig):
    def _shim(ledger, H_excess, V, T_, T_K, *a, **k):
        STATS["calls"] += 1
        val = orig(ledger, H_excess, V, T_, T_K, *a, **k)
        # (a) 无固相
        if any(m > spec.X_MIN and sp in SOLIDS for sp, m in ledger.items()):
            return val
        # (b) 在场 Ksp-OH 阳离子
        if not any(m > spec.X_MIN and sp in KSP_CATS
                   for sp, m in ledger.items()):
            return val
        # (c) 无未解离酸/碱储备
        if any(ledger.get(sp, 0.0) > 1e-6 for sp in UNDISSOC):
            return val
        ph_ch = _charge_consistent(ledger, V, T_K)
        if ph_ch is None:
            return val
        # (d)(e) 账本与 He 矛盾 且 两口径差 > 2
        net = 0.0
        for sp, m in ledger.items():
            if sp == "H_2O" or m <= 0 or sp.startswith("__"):
                continue
            try:
                net += spec.charge_of(sp) * m
            except Exception:                               # noqa: BLE001
                return val
        if abs(net + H_excess) <= 1e-6:
            return val
        if abs(ph_ch - val) <= 2.0:
            return val
        STATS["bypass"] += 1
        if len(_detail) < 12:
            _detail.append((round(val, 4), round(ph_ch, 4)))
        return ph_ch
    return _shim


_orig = eng.estimate_pH
eng.estimate_pH = _make_shim(_orig)
print(f"已装 shim（patch engine.estimate_pH）")
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
    print(f"  estimate_pH 调用 {STATS['calls']} 次，其中旁路触发 "
          f"{STATS['bypass']} 次")
    print(f"  旁路取值样例 (原值 -> 电荷自洽值): {_detail[:6]}")
    print(f"  终态 pH = {probe.get('pH')}  pH_solver = "
          f"{probe.get('pH_solver')}")
    print(f"  max|S| = {probe.get('max_abs_S')}   iters = {probe.get('iters')}"
          f"  exit = {probe.get('exit')}")
    led = dict(probe.get("ledger") or {})
    print(f"  终态账本（>1e-6）:")
    for sp, m in sorted(led.items(), key=lambda kv: -abs(kv[1])):
        if m > 1e-6 and not sp.startswith("__"):
            print(f"    {sp:24s} {m:12.6g}")
    # 断言细节
    r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
    amt = {}
    for e in r["production"] + r["final"]:
        amt[e["name"]] = max(amt.get(e["name"], 0.0), e["mol"])
    print(f"  Ga(OH)_3 产量 = {amt.get('Ga(OH)_3', 0.0):.6g} "
          f"（断言要 ≥ 0.9）")
finally:
    eng.estimate_pH = _orig
    assert eng.estimate_pH is _orig, "还原失败！"
    print("\n已还原 engine.estimate_pH（并断言还原成功）")

print("\n=== 判读 ===")
print("  若断言 PASS 且 Ga(OH)_3 ≥ 0.9 ⟹ **因果链成立**，按此判据改引擎源码。")
print("  若仍 FAIL ⟹ 修 pH 不足以修好 F31，需重新定位（别白改）。")
