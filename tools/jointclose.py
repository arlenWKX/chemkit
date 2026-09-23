# -*- coding: utf-8 -*-
"""第 203 轮 · 闭合行替换的**爆破半径测量**（改 joint.py 之前必做）。

计划改动：`joint.py::_solve_ph` 的闭合行由
    F_{m+1} = estimate_pH(led, He) − pH        （不连续启发式）
换成
    F_{m+1} = charge_pH(led) − pH              （精确电荷平衡，单调无分支）
这是 PHREEQC 的做法（`ph_unknown == charge_balance_unknown`）。

风险：闭合行把联立 pH 锚在 pH 机器读数上；换锚会移动**所有**走联立的
用例的解。必须先量：
  ① 引擎实际调用闭合行时，`charge_pH` 是否成功（返回 float）？
  ② 成功时 |charge_pH − estimate_pH| 有多大？只在"两口径本就不一致"的
     用例上才产生位移（那些正是目标病灶）；
  ③ 全库有多少例会移动、移动多少。

做法（只读）：给 `joint._solve_ph` 打一个**记录型**补丁——把 `_F` 里
闭合行那一步的 (estimate_pH, charge_pH) 采集下来，跑全量套件，然后汇总。
不改任何返回值，纯观测。

用法： python tools/jointclose.py [样本步长，默认 1=全量]
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
import chemkit.joint as joint                               # noqa: E402
import chemkit.acidbase as ab                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 1
T = load_tables()
cases = load_cases(None)
V_DEF, T_K = 1.0, 298.15

# ---- 观测钩子：包住 _solve_ph，逐次记录闭合行两口径 ----
# ⚠️ **必须 patch `eng._solve_ph`**：`engine.py` L50 是 `from .joint import
# _solve_ph`（**按名导入**），改 `joint._solve_ph` 对引擎热路径**无效**
# （lessons.md 已记这个坑两次）。同时 patch 两处以防万一。
REC = []
_orig = joint._solve_ph
_orig_eng = eng._solve_ph
CUR = {"case": "?"}


def _patched(ledger, H_excess, actives, V, T_K_, TT, gsup, p_ext_kpa,
             gas_escape, S_of, H_ION, WATER):
    rec = {"case": CUR["case"], "V": V, "n_act": len(actives)}
    try:
        rec["est"] = joint.estimate_pH(ledger, H_excess, V, TT, T_K_)
    except Exception as exc:                                # noqa: BLE001
        rec["est"] = f"EXC {type(exc).__name__}"
    try:
        rec["exact"] = ab.charge_pH(ledger, V, TT, T_K_)
    except Exception as exc:                                # noqa: BLE001
        rec["exact"] = f"EXC {type(exc).__name__}"
    try:
        rec["st"] = _orig(ledger, H_excess, actives, V, T_K_, TT, gsup,
                          p_ext_kpa, gas_escape, S_of, H_ION, WATER)[0]
    except Exception as exc:                                # noqa: BLE001
        rec["st"] = f"EXC {type(exc).__name__}"
    REC.append(rec)
    return _orig(ledger, H_excess, actives, V, T_K_, TT, gsup,
                 p_ext_kpa, gas_escape, S_of, H_ION, WATER)


joint._solve_ph = _patched
eng._solve_ph = _patched

print(f"全量 {len(cases)} 例，步长 {STEP}；采集联立闭合行两口径…")
n_ok = n_fail = 0
moved = []
for i, c in enumerate(cases):
    if i % STEP:
        continue
    CUR["case"] = c["name"]
    try:
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": V_DEF}, T)
    except Exception:                                       # noqa: BLE001
        pass
joint._solve_ph = _orig
eng._solve_ph = _orig_eng

print(f"\n=== ① 引擎实际调用闭合行的次数 ===")
print(f"  _solve_ph 调用 {len(REC)} 次")
for r in REC[:1]:
    print(f"  样例: {r}")

ok = [r for r in REC if isinstance(r.get("exact"), float)
      and isinstance(r.get("est"), float)]
bad = [r for r in REC if not isinstance(r.get("exact"), float)]
print(f"  charge_pH 成功 {len(ok)} 次 / 失败(含 None/异常) {len(bad)} 次")
if bad:
    kinds = {}
    for r in bad:
        k = str(r.get("exact"))[:40]
        kinds[k] = kinds.get(k, 0) + 1
    print("  失败形态:")
    for k, v in sorted(kinds.items(), key=lambda kv: -kv[1])[:8]:
        print(f"    {k!r} × {v}")

print(f"\n=== ② 两口径差值分布（仅成功样本）===")
if ok:
    ds = sorted(abs(r["exact"] - r["est"]) for r in ok)
    n = len(ds)
    for q, tag in ((0.5, "p50"), (0.9, "p90"), (0.99, "p99")):
        print(f"  |ΔpH| {tag} = {ds[min(int(n * q), n - 1)]:.4f}")
    print(f"  |ΔpH| max = {ds[-1]:.4f}")
    print(f"  差 <0.01 的占比: "
          f"{sum(1 for d in ds if d < 0.01) / n * 100:.1f}%")
    print(f"  差 >1.0  的占比: "
          f"{sum(1 for d in ds if d > 1.0) / n * 100:.1f}%")

print(f"\n=== ③ 位移最大的用例（这些是换锚会**改变解**的目标）===")
rows = sorted(ok, key=lambda r: -abs(r["exact"] - r["est"]))[:20]
print(f"  {'用例':38s} {'estimate':>10} {'charge_pH':>10} {'Δ':>8} {'n_act':>6} 状态")
for r in rows:
    print(f"  {r['case'][:38]:38s} {r['est']:>10.4f} {r['exact']:>10.4f} "
          f"{r['exact'] - r['est']:>+8.3f} {r['n_act']:>6} {r.get('st')}")

print(f"\n=== ④ 判读 ===")
print("  · 若绝大多数样本 |Δ|<0.01 ⟹ 换锚对现行通过用例**几乎无扰动**，")
print("    改动安全，且正是把'两口径本来就不一致'的那些例纠正过来。")
print("  · charge_pH 失败率高 ⟹ 必须保留启发式兜底（否则联立解会凭空失败）。")
print("  · charge_pH 成功且 |Δ| 大 ⟹ 那些用例的联立解**本来就锚错了**。")

out = os.path.join(ROOT, "logs", "jointclose.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"n_calls": len(REC), "n_ok": len(ok), "n_fail": len(bad),
               "worst": rows, "all": REC}, f, ensure_ascii=False, indent=1,
              default=str)
print(f"\n已写 {out}")
