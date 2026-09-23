# -*- coding: utf-8 -*-
"""第 204 轮 · 自洽性检测器**有效性验证**（先验尺子，再量东西）。

第 204 轮踩坑：第一版检测器直接用表内 `logb`（I→0）比较，报出 193 例
"自相矛盾"，**连已知正确的 `16 AlCl3+3NaOH`（pH 6.269、Al(OH)₃ 0.9978）
也被报 Δ−63.4** ⟹ 全是 SIT 项造成的假阳性。
lessons 明写"先验证尺子，再量东西"——故本脚本用两个**判据已知**的用例
校准检测器：

  · 正对照 `16 AlCl3+3NaOH`：套件 PASS、pH 6.269、固相 0.9978 ⟹
    检测器**必须报小偏差**（否则是假阳性，不可用）。
  · 负对照 `F31 GaCl3+3NaOH`：套件 FAIL、pH 1.62、固相 0 ⟹
    检测器**应当报大偏差**（否则漏检）。

用法： python tools/selfcheck_val.py
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
from chemkit.testsuit import load_cases, run_case           # noqa: E402

T = load_tables()
cases = {c["name"]: c for c in load_cases(None)}

LADDER = {}
for b in T.beta:
    if b.get("ligand") == "OH^-" and b.get("m", 1) == 1:
        LADDER.setdefault(b["center"], []).append(
            (b["nu"], b["logb"], b["complex"]))


def measure(c):
    V = float((c.get("cond") or {}).get("V_L", 1.0))
    T_K = float((c.get("cond") or {}).get("T_K", 298.15))
    probe = {}
    eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
              c.get("cond") or {"V_L": V}, T, _probe=probe)
    led = dict(probe.get("ledger") or {})
    pH = probe.get("pH")
    if not led or not isinstance(pH, float):
        return None
    pOH = spec.pKw_of(T_K) - pH
    I = round(spec.ionic_strength(led, V), 1)
    out = []
    for center, lad in LADDER.items():
        m_free = led.get(center, 0.0)
        if m_free <= spec.X_MIN:
            continue
        for nu, lb, cx in lad:
            m_cx = led.get(cx, 0.0)
            if m_cx <= spec.X_MIN:
                continue
            dz2 = spec.charge_of(cx) ** 2 - spec.charge_of(center) ** 2 - nu
            lb_eff = spec.sit_logK(lb, dz2, None, I)
            dev = math.log10(m_cx / m_free) - (lb_eff + nu * pOH)
            out.append((abs(dev), center, cx, dev, lb, lb_eff, I))
    out.sort(reverse=True)
    return probe, out


TARGETS = [("16 AlCl3+3NaOH", "正对照（套件 PASS）"),
           ("F31 GaCl3+3NaOH", "负对照（套件 FAIL）")]
print(f"{'用例':26s} {'套件':>6} {'pH':>8} {'I':>6} {'最大偏差':>10}  形态")
for pre, tag in TARGETS:
    hit = [v for n, v in cases.items() if n.startswith(pre)]
    if not hit:
        print(f"  {pre}: 找不到用例")
        continue
    c = hit[0]
    ok = run_case(c, T, verbose=False)
    m = measure(c)
    if m is None:
        print(f"  {c['name'][:26]:26s} 无探针")
        continue
    probe, out = m
    if not out:
        print(f"  {c['name'][:26]:26s} {'PASS' if ok else 'FAIL':>6} "
              f"{probe.get('pH'):8.3f} {'-':>6} {'(无梯成员)':>10}")
        continue
    a, center, cx, dev, lb, lbe, I = out[0]
    print(f"  {c['name'][:26]:26s} {'PASS' if ok else 'FAIL':>6} "
          f"{probe.get('pH'):8.3f} {I:6.2f} {dev:+10.2f}  "
          f"{cx} (logb {lb} → 条件 {lbe:.2f})  [{tag}]")
    for a2, _ct, cx2, dev2, lb2, lbe2, _I in out[1:4]:
        print(f"  {'':26s} {'':>6} {'':>8} {'':>6} {dev2:+10.2f}  {cx2}")

print("\n=== 判读 ===")
print("  正对照偏差应 **小**（< ~3）⟹ 检测器不误报；")
print("  负对照偏差应 **大** ⟹ 检测器能抓真矛盾态。")
print("  两者都满足才可以把检测器用于全库普查。")
