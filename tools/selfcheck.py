# -*- coding: utf-8 -*-
"""第 204 轮 · **终态自洽性检测器**（本轮的核心诊断）。

第 204 轮把 F31 的病根又推进了一层，且推翻了第 203 轮"pH 算错"的说法：

  引擎终态账本 `{[Ga(OH)₄]⁻ 0.7431, Ga³⁺ 0.2534, [Ga(OH)]²⁺ 0.0034}`，
  自报 pH = 1.62。用**引擎自己的 β 值**反推该账本在 pH 1.62 上的表观
  形成常数：
        logK_app = log([Ga(OH)₄⁻]/[Ga³⁺]) + 4·pH = 13.30
  而库内 `logβ₄ = 37.6` ⟹ **该账本在自己的 pH 上距平衡 25 个 log 单位**。
  没有任何 pH 能让这个账本成立（pH 1.62 要求 β₄=10^13.3；β₄=10^37.6
  要求 pH 12.38）⟹ **账本本身就自相矛盾**，不是"pH 读错了"。

  机制：走步用含 H⁺ 配体的配位腿（`Ga³⁺+4H₂O→[Ga(OH)₄]⁻+4H⁺`，每 mol
  耗 4 mol 碱）吃掉了强碱储备，使体系落到"碱已耗尽而配合物仍在"的态；
  此后正确性所需的**再平衡步**（解配位/再沉淀）受限于账本里只剩的
  游离质子，走不动 ⟹ 停在矛盾态。

本检测器把这条判据**做成可全库跑的量化口径**：
  对每个用例的终态，取账本里该金属的羟合梯（beta, ligand=OH⁻），
  用**引擎自己的 logβ** 与自报 pH 重算各形态比例，比较账本实际比例。
  偏差 > 阈（log 单位）即为"自相矛盾态"。

用法： python tools/selfcheck.py [阈值，默认 3.0] [步长，默认 1]
"""
import io
import json
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
from chemkit.testsuit import load_cases                     # noqa: E402

THR = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
STEP = int(sys.argv[2]) if len(sys.argv) > 2 else 1
T = load_tables()

# 羟合梯表：center -> [(nu, logb, complex)]
LADDER = {}
for b in T.beta:
    if b.get("ligand") == "OH^-" and b.get("m", 1) == 1:
        LADDER.setdefault(b["center"], []).append(
            (b["nu"], b["logb"], b["complex"]))
for k in LADDER:
    LADDER[k].sort()
print(f"羟合梯金属 {len(LADDER)} 种；阈值 {THR} log 单位；步长 {STEP}")

cases = load_cases(None)
rows = []
n_scanned = 0
for i, c in enumerate(cases):
    if i % STEP:
        continue
    V = float((c.get("cond") or {}).get("V_L", 1.0))
    T_K = float((c.get("cond") or {}).get("T_K", 298.15))
    probe = {}
    try:
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": V}, T, _probe=probe)
    except Exception:                                       # noqa: BLE001
        continue
    led = dict(probe.get("ledger") or {})
    pH = probe.get("pH")
    if not led or not isinstance(pH, float):
        continue
    n_scanned += 1
    pKw = spec.pKw_of(T_K)
    pOH = pKw - pH
    for center, lad in LADDER.items():
        m_free = led.get(center, 0.0)
        if m_free <= spec.X_MIN:
            continue
        # 账本里实际的梯成员
        act = [(nu, lb, cx) for nu, lb, cx in lad if led.get(cx, 0.0) > 0.0]
        if not act:
            continue
        # 用 logβ 与自报 pH 预测各成员与游离中心的比例
        # ⚠️ **必须用引擎的条件常数**（SIT 已开）：`logK_T` 吃 Cand，
        # 这里没有 Cand，故按 `speciation.sit_logK` 的同一公式手算，
        # Δz² 用 `center + nu·OH⁻ → complex` 的电荷平方差（与 `_beta_dz2`
        # 同一构造，只是配体换成 OH⁻ 且 z=−1）。
        worst = 0.0
        detail = []
        _chg = spec.charge_of
        I = round(spec.ionic_strength(led, V), 1)   # 与引擎同口径量化
        for nu, lb, cx in act:
            m_cx = led.get(cx, 0.0)
            if m_cx <= spec.X_MIN:
                continue
            lb_eff = lb
            try:
                dz2 = _chg(cx) ** 2 - _chg(center) ** 2 - nu * 1.0
                lb_eff = spec.sit_logK(lb, dz2, None, I)
            except Exception:                               # noqa: BLE001
                lb_eff = lb
            pred = lb_eff + nu * pOH
            actual = math.log10(m_cx / m_free)
            dev = actual - pred
            detail.append((abs(dev), cx, dev, m_cx))
            worst = max(worst, abs(dev))
        if worst > THR:
            rows.append({"case": c["name"], "center": center, "pH": pH,
                         "worst_dev": worst,
                         "detail": [[cx, round(dv, 2), mc]
                                    for _a, cx, dv, mc in sorted(
                                        detail, reverse=True)]})

rows.sort(key=lambda r: -r["worst_dev"])
print(f"\n扫描 {n_scanned} 例；**自相矛盾态 {len(rows)} 例**（偏差 > {THR} log）")
print(f"\n{'用例':34s} {'pH':>7} {'中心':10s} {'最大偏差':>9}  形态(账本 vs 预测)")
for r in rows[:30]:
    d = "; ".join(f"{cx} Δ{dv:+.1f}" for cx, dv, _m in r["detail"][:3])
    print(f"{r['case'][:34]:34s} {r['pH']:7.3f} {r['center'][:10]:10s} "
          f"{r['worst_dev']:9.2f}  {d[:56]}")

f31 = [r for r in rows if str(r["case"]).startswith("F31")]
print(f"\nF31: {'命中 ✓ ' + str(round(f31[0]['worst_dev'], 2)) + ' log' if f31 else '未命中'}")

print("\n=== 判读 ===")
print("  命中 = 该用例终态账本的物种分布与**它自己声称的 pH** 相差 > 阈值")
print("  log 单位（用引擎自己的 β 算）⟹ 该终态在热力学上不可能存在，")
print("  属于走步把体系带进了矛盾态。")
print("  ⚠️ 注意双向：也可能是**系统确实还没平衡**（欠收敛），此时是残差问题；")
print("  需与 `resid_live` 对照区分（残差小却偏差大 = 真矛盾态）。")

out = os.path.join(ROOT, "logs", "selfcheck.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"thr": THR, "n_scanned": n_scanned, "rows": rows},
              f, ensure_ascii=False, indent=1, default=str)
print(f"\n已写 {os.path.relpath(out, ROOT)}")
