# -*- coding: utf-8 -*-
"""第 271 轮 · **「可再分配族总量 / 账本固定电荷」分布普查**。

## 为什么量这个

第 270 轮查明 `H45` 的呈现 pH（5.89）与账本自身的电荷平衡（**11.885**，手算）
矛盾。两道门串联：

* `_deg`（h_c/o_c ∈ [0.1, 10]）—— 比值 12.9，不成立；
* `_exact_ok` ②（族成员 ≥2）—— 账本只有 1 个族成员。

撤 `_deg` **逐字节无变化**；把②从 `<2` 放宽到 `<1` ⟹ `H45` 修好（5.89→11.88 ✓）
但 `I31`/`IN1`（In 同一化学）残差 **0 → 1.353**、总量 **+3.608** ⟹ 已回退。

**关键猜想**：`_nf = 1` 时 `charge_pH` 的族再分配量取决于**该族总量**——
`H45` 的族总量是**痕量**（`Al³⁺ 1.46e-6`），相对账本固定电荷（~0.5）可忽略
⟹ `charge_pH` 退化为"账本电荷平衡"，**良定义**；而 In 那一类族量级不同。
本工具把这条猜想**量出来**，据此定条件（而不是拍阈值）。

## 口径

包 `speciation._exact_ok`（它拿得到 `ledger/He_res/T/V`，且**只在两道门
相关处被调用**）。对每次调用记录：

    M_fam   = Σ m(在 build_families 里的物种)          可再分配的总量
    Q_fix   = Σ |z·m|(不在族里的物种)                  固定电荷的**量级**
    E_fix   = |Σ z·m|(不在族里的物种)                  固定电荷的**净值**
    n_fam   = 族成员数（当前闸的判据）
    ratio   = M_fam / max(Q_fix, E_fix)               本工具要看的量

逐例聚合（max ratio、出现过的 n_fam 集合），并单独列出关注用例。

零侵入（只包模块级函数），串行跑（~2 min），运行期间**冻结 chemkit/**。
输出 `logs/famcharge_census.json`。

用法：python tools/famcharge_census.py
      python tools/famcharge_census.py H45 I31 La1
"""
from __future__ import annotations

import collections
import io
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.speciation as _sp                            # noqa: E402
from chemkit.acidbase import build_families                 # noqa: E402
from chemkit.candidates import WATER                        # noqa: E402
from chemkit.core import charge_of                          # noqa: E402
from chemkit.data import load_tables                        # noqa: E402

_ORIG = _sp._exact_ok
FAMS: frozenset = frozenset()
PER: dict = collections.defaultdict(
    lambda: {"n": 0, "ratio_max": 0.0, "nf": set(), "mfam_max": 0.0})
CASE = [""]


def _hook(ledger, He_res, T, V):
    r = _ORIG(ledger, He_res, T, V)
    m_fam = 0.0
    q_fix = e_fix = 0.0
    n_fam = 0
    for s, m in ledger.items():
        if m <= 0.0 or s == WATER or s.startswith("__"):
            continue
        if s in FAMS:
            m_fam += m
            n_fam += 1
        else:
            z = charge_of(s)
            q_fix += abs(z) * m
            e_fix += z * m
    den = max(q_fix, abs(e_fix))
    ratio = (m_fam / den) if den > 0 else float("inf")
    d = PER[CASE[0]]
    d["n"] += 1
    d["ratio_max"] = max(d["ratio_max"], ratio)
    d["mfam_max"] = max(d["mfam_max"], m_fam)
    d["nf"].add(n_fam)
    return r


def main(argv: list[str]) -> int:
    global FAMS
    pres = tuple(a for a in argv if not a.startswith("-"))
    T = load_tables()
    FAMS = frozenset(build_families(T))
    _sp._exact_ok = _hook
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    from chemkit.engine import judge
    t0 = time.time()
    for i, c in enumerate(cases, 1):
        nm = c["name"]
        if pres and not nm.startswith(pres):
            continue
        CASE[0] = nm
        try:
            judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                  c.get("cond") or {}, T)
        except Exception:                                   # noqa: BLE001
            pass
        if i % 300 == 0:
            print(f"  {i}/{len(cases)} … {time.time() - t0:.0f}s", flush=True)
    _sp._exact_ok = _ORIG

    out = {k: {"n": v["n"], "ratio_max": v["ratio_max"],
               "mfam_max": v["mfam_max"], "nf": sorted(v["nf"])}
           for k, v in PER.items()}
    with io.open(os.path.join(ROOT, "logs", "famcharge_census.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)

    # 分层：按 max ratio 分桶
    buckets = [(0, 1e-9), (1e-9, 1e-6), (1e-6, 1e-3), (1e-3, 1e-1),
               (1e-1, 1.0), (1.0, float("inf"))]
    print(f"\n共 {len(out)} 例有 `_exact_ok` 调用；按 **max(M_fam/Q_fix)** 分桶：")
    tot = 0
    for lo, hi in buckets:
        ks = [k for k, v in out.items() if lo <= v["ratio_max"] < hi]
        tot += len(ks)
        lbl = (f"[{lo:g}, {hi:g})" if hi != float("inf")
               else f"[{lo:g}, ∞)")
        print(f"   {lbl:<16} {len(ks):>5} 例")
    tgt = ("H45", "M03", "ALU2", "I31", "IN1", "La1", "Ce1", "Nd1",
           "MgN1", "NiE1", "F31")
    print("\n关注用例（`nf` = 出现过的族成员数；`ratio_max` = M_fam/Q_fix 最大）：")
    for pre in tgt:
        hit = [k for k in out if k.startswith(pre)]
        for k in hit:
            v = out[k]
            print(f"   {k[:40]:<40} n={v['n']:>5} nf={v['nf']} "
                  f"M_fam_max={v['mfam_max']:.3g} ratio_max={v['ratio_max']:.3g}")
    print("\n-> logs/famcharge_census.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
