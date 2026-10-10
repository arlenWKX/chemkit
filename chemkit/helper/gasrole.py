# -*- coding: utf-8 -*-
"""气体角色探针：逐例打印**气体物种在账本里的真实状态**，以及 `S_of` 气体
分支在**旧式/新式**两种口径下各自采用的活度。

## 为什么需要这个工具

第 286 轮查明：`S_of` 的气体分支把「自产气体逸出后的残余分压约定」`A_GAS`
（`P_RES/P_STD = 1 kPa / 101.325 kPa ≈ 1e-2`）当成了**溶质活度的下限**——

    a = max(min(c_g, H·p_ext), A_GAS)          # 旧式（有缺陷）

后果是**同一个物种、同一个状态**，只因为写在箭头哪一边就拿到两个标准态：

* 写在**产物**侧 ⟹ 命中 `typ=4` 气体分支 ⟹ `a = max(·, A_GAS)`（被**抬高**）
* 写在**反应物**侧 ⟹ 命中 `typ=2` 走 `_logc_of` 溶质口径 ⟹ `a = c_g`

⟹ `S_fwd ≠ −S_rev`（热力学硬性违反）。第 286 轮改成

    a = min(c_g, H·p_ext)  if H is not None and c_g > X_MIN  else A_GAS

即 `A_GAS` 只保留给"**气体不在账**（已逸出/痕量/无 Henry 数据）"的场合。

本工具对每个气体物种同时打印两种取值与 Δlog₁₀a，用来：
① 定位受影响的用例；② 核对新口径是否落在化学事实上；③ 事后复查（跑一次读多次）。

用法：`python chemkit/helper/gasrole.py Y09 Sn11 [更多前缀…]`
"""
from __future__ import annotations

import math
import sys

sys.path.insert(0, ".")

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.candidates import A_GAS, P_EXT_KPA, X_MIN   # noqa: E402
from chemkit.data import henry_of, load_tables           # noqa: E402
from chemkit.engine import judge                         # noqa: E402
from chemkit.system import Reaction                      # noqa: E402
from chemkit.testsuit import load_cases                  # noqa: E402


def _lg(x: float) -> float:
    return math.log10(x) if x > 0.0 else float("-inf")


def main(keys: list[str]) -> int:
    T = load_tables()
    cases = load_cases(None)
    picks = [c for c in cases if c["name"].startswith(tuple(keys))]
    if not picks:
        print("未匹配到用例")
        return 2
    for c in picks:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        V = float(cond.get("V_L", 1.0) or 1.0)
        T_K = float(cond.get("T_K", 298.15) or 298.15)
        probe: dict = {}
        r = judge(subs, cond, T, _probe=probe)
        led = dict(probe.get("ledger") or {})
        esc: dict = {}
        for e in (r.get("escaped") or []):
            if isinstance(e, dict):
                esc[e.get("name")] = e.get("mol", 0.0)
        print("=" * 78)
        print(f"{c['name']}")
        print(f"  pH={r.get('final_pH')}  degree={r.get('degree')}  "
              f"changed={r.get('changed')}  V={V}L  T={T_K}K")
        ne = Reaction(r).net_equation
        print(f"  净方程: {'(无)' if ne is None else ne.plain()}")
        st = [(s.get("equation"), round(s.get("extent") or 0.0, 6))
              for s in (r.get("steps") or [])]
        print(f"  步数={len(st)}  步骤: {st[:8]}")
        print("  产物: " + ", ".join(
            f"{p['name']}={p['mol']:.4g}" for p in (r.get("production") or [])))
        if esc:
            print("  逸出: " + ", ".join(f"{k}={v:.4g}" for k, v in esc.items()))
        if c.get("note"):
            print(f"  注: {c['note']}")
        gs = sorted({s for s in list(led) + list(esc) if s in T.gases})
        if not gs:
            print("  [气] 本用例无气体物种在账/逸出")
        for s in gs:
            n_led = led.get(s, 0.0)
            c_g = n_led / V
            H = henry_of(T, s, T_K)
            if H is None:
                print(f"  [气] {s:9s} 账n={n_led:.4g} c={c_g:.4g}  H=None"
                      f"  → 新旧同取 A_GAS={A_GAS:.4g} (log {_lg(A_GAS):+.3f})")
                continue
            csat = H * P_EXT_KPA
            a_old = max(min(c_g, csat), A_GAS)
            a_new = min(c_g, csat) if c_g > X_MIN else A_GAS
            print(f"  [气] {s:9s} 账n={n_led:.4g} c={c_g:.4g} H={H:.4g} "
                  f"c_sat={csat:.4g} | 旧 a={a_old:.4g}(log {_lg(a_old):+.3f})  "
                  f"新 a={a_new:.4g}(log {_lg(a_new):+.3f})  "
                  f"Δlog={_lg(a_new) - _lg(a_old):+.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["Y09"]))
