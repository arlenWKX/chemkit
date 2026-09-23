# -*- coding: utf-8 -*-
"""第 266 轮 · 对 `entry ∧ ¬hi` 的步，问 **`_pin_hi` 为什么是 False**。

第 266 轮实测：把 `_pin_mode` 从 `entry and hi` 改成 `entry`（冻结到入口口径），
残差质量 **302.8 → 268.6（−34.2）**、`resid_max` **25.415 → 24.193**（本会话
首次移动），但 **3 例 PASS->FAIL**，其中 `B26 FeCl3+NaOH(1:2)` 化学变错
（`Fe(OH)₃ = 0`，产物变成 `Fe₂O₃`）—— 而 B26 **正是 D14 当年引入 AND 冻结的
动因**（"B26 rec5 入口有固相、x_max 固相溶完"）。

⟹ 关键问题：**Al 族与 B26 的 `hi` 侧为什么不同？** 若能分开，"冻结到入口"
就只该用在其中一类上。

## 判据来源

`estimate_state` 的 pinned 块有三条出口（`PH_TAGS` 已标注）：

* `_pin` 为空（无"Ksp-OH 阳离子 + 其固相在场"）⟹ **无标签**；
* `charge_pH(pinned)` 返回 `None` ⟹ **无标签**；
* **储库耗尽回检失败** ⟹ 标签 `pinned耗尽回绝`。

故用 `PH_TAGS` 就能把"为什么 hi 侧不采用"分成**两类**（有无 `pinned耗尽回绝`），
不需要重写判据（lessons D2）。

用法：python tools/pin_hi_why.py I31 B26 Amp14 M01 E55
"""
from __future__ import annotations

import collections
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as _sp                            # noqa: E402
from chemkit.candidates import H_ION, WATER                 # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

_ORIG_EST = eng.estimate_pH
_ORIG_SE = eng.solve_extent


def _tags_of(led, he, V, T, T_K):
    save = _sp.PH_TAGS
    _sp.PH_TAGS = []
    try:
        _ORIG_EST(led, he, V, T, T_K)
        return list(_sp.PH_TAGS)
    finally:
        _sp.PH_TAGS = save


def _adopted(led, he, V, T, T_K):
    t = _tags_of(led, he, V, T, T_K)
    return bool(t) and t[-1] == "电荷平衡精确解(pinned)"


def main(argv: list[str]) -> int:
    pres = argv or ["I31"]
    T = load_tables()
    cases = {c["name"]: c for c in load_cases(None)}
    for pre in pres:
        hit = [v for k, v in cases.items() if k.startswith(pre)]
        if not hit:
            print(f"[跳过] 无 {pre}")
            continue
        c = hit[0]
        cond = c.get("cond") or {}
        V = float(cond.get("V_L", 1.0))
        T_K = float(cond.get("T_K", 298.15))
        REASON = collections.Counter()
        HOT: dict = {}

        def _se(cd, d, led, He, V_, T_K_, T_, *a, **k):
            out = _ORIG_SE(cd, d, led, He, V_, T_K_, T_, *a, **k)
            try:
                rr = cd.r if d > 0 else cd.pr
                pp = cd.pr if d > 0 else cd.r
                lim = [led.get(s, 0.0) / nu for s, nu in rr.items()
                       if s not in (WATER, H_ION) and nu > 0]
                xm = min(lim) if lim else 0.0
                if xm <= 0.0:
                    return out
                if not _adopted(led, He, V_, T_, T_K_):
                    return out
                led_hi = dict(led)
                for s in set(rr) | set(pp):
                    if s in (WATER, H_ION):
                        continue
                    led_hi[s] = max(0.0, led.get(s, 0.0)
                                    + (pp.get(s, 0) - rr.get(s, 0)) * xm)
                nu_H = pp.get(H_ION, 0) - rr.get(H_ION, 0)
                hi_t = _tags_of(led_hi, He + nu_H * xm, V_, T_, T_K_)
                if hi_t and hi_t[-1] == "电荷平衡精确解(pinned)":
                    return out
                # **判别假设**（B26 的 D14 注记说它"x_max 固相溶完"）：
                # hi 端该阳离子的 Ksp-OH 固相还在不在？在 ⟹ 钉住前提仍成立，
                # 只是 `charge_pH` 那边出了别的事；不在 ⟹ 前提真的没了。
                _sol = [e["solid"] for e in T.ksp
                        if e["pair"][1] == "OH^-" and e["pair"][0] in led_hi]
                _pres = [s for s in _sol if led_hi.get(s, 0.0) > 1e-6]
                solid_ok = bool(_pres)
                why = ("储库耗尽回绝" if "pinned耗尽回绝" in hi_t
                       else ("无标签" if not hi_t else "/".join(hi_t)))
                why = f"{why} | hi端固相{'在' if solid_ok else '**没了**'}"
                REASON[why] += 1
                HOT.setdefault(why, []).append(
                    (round(He, 6), round(xm, 8), round(out[0], 10), d))
            except Exception:                               # noqa: BLE001
                pass
            return out

        eng.estimate_pH = _ORIG_EST
        eng.solve_extent = _se
        try:
            eng.judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T)
        finally:
            eng.solve_extent = _ORIG_SE
        n = sum(REASON.values())
        print(f"\n=== {c['name']} ===  entry∧¬hi 共 {n} 段")
        for r, v in REASON.most_common():
            print(f"   {v:>6}  {r}")
        for r, lst in HOT.items():
            z = sum(1 for _he, _xm, x0, _d in lst if x0 <= 1e-4 * _xm)
            print(f"   [{r[:22]}] 其中 x*≈0（零推进）{z}/{len(lst)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
