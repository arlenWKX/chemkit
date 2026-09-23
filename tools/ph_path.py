# -*- coding: utf-8 -*-
"""第 276 轮 · **"呈现 pH 是哪条通路给的"三路并排**。

## 为什么需要

`tools/branch4_src.py H45` 实测：分支 4 的**最终选支是 `h_c = 0.374729`**
（⟹ pH **0.426**），而引擎**呈现 pH = 12.19**。也就是说分支 4 的意见被
**另一条通路覆盖**了。§1.12 是"两个 pH 口径打架"，但**到底是哪两条**、
各自给多少，必须逐路打出来才能定位。

本工具对**同一个终态账本**并排调用：

    ① estimate_pH（走步/呈现用的那条）+ 它的 PH_TAGS 标注（分支名）
    ② charge_pH（账本电荷平衡的精确解，`pinned=None` 与 `pinned=True` 各一次）
    ③ 直读 −log10(He/V)（残余游离强酸池，**精确记账**，不是估计）
    ④ 账本电荷平衡的**手算式**：Σz·n + He = 0 的余项（自洽性判据）

用法：python tools/ph_path.py H45
"""

from __future__ import annotations

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

import math                                                 # noqa: E402

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit. acidbase import charge_pH                     # noqa: E402
from chemkit.candidates import WATER                        # noqa: E402
from chemkit.core import charge_of                          # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    pre = args[0]
    T = load_tables()
    hit = [c for c in load_cases(None) if c["name"].split()[0] == pre]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    T_K = float(cond.get("T_K", 298.15))
    pr: dict = {}
    r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T,
                  _probe=pr)
    led = dict(pr.get("ledger") or {})
    He = pr.get("H_excess") or 0.0
    print(f"=== {c['name']} ===")
    print(f"  引擎呈现 pH = {r.get('final_pH')}   solver = {pr.get('pH_solver')}"
          f"   He = {He:+.8g}")
    print("  账本: " + "  ".join(f"{k}={v:.6g}" for k, v in sorted(led.items())
                                if k != WATER and v > 0))

    # ③ 直读：残余游离强酸池（精确记账）
    if He > 0.0:
        print(f"\n③ 直读 −log10(He/V) = {max(-1.0, -math.log10(He / V)):.4f}"
              f"   （游离强酸池，精确记账）")
    elif He < 0.0:
        print(f"\n③ 直读 pKw + log10(−He/V) = "
              f"{spec.pKw_of(T_K) + math.log10(-He / V):.4f}")

    # ④ 手算式：账本电荷平衡余项
    z = sum(charge_of(sp) * m for sp, m in led.items() if sp != WATER)
    print(f"④ Σz·n + He = {z + He:+.6g}   "
          f"（|·| ≈ 0 ⟹ 该 He 与账本物种**自洽**）")

    # ① estimate_pH + PH_TAGS
    spec.PH_TAGS = []
    try:
        p_est = spec.estimate_pH(dict(led), He, V, T, T_K)
        tags = list(spec.PH_TAGS)
    finally:
        spec.PH_TAGS = None
    print(f"\n① estimate_pH = {p_est:.6f}   tags = {tags}")

    # ①b 「直读早退」闸：哪些在账物种被判为"角色物种"（决定 _role_free）
    ru = spec._role_union(T)
    hits = [(sp, m) for sp, m in sorted(led.items(), key=lambda kv: -kv[1])
            if sp in ru]
    print(f"①b _role_union 命中 {len(hits)} 个在账物种"
          f"（阈值 ANN_MIN_EXTENT={spec.ANN_MIN_EXTENT:g} mol）：")
    for sp, m in hits:
        print(f"     {sp:<16}{m:>12.6g}  "
              f"{'**显著 ⟹ _role_free=False**' if m > spec.ANN_MIN_EXTENT else '痕量'}")

    # ①c 滴定后残余：直读早退用的是 `He_res / V`，不是传入的 He
    tit, He_res, led_t = spec._buffer_titration(
        dict(led), He, V, T, spec.pKw_of(T_K), T_K=T_K)
    print(f"①c _buffer_titration: tit={tit}   He_res={He_res:.8g}   "
          f"He_res/V={He_res / V:.6g}")
    if led_t is not led:
        ch = {k: round(v, 8) for k, v in led_t.items()
              if abs(v - led.get(k, 0.0)) > 1e-12 and k != WATER}
        print(f"     虚拟账本变化: {ch}")

    # ② charge_pH 两档
    for pin in (None, True):
        try:
            p_c = charge_pH(dict(led), He, V, T, T_K, pinned=pin)
        except Exception as exc:                            # noqa: BLE001
            p_c = f"EXC {exc!r}"
        print(f"② charge_pH(pinned={pin}) = {p_c}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
