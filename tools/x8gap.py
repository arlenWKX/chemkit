# -*- coding: utf-8 -*-
"""第 272 轮 · **X-8 口径缺口直测**：同一退出态上，"残差口径"与"求解器口径"
看到的 S 差 22.4 个数量级，差在哪个账本、哪个物种。

## 背景（本轮实测）

`tools/tension.py --case BR2` 对 `BR2 Cl2 过量+KBr 0.5 半量（卤素计量锚）`：

    live=22.433  pH=-0.042  exit=no-cands  iters=25
      |S|=22.433 kind=redox ext_max=0.05
      10HBrO + Cl_2 -> 5Br_2 + 2H^+ + 2ClO_3^-
      → 已达平衡（求解时 |S|=0.0 ext=0 x_max=0.05）
      走步判据：x*=6.402e-11  x_max=0.05

即 **同一退出态**：残差口径报 22.433（远未平衡）、求解器口径报 ≈0（已平衡）。
而 `tools/fscan.py` 读引擎自己的 `f` 闭包证实：该步 `f(0)=+7.6e-05`、`x_bis=0`，
**求根器没有错**——它看到的曲线确实在 x=0 就归零了。

engine.py L704-712 的注释已把成因写明：redox 分支的 `f` 用 `led_v`
（`estimate_state` 滴定后的**虚拟账本**），而残差/走步的**候选评估**用真实账本：

    return direction * S_of(c, led_v, …)          # redox f
    return direction * S_of(c, led_work, …)       # 非 redox f

本工具把两个账本**逐物种摆出来**，并按反应的 Q 对数贡献排序——要修的是
"哪个账本对 redox 有权威"（§7 X-8），不是求根器。

用法：python tools/x8gap.py BR2 10HBrO
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
from chemkit.candidates import H_ION, WATER                 # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def _txt(dd: dict) -> str:
    return " + ".join(f"{'' if n == 1 else n}{s}"
                      for s, n in dd.items() if s != WATER)


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    pre = args[0] if args else "BR2"
    needle = args[1] if len(args) > 1 else "10HBrO"
    T = load_tables()
    hit = [x for x in load_cases(None) if x["name"].split()[0] == pre]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))
    T_K = float(cond.get("T_K", 298.15))
    pr: dict = {}
    eng.judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T, _probe=pr)
    led = dict(pr.get("ledger") or {})
    He = pr.get("H_excess") or 0.0
    pH0 = pr.get("pH") or 7.0
    print(f"=== {c['name']} ===")
    print(f"退出态：pH={pH0:.4f}  He={He:+.8g}  S_max={pr.get('max_abs_S')}")
    print("真实账本（mol）: " + "  ".join(
        f"{k}={v:.6g}" for k, v in sorted(led.items())
        if k != WATER and v > 0))

    # 定位候选：用主循环的枚举口径（真账本 + 退出态 pH）
    cands = eng.enumerate_candidates(led, He, pH0, V, T_K, T, True)
    pick = None
    for cd in cands:
        for d in (1, -1):
            rr = cd.r if d > 0 else cd.pr
            pp = cd.pr if d > 0 else cd.r
            if needle not in f"{_txt(rr)} -> {_txt(pp)}":
                continue
            S = eng.S_of(cd, led, V, pH0, T_K, T, frozenset(),
                         eng.P_EXT_KPA, True, None)
            pick = (cd, d, S, rr, pp)
    if pick is None:
        print("[未匹配候选]")
        return 1
    cd, d, S_raw, rr, pp = pick
    print(f"\n候选：{_txt(rr)} -> {_txt(pp)}   d={d:+d}")
    print(f"残差口径 S_of(真账本, pH={pH0:.4f})   = {S_raw:+.6f}")

    # 求解器口径：estimate_state 的 led_v
    ph_v, led_v, _ = spec.estimate_state(
        dict(led), He, V, T, T_K, {}, frozenset(), None, None, pin_mode=False)
    S_v = eng.S_of(cd, led_v, V, ph_v, T_K, T, frozenset(),
                   eng.P_EXT_KPA, True, None)
    print(f"求解器口径 S_of(led_v, pH={ph_v:.4f}) = {S_v:+.6f}"
          f"   （注：f 还乘了 direction，此处按写向比较大小）")
    print(f"口径差 ΔlogQ = {S_v - S_raw:+.4f}")

    print(f"\n{'物种':<16}{'真账本':>13}{'led_v':>13}{'Δ':>13}"
          f"{'ν(写向)':>9}{'Δ对logQ':>12}")
    keys = sorted((set(led) | set(led_v) | set(rr) | set(pp)) - {WATER})
    nu = {}
    for s in set(rr) | set(pp):
        nu[s] = float(pp.get(s, 0.0)) - float(rr.get(s, 0.0))
    rows = []
    for s in keys:
        a = led.get(s, 0.0)
        b = led_v.get(s, 0.0)
        n = nu.get(s, 0.0)
        contrib = n * (math.log10(max(b, 1e-300) / V)
                       - math.log10(max(a, 1e-300) / V)) if a > 0 and b > 0 else 0.0
        if abs(b - a) > 1e-12 or n:
            rows.append((abs(contrib), s, a, b, b - a, n, contrib))
    for _ab, s, a, b, dd, n, contrib in sorted(rows, reverse=True)[:24]:
        print(f"{s:<16}{a:>13.6g}{b:>13.6g}{dd:>13.6g}{n:>9.4g}{contrib:>12.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
