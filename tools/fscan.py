# -*- coding: utf-8 -*-
"""第 272 轮 · 直接扫**引擎自己的 `f(x)`**（求根审计句柄），不复制任何式子。

## 为什么必须用引擎的句柄

`tools/xscan_ch.py` 对 `BR2` 的通道 `10HBrO + Cl₂ -> 5Br₂ + 2H⁺ + 2ClO₃⁻`
扫出"S 从 +22.43 掉到 −10.71"（Δx 仅 `x_max` 的 7.7%），据此**看似**存在
过零点；同一终态上 `solve_extent` 实测 `ext = 2.3e-07`。两者矛盾。

第 272 轮查明矛盾来源：`solve_extent` 的 `f` 对 **redox 候选**走的是
另一条式子（engine.py L689-715）——

    pH_x, led_v, _ = estimate_state(led_work, He + nu_H·x, …, pin_mode=…)
    return direction * S_of(c, led_v, …)      # ← 用 led_v（滴定后虚拟账本）

而**非** redox 候选用 `estimate_pH` + `led_work`。我此前的复刻工具
（`estimate_pH` + `led2`）因此量的是**另一条曲线**，结论不可信。

本工具改用 `ROOT_AUDIT["rec"][i]["f"]`——`_audit_bracket` 把 `f` **闭包
本身**存进了审计记录（engine.py L430）⟹ 扫的就是二分用的那个函数，
口径不可能不一致（"探头与规则同一判据"的机制化，第 267 轮教训）。

用法：python tools/fscan.py BR2 10HBrO [npts]
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

import chemkit.engine as eng                                # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def _scan(f, x_max: float, npts: int) -> None:
    """对数+线性混合网格扫 f，打印并报出符号翻转区间。"""
    xs = []
    for k in range(npts):
        # 前一半对数（近 0 处 f 最陡），后一半线性
        if k < npts // 2:
            lo, hi = 1e-9, 1e-1
            t = k / max(1, npts // 2 - 1)
            xs.append(x_max * 10.0 ** (__import__("math").log10(lo)
                                       + t * (__import__("math").log10(hi)
                                              - __import__("math").log10(lo))))
        else:
            t = (k - npts // 2) / max(1, npts - npts // 2 - 1)
            xs.append(x_max * (0.1 + 0.9 * t))
    xs = sorted(set(xs))
    print(f"\n{'x':>14} {'x/x_max':>10} {'f(x)':>14}")
    prev = None
    flips = []
    for x in xs:
        try:
            v = f(x)
        except Exception as exc:                            # noqa: BLE001
            print(f"{x:>14.6g} {x / x_max:>10.4g}   EXC {exc!r}")
            continue
        print(f"{x:>14.6g} {x / x_max:>10.4g} {v:>14.6g}")
        if prev is not None and (prev[1] > 0) != (v > 0):
            flips.append((prev[0], x))
        prev = (x, v)
    print(f"\n符号翻转区间数 = {len(flips)}")
    for a, b in flips[:12]:
        print(f"   ({a:.8g}, {b:.8g})   f={f(a):+.6g} -> {f(b):+.6g}")


def _bisect(f, x_max: float, tol: float = 1e-11) -> float:
    lo, hi = 0.0, x_max
    if f(0.0) <= 0:
        return 0.0
    if f(hi) > 0:
        return hi
    while hi - lo > tol * max(1.0, x_max):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return lo


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    pre = args[0] if args else "BR2"
    needle = args[1] if len(args) > 1 else "10HBrO"
    npts = int(args[2]) if len(args) > 2 else 32
    every = 1
    for a in argv:
        if a.startswith("--every="):
            every = int(a[8:])
    T = load_tables()
    # 用例名形如 "BR2 Cl2 饱和+KBr 0.5 摩尔过量卤素归中锚点" ⟹ 比首 token
    # （`startswith("BR2")` 会让 H1 命中 H10 等同族坑）
    hit = [x for x in load_cases(None) if x["name"].split()[0] == pre]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    cond = c.get("cond") or {}

    eng.ROOT_AUDIT = {"n": 0, "every": every, "rec": [], "case": c["name"]}
    try:
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T)
    finally:
        rec = eng.ROOT_AUDIT["rec"]
        eng.ROOT_AUDIT = None
    print(f"=== {c['name']} ===  审计记录 {len(rec)} 条")

    hits = [r for r in rec if needle in r["eq"]]
    if not hits:
        seen = sorted({r["eq"] for r in rec})
        print(f"[未匹配 {needle!r}] 本用例候选共 {len(seen)} 条，含 "
              f"'HBrO'/'Cl' 的：")
        for e in seen:
            if "HBrO" in e or "Cl" in e:
                print("   ", e)
        return 1
    print(f"匹配 {len(hits)} 条：{hits[0]['eq']}   kind={hits[0]['kind']}")
    xs = [r["x_bis"] for r in hits]
    print(f"  x_bis 范围 [{min(xs):.6g}, {max(xs):.6g}]  x_max="
          f"{hits[0]['x_max']:.6g}  flips 取值 {sorted({r['flips'] for r in hits})}"
          f"  f0 取值 {sorted({round(r['f0'], 4) for r in hits})}")
    for r in hits[:6]:
        print(f"  dir={r['dir']:+d} x_max={r['x_max']:.6g} f0={r['f0']:+.6g} "
              f"x_bis={r['x_bis']:.8g} x_ill={r['x_ill']:.8g} "
              f"n_bis={r['n_bis']} flips={r['flips']} "
              f"ph_closed={r['ph_closed']} cls_ok={r['cls_ok']}")

    # 取"驱动为正方向"的那条（f0>0）做细扫——没有根的那一侧才需要解释
    pick = None
    for r in hits:
        if r["f0"] > 0:
            pick = r
            break
    if pick is None:
        print("\n[全部 f0<=0] 该候选在所有迭代里探头侧都不驱动 ⟹ 零推进"
              "由**口径**造成（评估用虚拟账本、求解器用真实账本），"
              "非求根器缺陷。")
        pick = hits[0]
    print(f"\n--- 细扫 dir={pick['dir']:+d} x_max={pick['x_max']:.8g} "
          f"x_bis={pick['x_bis']:.8g} ---")
    _scan(pick["f"], pick["x_max"], npts)
    xb = _bisect(pick["f"], pick["x_max"])
    print(f"\n本工具二分根 = {xb:.8g}   引擎 x_bis = {pick['x_bis']:.8g}   "
          f"引擎 x_ill = {pick['x_ill']:.8g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
