# -*- coding: utf-8 -*-
"""第 273 轮 · **共轭形态重建起火普查**：某例的残差变化，究竟由**哪条候选的
哪个物种**的重建引起。

## 为什么需要

第 272 轮落地"共轭形态重建"后，全库是净赚（残差质量 −17.5%、套件 +7），
但 `H35 Ca(HCO3)2 加热 353K` 从 **0.363 → 10.055**、`Ni41 NiCl2+NaHCO3`
从 0.859 → 11.964。全局指标只能说"与重建有关"，**说不出是哪一条把走步
点燃的** —— 而这两例的走步都是**极限环**（`dev.py case H35 --all` 里
同一组步骤逐字重复），要定位就必须知道每一步的 `S` 是被谁抬起来的。

本工具开 `engine.FORM_AUDIT` 跑一个用例，把**每次重建**按
`(候选式, 被重建物种)` 聚合，并给出该次重建把 `log c` 从
`ACT_FLOOR = 1e-12`（`floor_logc`）抬到了多少 —— 抬升量正比于它对该反应
`logQ` 项的贡献。

用法：
    python tools/formfire.py H35            # 单例聚合
    python tools/formfire.py Ni41 --top 15
    python tools/formfire.py --list         # 只列用例名
"""

from __future__ import annotations

import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
from chemkit.candidates import WATER                        # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def _txt(dd) -> str:
    return " + ".join(f"{'' if n == 1 else n}{s}"
                      for s, n in dd.items() if s != WATER)


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    top = 12
    for a in argv:
        if a.startswith("--top="):
            top = int(a[6:])
    T = load_tables()
    cases = load_cases(None)
    if not args or "--list" in argv:
        print("用例名（首 token）：")
        print("  " + "  ".join(sorted({c["name"].split()[0] for c in cases})))
        return 0
    pre = args[0]
    hit = [c for c in cases if c["name"].split()[0] == pre]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    eng.FORM_AUDIT = []
    try:
        r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                      c.get("cond") or {}, T)
    finally:
        fires = eng.FORM_AUDIT
        eng.FORM_AUDIT = None
    print(f"=== {c['name']} ===  pH={r.get('final_pH')}  "
          f"重建起火 **{len(fires)}** 次")

    agg: dict = defaultdict(lambda: [0, 0.0, 0.0])
    for f in fires:
        k = (f["eq"] and f["eq"].kind, f["s"], f["p"], round(f["base"]))
        g = agg[k]
        g[0] += 1
        g[1] = max(g[1], f["c"])
        g[2] = max(g[2], f["pH"])
    print(f"\n{'候选 kind':<10}{'被重建':<14}{'伙伴':<14}{'碱侧':>5}"
          f"{'次数':>6}{'最大 c':>12}{'最大 pH':>9}")
    for (kind, s, p, base), (n, cmax, ph) in sorted(
            agg.items(), key=lambda kv: -kv[1][0])[:top]:
        print(f"{str(kind):<10}{s:<14}{p:<14}{base:>5}{n:>6}"
              f"{cmax:>12.4g}{ph:>9.3f}")

    # 抬升量最大的若干次（floor=1e-12 ⟹ 抬升 = log10(c) + 12）
    print(f"\n单次抬升最大的 {min(top, len(fires))} 次"
          f"（相对 ACT_FLOOR = 1e-12）：")
    import math
    seen = set()
    rank = []
    for f in fires:
        import math as _m
        lift = _m.log10(max(f["c"], 1e-300)) + 12.0
        key = (f["eq"] and _txt(f["eq"].r) + " -> " + _txt(f["eq"].pr),
               f["s"])
        if key in seen:
            continue
        seen.add(key)
        rank.append((lift, key, f))
    rank.sort(key=lambda t: -t[0])
    for lift, (eq, s), f in rank[:top]:
        print(f"  ΔlogQ 项 {lift:>7.3f}  {s:<14} ← {f['p']:<14} "
              f"c={f['c']:.4g}  pH={f['pH']:.3f}")
        print(f"        {eq}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
