# -*- coding: utf-8 -*-
"""第 282 轮 · **同投料跨温度对比**（`TE1`–`TE4` 那条"非 298 K 才错"的线）。

## 为什么需要

第 281 轮把 `TE1`–`TE4`（`NH4Cl + 半量 NaOH`，273.15/323.15/348.15/373.15 K）
的残差收窄成"**温度相关**"：新增的 `HN1`/`HN2` 锚在 **298.15 K** 上
`resid = 0.0000`、pH 逐位等于 pKa，而同化学的非 298 K 四例残差**都是 1.878**。

要定位"温度到底改了哪一步"，最快的办法是**把同一投料在一串温度上跑一遍**，
并排看：残差 / 步数 / 每步的 `logK`、`S`、`extent`，以及 `He` 的终值。

判读提纲：
* 若 **298.15 K 干净、两侧立刻坏** ⟹ 某个按 `T_K` 缓存/短路的东西只在 298 命中；
* 若**全温度都坏** ⟹ 与温度无关，是体系本身的走步问题；
* 若 `logK` 随 T 变而 `S` 不随之归零 ⟹ `logK` 与 `pH` 两条温度通道口径不一致。

用法：
    python tools/tcmp.py TE1                       # 用该例的投料、扫温度
    python tools/tcmp.py TE1 --t=273.15,298.15,323.15
    python tools/tcmp.py TE1 --full                # 打印完整步表
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
from chemkit.converg import _live                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    pre = args[0]
    temps = [273.15, 278.15, 288.15, 298.15, 308.15, 323.15, 348.15, 373.15]
    for a in argv:
        if a.startswith("--t="):
            temps = [float(x) for x in a[4:].split(",")]
    # **浓度扫描**（第 282 轮加）：TE1(0.01/0.005) 坏而 HN1(0.02/0.01)、
    # HN2(0.002/0.001) 好 —— 非单调 ⟹ 不是"稀/浓"而是某个**阈值/分支**。
    scale = None
    for a in argv:
        if a.startswith("--scale="):
            scale = [float(x) for x in a[8:].split(",")]
    full = "--full" in argv
    T = load_tables()
    hit = [c for c in load_cases(None) if c["name"].split()[0] == pre]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    print(f"=== {c['name']} ===")
    print(f"  投料 {c['subs']}（cond 会被本工具的 --t 覆盖）")
    if scale:
        print("\n=== 浓度扫描（投料 × scale，298.15 K 若未另给 --t）===")
        tk0 = temps[0] if len(temps) == 1 else 298.15
        print(f"{'scale':>9}{'c(NH4Cl)':>12}{'pH':>9}{'He':>12}{'resid':>9}"
              f"{'步数':>6}")
        for sc in scale:
            cond = dict(c.get("cond") or {})
            cond["T_K"] = tk0
            pr2: dict = {}
            r2 = eng.judge([{"name": s[0], "mol": float(s[1]) * sc}
                            for s in c["subs"]], cond, T, _probe=pr2)
            st2 = r2.get("steps") or []
            print(f"{sc:>9g}{0.01 * sc:>12.5g}"
                  f"{(r2.get('final_pH') or 0):>9.3f}"
                  f"{(pr2.get('H_excess') or 0):>12.6g}"
                  f"{abs(_live(pr2.get('active') or [])):>9.4f}{len(st2):>6}")
        return 0
    print(f"\n{'T_K':>9}{'pH':>9}{'He':>12}{'resid':>9}{'iters':>7}"
          f"{'步数':>6}  {'最慢/最怪的一步'}")
    for tk in temps:
        cond = dict(c.get("cond") or {})
        cond["T_K"] = tk
        pr: dict = {}
        r = eng.judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      cond, T, _probe=pr)
        steps = r.get("steps") or []
        rv = abs(_live(pr.get("active") or []))
        top = None
        for a in (pr.get("active") or []):
            if not _live([a]):
                continue
            if top is None or abs(a.get("S", 0)) > abs(top.get("S", 0)):
                top = a
        note = ""
        if top is not None:
            note = (f"|S|={abs(top['S']):.3f} {top.get('kind')} "
                    f"{str(top.get('eq'))[:44]}")
        print(f"{tk:>9.2f}{r.get('final_pH') or 0:>9.3f}"
              f"{pr.get('H_excess') or 0:>12.6g}{rv:>9.4f}"
              f"{pr.get('iters') or 0:>7}{len(steps):>6}  {note}")
        if full:
            for s in steps:
                print(f"      [{(s.get('kind') or ''):<9}] "
                      f"x={s.get('extent'):<12.6g} "
                      f"logK={s.get('logK'):<8.4g} S={s.get('S'):<8.4g} "
                      f"{s.get('equation')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
