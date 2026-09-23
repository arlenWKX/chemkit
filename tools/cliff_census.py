# -*- coding: utf-8 -*-
"""第 265 轮 · **单次 `solve_extent` 探针内的 pH 跳变普查**（给 §1.1 决策用的证据）。

## 为什么量这个

第 265 轮把 `I31` 的残差查到根：**同一次迭代里有两个 pH 口径**。
`tools/solve_trace.py` 实测（最后一次匹配探针）：

    x=0（两端冻结判定）   pH 6.1680   分支 `电荷平衡精确解(pinned)`
    x=x_max（判定）       pH 11.3174  分支 `碱侧max`
    ⟹ `_pin_mode = entry and hi = False`（整步弃用 pinned）
    其后所有 `f(x)` 求值      pH 8.98→11.32  分支 `碱侧max`

即：候选评估/残差探针报的是 **6.168（S=+2.816，看着有驱动）**，
而 `solve_extent` 自己求根用的是 **8.984（S≈0，本就在平衡）** ⟹ `x*≈4e-10`
⟹ "挑得起走不动"（lessons **D2** 口径不一致）。

修法试过一条并**被实测否掉**（见 log 第 265 轮）：给 pinned 解加"质子不变量回检"
⟹ `I31` 残差 2.816 → **15.80**（更差：走步停在铟相对 In(OH)₃ **过饱和**的态），
因为 `charge_pH` 的 docstring 本就写明 **pinned 是为"账本溶解量本身与电荷平衡
不自洽"而生**（L08 型）。⟹ 这是 §1.1 表示层族的成员，不是局部可修的。

本工具量化这一族的**规模**：单次 `solve_extent` 内 pH 摆动 > 1 的次数与用例数，
以及其中有多少以 `x*≈0` 收场（= 零推进的真正来源）。

## 口径

* 包 `engine.estimate_pH`（模块全局名调用，补丁有效）记录每次返回的 pH；
  包 `engine.solve_extent` 按调用切分记录段。
* `span = max(pH) − min(pH)`；`x*` 取 `solve_extent` 的返回。
* **不改 chemkit 一行**。串行跑（~2 min），运行期间**冻结 chemkit/**。

用法：python tools/cliff_census.py            # 全库
      python tools/cliff_census.py N15 I31    # 只看这些前缀
输出 `logs/cliff_census.json`
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

import chemkit.engine as eng                                # noqa: E402
from chemkit.data import load_tables                        # noqa: E402

_ORIG_EST = eng.estimate_pH
_ORIG_SE = eng.solve_extent

CUR: list = []
ON = [False]
CALLS = [0]
BIG = [0]            # span > 1
BIG0 = [0]           # span > 1 且 x* ≈ 0（微步）
PER = collections.Counter()


def _est(led, he, V, T, T_K, *a, **k):
    r = _ORIG_EST(led, he, V, T, T_K, *a, **k)
    if ON[0]:
        CUR.append(r)
    return r


def _se(cd, d, led, He, V, T, T_K, *a, **k):
    CUR.clear()
    ON[0] = True
    try:
        out = _ORIG_SE(cd, d, led, He, V, T, T_K, *a, **k)
    finally:
        ON[0] = False
    CALLS[0] += 1
    if len(CUR) >= 2:
        span = max(CUR) - min(CUR)
        if span > 1.0:
            BIG[0] += 1
            PER[CASE[0]] += 1
            if out[0] <= max(1e-9, 1e-4 * out[1]):
                BIG0[0] += 1
    return out


CASE = [""]


def main(argv: list[str]) -> int:
    pres = tuple(a for a in argv if not a.startswith("-"))
    eng.estimate_pH = _est
    eng.solve_extent = _se
    T = load_tables()
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    t0 = time.time()
    rows = []
    for i, c in enumerate(cases, 1):
        nm = c["name"]
        if pres and not nm.startswith(pres):
            continue
        CASE[0] = nm
        before = PER[nm]
        try:
            eng.judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T)
        except Exception:                                       # noqa: BLE001
            pass
        rows.append({"name": nm, "n_big": PER[nm] - before})
        if i % 150 == 0:
            print(f"  {i}/{len(cases)} … {time.time() - t0:.0f}s", flush=True)
    eng.estimate_pH = _ORIG_EST
    eng.solve_extent = _ORIG_SE
    hot = sorted((r for r in rows if r["n_big"]), key=lambda r: -r["n_big"])
    out = {"solve_extent_calls": CALLS[0], "span_gt1": BIG[0],
           "span_gt1_zero_ext": BIG0[0],
           "rate": BIG[0] / max(1, CALLS[0]),
           "n_cases_big": len(hot), "n_cases": len(rows),
           "hot": hot, "wall_s": round(time.time() - t0, 1)}
    with io.open(os.path.join(ROOT, "logs", "cliff_census.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(f"\n`solve_extent` 调用 {CALLS[0]} 次")
    print(f"单次探针内 pH 摆动 > 1 的：**{BIG[0]}** 次 = "
          f"**{out['rate'] * 100:.3f}%**")
    print(f"其中以 `x*≈0` 收场（= 真零推进）：**{BIG0[0]}** 次")
    print(f"落在 {out['n_cases_big']} / {out['n_cases']} 个用例上；"
          f"串行墙钟 {out['wall_s']}s")
    print("\n摆动最多的 25 例：")
    for r in hot[:25]:
        print(f"   {r['name'][:44]:<44} {r['n_big']}")
    print("\n-> logs/cliff_census.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
