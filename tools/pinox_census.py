# -*- coding: utf-8 -*-
"""第 290 轮 · **"固相钉住含氧酸根"前提的全库触发面普查**（零侵入）。

目的：在给 `estimate_state` 加"固相 + 其含氧酸根同时在账 ⟹ 酸根游离量被
`[A] = β_k·Ksp·oh^(k−z)` 钉住、与电荷平衡联立单未知量求解"这一分支之前，
先量出**这个前提在全部估计调用里命中哪些用例、多少次**（爆破半径枚举，
lessons §2"改闸前先枚举新放行集合"）。

前提（与拟议实现同一口径）：
  · 含氧酸根：beta 表 OH⁻ 配离子、m=1、净负电、ν > y/x（其中心的 OH⁻ 型
    Ksp 固相的 y/x），账本量 > ANN_MIN_EXTENT；
  · 该中心的 Ksp 固相同时在账，量 > ANN_MIN_EXTENT。

做法：`estimate_state` 以模块全局名调用 `_buffer_titration`，把
`speciation._buffer_titration` 包一层即可拿到**滴定后的虚拟账本**
（前提正是在它上面评估），不改 chemkit/ 一行。

输出 `logs/pinox_census.json`：
  * `pairs`     —— (solid, oxyanion) 命中次数直方图
  * `per_case`  —— 逐例：命中的 (solid, oxyanion, 次数, 最大酸根量, 最大固相量)
  * `n_calls`   —— `_buffer_titration` 总调用次数

⚠️ 串行运行（分钟级）。运行期间**冻结 chemkit/**。
"""
from __future__ import annotations

import collections
import io
import json
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

import chemkit.speciation as _sp                  # noqa: E402
from chemkit.candidates import ANN_MIN_EXTENT     # noqa: E402
from chemkit.core import charge_of                # noqa: E402
from chemkit.data import load_tables              # noqa: E402

THR = ANN_MIN_EXTENT    # 显著性用引擎同一把尺子（mol）

CUR_CASE = ["<init>"]
PAIRS = collections.Counter()
PER_CASE = collections.defaultdict(collections.Counter)
N_CALLS = [0]


def _pin_ox_table(T) -> dict:
    """含氧酸根 -> (solid, n)（与拟议实现同一推导；此处只需识别，不算系数）。"""
    ksp_oh = {}
    for e in T.ksp:
        c0, a0 = e["pair"]
        if a0 == "OH^-" and c0 not in ksp_oh:
            ksp_oh[c0] = e
    out = {}
    for b in T.beta:
        if b["ligand"] != "OH^-" or b.get("m", 1) != 1:
            continue
        cx = b.get("complex")
        if not cx or cx in T.solids or charge_of(cx) >= 0:
            continue
        e = ksp_oh.get(b["center"])
        if e is None:
            continue
        x, y = _sp._ksp_xy(e)
        n = b.get("nu", 1) - y / x
        if n <= 0:
            continue
        out[cx] = (e["solid"], n)
    return out


def main() -> int:
    T = load_tables()
    pin_ox = _pin_ox_table(T)
    print(f"含氧酸根钉住表（{len(pin_ox)} 条）：")
    for cx, (solid, n) in sorted(pin_ox.items()):
        print(f"  {cx:24s} -> {solid:16s} n={n:g}")

    orig = _sp._buffer_titration

    def wrapped(ledger, H_excess, V, T_, pKw, **kw):
        N_CALLS[0] += 1
        tit, he, led = orig(ledger, H_excess, V, T_, pKw, **kw)
        for cx, (solid, _n) in pin_ox.items():
            m_cx = led.get(cx, 0.0)
            if m_cx <= THR:
                continue
            m_s = led.get(solid, 0.0)
            if m_s <= THR:
                continue
            PAIRS[(solid, cx)] += 1
            PER_CASE[CUR_CASE[0]][(solid, cx)] += 1
        return tit, he, led

    _sp._buffer_titration = wrapped
    try:
        with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                     encoding="utf-8") as fh:
            cases = json.load(fh)
        from chemkit.engine import judge
        for c in cases:
            CUR_CASE[0] = c["name"].split()[0]
            try:
                judge([{"name": n, "mol": m} for n, m in c["subs"]],
                      c.get("cond") or {}, T)
            except Exception as exc:                    # noqa: BLE001
                print(f"[warn] {CUR_CASE[0]}: {exc!r}")
    finally:
        _sp._buffer_titration = orig

    out = {
        "n_calls": N_CALLS[0],
        "thr_mol": THR,
        "pairs": {f"{s} | {cx}": n for (s, cx), n in PAIRS.most_common()},
        "per_case": {case: {f"{s} | {cx}": n for (s, cx), n in cnt.most_common()}
                     for case, cnt in sorted(PER_CASE.items())},
    }
    dst = os.path.join(ROOT, "logs", "pinox_census.json")
    with io.open(dst, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(f"\n总调用 {N_CALLS[0]} 次；命中用例 {len(PER_CASE)} 个：")
    for case, cnt in sorted(PER_CASE.items()):
        tot = sum(cnt.values())
        pairs = ", ".join(f"{s}|{cx}×{n}" for (s, cx), n in cnt.most_common())
        print(f"  {case:10s} {tot:8d} 次   {pairs}")
    print(f"已写 {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
