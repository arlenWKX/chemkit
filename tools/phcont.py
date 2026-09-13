"""原型实验：把"固相储库"从**分支选择**改成**同一个电荷平衡方程的一项**，
看 pH(x) 的跳变是否消失（§7 X 的验收判据），以及断言是否仍全绿。

**背景**（§7 X-3/X-4，已实测证实）：`estimate_pH` 是四分支启发式，x 推进
到固相出现/消失的点时**换分支**，pH 在一个网格步内跳 3.7–9.1 个单位
（全库 23.2% 的求根括号跳变 >1 单位）。化学上 pH(x) 是连续的滴定曲线，
跳变纯属模型假象；它让 S(x) 长出"−/+ 口袋"，既挡住更快的求根法
（Illinois 2.04× 被判否），也解释了走步的螺旋式慢收敛（§7 T）。

**本实验的做法**（不改生产路径，只做猴补 + 度量）：
  1. 找出账本里**在场的氢氧化物固相** `M_x(OH)_y`（量 > X_MIN）；
  2. 该固相在场时，自由阳离子浓度被 Ksp 钉住：
     `[M] = 10^((y·(pKw−pH) − pKsp)/x)`——作为 `charge_pH(pinned=...)` 的
     一项，同时从固定电荷里删掉账本那条自由离子量（避免双重记账）；
  3. 其余一切交给**同一个电荷平衡方程**（族分布 + 水 + 固定离子），
     不再问"哪个分支"。
  4. 无固相在场、或无根时**原样交回** estimate_pH（对照不变）。

**判据**：① pH 跳变普查（应显著下降）；② 全量断言（化学不能变坏）。

用法：
    python tools/phcont.py              # 全量：断言 + 普查
    python tools/phcont.py -e 40        # 普查抽样间隔
    python tools/phcont.py --orig       # 同一把尺子量基线
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from chemkit import engine                           # noqa: E402
from chemkit import testsuit as ts                   # noqa: E402
from chemkit.acidbase import charge_pH               # noqa: E402
from chemkit.candidates import X_MIN, _ksp_xy        # noqa: E402
from chemkit.acidbase import ledger_charge as _ledger_charge   # noqa: E402
from chemkit.core import charge_of, pKw_of           # noqa: E402
from chemkit.data import load_tables                 # noqa: E402
from chemkit.speciation import _pksp                 # noqa: E402
from roots import census                             # noqa: E402

_ORIG = engine.estimate_pH
_STAT = {"hit": 0, "fall": 0, "none": 0}
_PIN: dict = {}


def _pins_for(T) -> dict:
    """{固相名: (z_cat, ksp条目, x, y, cat)}——只取阴离子为 OH⁻ 的条目。"""
    got = _PIN.get(id(T))
    if got is None:
        got = {}
        for e in T.ksp:
            cat, an = e["pair"]
            if an != "OH^-":
                continue
            x, y = _ksp_xy(e)
            got[e["solid"]] = (charge_of(cat), e, x, y, cat)
        _PIN[id(T)] = got
    return got


def _ph_pinned(ledger, H_excess, V, T, T_K, bt_cache=None, touch=None):
    """estimate_pH 的替身：固相在场 → 钉住式电荷平衡；否则原样交回。"""
    pins = []
    led2 = None
    for sp, m in ledger.items():
        if m <= X_MIN or sp.startswith("__"):
            continue
        info = _pins_for(T).get(sp)
        if info is None:
            continue
        z, e, x, y, cat = info
        pins.append((z, _pksp(e, T_K), x, y))
        if led2 is None:
            led2 = dict(ledger)
        led2.pop(cat, None)          # 钉住值取代账本条目（避免双重记账）
    if not pins:
        _STAT["fall"] += 1
        return _ORIG(ledger, H_excess, V, T, T_K, bt_cache, touch)
    ph = charge_pH(led2, V, T, T_K, tol=1e-9, fast=True, pinned=tuple(pins))
    if ph is None:
        _STAT["none"] += 1
        return _ORIG(ledger, H_excess, V, T, T_K, bt_cache, touch)
    _STAT["hit"] += 1
    return ph


def _ph_exact(ledger, H_excess, V, T, T_K, bt_cache=None, touch=None):
    """estimate_pH 的替身：**账本自洽就直接解电荷平衡**（统一方程），否则回退机器。

    前提是本轮实测的：97.4% 的探头态满足 `|Σz·n + He| ≤ 1e-6`（p50 5.6e-17），
    此时账本内部的电荷平衡就是真实约束（§7 F/N 的 L08 反例正是不满足的那 2.6%）。
    好处是 **pH 不再有分支阈值**（酸侧 max / OH⁻ 直读 / 两性中点 / 缓冲对
    全是同一个方程的近似），pH(x) 因此天然连续。
    """
    if abs(_ledger_charge(ledger) + H_excess) <= 1e-6:
        ph = charge_pH(ledger, V, T, T_K, tol=1e-9, fast=True)
        if ph is not None:
            _STAT["hit"] += 1
            return ph
        _STAT["none"] += 1
    _STAT["fall"] += 1
    return _ORIG(ledger, H_excess, V, T, T_K, bt_cache, touch)


def main() -> None:
    every = 40
    if "-e" in sys.argv:
        every = int(sys.argv[sys.argv.index("-e") + 1])
    want: tuple[str, ...] = ()
    if "--cases" in sys.argv:
        want = tuple(sys.argv[sys.argv.index("--cases") + 1:])
    patch = "--orig" not in sys.argv
    mode = ("统一方程（账本自洽即解电荷平衡）" if "--exact" in sys.argv
            else "钉住式电荷平衡（原型）")
    detail = bool(want)
    if patch:
        engine.estimate_pH = _ph_exact if "--exact" in sys.argv else _ph_pinned
    T = load_tables()
    cases = ts.load_cases(None)
    if want:
        cases = [c for c in cases if c["name"].startswith(want)]
    engine.ROOT_AUDIT = {"n": 0, "every": every, "rec": [], "case": None}
    ts.FAILS.clear()
    ts.RESULTS.clear()
    bad = []
    for c in cases:
        n0 = len(ts.RESULTS)
        engine.ROOT_AUDIT["case"] = c["name"]
        ts.run_case(c, T, verbose=False)
        if len(ts.RESULTS) > n0 and not ts.RESULTS[-1]["ok"]:
            bad.append(c["name"])
            if detail:
                print(f"\n## {c['name']}")
                for e in ts.RESULTS[-1]["errors"]:
                    print("   ", e)
    rec = engine.ROOT_AUDIT["rec"]
    engine.ROOT_AUDIT = None
    engine.estimate_pH = _ORIG
    print(f"\n== {'基线（原 pH 机器）' if not patch else mode} ==")
    print(f"断言：{len(cases) - len(bad)}/{len(cases)} 例通过（红 {len(bad)}）")
    if patch:
        print(f"钉住档命中 {_STAT['hit']} 次 / 回退原机器 {_STAT['fall']} 次"
              f" / 无根再回退 {_STAT['none']} 次")
        if bad:
            print("翻红用例：", bad[:20])
    census(rec, label="pH 连续性与换根普查")


if __name__ == "__main__":
    main()
