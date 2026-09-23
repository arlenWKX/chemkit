# -*- coding: utf-8 -*-
"""第 266 轮 · **锋利的"口径不一致"判据**（第 265 轮 `cliff_census.py` 的正名版）。

## 为什么换判据

第 265 轮用「单次 `solve_extent` 探针内 pH 的 `max − min`」普查，得 36.8%。
**那是上界不是跳变率**：探针的 x 区间本可从 0 到 `x_max`（可达 ~1 mol），
宽区间本来就会给出宽 pH 跨度。

第 265 轮把 `I31` 查到底后拿到了**真正的机制**：
同一次迭代里，**候选评估**用主循环的 pH（`pin_mode=None`），
**求解器**用两端冻结后的 `pin_mode`（可能是 `False`）——两个口径在 `x=0` 处
就给出不同的 pH（I31：6.168 vs 8.984），于是"挑得起走不动"（lessons **D2**）。

## 本工具的判据（**精确**，不是跨度上界）

对每一次 `solve_extent` 调用：

1. 记下**进入前**主循环最后一次 `estimate_pH` 的返回值 = `pH_eval`
   （那就是候选评估与残差探针用的 pH）；
2. 记录 `solve_extent` 内部的全部 `estimate_pH` 调用。前两次是两端冻结判定
   （`pin_mode` 未传），**第三次起的 `pin_mode` 就是 `_pin_mode`**（走步口径）；
3. 用**同一口径**在 `x=0` 重算 `pH_solver0 = estimate_pH(led, He, ..., pin_mode=_pin_mode)`；
4. **口径差** `Δ = |pH_eval − pH_solver0|`。

`Δ` 大 = 引擎"挑的时候看到的状态"与"走的时候用的状态"不是同一个 ⟹
lessons D2 的**直接**量化。**这是准确值，不是上界。**

**内置正对照**（判据必须自校准，第 257/261 轮教训）：
当 `_pin_mode is None`（两端判定相同口径）时，`pH_solver0` 应当**逐位等于**
段内第一次 `estimate_pH` 的返回值（同一 led/He/口径）。工具报告该吻合率；
**吻合率不是 100% 就说明探针本身有问题，结论作废**。

零侵入（只包 `engine` 的模块级函数），不改 `chemkit/` 一行。
串行跑（~2 min），运行期间**冻结 chemkit/**。
输出 `logs/cliff_census2.json`。

用法：python tools/cliff_census.py            # 全库
      python tools/cliff_census.py I31 N15    # 只看这些前缀
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

CUR: list = []          # [(he, pH, pin_mode)] 本次 solve_extent 内
ON = [False]
LAST = {"pH": None}     # 主循环（求解之外）最后一次 estimate_pH 的返回值
CALLS = [0]
MIS = 0                 # 口径差 > 0.05 的次数
MIS0 = 0                # 且 x* ≈ 0
CAL_OK = 0              # 正对照：pm is None 时 ph0 与段内首点吻合
CAL_N = 0
PER = collections.Counter()
PER0 = collections.Counter()
CASE = [""]
LENS: list = []
ERRS: list = []
MODEL_N = [0]
MODEL_OK = [0]
PRED_N = [0]
MODEL_BAD: list = []


def _est(led, he, V, T, T_K, *a, **k):
    r = _ORIG_EST(led, he, V, T, T_K, *a, **k)
    if ON[0]:
        CUR.append((he, r, k.get("pin_mode", "unset")))
    else:
        LAST["pH"] = r
    return r


def _adopted(led, he, V, T, T_K):
    """复刻 `solve_extent` 的 `_pin_adopted`：该 (ledger, He) 上会不会采用 pinned。

    判据**读回引擎自己的 `PH_TAGS` 标注**（D2：不复制判据）。
    """
    import chemkit.speciation as _sp
    save = _sp.PH_TAGS
    _sp.PH_TAGS = []
    try:
        _ORIG_EST(led, he, V, T, T_K)
        tags = _sp.PH_TAGS
    finally:
        _sp.PH_TAGS = save
    return bool(tags) and tags[-1] == "电荷平衡精确解(pinned)"


def _se(cd, d, led, He, V, T_K, T, *a, **k):
    """⚠️ `solve_extent` 的签名是 `(c, d, ledger, He, V, T_K, T, …)` ——
    **`T_K` 在 `T` 之前**（与 `estimate_pH(ledger, He, V, T, T_K, …)` 相反）。
    第 266 轮初版按 `estimate_pH` 的顺序写形参名，导致回算时把 `Tables`
    当成了 `T_K` ⟹ `TypeError: unhashable type: 'Tables'`，
    且异常被外层 `except` 吞掉后表现为"一次都不命中"。
    **教训：包一个函数前先把它的签名抄对，别按"另一个同族函数"的顺序猜。**
    """
    global MIS, MIS0, CAL_OK, CAL_N
    ph_eval = LAST["pH"]
    CUR.clear()
    ON[0] = True
    try:
        out = _ORIG_SE(cd, d, led, He, V, T_K, T, *a, **k)
    finally:
        ON[0] = False
    CALLS[0] += 1
    LENS.append(len(CUR))
    if len(CUR) < 3:
        return out
    pm = None
    for _he, _ph, _p in CUR:
        if _p != "unset":
            pm = _p
            break
    try:
        ph0 = _ORIG_EST(led, He, V, T, T_K, pin_mode=pm)
    except Exception as _exc:                               # noqa: BLE001
        if not ERRS:
            ERRS.append(f"{type(_exc).__name__}: {_exc}  pm={pm!r}")
        return out
    if pm is None:
        CAL_N += 1
        if abs(ph0 - CUR[0][1]) <= 1e-9:
            CAL_OK += 1
    if ph_eval is None:
        return out
    # **机制模型的验证**（第 266 轮）：预测 `_pin_mode == entry and hi`。
    #   模型若对，则"口径不一致" ⟺ `entry=True 且 hi=False`
    #   （entry=False ⟹ pin_mode=False ⟹ x=0 处两条口径都跳过 pinned ⟹ 相同；
    #    entry=True,hi=True ⟹ pin_mode=True ⟹ x=0 处都是 pinned ⟹ 相同）。
    # 这不只是解释：**它同时是一个对引擎行为的可检验断言**。
    try:
        entry = _adopted(led, He, V, T, T_K)
        rr = cd.r if d > 0 else cd.pr
        pp = cd.pr if d > 0 else cd.r
        from chemkit.candidates import H_ION as _H, WATER as _W
        lim = [led.get(s, 0.0) / nu for s, nu in rr.items()
               if s not in (_W, _H) and nu > 0]
        xm = min(lim) if lim else 0.0
        hi = False
        if xm > 0.0:
            led_hi = dict(led)
            for s in set(rr) | set(pp):
                if s in (_W, _H):
                    continue
                led_hi[s] = max(0.0, led.get(s, 0.0)
                                + (pp.get(s, 0) - rr.get(s, 0)) * xm)
            nu_H = pp.get(_H, 0) - rr.get(_H, 0)
            hi = _adopted(led_hi, He + nu_H * xm, V, T, T_K)
        MODEL_N[0] += 1
        if (entry and hi) == bool(pm):
            MODEL_OK[0] += 1
        else:
            MODEL_BAD.append((CASE[0], entry, hi, pm))
        if entry and not hi:
            PRED_N[0] += 1
    except Exception:                                       # noqa: BLE001
        pass
    delta = abs(ph_eval - ph0)
    if delta > 0.05:
        MIS += 1
        PER[CASE[0]] += 1
        if out[0] <= max(1e-9, 1e-4 * out[1]):
            MIS0 += 1
            PER0[CASE[0]] += 1
    return out


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
        b, b0 = PER[nm], PER0[nm]
        try:
            eng.judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T)
        except Exception:                                   # noqa: BLE001
            pass
        rows.append({"name": nm, "n_mis": PER[nm] - b, "n_mis0": PER0[nm] - b0})
        if i % 150 == 0:
            print(f"  {i}/{len(cases)} … {time.time() - t0:.0f}s", flush=True)
    eng.estimate_pH = _ORIG_EST
    eng.solve_extent = _ORIG_SE
    hot = sorted((r for r in rows if r["n_mis"]), key=lambda r: -r["n_mis"])
    hot0 = sorted((r for r in rows if r["n_mis0"]), key=lambda r: -r["n_mis0"])
    out = {"solve_extent_calls": CALLS[0], "mismatch_gt05": MIS,
           "mismatch_gt05_zero_ext": MIS0,
           "rate": MIS / max(1, CALLS[0]),
           "n_cases_mis": len(hot), "n_cases_mis0": len(hot0),
           "n_cases": len(rows),
           "calib_ok": CAL_OK, "calib_n": CAL_N,
           "hot": hot, "hot0": hot0, "wall_s": round(time.time() - t0, 1)}
    with io.open(os.path.join(ROOT, "logs", "cliff_census2.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(f"\n=== 判据自校准（正对照）===")
    print(f"  段长分布（前 8）：{collections.Counter(LENS).most_common(8)}")
    print(f"  探针异常：{ERRS[:1]}")
    print(f"\n=== 机制模型验证：`_pin_mode == (_pin_entry and _pin_hi)` ===")
    print(f"  可验证段 {MODEL_N[0]}；吻合 **{MODEL_OK[0]}** "
          f"= {(MODEL_OK[0] / MODEL_N[0] * 100 if MODEL_N[0] else 0):.3f}%")
    print(f"  其中 `entry and not hi` 的段：**{PRED_N[0]}**"
          f"（模型预测**只有**这些段会出现口径不一致）")
    if MODEL_BAD:
        print(f"  ⚠️ 反例 {len(MODEL_BAD)} 条，前 3：{MODEL_BAD[:3]}")
    print(f"  `_pin_mode is None` 的段：{CAL_N} 次；"
          f"`pH_solver0` 与段内首点**逐位吻合** {CAL_OK} 次 "
          f"= {(CAL_OK / CAL_N * 100 if CAL_N else 0):.2f}%")
    if CAL_N and CAL_OK != CAL_N:
        print("  ⚠️ **吻合率不是 100% ⟹ 探针口径有偏差，下面的数字不可用**")
    print(f"\n=== 口径差（评估 pH vs 求解 pH@x=0）> 0.05 ===")
    print(f"  `solve_extent` 调用 {CALLS[0]} 次")
    print(f"  口径不一致：**{MIS}** 次 = **{out['rate'] * 100:.3f}%**"
          f"（第 265 轮的跨度上界是 36.756%）")
    print(f"  其中以 `x*≈0` 收场（真「挑得起走不动」）：**{MIS0}** 次")
    print(f"  落在 {out['n_cases_mis']} / {out['n_cases']} 例上"
          f"（其中零推进 {out['n_cases_mis0']} 例）")
    print("\n口径差最多的 25 例：")
    for r in hot[:25]:
        print(f"   {r['name'][:44]:<44} 不一致 {r['n_mis']:>5}"
              f"（零推进 {r['n_mis0']}）")
    print("\n-> logs/cliff_census2.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
