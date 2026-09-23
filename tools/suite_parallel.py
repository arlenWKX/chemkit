# -*- coding: utf-8 -*-
"""第 206 轮 · 并行跑套件（多进程）—— 回答"能否用多进程提高效率"。

为什么多进程可行（已实测）：
  · `Tables` 载入仅 **182 ms**、pickle 后 0.16 MB ⟹ 每个 worker 自己载一份
    的启动成本远小于单例耗时；
  · 本机 `cpu_count = 8`，且**只有 `spawn`** 启动方式（无 fork）⟹
    worker 不继承父进程状态，**必须自己载表**（这反而是干净的）；
  · 用例之间**互相独立**（各自 `judge()` 自己的 subs/cond）。

为什么不改 `testsuit.run_case`：它写模块级全局（`PASS_N`/`FAILS`/`RESULTS`），
改它有回归风险。本脚本用**进程隔离**绕开：每个 worker 只回传
`(name, ok, errors, ms)`，父进程汇总。

用法：
  python tools/suite_parallel.py            # 自动选 worker 数
  python tools/suite_parallel.py 6          # 指定 worker 数
  python tools/suite_parallel.py 6 --check  # 与串行结果逐例对拍（慢但严谨）
"""
import io
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import multiprocessing as mp                               # noqa: E402
from concurrent.futures import ProcessPoolExecutor         # noqa: E402

_T = None                       # worker 内缓存的 Tables


def _worker_init():
    """每个 worker 进程启动时载一次表（spawn 下无继承，必须自己载），
    并做一次**预热**。

    ⚠️ 第 206 轮实测教训：不预热时，每个 worker 的**首个用例**要独自承担
    全部惰性缓存构建（`_redox_pair_static` / `logk_poly` / `build_derived`…），
    于是套件里位置最靠前的简单用例（`1 HCl+NaOH`、`3 HCl+NaOH 碱过量`…）
    显示成 12~24 s —— 看着像算法病态，其实**纯冷启动**：同一用例在已预热的
    进程里只要 **9~270 ms**。`testsuit.main` 本就有预热，并行器必须照做，
    否则 8 个 worker 各付一次，还污染"最慢用例"榜。
    """
    global _T
    from chemkit.data import load_tables
    from chemkit.engine import judge
    _T = load_tables()
    judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, _T)   # 预热


def _worker_run(case):
    """跑单个用例，返回 (name, ok, errors, ms)。"""
    from chemkit import testsuit as ts
    name = case["name"]
    # 清掉该 worker 内的全局残余，只取本次追加的失败项
    n0 = len(ts.FAILS)
    _t0 = time.perf_counter()
    try:
        ok = ts.run_case(case, _T, verbose=False)
        errs = list(ts.FAILS[n0:])
    except Exception as exc:                                # noqa: BLE001
        ok = False
        errs = [f"运行异常 {type(exc).__name__}: {exc}"]
    del ts.FAILS[n0:]
    return (name, bool(ok), errs, round((time.perf_counter() - _t0) * 1000, 1))


def run_parallel(n_workers, cases, chunksize=1):
    """chunksize 取 1：单例耗时差异极大（有 100s+ 的慢例），块太大会造成
    **尾部长尾**（实测 chunksize=4 时 800→1000 例期间墙钟几乎不前进）。"""
    t0 = time.perf_counter()
    out = []
    with ProcessPoolExecutor(max_workers=n_workers,
                             initializer=_worker_init) as ex:
        for i, res in enumerate(ex.map(_worker_run, cases,
                                       chunksize=chunksize), 1):
            out.append(res)
            if i % 200 == 0:
                el = time.perf_counter() - t0
                print(f"  {i}/{len(cases)}  elapsed {el:.1f}s", flush=True)
    return out, time.perf_counter() - t0


def run_serial(cases):
    from chemkit.data import load_tables
    from chemkit import testsuit as ts
    T = load_tables()
    t0 = time.perf_counter()
    out = []
    for c in cases:
        n0 = len(ts.FAILS)
        try:
            ok = ts.run_case(c, T, verbose=False)
            errs = list(ts.FAILS[n0:])
        except Exception as exc:                            # noqa: BLE001
            ok, errs = False, [f"{type(exc).__name__}: {exc}"]
        del ts.FAILS[n0:]
        out.append((c["name"], bool(ok), errs, None))
    return out, time.perf_counter() - t0


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    do_check = "--check" in sys.argv
    n = int(args[0]) if args else max(1, min(8, (os.cpu_count() or 4)))
    from chemkit.testsuit import load_cases
    cases = load_cases(None)
    print(f"用例 {len(cases)} 个；worker 数 = {n}（cpu_count="
          f"{os.cpu_count()}，spawn 模式）")

    res_p, t_p = run_parallel(n, cases)
    n_ok = sum(1 for _nm, ok, _e, _m in res_p if ok)
    fails = [nm for nm, ok, _e, _m in res_p if not ok]
    print(f"\n[并行 {n} workers] 用例 {n_ok}/{len(cases)}  墙钟 {t_p:.1f}s")

    # 逐例耗时 → 解释"为什么并行只有 ~2x"（长尾决定墙钟上限）
    timed = sorted(((m or 0.0, nm) for nm, _ok, _e, m in res_p), reverse=True)
    tot = sum(m for m, _nm in timed)
    if timed:
        print(f"\n串行 CPU 合计 {tot / 1000:.1f}s；"
              f"理论下限 = max(最慢单例, CPU合计/worker数)")
        print(f"最慢 8 例:")
        for m, nm in timed[:8]:
            print(f"  {m / 1000:7.2f}s  {nm}")
        print(f"最慢单例 = {timed[0][0] / 1000:.2f}s ⟹ 并行墙钟**不可能**低于它")

    # ===== 一次运行写全部产物：**不得为看另一面而重跑** =====
    import json
    import datetime
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    payload = {
        "n_workers": n, "wall_s": round(t_p, 2),
        "n_cases": len(cases), "n_ok": n_ok, "n_fail": len(fails),
        "cpu_total_s": round(tot / 1000, 1),
        "slowest_case_s": round(timed[0][0] / 1000, 2) if timed else None,
        "fails": fails,
        "cases": [{"name": nm, "ok": ok, "errors": e, "ms": m}
                  for nm, ok, e, m in res_p],
        "slowest": [{"name": nm, "ms": m} for m, nm in timed[:40]],
    }
    jp = os.path.join(ROOT, "logs", f"suite-parallel-{stamp}.json")
    with io.open(jp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=1)
    lp = os.path.join(ROOT, "logs", "suite-parallel-latest.json")
    with io.open(lp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=1)
    print(f"\n[留档] {os.path.relpath(jp, ROOT)}")
    print(f"[留档] {os.path.relpath(lp, ROOT)}（机读，含逐例 ok/errors/ms）")
    print("  ⟹ 后续要看慢例/失败清单/耗时，**读这个文件**，不要重跑"
          "（`python tools/suite_show.py`）")

    if do_check:
        res_s, t_s = run_serial(cases)
        n_ok_s = sum(1 for _nm, ok, _e, _m in res_s if ok)
        print(f"[串行]           用例 {n_ok_s}/{len(cases)}  墙钟 {t_s:.1f}s")
        print(f"加速比 = {t_s / t_p:.2f}x")
        # 逐例对拍
        d = {nm: (ok, tuple(e)) for nm, ok, e, _m in res_p}
        ds = {nm: (ok, tuple(e)) for nm, ok, e, _m in res_s}
        diff = [nm for nm in ds if d.get(nm) != ds[nm]]
        print(f"逐例对拍不一致: {len(diff)} 例")
        for nm in diff[:10]:
            print(f"  ✗ {nm}")
            print(f"     并行: {d.get(nm)}")
            print(f"     串行: {ds[nm]}")
        if not diff:
            print("  ✓ 逐例结果（通过性 + 错误文本）完全一致 ⟹ 并行安全")
    else:
        print(f"（加 --check 可串行对拍并测加速比）")

    print(f"\n失败 {len(fails)} 例")
    for nm in fails[:15]:
        print(f"  {nm}")


if __name__ == "__main__":
    main()
