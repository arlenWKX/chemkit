# -*- coding: utf-8 -*-
"""并行跑套件 —— **薄封装**，真正的实现在 `chemkit/testsuit.py`。

⚙️ 第 241 轮重构（根因修，不是补丁）：
本脚本原先自己建 `ProcessPoolExecutor` 并只回传 `(name, ok, errors, ms)`，
另写一份 `logs/suite-parallel-latest.json`。它有三个真实缺陷：

1. **错误文本丢失**：它从 `ts.FAILS[n0:]` 取本次失败项，而 `testsuit` 在
   失败时往 `FAILS` 里 append 的是**用例名、不是错误文本**（错误文本进
   `RESULTS[i]["errors"]`）⟹ 留档里 `errors` 只有 `["9 BaSO4+Na2CO3(浓)"]`
   这种"名字当错误"的垃圾，看失败原因必须回去翻控制台。
2. **质量口径字段全缺**：没有 `resid_live` / `resid_max` / `iters` /
   `exit` / `pH_solver` / `resid_src_eq`，于是"残差"还得另跑 `converg`。
3. **两份实现必然漂移**：`testsuit` 的 `_par_init` 预热、`_SLOW_FIRST` 慢例
   优先调度、`_enrich_results_with_probe` 富化，本脚本一个都没有——包括那条
   代价高昂的预热教训（不预热则每 worker 首例假性 12~24 s）。

现在它只做一件事：调 `testsuit.main(jobs=N)`（**唯一实现**），把并行统计
写进 `logs/suite-parallel-latest.json`（只含墙钟/worker/CPU 合计/最慢榜，
逐例结果不在此处复制 —— 那在 `logs/suite-latest.json`）。

用法：
  python tools/suite_parallel.py            # 自动选 worker 数
  python tools/suite_parallel.py 6          # 指定 worker 数
  python tools/suite_parallel.py 6 --check  # 额外串行跑一遍对拍（慢但严谨）
"""
import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit import testsuit as ts                            # noqa: E402


def _default_jobs():
    """按**物理核**取，不按逻辑核取满。

    实测（i7-8650U = 4 物理 / 8 逻辑）：2w 74.8s / 4w 55.5s / 6w 50.0s /
    8w 51.5s ⟹ 超过物理核只增争用（单例 CPU 成本 111→220 ms）。
    """
    cpu = os.cpu_count() or 2
    return max(1, min(4, cpu // 2))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    do_check = "--check" in sys.argv
    jobs = int(args[0]) if args else _default_jobs()
    if jobs <= 1:
        print("!! jobs<=1 请直接跑 `python tools/dev.py suite --jobs=1`")
        return 2

    # 捕获 canonical 并行路径的墙钟（不改 testsuit 接口：临时换引用再还原）
    wall = {}

    def _spy(cases, n):
        r = _orig(cases, n)
        wall["s"] = r[1]
        return r

    _orig = ts.run_cases_parallel
    ts.run_cases_parallel = _spy
    try:
        print(f"[并行] worker={jobs}（cpu_count={os.cpu_count()}，spawn 模式）")
        rc = ts.main(jobs=jobs)
    finally:
        ts.run_cases_parallel = _orig

    if do_check:
        print("\n[对拍] 串行再跑一遍（严谨但慢）…")
        keep = {r["name"]: (r["ok"], tuple(r["errors"])) for r in ts.RESULTS}
        ts.main(jobs=1)
        now = {r["name"]: (r["ok"], tuple(r["errors"])) for r in ts.RESULTS}
        diffs = [n for n in keep if now.get(n) != keep[n]]
        print(f"[对拍] 逐例不一致 {len(diffs)} 例"
              + ("（通过性 + 错误文本完全一致 ⟹ 并行安全）" if not diffs else ""))
        for n in diffs[:10]:
            print(f"  ✗ {n}\n     并行: {keep[n]}\n     串行: {now[n]}")

    # ---- 并行统计侧车（逐例结果不在这里复制）----
    times = sorted((t, n) for t, n in ts.TIMES)
    tot = sum(t for t, _ in times)
    top = [{"name": n, "ms": round(t * 1000, 1)}
           for t, n in reversed(times[-40:])]
    n_ok = sum(1 for r in ts.RESULTS if r["ok"])
    payload = {
        "n_workers": jobs,
        "wall_s": round(wall.get("s", 0.0), 2),
        "n_cases": len(ts.RESULTS),
        "n_ok": n_ok,
        "n_fail": len(ts.RESULTS) - n_ok,
        "cpu_total_s": round(tot, 1),
        "slowest_case_s": round(times[-1][0], 2) if times else None,
        "slowest": top,
        "cases_json": "logs/suite-latest.json",
        "note": "逐例 ok/errors/ms/pH/resid_* 在 logs/suite-latest.json"
                "（dev.py suite 写，本文件不复制）",
    }
    lp = os.path.join(ROOT, "logs", "suite-parallel-latest.json")
    with io.open(lp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=1)
    print(f"\n[留档] {os.path.relpath(lp, ROOT)}"
          f"（worker 数 / 墙钟 / CPU 合计 / 最慢榜）")
    print("[留档] logs/suite-latest.json（逐例 ok/errors/ms/pH/resid_*）")
    print("  ⟹ 看任何一面都用 `python tools/suite_show.py`，**不要重跑**")
    return rc


if __name__ == "__main__":
    sys.exit(main())
