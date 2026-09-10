#!/usr/bin/env python3
"""性能剖析：全用例 3 遍计时取均值，输出 >50ms 热点排名与分布统计。

用法：
    python bench.py                 # 全量剖析
    python bench.py --top 30        # 只显示前 30 慢
    python bench.py --only "名称"   # 只跑匹配的用例（子串匹配）
    python bench.py --profile "名称"  # cProfile 单例
"""
import argparse
import gc
import json
import os
import sys
import time
import statistics
import cProfile
import pstats
import io

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chemkit.data import load_tables
from chemkit.engine import judge

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "chemkit/data/tests.json")


def load_cases():
    with open(DATA, encoding="utf-8") as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--only", type=str, default=None)
    ap.add_argument("--profile", type=str, default=None)
    ap.add_argument("--repeats", type=int, default=3)
    args = ap.parse_args()

    T = load_tables()
    judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, T)  # 预热

    cases = load_cases()
    if args.profile:
        pat = args.profile
        c = next((x for x in cases if pat in x["name"]), None)
        if c is None:
            print(f"未找到用例 {pat}")
            return
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        pr = cProfile.Profile()
        pr.enable()
        judge(subs, c.get("cond") or {"V_L": 1.0}, T)
        pr.disable()
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
        ps.print_stats(35)
        print(s.getvalue())
        return

    if args.only:
        cases = [c for c in cases if args.only in c["name"]]

    results = []   # (mean_t, name, cond_kind, n_steps)
    t_all0 = time.time()
    for ci, c in enumerate(cases):
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        # 用 sys 记录步数（从结果）
        times = []
        r = None
        for _ in range(args.repeats):
            gc.disable()           # 计时区禁 GC（模板缓存 ~70 万存活对象，
            t0 = time.perf_counter()   # gen-2 全量回收可达数百 ms，噪声不可控）
            r = judge(subs, dict(cond), T)
            times.append(time.perf_counter() - t0)
            gc.enable()
        if ci % 50 == 49:          # 周期性回收周期垃圾（多为引用计数已回收，
            gc.collect()           # 此处兜底），频率足够低不干扰计时
        mean_t = statistics.mean(times)
        iso = "iso" if cond.get("isothermal", True) else "non-iso"
        tk = float(cond.get("T_K", cond.get("T_C", 25.0) + 273.15))
        nonstd = "T%+d" % round(float(tk) - 298.15) if float(tk) != 298.15 else "T0"
        results.append((mean_t, c["name"], f"{iso}/{nonstd}",
                        len(r.get("steps", []))))
    t_all = time.time() - t_all0

    results.sort(reverse=True)
    print(f"共 {len(results)} 例 × {args.repeats} 遍，墙钟 {t_all:.1f}s\n")
    print("---- 最慢前 %d ----" % args.top)
    for t, name, kind, ns in results[:args.top]:
        print(f"  {t*1000:8.1f}ms  [{kind}]  steps={ns:<3d} {name}")

    ts = [t for t, *_ in results]
    print("\n---- 分布 ----")
    for thr in (0.01, 0.02, 0.05, 0.1, 0.5):
        n = sum(1 for t in ts if t > thr)
        print(f"  >{thr*1000:6.0f}ms: {n:4d} 例 ({100.0*n/len(ts):5.1f}%)")
    print(f"\n  均值 {statistics.mean(ts)*1000:.1f}ms  中位 {statistics.median(ts)*1000:.1f}ms  "
          f"最值 {max(ts)*1000:.1f}ms  P90 {sorted(ts)[int(0.9*len(ts))]*1000:.1f}ms")

    # 恒温/非恒温分组
    groups = {}
    for t, name, kind, ns in results:
        g = kind.split("/")[0]
        groups.setdefault(g, []).append(t)
    for g, ts2 in groups.items():
        if len(ts2) > 3:
            print(f"\n  [{g}] n={len(ts2)} 均值 {statistics.mean(ts2)*1000:.1f}ms "
                  f"中位 {statistics.median(ts2)*1000:.1f}ms 最值 {max(ts2)*1000:.1f}ms")
    # 温度分组
    groups2 = {}
    for t, name, kind, ns in results:
        g = kind.split("/")[1]
        groups2.setdefault(g, []).append(t)
    for g in sorted(groups2):
        ts2 = groups2[g]
        if len(ts2) > 3:
            print(f"  [T {g}] n={len(ts2)} 均值 {statistics.mean(ts2)*1000:.1f}ms "
                  f"中位 {statistics.median(ts2)*1000:.1f}ms 最值 {max(ts2)*1000:.1f}ms")


if __name__ == "__main__":
    main()
