# -*- coding: utf-8 -*-
"""第 207 轮 · worker 数 A/B 实测（使用者建议按物理核数取 4）。

本机实测：`i7-8650U` = **4 物理核 / 8 逻辑核（HT）**。
Python 进程 CPU 密集且各有独立解释器 ⟹ 8 个 worker **超订** 4 个物理核，
HT 只共享执行单元、不增吞吐 ⟹ 8 未必快于 4，还可能因调度/缓存争用更慢。

⚠️ **spawn 陷阱（本轮踩到）**：Windows 只有 `spawn`，子进程会**重新 import
本模块**。若脚本主体在模块顶层执行，子进程就会再次建池 ⟹ 递归派生 ⟹
`BrokenProcessPool`。**凡是自己建进程池的脚本，主体必须放进
`if __name__ == "__main__":`**（库函数 `testsuit` 本身没问题）。

一次运行测 2/4/6/8，并一次写全部产物（不为看别的对照而重跑）。

用法： python tools/workers_ab.py [重复次数，默认 1]
"""
import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def main():
    from chemkit import testsuit
    from chemkit.testsuit import load_cases

    rep = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    counts = (2, 4, 6, 8)
    cases = load_cases(None)
    print(f"用例 {len(cases)}；物理核 4 / 逻辑核 8；重复 {rep} 次")
    print(f"{'workers':>8} {'墙钟s':>9} {'用例':>12} {'CPU合计s':>10} "
          f"{'最慢单例s':>10}")

    out = {}
    for n in counts:
        walls, oks, cpus, slow = [], [], [], []
        for _r in range(rep):
            testsuit.TIMES.clear()
            testsuit.RESULTS.clear()
            testsuit.FAILS.clear()
            n_ok, wall = testsuit.run_cases_parallel(cases, n)
            walls.append(wall)
            oks.append(n_ok)
            ms = [t for t, _nm in testsuit.TIMES]
            cpus.append(sum(ms))
            slow.append(max(ms) if ms else 0.0)
        w = min(walls)
        print(f"{n:>8} {w:>9.1f} {oks[0]:>9}/{len(cases)} "
              f"{cpus[0]:>10.1f} {slow[0]:>10.2f}", flush=True)
        out[n] = {"wall_s": round(w, 1), "n_ok": oks[0],
                  "cpu_total_s": round(cpus[0], 1),
                  "slowest_case_s": round(slow[0], 2),
                  "walls_all": [round(x, 1) for x in walls]}

    base = out[max(counts)]["wall_s"]
    print(f"\n相对 {max(counts)} workers 的加速：")
    for n in counts:
        print(f"  {n:>2} workers: {base / out[n]['wall_s']:.3f}x  "
              f"（墙钟 {out[n]['wall_s']}s）")
    best = min(out, key=lambda k: out[k]["wall_s"])
    print(f"\n最快 = **{best} workers**（{out[best]['wall_s']}s）")

    jp = os.path.join(ROOT, "logs", "workers_ab.json")
    with io.open(jp, "w", encoding="utf-8") as f:
        json.dump({"physical_cores": 4, "logical": 8, "results": out},
                  f, ensure_ascii=False, indent=1)
    print(f"已写 {os.path.relpath(jp, ROOT)}")
    print("\n判读：物理核=4 ⟹ 4 worker 通常最优；>4 只增加 HT 争用。")


if __name__ == "__main__":
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()
