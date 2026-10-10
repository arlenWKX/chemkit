# -*- coding: utf-8 -*-
"""tools/order_audit.py —— **执行顺序依赖审计**（§7 X-39 第 7 项）。

问题：结果是否依赖**用例执行顺序**（静态层/缓存/预热路径）？这类缺陷不会在
单例复跑里显形，只在"换一个顺序"时冒出来。历史锚点：J14 `PbCl₂@363K` 曾冷启动
给 pH 3.0 / Pb²⁺ 0.015631，而先判过任一其它温度后给 pH 6.21 / 0.027474
（§7 X-29 顺带发现）——后经 X-31（OH⁻ 模板 pKw 双折算根治）消失。

指纹 = 每例的 (pH, degree, changed, reacted, 净方程, 步数, H_excess, 终态物种表)。
三个顺序各在**独立进程**里跑（同一进程内逐例独立，跨进程换序才暴露路径依赖）：

    python tools/order_audit.py                # 自跑 normal/reverse/shuffle 并比较
    python tools/order_audit.py --dump ORD F   # 只落一份指纹（ORD: normal/reverse/shuffle）
    python tools/order_audit.py --cmp A B      # 比较两份指纹

退出码：0 = 无差异；1 = 有差异（明细打到 stdout）。
"""
from __future__ import annotations

import argparse
import io
import json
import os
import random
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                        # pragma: no cover
        pass

LAB = ("pH", "deg", "changed", "reacted", "eq", "steps", "He", "final")
SEED = 20260101


def dump(order: str, path: str) -> int:
    from chemkit import Reaction
    from chemkit.data import load_tables
    from chemkit.engine import judge
    from chemkit.testsuit import load_cases

    T = load_tables()
    cases = load_cases(None)
    if order == "reverse":
        cases = list(reversed(cases))
    elif order == "shuffle":
        random.Random(SEED).shuffle(cases)
    data = {}
    for c in cases:
        r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
        R = Reaction(r)
        data[c["name"]] = [
            None if r["final_pH"] is None else round(r["final_pH"], 6),
            r["degree"], bool(r["changed"]), bool(r["reacted"]),
            None if R.net_equation is None else R.net_equation.plain(),
            len(r["steps"]), round(r.get("H_excess", 0.0), 9),
            sorted((e["name"], round(e["mol"], 9)) for e in r["final"]
                   if e["mol"] > 1e-7),
        ]
    with io.open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, sort_keys=True)
    print(f"[{order}] {len(data)} 例 -> {path}")
    return 0


def cmp_files(pa: str, pb: str) -> int:
    a = json.load(io.open(pa, encoding="utf-8"))
    b = json.load(io.open(pb, encoding="utf-8"))
    diff = [k for k in a if a[k] != b.get(k)]
    print(f"[cmp] {os.path.basename(pa)} n={len(a)}  vs  "
          f"{os.path.basename(pb)} n={len(b)}   diff={len(diff)}")
    for k in diff[:15]:
        x, y = a[k], b[k]
        print(f" * {k}")
        for i, lab in enumerate(LAB):
            if x[i] != y[i]:
                print(f"    {lab}: {str(x[i])[:70]}  ->  {str(y[i])[:70]}")
    return 1 if diff else 0


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("--dump", nargs=2, metavar=("ORDER", "FILE"))
    ap.add_argument("--cmp", nargs=2, metavar=("A", "B"))
    args = ap.parse_args(argv)
    if args.dump:
        return dump(args.dump[0], args.dump[1])
    if args.cmp:
        return cmp_files(args.cmp[0], args.cmp[1])
    # 默认：三个顺序各起一个**子进程**（路径依赖只在跨进程换序时显形）
    files = {}
    for order in ("normal", "reverse", "shuffle"):
        path = os.path.join(ROOT, f".tmp_ord_{order}.json")
        r = subprocess.run([sys.executable, os.path.abspath(__file__),
                            "--dump", order, path],
                           cwd=ROOT, capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        print((r.stdout or "").strip() or (r.stderr or "")[-200:])
        if r.returncode != 0:
            return 2
        files[order] = path
    rc = 0
    rc |= cmp_files(files["normal"], files["reverse"])
    rc |= cmp_files(files["normal"], files["shuffle"])
    print("[判定] 无顺序依赖" if rc == 0
          else "[判定] **存在顺序依赖**（见上方明细）")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
