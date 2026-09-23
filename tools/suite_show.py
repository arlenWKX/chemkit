# -*- coding: utf-8 -*-
"""第 206 轮 · 读套件留档（**替代重跑**）。

铁律（使用者指出）：**跑一次，写全部产物；之后只读，不为看另一面而重跑。**
`tools/suite_parallel.py` 一次运行写 `logs/suite-parallel-latest.json`
（含逐例 ok / errors / ms，以及慢例榜与失败清单）。本脚本负责展示。

用法：
  python tools/suite_show.py                 # 概要 + 失败清单
  python tools/suite_show.py fails           # 只看失败
  python tools/suite_show.py slow [N]        # 最慢 N 例（默认 20）
  python tools/suite_show.py case <前缀>      # 某例的 ok/errors/ms
  python tools/suite_show.py diff <另一份>    # 两份留档的差异（改前/改后对比）
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

LATEST = os.path.join(ROOT, "logs", "suite-parallel-latest.json")


def _load(path):
    with io.open(path, encoding="utf-8") as f:
        return json.load(f)


def main():
    if not os.path.exists(LATEST):
        print("没有留档。先跑一次： python tools/suite_parallel.py 8")
        return 1
    d = _load(LATEST)
    mode = sys.argv[1] if len(sys.argv) > 1 else "summary"

    if mode == "summary":
        print(f"留档: {os.path.relpath(LATEST, ROOT)}")
        print(f"  worker={d['n_workers']}  墙钟={d['wall_s']}s  "
              f"串行CPU合计={d['cpu_total_s']}s  最慢单例="
              f"{d['slowest_case_s']}s")
        print(f"  用例 {d['n_ok']}/{d['n_cases']} PASS，失败 {d['n_fail']} 例")
        print(f"\n失败清单（{d['n_fail']}）:")
        for nm in d["fails"]:
            print(f"  {nm}")
        print(f"\n（改用 `fails` / `slow [N]` / `case <前缀>` / "
              f"`diff <另一份>` 看其它面——**不需要重跑**）")

    elif mode == "fails":
        for c in d["cases"]:
            if not c["ok"]:
                print(f"{c['name']}")
                for e in c["errors"]:
                    print(f"    {e}")

    elif mode == "slow":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 20
        for r in d["slowest"][:n]:
            print(f"  {r['ms'] / 1000:7.2f}s  {r['name']}")

    elif mode == "case":
        pre = sys.argv[2] if len(sys.argv) > 2 else ""
        hit = [c for c in d["cases"]
               if c["name"].startswith(pre) or c["name"] == pre]
        if not hit:
            print(f"无匹配 {pre!r}")
            return 1
        for c in hit:
            print(f"{c['name']}  ok={c['ok']}  {c['ms']}ms")
            for e in c["errors"]:
                print(f"    {e}")

    elif mode == "diff":
        other = sys.argv[2]
        o = _load(other if os.path.isabs(other)
                  else os.path.join(ROOT, other))
        a = {c["name"]: c for c in d["cases"]}
        b = {c["name"]: c for c in o["cases"]}
        only_a = [n for n in a if n not in b]
        only_b = [n for n in b if n not in a]
        flipped = []
        for n in a:
            if n in b and a[n]["ok"] != b[n]["ok"]:
                flipped.append((n, b[n]["ok"], a[n]["ok"]))
        print(f"A(latest) {d['n_ok']}/{d['n_cases']}  "
              f"B({os.path.basename(other)}) {o['n_ok']}/{o['n_cases']}")
        print(f"仅在 A: {len(only_a)}   仅在 B: {len(only_b)}")
        print(f"**通过性翻转: {len(flipped)} 例**")
        for n, ob, oa in flipped:
            print(f"  {n}: {'PASS' if ob else 'FAIL'} -> "
                  f"{'PASS' if oa else 'FAIL'}")
            if not oa:
                for e in a[n]["errors"][:3]:
                    print(f"      A: {e}")
            if not ob:
                for e in b[n]["errors"][:3]:
                    print(f"      B: {e}")
        # 错误文本变化（同为 FAIL 但原因不同）
        chg = [(n, b[n]["errors"], a[n]["errors"]) for n in a
               if n in b and not a[n]["ok"] and not b[n]["ok"]
               and a[n]["errors"] != b[n]["errors"]]
        if chg:
            print(f"\n失败原因变化（仍 FAIL）: {len(chg)} 例")
            for n, eb, ea in chg[:10]:
                print(f"  {n}")
                print(f"      B: {eb[:2]}")
                print(f"      A: {ea[:2]}")
    else:
        print(f"未知模式 {mode!r}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
