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
FLAT = os.path.join(ROOT, "logs", "suite-latest.json")


def _load(path):
    with io.open(path, encoding="utf-8") as f:
        return json.load(f)


# --- 两种留档形态归一 ---------------------------------------------------
# 形态 P（并行器 `tools/suite_parallel.py`）：dict，含 n_workers/wall_s/
#        fails/slowest/cases。
# 形态 F（`dev.py suite`）：**顶层 list**，逐例 dict。
# 两者都含逐例 ok/errors/ms。**读档器必须两种都吃**，否则会出现
# "留档明明在、工具说没留档"这种最浪费时间的假故障。
def _norm(o, path):
    if isinstance(o, list):
        cases = o
        n_ok = sum(1 for c in cases if c["ok"])
        return {
            "cases": cases,
            "n_cases": len(cases),
            "n_ok": n_ok,
            "n_fail": len(cases) - n_ok,
            "fails": [c["name"] for c in cases if not c["ok"]],
            "slowest": sorted(cases, key=lambda c: -c.get("ms", 0)),
            "n_workers": None,
            "wall_s": None,
            "cpu_total_s": None,
            "slowest_case_s": round(max((c.get("ms", 0)
                                         for c in cases), default=0) / 1000.0, 1),
            "path": path,
        }
    d = dict(o)
    d["path"] = path
    return d


def _resolve(tag):
    """tag 为 None 时用默认留档（并行档优先，缺失则退到 flat 档）。"""
    if tag is None:
        for p in (LATEST, FLAT):
            if os.path.exists(p):
                return p
        return None
    return tag if os.path.isabs(tag) else os.path.join(ROOT, tag)


def main():
    p0 = _resolve(None)
    if p0 is None:
        print("没有留档。先跑一次： python tools/dev.py suite")
        return 1
    d = _norm(_load(p0), p0)
    mode = sys.argv[1] if len(sys.argv) > 1 else "summary"

    if mode == "summary":
        print(f"留档: {os.path.relpath(d['path'], ROOT)}")
        if d["n_workers"] is not None:
            print(f"  worker={d['n_workers']}  墙钟={d['wall_s']}s  "
                  f"串行CPU合计={d['cpu_total_s']}s  最慢单例="
                  f"{d['slowest_case_s']}s")
        else:
            print(f"  （dev.py suite 形态：无并行统计；"
                  f"最慢单例={d['slowest_case_s']}s）")
        print(f"  用例 {d['n_ok']}/{d['n_cases']} PASS，失败 {d['n_fail']} 例")
        print(f"\n失败清单（{d['n_fail']}）:")
        for nm in d["fails"]:
            print(f"  {nm}")
        print("\n（改用 `fails` / `slow [N]` / `case <前缀>` / "
              "`diff <旧> [新]>` 看其它面——**不需要重跑**）")

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
        if len(sys.argv) < 3:
            print("用法: suite_show.py diff <旧留档> [新留档]\n"
                  "  省略第二个路径时 A 取默认留档（并行档优先）。\n"
                  "  ⚠️ **两份都应是逐例形态**（`logs/suite-latest.json`）——\n"
                  "  并行侧车档只有墙钟/最慢榜，没有逐例字段，不能当 A/B。")
            return 1
        if len(sys.argv) >= 4:
            p0 = _resolve(sys.argv[3])
            d = _norm(_load(p0), p0)
        p0 = d["path"]
        p1 = _resolve(sys.argv[2])
        ob = _norm(_load(p1), p1)
        a = {c["name"]: c for c in d["cases"]}
        b = {c["name"]: c for c in ob["cases"]}
        missing = [k for k in ("pH", "resid_live", "errors")
                   if not any(c.get(k) is not None for c in d["cases"])]
        if missing:
            print(f"!! A 缺少逐例字段 {missing} ⟹ 它不是**逐例档**；"
                  "数值漂移/失败原因两节不可用。\n"
                  "   请用 `logs/suite-latest.json`（`dev.py suite` 写）。\n")
        only_a = [n for n in a if n not in b]
        only_b = [n for n in b if n not in a]
        flipped = []
        for n in a:
            if n in b and a[n]["ok"] != b[n]["ok"]:
                flipped.append((n, b[n]["ok"], a[n]["ok"]))
        print(f"A = {os.path.relpath(d['path'], ROOT)}  "
              f"{d['n_ok']}/{d['n_cases']}")
        print(f"B = {os.path.relpath(ob['path'], ROOT)}  "
              f"{ob['n_ok']}/{ob['n_cases']}")
        print(f"仅在 A: {len(only_a)}   仅在 B: {len(only_b)}")
        print(f"**通过性翻转: {len(flipped)} 例**")
        for n, ob_ok, oa_ok in flipped:
            print(f"  {n}: {'PASS' if ob_ok else 'FAIL'} -> "
                  f"{'PASS' if oa_ok else 'FAIL'}")
            if not oa_ok:
                for e in a[n]["errors"][:3]:
                    print(f"      A: {e}")
            if not ob_ok:
                for e in b[n]["errors"][:3]:
                    print(f"      B: {e}")
        if missing:
            return 0
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
        # 数值漂移（通过性不变但 pH/残差动了）——发版前必看
        drift = []
        for n in a:
            if n not in b:
                continue
            x, y = a[n], b[n]
            dp = abs((x.get("pH") or 0.0) - (y.get("pH") or 0.0))
            dr = abs((x.get("resid_live") or 0.0)
                     - (y.get("resid_live") or 0.0))
            if dp > 1e-6 or dr > 1e-6:
                drift.append((n, y.get("pH"), x.get("pH"),
                              y.get("resid_live"), x.get("resid_live")))
        print(f"\n数值漂移（pH 或 resid_live 变化 >1e-6）: {len(drift)} 例")
        for n, pb, pa, rb, ra in drift[:15]:
            print(f"  {n}: pH {pb} -> {pa}   resid {rb} -> {ra}")
    else:
        print(f"未知模式 {mode!r}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
