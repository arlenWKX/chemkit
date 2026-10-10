# -*- coding: utf-8 -*-
"""`tools/inspect.py` —— **结果读取器**（第 146 轮）：把套件/指标的结构化留档读成
人可读摘要，避免"为了看细节而重跑套件"。

## 为什么存在（用户既定纪律）
一次运行就把原始信息全量落盘（`dev.py suite` → `.tmp_dev_results.json`、
`dev.py perf` → `.tmp_dev_converg.json`），后续**用本工具反复读**，
不要在生成端截断、也不要靠重跑换取信息。

## 子命令
    python tools/inspect.py fails [--grep 关键词] [--full]
        列出套件留档里的失败用例（名称 + 全部 errors + note 摘要）
    python tools/inspect.py case <前缀>
        某例的留档详情（ok/errors/pH/degree/changed/annotations/净方程）
    python tools/inspect.py stats
        留档概览：总例数、失败数、耗时分布、按前缀分组
    python tools/inspect.py converg [--top N]
        读 `.tmp_dev_converg.json`（perf 留档），列残差最高的 N 例

留档路径可用 `--results` / `--converg` 覆盖（默认 `.tmp_dev_results.json` /
`.tmp_dev_converg.json`）。本工具**只读**，不改任何状态、不跑套件。
"""
import argparse
import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))   # 仓库根（chemkit/helper/x.py -> 上三级）
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def _load(path: str):
    if not os.path.exists(path):
        print(f"!! 留档不存在：{path}（先跑一次 `dev.py suite` / `dev.py perf`）")
        return None
    with io.open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _fmt_case(r: dict, full: bool = False) -> str:
    out = [f"  [{'ok ' if r.get('ok') else 'FAIL'}] {r.get('name')}"
           f"  pH={r.get('pH')} degree={r.get('degree')}"
           f" changed={r.get('changed')} {r.get('ms')}ms"]
    for e in (r.get("errors") or []):
        out.append(f"      · {e if full else e[:220]}")
    if r.get("note"):
        out.append(f"      note: {r['note'][:200] if not full else r['note']}")
    if r.get("net_equation"):
        out.append(f"      eq: {r['net_equation']}")
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("cmd", choices=["fails", "case", "stats", "converg", "diff"])
    ap.add_argument("arg", nargs="?", default=None)
    ap.add_argument("arg2", nargs="?", default=None)
    ap.add_argument("--grep", default=None, help="只列名字含该子串的失败")
    ap.add_argument("--full", action="store_true", help="不截断 errors/note")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--results", default=os.path.join(ROOT, "logs", "suite-latest.json"))
    ap.add_argument("--converg", default=os.path.join(ROOT, "logs", "perf-latest.json"))
    a = ap.parse_args()

    if a.cmd == "diff":
        # 两份 converg 留档逐例对比——"哪一例变慢了/残差涨了"必须能在**不重跑**
        # 的前提下回答（历史基线就留在 git 里，`git show <rev>:converg-baseline.json`）。
        da, db = _load(a.arg), _load(a.arg2)
        if da is None or db is None:
            return 1
        ma = {c["name"]: c for c in da["cases"]}
        mb = {c["name"]: c for c in db["cases"]}
        rows = []
        for nm in sorted(set(ma) & set(mb)):
            x, y = ma[nm], mb[nm]
            rows.append((y.get("iters", 0) - x.get("iters", 0),
                         (y.get("resid_live") or 0.0) - (x.get("resid_live") or 0.0),
                         nm, x, y))
        print(f"逐例对比：{len(rows)} 例（Δiters 降序）")
        print(f"{'Δiters':>8} {'Δresid':>8}  {'iters 前->后':>17}"
              f"  {'resid 前->后':>19}  用例")
        for di, dr, nm, x, y in sorted(rows, reverse=True)[:a.top]:
            print(f"{di:8d} {dr:8.3f}  {x.get('iters', 0):7d} -> {y.get('iters', 0):<7d}"
                  f"  {x.get('resid_live') or 0.0:8.3f} -> {y.get('resid_live') or 0.0:<8.3f}"
                  f"  {nm}")
        print(f"\n合计 Δiters = {sum(r[0] for r in rows):+d}；"
              f"变慢 {sum(1 for r in rows if r[0] > 0)} 例，"
              f"变快 {sum(1 for r in rows if r[0] < 0)} 例")
        wors = sorted((r for r in rows if r[1] > 0.05), key=lambda r: -r[1])[:8]
        if wors:
            print("残差变差最多（Δresid_live）：")
            for _di, dd, nm, x, y in wors:
                print(f"   {dd:7.3f}  {nm}  ({(x.get('resid_live') or 0):.3f} -> "
                      f"{(y.get('resid_live') or 0):.3f})")
        return 0

    if a.cmd == "converg":
        d = _load(a.converg)
        if d is None:
            return 1
        rows = d if isinstance(d, list) else d.get("cases", [])
        def resid(r):
            return max(abs(r.get("resid_live") or 0.0), abs(r.get("max_abs_S") or 0.0))
        rows = sorted(rows, key=lambda r: -resid(r))[:a.top]
        print(f"{'resid':>9}  {'iters':>6}  用例")
        for r in rows:
            print(f"{resid(r):9.3f}  {r.get('iters', 0):6d}  {r.get('name')}")
        return 0

    d = _load(a.results)
    if d is None:
        return 1
    if a.cmd == "stats":
        bad = [r for r in d if not r.get("ok")]
        ms = sorted(r.get("ms", 0.0) for r in d)
        pre = {}
        for r in d:
            nm = r.get("name", "")
            key = nm.split(" ")[0] if nm[:1].isalpha() else "数字"
            pre[key] = pre.get(key, 0) + 1
        print(f"总例数 {len(d)}  失败 {len(bad)}  通过 {len(d) - len(bad)}")
        if ms:
            print(f"耗时 ms：p50={ms[len(ms) // 2]:.1f} max={ms[-1]:.1f}")
        top = sorted(pre.items(), key=lambda kv: -kv[1])[:12]
        print("前缀分布（前 12）：", ", ".join(f"{k}:{v}" for k, v in top))
        return 0
    if a.cmd == "case":
        hit = [r for r in d if str(r.get("name", "")).startswith(a.arg or "")]
        if not hit:
            print(f"!! 留档里没有以 {a.arg!r} 开头的用例")
            return 1
        for r in hit:
            print(_fmt_case(r, full=True))
        return 0
    bad = [r for r in d if not r.get("ok")]
    if a.grep:
        bad = [r for r in bad if a.grep in str(r.get("name", ""))]
    print(f"失败 {len(bad)} 例")
    for r in bad:
        print(_fmt_case(r, full=a.full))
    return 0


if __name__ == "__main__":
    sys.exit(main())
