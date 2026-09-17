# -*- coding: utf-8 -*-
"""tools/frozen_audit.py —— **冻结集审计**（§7 X-38）。

问题：`engine.resid_live_ok()`（质量口径的单一定义）显式排除 `frozen` 候选，
于是"冻结"成了残差指标看不见的盲区。实测反例 N30/D32：走步把
`Fe^{3+} + 3SCN^- -> [Fe(SCN)_3]` **两侧都在场、|S| = 4.02** 的通道永久冻结，
再以 `no-cands` 退出，而 `resid_src` 报的是另一条 `|S| = 0.0 / 已达平衡`
的候选 ⟹ 指标全绿、化学错在主产物上。

本工具做全库普查，回答三件事：
  1. **规模**：多少例的出口状态里存在 `two_sided ∧ frozen ∧ |S| > 阈值` 的通道？
  2. **盲区**：其中多少例的"活残差"（resid_src）反而很小（指标看不见）？
  3. **来源**：按 kind 分布 + |S| 最大的若干条（含 ext_max，用于判断是否
     "有实在的量可动"而不是痕量伪驱动）。

用法（仓库根）：
    python tools/frozen_audit.py [--thr 0.5] [--top 25] [--prefix N30 D32]
                                 [--cases 路径.json]

只读诊断：不改任何求解行为，只读 `judge(..., _probe={})` 的出口画像。
"""
from __future__ import annotations

import argparse
import io
import json
import os
import statistics
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                        # pragma: no cover
        pass

from chemkit.data import load_tables                          # noqa: E402
from chemkit.engine import judge                              # noqa: E402
from chemkit.testsuit import load_cases                       # noqa: E402


def _rows(probe: dict, thr: float, min_ext: float) -> list[dict]:
    """该例出口画像里 `two_sided ∧ frozen ∧ |S| > thr ∧ ext_max ≥ min_ext`。

    `ext_max` 是驱动方向反应物的**化学计量上限**（引擎探针给的"最多能走多少"）：
    |S| 大而 ext_max ≈ 0 是**痕量伪驱动**（引擎自己的注释即此，D38 型），
    冻结它是对的；要抓的是"有实在的量可动却被冻住"。默认阈值取引擎自身的
    显著度判据 `ANN_MIN_EXTENT = 1e-3`（与 `resid_live_ok` 同一把尺子）。
    """
    out = []
    for a in probe.get("active") or ():
        if not a.get("two_sided") or not a.get("frozen"):
            continue
        if abs(a.get("S") or 0.0) <= thr:
            continue
        if (a.get("ext_max") or 0.0) < min_ext:
            continue
        out.append(a)
    return out


def scan(cases: list[dict], T, thr: float, min_ext: float) -> dict:
    hit: list[tuple[str, dict, dict]] = []     # (用例名, 探针, 候选)
    stat = {"n": 0, "n_hit": 0, "n_blind": 0, "ch": 0, "kinds": Counter(),
            "absS": [], "live_max": [], "exts": [], "pairs": Counter()}
    for c in cases:
        stat["n"] += 1
        pr: dict = {}
        judge([{"name": n, "mol": m} for n, m in c["subs"]],
              c.get("cond") or {"V_L": 1.0}, T, _probe=pr)
        rs = _rows(pr, thr, min_ext)
        if not rs:
            continue
        stat["n_hit"] += 1
        stat["ch"] += len(rs)
        live = abs((pr.get("resid_src") or {}).get("S") or 0.0)
        stat["live_max"].append(live)
        if live <= 0.1:
            stat["n_blind"] += 1
        for a in rs:
            stat["kinds"][a.get("kind") or "?"] += 1
            stat["absS"].append(abs(a["S"]))
            stat["exts"].append(a.get("ext_max") or 0.0)
            hit.append((c["name"], pr, a))
    hit.sort(key=lambda t: -abs(t[2]["S"]))
    by_ext = sorted(hit, key=lambda t: -(t[2].get("ext_max") or 0.0))
    return {"stat": stat, "hit": hit, "by_ext": by_ext}


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("--thr", type=float, default=0.5,
                    help="|S| 阈值（默认 0.5 个 log 单位）")
    ap.add_argument("--min-ext", type=float, default=1e-3,
                    help="可动上限下限（默认 1e-3 = 引擎 ANN_MIN_EXTENT）")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--prefix", nargs="*", default=None)
    ap.add_argument("--cases", default=None)
    args = ap.parse_args(argv)

    T = load_tables()
    cases = load_cases(args.cases)
    if args.prefix:
        cases = [c for c in cases if c["name"].startswith(tuple(args.prefix))]
    r = scan(cases, T, args.thr, args.min_ext)
    st, hit = r["stat"], r["hit"]

    print(f"== frozen_audit ==  判据: two_sided ∧ frozen ∧ |S| > {args.thr} "
          f"∧ ext_max ≥ {args.min_ext:g}")
    print(f"扫描用例 {st['n']} 例")
    print(f"命中用例 {st['n_hit']} 例（{st['n_hit'] / max(st['n'], 1) * 100:.1f}%），"
          f"通道 {st['ch']} 条")
    if st["absS"]:
        s = sorted(st["absS"])
        print(f"|S| 分布: p50={s[len(s) // 2]:.3f} p90={s[int(len(s) * 0.9)]:.3f} "
              f"max={s[-1]:.3f}   ≥1 的通道 {sum(1 for x in s if x >= 1)} 条")
        e = sorted(st["exts"])
        print(f"ext_max 分布: p50={e[len(e) // 2]:.3g} p90={e[int(len(e) * 0.9)]:.3g} "
              f"max={e[-1]:.3g}   ≥0.01 的通道 {sum(1 for x in e if x >= 0.01)} 条")
    print(f"**指标盲区**（活残差 resid_src ≤ 0.1 却命中）：{st['n_blind']} 例"
          f"（占命中 {st['n_blind'] / max(st['n_hit'], 1) * 100:.0f}%）")
    print(f"按 kind 分布: {dict(st['kinds'].most_common())}")

    print(f"\n-- |S| 最大的 {args.top} 条 --")
    seen: set[str] = set()
    shown = 0
    for name, pr, a in hit:
        if shown >= args.top:
            break
        shown += 1
        dup = "  [同例]" if name in seen else ""
        seen.add(name)
        live = (pr.get("resid_src") or {})
        print(f"\n{shown:>3}. |S|={abs(a['S']):.3f}  {name}{dup}")
        print(f"     {a['eq'][:96]}")
        print(f"     kind={a.get('kind')} ext_max={a.get('ext_max')} "
              f"slow={int(bool(a.get('slow')))} blocked={int(bool(a.get('blocked')))}")
        print(f"     该例活残差: |S|={abs(live.get('S') or 0.0):.3f} "
              f"（{str(live.get('eq') or '无')[:60]}）")

    print(f"\n-- ext_max 最大的 {args.top} 条（有实在的量可动却被冻住）--")
    for i, (name, pr, a) in enumerate(r["by_ext"][:args.top], 1):
        live = (pr.get("resid_src") or {})
        print(f"{i:>3}. ext_max={a.get('ext_max'):<10.4g} |S|={abs(a['S']):<8.3f} "
              f"{name}")
        print(f"     {a['eq'][:96]}")
        print(f"     该例活残差: |S|={abs(live.get('S') or 0.0):.3f}")

    out = os.path.join(ROOT, ".tmp_frozen_audit.json")
    with io.open(out, "w", encoding="utf-8") as f:
        json.dump({"thr": args.thr, "stat": {
            k: (dict(v) if isinstance(v, Counter) else v)
            for k, v in st.items()},
            "hit": [{"case": n, "S": a["S"], "eq": a["eq"], "kind": a.get("kind"),
                     "ext_max": a.get("ext_max")} for n, _p, a in hit]},
            f, ensure_ascii=False, indent=1)
    print(f"\n[落盘] {os.path.basename(out)}（全量命中明细）")
    return 1 if st["n_hit"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
