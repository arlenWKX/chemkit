"""自缓冲步普查（§7 X-33）：全库有多少步"本步产物替本步吃掉质子"。

判据（与修法同一判据，见 engine.SELFBUF_AUDIT 与 speciation.complex_capacity）：
候选产物里 β_pka 配离子的总吸收容量 Σ(d·νH⁺) ≥ 本步净释出 H⁺ ⟹ 滴定把这批
质子原地退回给产物自己，`He_res ≡ 0` 与步长无关 ⟹ pH 与步长解耦，该步会跑到
"产物/反应物 = 10^(pH−pK)"的幻影量。

本工具**不改任何行为**，只把规模量出来（哪些配离子、哪些用例、幻影放大倍数），
为"走步全局化 + 落地禁令"排工期。

用法：
    python tools/selfbuf.py [-n 25] [--only 前缀,前缀] [--case 用例前缀]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter

sys.path.insert(0, ".")

import chemkit.engine as E                                  # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.engine import judge                            # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=25)
    ap.add_argument("--only", default=None, help="逗号分隔的用例名前缀")
    ap.add_argument("--case", default=None, help="单用例前缀（打印其全部命中步）")
    args = ap.parse_args()

    T = load_tables()
    cases = load_cases(None)
    if args.only:
        pre = tuple(x.strip() for x in args.only.split(",") if x.strip())
        cases = [c for c in cases if c["name"].startswith(pre)]
    if args.case:
        cases = [c for c in cases if c["name"].startswith(args.case)]

    sp_hits: Counter = Counter()
    kind_hits: Counter = Counter()
    case_hits: Counter = Counter()
    n_cases = 0
    for c in cases:
        E.SELFBUF_AUDIT = []
        r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
        hits = E.SELFBUF_AUDIT
        E.SELFBUF_AUDIT = None
        if not hits:
            continue
        n_cases += 1
        # 去重：同一步在二分内被反复求值 ⟹ 按 (产物, νH⁺) 记一次
        seen = set()
        for h in hits:
            key = (h["prod"], h["nu_h"])
            if key in seen:
                continue
            seen.add(key)
            for s in h["prod"]:
                sp_hits[s] += 1
            kind_hits[h["kind"]] += 1
            case_hits[c["name"]] += 1
            if args.case:
                print(f"  s={h['eq']:<34} νH⁺={h['nu_h']:.0f} "
                      f"容量={h['cap']:.0f} 产物={h['prod']}")
                print(f"     {h['r']} -> {h['pr']}")
        if args.case:
            print(f"  [用例 {c['name']}] pH={r.get('final_pH')} "
                  f"步数={len(r.get('steps') or [])} 命中 {len(seen)} 类")
    print(f"\n普查 {len(cases)} 例：**命中自缓冲步的用例 {n_cases} 例**"
          f"（{100.0 * n_cases / max(len(cases), 1):.1f}%）")
    print("  按步类型:", dict(kind_hits))
    print("  按配离子（命中用例数）:")
    for sp, k in sp_hits.most_common(args.n):
        print(f"    {sp:<24} {k}")
    print("  按用例（命中类数）:")
    for name, k in case_hits.most_common(args.n):
        print(f"    {name[:52]:<54} {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
