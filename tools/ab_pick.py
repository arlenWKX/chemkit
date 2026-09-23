# -*- coding: utf-8 -*-
"""第 264 轮 · 拾取规则的 A/B：新规则（同产物可动量回退）vs 旧规则（纯 S 贪心）。

`_pick_ev` 定义在 `engine.py` 内、且被 `engine.py` 内部按**模块全局名**调用，
故把 `engine._pick_ev` 换回旧实现即得"改前行为"（同一进程内可比）。

用法：
    python tools/ab_pick.py I31 T52 H32 K01
    python tools/ab_pick.py I31 --json logs/abpick_i31.json
"""
from __future__ import annotations

import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as _eng                     # noqa: E402
from chemkit.converg import _live                 # noqa: E402
from chemkit.data import load_tables              # noqa: E402

_NEW = _eng._pick_ev


def _old(evals, ledger):
    """第 264 轮之前的口径：纯 S 贪心。"""
    return max(evals, key=lambda e: e[2])


def run(case, T, new: bool):
    _eng._pick_ev = _NEW if new else _old
    try:
        from chemkit.engine import judge
        pr = {}
        subs = [{"name": s[0], "mol": float(s[1])} for s in case["subs"]]
        r = judge(subs, case.get("cond") or {}, T, _probe=pr)
    finally:
        _eng._pick_ev = _NEW
    return {
        "pH": r.get("final_pH"),
        "resid": abs(_live(pr.get("active") or [])),
        "iters": pr.get("iters"), "exit": pr.get("exit"),
        "steps": [(s["equation"], round(s["extent"], 8))
                  for s in (r.get("steps") or [])],
        "final": {e["name"]: e["mol"] for e in (r.get("final") or [])
                  if e["mol"] > 1e-9},
        "net": r.get("net_equation"),
    }


def main(argv: list[str]) -> int:
    names = [a for a in argv if not a.startswith("-")]
    out_path = None
    if "--json" in argv:
        out_path = argv[argv.index("--json") + 1]
        names = [n for n in names if n != out_path]
    T = load_tables()
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = {c["name"]: c for c in json.load(fh)}
    res = {}
    for pre in names:
        hit = [c for k, c in cases.items() if k.startswith(pre)]
        if not hit:
            print(f"[跳过] 无匹配前缀 {pre}")
            continue
        c = hit[0]
        o, n = run(c, T, False), run(c, T, True)
        res[c["name"]] = {"old": o, "new": n}
        print(f"\n===== {c['name']} =====")
        for tag, x in (("旧(纯S贪心)", o), ("新(同产物回退)", n)):
            print(f"  {tag}: pH={x['pH']} resid={x['resid']:.4f} "
                  f"iters={x['iters']} exit={x['exit']} steps={len(x['steps'])}")
        sk = sorted(set(o["final"]) | set(n["final"]))
        d = [(k, o["final"].get(k, 0.0), n["final"].get(k, 0.0)) for k in sk
             if abs(o["final"].get(k, 0.0) - n["final"].get(k, 0.0)) > 1e-9]
        print("  末态差异:")
        for k, a, b in d:
            print(f"     {k:<22} {a:>12.5g} -> {b:>12.5g}")
        print(f"  旧净方程: {str(o['net'])[:100]}")
        print(f"  新净方程: {str(n['net'])[:100]}")
    if out_path:
        with io.open(out_path, "w", encoding="utf-8") as fh:
            json.dump(res, fh, ensure_ascii=False, indent=1)
        print(f"\n-> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
