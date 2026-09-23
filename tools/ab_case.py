# -*- coding: utf-8 -*-
"""第 263 轮 · 逐例 A/B 诊断：新路径开/关的走步对比。

用 monkey-patch 关掉**新增的那一条**（两性支路的精确解，`min_fams=1`），
保留 `_deg` 原路径 ⟹ 同一进程内得到"改前/改后"，可逐例对拍。

用法：
    python tools/ab_case.py N15
    python tools/ab_case.py N15 Z15 Q04 TC1 --json logs/ab263.json
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

from chemkit import speciation as _sp            # noqa: E402
from chemkit.data import load_tables             # noqa: E402

_ORIG = _sp._exact_ok


def _off(ledger, He_res, T, V):
    """关掉精确解接管（⚠️ 第 263 轮收窄后 `_exact_ok` 被**两处**共用，
    故这是"两个接管点一起关"，只作粗探针；权威 A/B 用两份套件档对拍
    `tools/acc_metrics.py` / `tools/suite_show.py diff`）。"""
    return False


def run(case: dict, T, on: bool):
    _sp._exact_ok = _ORIG if on else _off
    from chemkit.converg import _live
    from chemkit.engine import judge
    pr = {}
    subs = [{"name": s[0], "mol": float(s[1])} for s in case["subs"]]
    r = judge(subs, case.get("cond") or {}, T, _probe=pr)
    led = {e["name"]: e["mol"] for e in (r.get("final") or [])
           if e["mol"] > 1e-9}
    # ⚠️ `resid_live` **不在** judge 的 probe 里：它是
    # `testsuit._enrich_results_with_probe` 用 `converg._live(active)` 算的
    # （质量口径，"值得解且 walk 会解的两侧平衡最大 |S|"）。第 263 轮初版
    # 本工具读 `pr["resid_live"]` 恒得 0，白跑一次排查。
    return {
        "pH": r.get("final_pH"),
        "pH_solver": pr.get("pH_solver"),
        "resid": abs(_live(pr.get("active") or [])),
        "iters": pr.get("iters"), "exit": pr.get("exit"),
        "steps": [(s["equation"], round(s["extent"], 6))
                  for s in (r.get("steps") or [])],
        "final": {k: round(v, 8) for k, v in sorted(led.items())},
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
        off = run(c, T, on=False)
        on = run(c, T, on=True)
        res[c["name"]] = {"off": off, "on": on}
        print(f"\n===== {c['name']} =====")
        print(f"  改前: pH={off['pH']} solver={off['pH_solver']} "
              f"resid={off['resid']:.4f} iters={off['iters']} exit={off['exit']}")
        print(f"  改后: pH={on['pH']} solver={on['pH_solver']} "
              f"resid={on['resid']:.4f} iters={on['iters']} exit={on['exit']}")
        sk = {k for k in off["final"]} | {k for k in on["final"]}
        diff = [(k, off["final"].get(k, 0.0), on["final"].get(k, 0.0))
                for k in sorted(sk)
                if abs(off["final"].get(k, 0.0) - on["final"].get(k, 0.0)) > 1e-9]
        print("  末态差异:")
        for k, a, b in diff:
            print(f"     {k:<20} {a:>12.3e} -> {b:>12.3e}")
        print("  改前步骤:")
        for e, x in off["steps"]:
            print(f"     {e[:78]:<78} x={x:.6g}")
        print("  改后步骤:")
        for e, x in on["steps"]:
            print(f"     {e[:78]:<78} x={x:.6g}")
    _sp._exact_ok = _ORIG
    if out_path:
        with io.open(out_path, "w", encoding="utf-8") as fh:
            json.dump(res, fh, ensure_ascii=False, indent=1)
        print(f"\n-> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
