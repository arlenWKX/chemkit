# -*- coding: utf-8 -*-
"""第 265 轮 · 记录**某条候选**在一次 `solve_extent` 里被探到的每个 (pH, 分支) 点。

问题：`tools/xscan_ch.py` 显示 `I31` 的残差通道
（`[In(OH)₄]⁻ + H⁺ -> In(OH)₃`，d=−1）在 `x ∈ [0, 8.3e-5]` 上
`estimate_pH` **恒为 6.1679**、随后**跳到 9.9979**（跨 +3.83）。
但那是**默认 `pin_mode=None`（逐点自判定）**下的读数；
走步实际传的是 `_pin_mode`（两端一致才整步 ON，D14 第 201 轮）。
必须看**走步真正看到的那串点**，否则判据不同口径（第 253/263 轮教训）。

做法：包 `engine.solve_extent`（模块全局名调用，补丁有效），
只对匹配的候选记录其内部 `estimate_pH` 的每次读数与分支标签。

用法：python tools/solve_trace.py I31 "[In(OH)_4]^- + H^+" [最多记录点数]
"""
from __future__ import annotations

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

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.candidates import H_ION, WATER                 # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def main(argv: list[str]) -> int:
    pre = argv[0] if argv else "I31"
    needle = argv[1] if len(argv) > 1 else "[In(OH)_4]^- + H^+"
    cap = int(argv[2]) if len(argv) > 2 else 60

    T = load_tables()
    hit = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
           if n.startswith(pre)]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]

    _orig_est = eng.estimate_pH
    _orig_se = eng.solve_extent
    REC = []          # 本次 solve_extent 内的 (he, pH, tag, pin_mode)
    STATE = {"active": False}

    def logged_est(led, he, V, T_, T_K, *a, **k):
        tags = []
        if STATE["active"]:
            save = spec.PH_TAGS
            spec.PH_TAGS = tags
        try:
            r = _orig_est(led, he, V, T_, T_K, *a, **k)
        finally:
            if STATE["active"]:
                spec.PH_TAGS = save
        if STATE["active"]:
            REC.append((he, r, "/".join(tags),
                        k.get("pin_mode", "—")))
        return r

    OUT: list[str] = []
    LAST = {"rec": None}

    def looked_se(cd, d, led, He, V, T_, T_K, *a, **k):
        rr = cd.r if d > 0 else cd.pr
        txt = " + ".join(f"{'' if n == 1 else n}{s}" for s, n in rr.items())
        if needle and needle not in txt:
            return _orig_se(cd, d, led, He, V, T_, T_K, *a, **k)
        REC.clear()
        STATE["active"] = True
        try:
            out = _orig_se(cd, d, led, He, V, T_, T_K, *a, **k)
        finally:
            STATE["active"] = False
        LAST["rec"] = (cd, d, He, out, list(REC))
        return out

    eng.estimate_pH = logged_est
    eng.solve_extent = looked_se
    try:
        pr = {}
        r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
    finally:
        eng.estimate_pH = _orig_est
        eng.solve_extent = _orig_se

    cd, d, He0, out, rec = LAST["rec"]
    OUT.append(f"用例 {c['name']}: pH={r.get('final_pH')} "
               f"resid={pr.get('max_abs_S')}  He_final={pr.get('H_excess')}")
    OUT.append(f"\n=== **最后一次**匹配探针：d={d:+d}  He0={He0:.8g}  "
               f"x*={out[0]:.6g}  x_max={out[1]:.6g} ===")
    OUT.append(f"    候选 r={dict(cd.r)}  pr={dict(cd.pr)}  kind={cd.kind}")
    OUT.append(f"{'He':>15} {'pH':>9} {'pin_mode':>9}  {'分支':<22} S")
    for i, (he, ph, tag, pm) in enumerate(rec):
        s = ""
        if i + 1 < len(rec):
            s = ""
        OUT.append(f"{he:>15.6g} {ph:>9.4f} {str(pm):>9}  {tag:<22}")
    jumps = [(i, rec[i - 1][1], rec[i][1], rec[i][0]) for i in range(1, len(rec))
             if abs(rec[i][1] - rec[i - 1][1]) > 0.05]
    OUT.append(f"\n    共探 {len(rec)} 点；pH 跳变 {len(jumps)} 处：")
    for i, a_, b_, he_ in jumps[:10]:
        OUT.append(f"      第 {i} 点 He={he_:.6g}  {a_:.4f} -> {b_:.4f}  "
                   f"({b_ - a_:+.4f})")
    path = os.path.join(ROOT, "logs", f"solvetrace-{pre}.txt")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(OUT) + "\n")
    print(f"最后一次匹配探针已写 -> {os.path.relpath(path, ROOT)}")
    print(f"  探 {len(rec)} 点，pH 跳变 {len(jumps)} 处；"
          f"x*={out[0]:.6g} x_max={out[1]:.6g}")
    for i, a_, b_, he_ in jumps[:5]:
        print(f"    第 {i} 点 He={he_:.6g}  {a_:.4f} -> {b_:.4f} ({b_ - a_:+.4f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
