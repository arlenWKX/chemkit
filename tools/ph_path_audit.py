# -*- coding: utf-8 -*-
"""第 252 轮 · `exact_proton_pH` vs `estimate_pH` 的**全库影响面**。

**背景（第 251 轮）**：1:1 缓冲点上 `pH` 必然等于 `pKa(T)`，实测
`estimate_pH` 正确、而 `exact_proton_pH` 偏酸（`H_2O_2` −0.335、
酚@373K −0.374）。`presentation_pH` 在 `exact_proton_pH` 非 None 时**优先用它**
⟹ 偏差会偶合进 `final_pH`。本轮先量影响面，**不先改代码**
（纪律：先测量、再动手；改闸前先枚举新放行/新影响集合）。

做法：对**每个用例**跑一次 `judge(..., _probe={})`，用探针的
`ledger`/`H_excess` 分别调两条通路，统计 `|Δ|` 分布与最大者。
**只读，不改任何状态。**

用法： python tools/ph_path_audit.py [--write]
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


def main():
    from chemkit.data import load_tables
    from chemkit.engine import judge, estimate_pH, exact_proton_pH
    import chemkit.testsuit as ts

    T = load_tables()
    cases = ts.load_cases(None)
    rows = []
    nerr = 0
    for c in cases:
        cond = c.get("cond") or {"V_L": 1.0}
        tk = cond.get("T_K") or 298.15
        V = cond.get("V_L") or 1.0
        pr = {}
        try:
            judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T,
                  _probe=pr)
        except Exception:                                   # noqa: BLE001
            nerr += 1
            continue
        led = pr.get("ledger")
        if not led:
            continue
        he = pr.get("H_excess", 0.0)
        est = estimate_pH(led, he, V, T, tk)
        ex = exact_proton_pH(led, he, V, T, tk)
        if est is None or ex is None:
            rows.append((c["name"], tk, est, ex, None))
            continue
        rows.append((c["name"], tk, est, ex, ex - est))

    ok = [r for r in rows if r[4] is not None]
    print(f"用例 {len(cases)}；拿到两路 pH 的 {len(ok)}（异常 {nerr}）；"
          f"exact 返回 None 的 {len(rows) - len(ok)}")
    if not ok:
        return 1
    ds = sorted((abs(r[4]), r) for r in ok)
    n = len(ok)
    print("\n=== |Δ| = |exact − estimate| 分布 ===")
    for thr in (0.001, 0.01, 0.05, 0.1, 0.5, 1.0):
        k = sum(1 for a, _r in ds if a > thr)
        print(f"   > {thr:<6} : {k:5d} 例  ({100.0 * k / n:.1f}%)")
    print(f"   max    : {ds[-1][0]:.4f}")
    print("\n=== |Δ| 最大的 20 例 ===")
    print("  %-52s %-8s %-9s %-9s %s"
          % ("用例", "T/K", "estimate", "exact", "Δ"))
    for a, (nm, tk, est, ex, d) in reversed(ds[-20:]):
        print("  %-52s %-8.2f %-9.4f %-9.4f %+.4f" % (nm[:52], tk, est, ex, d))
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "ph_path_audit.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"name": r[0], "T": r[1], "estimate": r[2],
                        "exact": r[3], "d": r[4]} for r in rows], f,
                      ensure_ascii=False, indent=1)
        print(f"\n[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
