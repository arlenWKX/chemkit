# -*- coding: utf-8 -*-
"""第 263 轮 · 验收指标对拍（**只读档**）。

按 handoff §1.3 定型的主线指标：`resid_max` + 残差质量总和；
并列出零推进家族的残差贡献。用法：

    python tools/acc_metrics.py logs/suite-r262-1293.json logs/suite-latest.json
"""
from __future__ import annotations

import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def load(p: str) -> list[dict]:
    with io.open(p, encoding="utf-8") as fh:
        return json.load(fh)


def metrics(rows: list[dict]) -> dict:
    vals = [(abs(r.get("resid_live") or 0.0), r) for r in rows]
    tot = sum(v for v, _ in vals)
    nz = [r for v, r in vals if v > 0.1]
    return {
        "n": len(rows),
        "resid_max": max(v for v, _ in vals),
        "resid_max_name": max(vals, key=lambda t: t[0])[1].get("name"),
        "mass": tot,
        "n_gt_0.1": len(nz),
        "mass_gt_0.1": sum(abs(r.get("resid_live") or 0.0) for r in nz),
        "n_gt_1": sum(1 for r in nz if abs(r.get("resid_live") or 0.0) > 1.0),
        "n_gt_5": sum(1 for r in nz if abs(r.get("resid_live") or 0.0) > 5.0),
        "fails": sum(1 for r in rows if not r.get("ok")),
    }


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("-")]
    if len(args) < 2:
        print(__doc__)
        return 2
    A, B = load(args[0]), load(args[1])
    ma, mb = metrics(A), metrics(B)
    print(f"{'指标':<16} {'A = ' + os.path.basename(args[0]):>28} "
          f"{'B = ' + os.path.basename(args[1]):>28}   Δ(B−A)")
    for k in ("n", "fails", "resid_max", "mass", "n_gt_0.1", "mass_gt_0.1",
              "n_gt_1", "n_gt_5"):
        va, vb = ma[k], mb[k]
        d = (vb - va) if isinstance(va, (int, float)) else ""
        ds = f"{d:+.4f}" if isinstance(d, float) else f"{d:+d}"
        print(f"{k:<16} {va:>28} {vb:>28}   {ds}")
    print(f"\nA resid_max 例 = {ma['resid_max_name']}")
    print(f"B resid_max 例 = {mb['resid_max_name']}")
    # 逐例残差变化（只看放大方向）
    ia = {r.get("name"): abs(r.get("resid_live") or 0.0) for r in A}
    ib = {r.get("name"): abs(r.get("resid_live") or 0.0) for r in B}
    delta = sorted(((ib[k] - ia.get(k, 0.0), k) for k in ib),
                   key=lambda t: t[0])
    print("\n残差**减少**最多 8 例：")
    for d, k in delta[:8]:
        print(f"   {str(k)[:40]:<40} {ia.get(k, 0):>8.4f} -> {ib[k]:>8.4f} ({d:+.4f})")
    print("\n残差**增加**最多 8 例：")
    for d, k in delta[::-1][:8]:
        if d <= 0:
            break
        print(f"   {str(k)[:40]:<40} {ia.get(k, 0):>8.4f} -> {ib[k]:>8.4f} ({d:+.4f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
