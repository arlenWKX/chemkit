# -*- coding: utf-8 -*-
"""第 269 轮 · **羟合梯缺口全库普查**。

## 背景

第 268 轮把 `F31 GaCl3+3NaOH`（`resid_max` 13.305）归因为**数据缺口**：
库内 `Ga³⁺` 的 OH 梯只有 ν=1 与 ν=4（缺 ν=2,3）⟹ 梯子**断成两段**
⟹ `build_families`（只读 `T.pka`）连不成一条可再分配的族
⟹ `_exact_ok` 第②条（族成员 ≥2）恒挡 ⟹ 精确解永不启用。

本轮先量**这个缺口的全库规模**：哪些阳离子的 OH 梯有 ν 空洞、
哪些是连通的、各出现在多少用例的账本里。

## 判据

* 数据源：`T.beta` 中 `ligand == 'OH^-'` 且 `m == 1`（单核）的条目。
* 「空洞」= 已有序号中缺失的整数（`1..max(ν)` 内）。
* 「连通」= ν 集合是 `1..max` 的连续段（子梯之间没有断口）。
* 用例关联：跑一遍全库，统计每个中心的梯物种出现在多少例的末态账本里。

**零侵入**（只读表 + 跑 `judge`），不改 `chemkit/` 一行。
输出 `logs/ladder_gap_census.json`。

用法：python tools/ladder_gap_census.py
"""
from __future__ import annotations

import collections
import io
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.data import load_tables                        # noqa: E402


def main() -> int:
    T = load_tables()
    # 中心 -> [(nu, logb, complex)]
    lad: dict[str, list] = {}
    for b in T.beta:
        if b.get("ligand") != "OH^-" or b.get("m", 1) != 1:
            continue
        lad.setdefault(b["center"], []).append(
            (b.get("nu"), b.get("logb"), b["complex"]))
    # Ksp-OH 阳离子（有固相储库的）
    ksp_cat = {e["pair"][0]: e["solid"] for e in T.ksp
               if e["pair"][1] == "OH^-"}

    rows = []
    for cat, ent in sorted(lad.items()):
        ent.sort(key=lambda t: (t[0] or 0))
        nus = sorted({n for n, _l, _c in ent if n})
        holes = [k for k in range(1, (max(nus) if nus else 0) + 1)
                 if k not in nus]
        rows.append({
            "center": cat, "nus": nus, "holes": holes,
            "n": len(nus),
            "connected": not holes,
            "solid": ksp_cat.get(cat),
            "complexes": [c for _n, _l, c in ent],
        })

    print(f"=== 单核 OH 梯：{len(rows)} 个中心 ===")
    disc = [r for r in rows if r["holes"]]
    conn = [r for r in rows if not r["holes"]]
    print(f"  连通（无空洞）：**{len(conn)}**；有空洞（断成多段）：**{len(disc)}**")
    print(f"\n{'中心':<16} {'ν 集合':<16} {'空洞':<10} {'固相':<12} 复杂度")
    for r in rows:
        nu = ",".join(str(x) for x in r["nus"])
        ho = ",".join(str(x) for x in r["holes"]) or "—"
        print(f"{r['center'][:16]:<16} {nu[:16]:<16} {ho[:10]:<10} "
              f"{str(r['solid'] or '—')[:12]:<12} "
              f"{' / '.join(c[:22] for c in r['complexes'])[:56]}")

    # 用例关联
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    from chemkit.engine import judge
    hit = collections.Counter()
    t0 = time.time()
    for i, c in enumerate(cases, 1):
        pr = {}
        try:
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
        except Exception:                                   # noqa: BLE001
            continue
        led = {e["name"]: e["mol"] for e in (r.get("final") or [])}
        for row in rows:
            if any(cx in led and led[cx] > 1e-9 for cx in row["complexes"]):
                hit[row["center"]] += 1
        if i % 300 == 0:
            print(f"  {i}/{len(cases)} … {time.time() - t0:.0f}s", flush=True)

    print("\n=== 梯物种出现在多少例的末态账本里 ===")
    for row in sorted(rows, key=lambda r: -hit[r["center"]]):
        h = hit[row["center"]]
        if not h:
            continue
        flag = "  ← **有空洞**" if row["holes"] else ""
        print(f"   {row['center'][:18]:<18} {h:>4} 例   ν={row['nus']}"
              f"{flag}")
    out = {"rows": rows, "case_hits": dict(hit),
           "n_centers": len(rows), "n_connected": len(conn),
           "n_disconnected": len(disc),
           "wall_s": round(time.time() - t0, 1)}
    with io.open(os.path.join(ROOT, "logs", "ladder_gap_census.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print("\n-> logs/ladder_gap_census.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
