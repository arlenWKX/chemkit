# -*- coding: utf-8 -*-
"""第 271 轮 · **卤素置换的方向与计量锚**（本轮定位的 `resid_max` 榜首体系）。

## 动因

第 271 轮定位：当前 `resid_max` 榜首是**三例同一化学**的氧化还原体系 ——
`36 Cl2+NaBr` / `N01 Cl2+KBr` / `H1 Cl2+NaBr 置换`（各 **9.434**），
残差通道 `3Cl₂ + Br⁻ -> 6Cl⁻ + 6H⁺ + BrO₃⁻`（S = −9.434，
反向驱动 +9.434，`x*` = 5.2e-7 而 `x_max` = 1.7e-3 ⟹ 零推进）。

本轮为该体系补**方向与计量**锚，断言只用两条化学事实：

* **氧化性强弱**：`E°(Cl₂/Cl⁻) = +1.36 V > E°(Br₂/Br⁻) = +1.08 V`
  ⟹ **Cl₂ 能把 Br⁻ 置换出来，Br₂ 不能把 Cl⁻ 置换出来**（方向是**单向**的）；
* **计量**：`Cl₂ + 2Br⁻ -> Br₂ + 2Cl⁻` ⟹ `n(Br₂) = n(Br⁻)/2`（质量守恒）。

两条都与引擎无关，且**可判**（定量下限 + 反向的"不得生成"）。

## 用法

    python tools/add_halogen_cases.py            # dry-run
    python tools/add_halogen_cases.py --check
    python tools/add_halogen_cases.py --write
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

E_CL2 = "+1.36 V"
E_BR2 = "+1.08 V"


def build() -> list[dict]:
    return [
        {
            "name": "BR1 Cl2+KBr 2 等计量置换（卤素方向锚）",
            "subs": [["Cl_2", 1], ["KBr", 2]],
            "has": {"Br_2": 0.9},
            "has_not": {"Cl_2": 0.01},
            "note": (
                f"第 271 轮扩充 · **卤素置换的方向锚**。"
                f"`E°(Cl₂/Cl⁻) = {E_CL2} > E°(Br₂/Br⁻) = {E_BR2}` ⟹ "
                f"Cl₂ 能把 Br⁻ **定量**置换为 Br₂（1:2 计量，"
                f"`Cl₂ + 2Br⁻ -> Br₂ + 2Cl⁻`）⟹ 断言取下限 0.9，"
                f"并断言 Cl₂ 基本耗尽（上限 0.01）。两条都只用氧化性强弱"
                f"与质量守恒，**不含引擎输出**。"
                f"本体系的既有用例（`N01`/`H1`/`36`）正是当前 `resid_max` 榜首"
                f"（9.434）——本条锁住**产物侧**：无论内部走步如何，"
                f"Br₂ 必须定量生成。"),
        },
        {
            "name": "BR2 Cl2 过量+KBr 0.5 半量（卤素计量锚）",
            "subs": [["Cl_2", 1], ["KBr", 0.5]],
            "has": {"Br_2": 0.225},
            "note": (
                f"第 271 轮扩充 · 同族的**计量侧**：KBr 是限量者（0.5 mol），"
                f"按 `Cl₂ + 2Br⁻ -> Br₂ + 2Cl⁻` ⟹ `n(Br₂) = 0.5/2 = 0.25`"
                f"（断言取下限 0.225 = 90%）。Cl₂ 过量 ⟹ 不断言 Cl₂ 耗尽。"
                f"**只锁 Br₂ 的产量**（质量守恒），与 `BR1`（等计量）配对"
                f"锁住「谁限量」这一维。"),
        },
        {
            "name": "BR3 Br2+KCl 2 逆方向不置换（卤素方向锚·负）",
            "subs": [["Br_2", 1], ["KCl", 2]],
            "has_not": {"Cl_2": 0.01},
            "has": {"Br_2": 0.9},
            "note": (
                f"第 271 轮扩充 · **卤素置换的负方向锚**："
                f"`E°(Br₂/Br⁻) = {E_BR2} < E°(Cl₂/Cl⁻) = {E_CL2}` ⟹ "
                f"**Br₂ 不能把 Cl⁻ 氧化成 Cl₂** ⟹ 断言 `Cl₂ < 0.01`（上限）。"
                f"同时 Br₂ 应基本保留（≥0.9）——"
                f"两条都是**化学事实**（氧化性强弱 + 质量守恒）。"
                f"与 `BR1` 配对把「方向」这一维在**两侧**都锁住："
                f"正向必须发生、反向必须不发生。"),
        },
    ]


def main(argv: list[str]) -> int:
    from chemkit.data import load_tables
    T = load_tables()
    new = build()
    path = os.path.join(ROOT, "chemkit", "data", "tests.json")
    with io.open(path, encoding="utf-8") as fh:
        db = json.load(fh)
    have = {c["name"] for c in db}
    for c in new:
        spec = c.get("has") or {}
        ban = c.get("has_not") or {}
        print(f"{c['name'][:48]:<48} has={str(spec)[:20]:>20} "
              f"not={str(ban)[:16]:>16}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<48} {'引擎末态':>26} {'resid':>8} {'判读':>6}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            fin = {e["name"]: e["mol"] for e in (r.get("final") or [])}
            want = c.get("has") or {}
            ban = c.get("has_not") or {}
            keys = list(want) + list(ban)
            got = " ".join(f"{k}={fin.get(k, 0.0):.4g}" for k in keys)
            ok = (all(fin.get(k, 0.0) >= v for k, v in want.items())
                  and all(fin.get(k, 0.0) < v for k, v in ban.items()))
            print(f"{c['name'][:48]:<48} {got[:26]:>26} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>6}  pH={r.get('final_pH')}")
        return 0

    if "--write" not in argv:
        print(f"\n[dry-run] 将新增 {sum(1 for c in new if c['name'] not in have)} 条"
              f"；加 --write 写入。")
        return 0
    add = [c for c in new if c["name"] not in have]
    db.extend(add)
    with io.open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(db, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print(f"\n[写入] 新增 {len(add)} 条 ⟹ 共 {len(db)} 条")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
