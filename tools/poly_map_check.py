# -*- coding: utf-8 -*-
"""第 250 轮 · `logk_poly` 映射层重造前的**决定性测量**：交集到底有没有。

第 249 轮已查明三件事：接线在（`_poly_logk` 优先于 van't Hoff）、
命中率 0、以及三个独立成因（`H+` 未归一化 / 固相键不同构 / 体系不相交）。
但当时"体系不相交"是**按物种名**估的（18.3%），不够决定性 ——
**真正该问的是：按运行时口径重造键之后，有没有候选能对上？**

本脚本：
  1. 跑**全部**用例，收集真实运行产生的候选键（`_poly_key` 的实参）；
  2. 对 54 条存储条目，按运行时 `_poly_norm` 重造键
     （**剔掉 `H+`/`H2O`** —— 这正是第 249 轮查出的漏项）；
  3. 精确回答"重造后有 n 条能命中、分别影响哪些用例"。

用法： python tools/poly_map_check.py [--write]
"""
import io
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# PHREEQC 写法 -> 我们 _poly_norm 之后的形态（用于剔除溶剂/质子）
SOLVENT_AFTER_NORM = {"h2o", "h", "e"}


def species_of(rxn: str) -> list[str]:
    out = []
    for side in rxn.split("="):
        for t in side.split("+"):
            t = re.sub(r"^\d+(\.\d+)?\s*", "", t.strip())
            if t:
                out.append(t)
    return out


def main():
    import chemkit.candidates as C
    from chemkit.data import load_tables
    from chemkit.engine import judge
    import chemkit.testsuit as ts

    d = json.load(io.open(os.path.join(ROOT, "chemkit/data/logk_poly.json"),
                          encoding="utf-8"))
    entries = d["entries"]

    # ① 重造键（剔溶剂/质子）
    rebuilt = {}
    for e in entries:
        sps = species_of(e["reaction"])
        ks = []
        for s in sps:
            n = C._poly_norm(s)
            if n and n not in SOLVENT_AFTER_NORM:
                ks.append(n)
        key = tuple(sorted(set(ks)))
        if key:
            rebuilt.setdefault(key, []).append(e)
    print(f"54 条里重造出**非空键**的：{len(rebuilt)} 种键 "
          f"（其余 {len(entries) - sum(len(v) for v in rebuilt.values())} 条"
          f"是纯固相/单物种，归不进 Cand 键）\n")

    # ② 跑全部用例，收集候选键 -> 用例名
    T = load_tables()
    seen = {}
    orig = C._poly_key

    def spy(r, pr):
        k = orig(r, pr)
        seen.setdefault(k, 0)
        seen[k] += 1
        return k

    C._poly_key = spy
    cases = ts.load_cases(None)
    nerr = 0
    for cse in cases:
        try:
            judge([{"name": n, "mol": m} for n, m in cse["subs"]],
                  cse.get("cond") or {"V_L": 1.0}, T)
        except Exception:                                   # noqa: BLE001
            nerr += 1
    C._poly_key = orig
    print(f"跑完全部 {len(cases)} 例（异常 {nerr}）；"
          f"真实运行产生的候选键 {len(seen)} 个\n")

    # ③ 精确求交
    hit = {k: v for k, v in rebuilt.items() if k in seen}
    print(f"**重造键 ∩ 真实候选键 = {len(hit)} 个**")
    for k, es in hit.items():
        print(f"  键 {k}  出现 {seen[k]} 次")
        for e in es:
            print(f"     <- {e['reaction']}  (logk298={e['logk298']}, "
                  f"src={e['src']})")
    if not hit:
        print("\n⟹ **结论（决定性）**：即使修好 `H+` 归一化漏项、按运行时口径"
              "重造键，仍然**一条都对不上**。")
        print("   说明这 54 条覆盖的体系确实不在套件的候选空间里 ⟹")
        print("   **不是映射层能修的**，映射层重造这条路线应当**关闭**。")
        print("   若将来要利用它们，正确做法是**扩大套件的化学体系**")
        print("   （HF/氟化物、硫化物、石膏/重晶石族），而不是修键。")

    # ④ 顺便报告：这些体系里的关键物种，库内到底有没有
    print("\n--- 附：那些体系的物种在库内是否存在（归一化后）---")
    uni = set()
    for k in T.thermo:
        uni.add(C._poly_norm(k))
    for e in T.beta:
        for k in ("complex", "center", "ligand"):
            uni.add(C._poly_norm(e[k]))
    for e in T.pka:
        uni.add(C._poly_norm(e["acid"]))
        uni.add(C._poly_norm(e["base"]))
    allsp = set()
    for e in entries:
        for s in species_of(e["reaction"]):
            n = C._poly_norm(s)
            if n and n not in SOLVENT_AFTER_NORM:
                allsp.add(n)
    inn = sorted(s for s in allsp if s in uni)
    out = sorted(s for s in allsp if s not in uni)
    print(f"  多项式条目涉及的物种 {len(allsp)} 个；库内有 {len(inn)}、"
          f"库内无 {len(out)}")
    print(f"  库内有的: {inn}")
    print(f"  库内无的（前 25）: {out[:25]}")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "poly_map_check.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump({"rebuilt_keys": [list(k) for k in rebuilt],
                       "n_candidate_keys": len(seen),
                       "hits": [list(k) for k in hit],
                       "species_in_lib": inn, "species_missing": out},
                      f, ensure_ascii=False, indent=1)
        print(f"\n[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
