# -*- coding: utf-8 -*-
"""第 249 轮 · `logk_poly` 为什么一条也用不上（键口径不匹配的定量诊断）。

**背景**：`data/logk_poly.json` 有 **54 条** PHREEQC `-analytical_expression`
温度多项式，`candidates._poly_logk` **优先于** van't Hoff 使用它
（`logK_T` L461-468）。但 150 例抽样产生 448 个候选键，**命中 0 个** ⟹
这 54 条是**死数据**。

**本脚本要回答**：是"键口径不一致"（可修）还是"体系根本不相交"（修键也没用）？
做法：对每条存储条目，
  ① 用**运行时同一个** `_poly_norm` 归一它的 `reaction` 里各物种，造出"应有的键"；
  ② 看这个键是否出现在**真实运行**产生的候选键集合里；
  ③ 分类统计。

用法： python tools/poly_key_audit.py [--write]
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


def species_of_reaction(rxn: str) -> list[str]:
    """把 PHREEQC 风格反应式拆成物种名列表（去系数）。"""
    out = []
    for side in re.split(r"[=]", rxn):
        for term in side.split("+"):
            t = term.strip()
            if not t:
                continue
            t = re.sub(r"^\d+(\.\d+)?\s*", "", t)      # 去系数
            if t and t not in ("H2O", "e-"):
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
    tab = C._load_poly()
    print(f"存储条目 {len(entries)} 条；载入表 {len(tab)} 键\n")

    # ① 收集真实运行产生的候选键
    T = load_tables()
    seen = {}
    orig = C._poly_key

    def spy(r, pr):
        k = orig(r, pr)
        seen[k] = seen.get(k, 0) + 1
        return k
    C._poly_key = spy
    import random
    random.seed(7)
    cases = ts.load_cases(None)
    for cse in random.sample(cases, 200):
        try:
            judge([{"name": n, "mol": m} for n, m in cse["subs"]],
                  cse.get("cond") or {"V_L": 1.0}, T)
        except Exception:                                   # noqa: BLE001
            pass
    C._poly_key = orig
    print(f"真实运行产生候选键 {len(seen)} 个\n")

    # ② 对每条存储条目，按运行时的归一化重造键，看能否命中
    ok = 0
    rows = []
    for e in entries:
        sps = species_of_reaction(e["reaction"])
        rebuilt = tuple(sorted({C._poly_norm(s) for s in sps} - {"", "e"}))
        stored = tuple(e["key"])
        hit = rebuilt in seen
        if hit:
            ok += 1
        rows.append((hit, stored, rebuilt, e["reaction"][:44], e["src"]))
    print("  %-6s %-34s %-34s %s" % ("命中", "存储键", "按运行时重造的键", "反应式"))
    for hit, stored, rebuilt, rxn, src in rows[:24]:
        print("  %-6s %-34s %-34s %s" % ("★" if hit else "", str(stored)[:34],
                                         str(rebuilt)[:34], rxn))
    print(f"\n按运行时归一化后可命中的条目：{ok} / {len(entries)}")
    if ok == 0:
        print("\n⟹ 结论：**不只是键口径问题**——即使按运行时口径重造键也一条都不命中。")
        print("   说明这批多项式覆盖的体系（Amm.dat 的 HF/NaF/MgF/硫化物/石膏…）")
        print("   与套件实际走到的候选**没有交集**，或重造键仍与候选键不同构。")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "poly_key_audit.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"hit": h, "stored": list(s), "rebuilt": list(r),
                        "reaction": rx, "src": sc}
                       for h, s, r, rx, sc in rows], f,
                      ensure_ascii=False, indent=1)
        print(f"[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
