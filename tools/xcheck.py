# -*- coding: utf-8 -*-
"""第 185 轮：**命名映射 + 方向归一**的跨库核对（pKa 与 Ksp）。

映射规则（chemkit TeX 形 → PHREEQC 形）：
    ^{2-} -> -2 ; ^+ -> + ; ^- -> - ; _n -> n ; 去 [ ] ; Fe^{2+} -> Fe+2
方向归一：库反应若写 A- + H+ = HA，其 logK = +pKa；写 HA = A- + H+，logK = −pKa。
"""
import importlib.util
import io
import os
import re
import sys

sys.path.insert(0, ".")
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

ROOT = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "_phq", os.path.join(ROOT, "_phqparse2.py"))
pp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pp)

DB = os.path.join(ROOT, "thirdparty", "phreeqcrm-3.8.6-17100", "database")


def to_phq(name: str) -> str:
    """chemkit → PHREEQC。"""
    s = name
    s = re.sub(r"\^\{([^}]*)\}", r"\1", s)      # ^{2-} -> 2- ; ^+ -> +
    s = re.sub(r"_(\d+)", r"\1", s)             # H_2O -> H2O
    s = s.replace("[", "").replace("]", "")
    s = re.sub(r"([+-])$", lambda m: m.group(1), s)
    return s


def key(name: str):
    """归一化键：去电荷、去下划线与括号 ⟹ 用于跨库匹配。"""
    s = to_phq(name)
    return re.sub(r"[\+\-]", "", s)


def rx_species(rx: str):
    """('lhs', 'rhs') -> (物种键集合, 该物种是否在产物侧, 系数)"""
    parts = rx.split("=")
    if len(parts) != 2:
        return None
    out = {}
    for side, is_pr in ((parts[0], False), (parts[1], True)):
        for tok in side.split("+"):
            tok = tok.strip()
            if not tok:
                continue
            m = re.match(r"^([\d\.]+)\s*(.*)$", tok)
            coef, nm = (float(m.group(1)), m.group(2)) if m else (1.0, tok)
            k = key(nm)
            if k and k not in ("H2O", "H", "e"):
                out[k] = (is_pr, coef, to_phq(nm))
    return out


IDX = {}
for f in sorted(x for x in os.listdir(DB) if x.endswith(".dat")):
    d = pp.parse(os.path.join(DB, f))
    for sec in ("SOLUTION_SPECIES", "PHASES"):
        for nm, v in d.get(sec, {}).items():
            if v.get("logk") is None or not v.get("reaction"):
                continue
            sp = rx_species(v["reaction"])
            if not sp:
                continue
            IDX.setdefault(frozenset(sp), []).append(
                (f, nm, v["logk"], v["reaction"], sp))

from chemkit.data import load_tables                        # noqa: E402
T = load_tables()

print(f"索引 {len(IDX)} 组反应\n=== pKa 核对（chemkit 有对应反应者）===")
print(f"{'酸/碱':<24}{'chemkit':>8}  库内（库:logK→归一后 pKa）")
hit = 0
for e in T.pka:
    ka, kb = key(e["acid"]), key(e["base"])
    if ka == kb:
        continue
    got = IDX.get(frozenset({ka, kb}))
    if not got:
        continue
    hit += 1
    cells = []
    for f, nm, lk, rx, sp in got[:3]:
        # 方向：若酸在反应物侧，则 logK = −pKa
        acid_is_reactant = not sp[ka][0]
        pka_equiv = -lk if acid_is_reactant else lk
        cells.append(f"{f.split('.')[0]}:{pka_equiv:.2f}")
    print(f"{(e['acid'] + '/' + e['base']):<24}{e['pka']:>8.2f}  " + "; ".join(cells))
print(f"\n（pKa 命中 {hit}/{len(T.pka)}）")

print("\n=== Ksp 核对（阴离子型：固体 ⇌ 离子）===")
hit = 0
for e in T.ksp:
    solid, cat, an = e["solid"], e["pair"][0], e["pair"][1]
    ks = key(solid)
    got = IDX.get(frozenset({ks})) if False else None
    # 用"含 solid 键的反应"反查
    cand = []
    for k, lst in IDX.items():
        if ks in k:
            cand.extend(lst)
    if not cand:
        continue
    hit += 1
    if hit <= 22:
        f, nm, lk, rx, sp = cand[0]
        print(f"{solid:<14}{e['pKsp']:>7.2f}  {f.split('.')[0]}:{nm} logK={lk:g}  "
              f"rx={rx[:44]}")
print(f"\n（Ksp 命中 {hit}/{len(T.ksp)}）")
