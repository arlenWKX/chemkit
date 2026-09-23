# -*- coding: utf-8 -*-
"""第 203 轮 · F31 数据核查（v2，修正键名）：Ga 羟合梯 / Ga(OH)_3 Ksp。

v1 教训：`T.beta` 的 `ligand` 字段对羟合物**不是** `"OH^-"`（v1 因此打印
"无 OH^- 配合物条目"，那是过滤器错，不是数据缺），`T.ksp` 用 `pKsp`
而非 `logk`。诊断工具自己出错时最危险——它会伪装成"数据缺失"。

用法： python tools/f31_data.py
"""
import io
import json
import math
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

from chemkit.data import load_tables                        # noqa: E402

T = load_tables()


def _dump(tag, rows):
    print(f"  {tag}: {len(rows)} 条")
    for r in rows:
        print("    " + json.dumps(r, ensure_ascii=False))


print("=== ⓪ ligand 字段取值普查（v1 过滤器错在哪）===")
ligs = {}
for b in T.beta:
    ligs[b.get("ligand")] = ligs.get(b.get("ligand"), 0) + 1
for k, v in sorted(ligs.items(), key=lambda kv: -kv[1]):
    print(f"    {str(k):14s} {v}")

print("\n=== ① Ga 相关 beta 条目 ===")
_dump("beta 含 Ga", [b for b in T.beta if "Ga" in str(b.get("complex"))])
print("\n=== ② Ga 相关 ksp 条目 ===")
_dump("ksp 含 Ga", [e for e in T.ksp if "Ga" in str(e.get("solid"))
                    or "Ga" in str(e.get("pair"))])

print("\n=== ③ 羟合梯普查（complex 名含 (OH)）===")
hydro = [b for b in T.beta if "(OH)" in str(b.get("complex"))]
print(f"  共 {len(hydro)} 条；中心体 = "
      f"{sorted({b.get('center') for b in hydro})}")
for center in ("Ga^{3+}", "Al^{3+}", "Cr^{3+}", "Fe^{3+}"):
    rows = [b for b in hydro if b.get("center") == center]
    if rows:
        print(f"  --- {center} ---")
        for b in sorted(rows, key=lambda r: (r.get("m", 1), r.get("nu", 0))):
            print(f"    m={b.get('m')} nu={b.get('nu')} "
                  f"{b.get('complex'):26s} logb={b.get('logb')}")

print("\n=== ④ 氢氧化物固相普查（pair 含 OH）===")
ohk = [e for e in T.ksp if "OH" in str(e.get("pair"))]
print(f"  共 {len(ohk)} 条; solids = {sorted({e.get('solid') for e in ohk})}")
for e in ohk:
    if "Ga" in str(e.get("solid")) or "Ga" in str(e.get("pair")):
        print("   " + json.dumps(e, ensure_ascii=False))

print("\n=== ⑤ Ga(OH)_3 两性窗口（由引擎数据推出）===")
_G3 = "Ga^{3+}"
b4 = [b for b in hydro if b.get("center") == _G3 and b.get("nu") == 4
      and b.get("m", 1) == 1]
ks = [e for e in ohk if _G3 in str(e.get("pair"))]
if b4 and ks:
    lb4 = b4[0]["logb"]
    lksp = -ks[0]["pKsp"]      # pKsp -> logK
    lk = lksp + lb4
    print(f"  β₄([Ga(OH)₄]⁻) = {lb4}   (文献 298K ≈ 34.2)")
    print(f"  pKsp(Ga(OH)₃)  = {ks[0]['pKsp']}  ⟹ logKsp = {lksp}"
          f"   (文献 ≈ -35.5 ~ -37)")
    print(f"  Ga(OH)₃(s) + OH⁻ -> [Ga(OH)₄]⁻ : logK = logKsp + logβ₄ "
          f"= {lk:.2f}")
    print("  文献两性溶解 logK ≈ +0.3 ~ +1.3 (1 M OH⁻ 下 ~2~20 M)")
else:
    print(f"  缺条目：β₄={len(b4)} ksp={len(ks)}")

print("\n=== ⑥ 独立解 F31 平衡（1 M Ga, 3 M NaOH, 假设 Ga(OH)₃ 饱和）===")
if b4 and ks:
    lksp = -ks[0]["pKsp"]
    lb4 = b4[0]["logb"]
    best = None
    poh = 0.0
    while poh <= 8.0:
        oh = 10 ** (-poh)
        h = 10 ** (-(14.0 - poh))
        ga3 = 10 ** lksp / oh ** 3
        ga4 = 10 ** (lksp + lb4) * oh
        # 电荷：Na+3 - Cl-3 = 0 ⟹ 3[Ga3+] + [H+] = [Ga4-] + [OH-]
        chg = 3 * ga3 + h - ga4 - oh
        if best is None or abs(chg) < abs(best[1]):
            best = (poh, chg, ga3, ga4, h, oh)
        poh += 0.002
    poh, chg, ga3, ga4, h, oh = best
    print(f"  饱和线上电荷平衡点: pH = {14.0 - poh:.3f} (pOH {poh:.3f})")
    print(f"    [Ga3+]={ga3:.4g} M  [Ga(OH)4-]={ga4:.4g} M  "
          f"[H+]={h:.3g}  [OH-]={oh:.3g}")
    print(f"    残余电荷 = {chg:+.4g} M")
    print(f"    总 Ga = {ga3 + ga4:.4g} M（投料 1 M）")
    print(f"    Ga(OH)3(s) 析出量 = {1.0 - (ga3 + ga4):.4f} mol "
          f"（标准要求 ≥ 0.9）")
    print("\n  ⟹ 化学正解：pH ≈ 12（碱侧饱和），析出 ≈1.0 mol 固相。")
else:
    print("  数据缺失，无法核算")
