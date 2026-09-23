# -*- coding: utf-8 -*-
"""第 223 轮 · Pitzer 数据完备性核查（决定实现是否可行）。

`tools/pitzer_survey.py` 已给出结构与 kind 分布，本脚本只回答
**实现可行性**的两个硬指标：

  ① **Cφ 有没有？** 标准 Pitzer 方程
        φ − 1 = (2/(Σm))·[ −A·I^1.5/(1+1.2√I) + ΣΣ m_i m_j (B_ij + Z·C_ij) + … ]
     中 `C_ij = Cφ_ij / (2√(|z_i z_j|))` 是**必需项**（高 I 下不可省）。
     `_kinds` 里列了 `cphi`，但需核实**条目数是否 > 0**。
  ② **离子名能否全部映射到库内物种？**（含畸形名如 `+` / `=`）

用法： python tools/pitzer_check.py
"""
import collections
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

from chemkit.data import load_tables                        # noqa: E402

T = load_tables()
pj = json.load(io.open(os.path.join(ROOT, "chemkit", "data",
                                    "pitzer.json"), encoding="utf-8"))
ents = pj["entries"]
kinds = collections.Counter(e.get("kind") for e in ents)

print("=== ① Cφ（cphi）完备性 —— 决定标准式能否成立 ===")
print(f"  kind 分布: {dict(kinds)}")
print(f"  **cphi 条目数 = {kinds.get('cphi', 0)}**")
if kinds.get("cphi", 0) == 0:
    print("  ⟹ **Cφ 完全缺失**。标准 Pitzer 的 `+Z·C_ij` 项无法计算，")
    print("     只能做**截断到 β⁰/β¹(+β²)** 的简化式：")
    print("     lnγ± ≈ −A√I/(1+1.2√I) + Σ_j m_j (2B_Mj)  （缺 C 项）")
    print("     该项在高 I（>3）下贡献显著，**省略会系统性低估活度系数**。")

print("\n=== ② 离子名映射 ===")
allions = collections.Counter()
for e in ents:
    for nm in (e.get("ions") or []):
        allions[nm] += 1
print(f"  不同离子名 {len(allions)} 个；出现最多 12: {allions.most_common(12)}")
bad = [n for n in allions if n in ("+", "-", "=", "") or len(n) <= 1]
print(f"  **畸形名**（长度≤1 或纯符号）: {bad}")

# chemkit 端物种名集合
cands = set()
for b in T.beta:
    cands.add(b["center"])
    cands.add(b["complex"])
cands |= set(T.pka_acid) | set(T.pka_base)
cands |= {e["pair"][0] for e in T.ksp} | {e["pair"][1] for e in T.ksp}
cands |= set(T.solids)

# 现成映射线索：tools/xcheck.py 是否已有 PHREEQC->chemkit 映射？
xp = os.path.join(ROOT, "tools", "xcheck.py")
txt = io.open(xp, encoding="utf-8").read() if os.path.exists(xp) else ""
print(f"  tools/xcheck.py {'存在' if txt else '不存在'}"
      f"{'，含映射表' if ('MAP' in txt or 'map' in txt) else '，未见映射表'}")
print(f"  直接同名命中: {sum(1 for n in allions if n in cands)}/{len(allions)}")

print("\n=== 判读与建议 ===")
if kinds.get("cphi", 0) == 0:
    print("  · Cφ 缺失 ⟹ 只能实现**简化 Pitzer**；必须：")
    print("    ① 在模块头与数据 note 里**显式写明这是截断式**及其高 I 偏差方向；")
    print("    ② 验收锚点不能要求高 I 精度（拿不到 Cφ 就达不到）；")
    print("    ③ 或在 pitzer.dat 之外补 Cφ（那是**数据工作**，且要权威源）。")
print("  · 名称映射需要一层（PHREEQC 记法 `Cl-`/`SO4-2` vs chemkit "
      "`Cl^-`/`SO_4^{2-}`），可复用 tools/xcheck.py 的既有映射思路。")
