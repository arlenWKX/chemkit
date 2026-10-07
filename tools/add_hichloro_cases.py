# -*- coding: utf-8 -*-
"""第 292 轮 · **浓氯中高阶氯合形态锚（Cd/Zn 对照）**。

## 动因

第 292 轮查明"高阶氯络合物族"10 例分两类：
① Cu 6 例（39/J12/O03/RX1/Y14/U10）——走步与账本**本就正确**（与库值
   独立不动点推导逐位吻合），失败只是弱池折叠呈现口径 vs 旧标总方程式
   期望（使用者裁定：产物复杂可不查总方程式 ⟹ 重裁为形态分布锚）；
② Cd/Pb/Sn(IV) 4 例——**数据缺口**：库内只有一氯条目。本轮把测试 note
   里**已在案**的两个值补入库（Cd β₄=2.0；Pb β₃=2.0 Luo 2007），
   I16/U06/DB11 重裁结案；D05（[SnCl₆]²⁻ 无在案出处）留作数据缺口。

本组锚锁**新入库的高阶氯合数据真的在驱动形态分布**（与引擎无关的推导）：

* `CD3`：Cd²⁺ 在 6 M Cl⁻ 中 **β₄ 主导**（[CdCl₄]²⁻ ≈ 99%）——
  logK₄_c=2.0−4D(I_eff≈3.1)≈1.02 ⟹ x₄/z=10.5×6⁴≈1.4e4；
* `ZN6`：Zn²⁺ 在 6 M Cl⁻ 中 **β₃ 主导**（[ZnCl₃]⁻ ≈ 72%）——
  Zn 没有 β₄ 条目（IUPAC 只到 ν=3），同浓度下分布顶点停在三氯。
  两条构成"同条件、不同金属、分布顶点不同"的对照，盯 β 阶梯数据的
  **整体形状**（任何一阶 logβ 漂移都会移动顶点）。

## 用法

    python tools/add_hichloro_cases.py            # dry-run
    python tools/add_hichloro_cases.py --check
    python tools/add_hichloro_cases.py --write
"""
from __future__ import annotations

import io
import json
import math
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

_TAG = "第 292 轮扩充 · 高阶氯合形态锚。"

# 库值阶梯（chemkit/data/beta.json）：(络合物, logβ°, Δz²累积, ε, z²(络合物), 结合Cl数)
CD_LADDER = [
    ("[CdCl]^+",      1.98, -4, 0.15, 1, 1),
    ("[CdCl_4]^{2-}", 2.00, -4, None, 4, 4),
]
ZN_LADDER = [
    ("[ZnCl]^+",   0.40, -4, -0.14, 1, 1),
    ("[ZnCl_2]",   0.69, -6, -0.20, 0, 2),
    ("[ZnCl_3]^-", 0.18, -6, -0.05, 1, 3),
]


def D(I: float) -> float:
    s = math.sqrt(I)
    return 0.51 * s / (1 + 1.5 * s)


def solve(m_tot: float, cl_tot: float, ladder):
    """累积 β 的 SIT-I_eff 不动点（与 tools/cucl_spec.py 同口径）。"""
    c = cl_tot * 0.7
    fr = [0.0] * len(ladder)
    for _ in range(500):
        z = m_tot / (1 + sum(fr))
        new = []
        for k, (name, lb, dz2, eps, z2, ncl) in enumerate(ladder):
            I = 0.5 * (4 * z + c + z2 * z * fr[k])
            lg = lb + dz2 * D(I) + (eps * I if eps is not None else 0.0)
            new.append((10 ** lg) * c ** ncl)
        fr = [0.5 * a + 0.5 * b for a, b in zip(fr, new)]
        z = m_tot / (1 + sum(fr))
        xs = [z * f for f in fr]
        c_new = cl_tot - sum(x * e[5] for x, e in zip(xs, ladder))
        if abs(c_new - c) < 1e-12:
            break
        c = 0.5 * c + 0.5 * c_new
    z = m_tot / (1 + sum(fr))
    return z, [z * f for f in fr], c


def build() -> list[dict]:
    z_cd, xs_cd, c_cd = solve(0.05, 6.1, CD_LADDER)   # NaCl 6 + CdCl2 自带 0.1
    z_zn, xs_zn, c_zn = solve(0.10, 6.2, ZN_LADDER)
    return [
        {
            "name": "CD3 CdCl2 0.05+NaCl 6 四氯主导锚（β4 数据驱动）",
            "subs": [["CdCl_2", 0.05], ["NaCl", 6.0]],
            "has": {"[CdCl_4]^{2-}": 0.04},
            "has_range": {"Cd^{2+}": [0.0, 0.0001], "[CdCl]^+": [0.0, 0.003]},
            "note": (
                _TAG + f"**库值独立推导**：logK₄_c=2.0−4·D(I_eff≈3.1)=1.02、"
                f"logK₁_c=1.98−4·D+0.15·I≈1.46 ⟹ x₄:x₁:z≈1.4e4:174:1 ⟹ "
                f"[CdCl₄]²⁻≈{xs_cd[1]:.4f}（~99%）、[CdCl]⁺≈{xs_cd[0]:.2e}、"
                f"游离 Cd²⁺≈{z_cd:.1e}（不动点求解见本文件 `solve`，"
                f"与引擎无关）。"),
        },
        {
            "name": "ZN6 ZnCl2 0.1+NaCl 6 三氯主导锚（与 CD3 对照）",
            "subs": [["ZnCl_2", 0.1], ["NaCl", 6.0]],
            "has": {"[ZnCl_3]^-": 0.05},
            "has_range": {"Zn^{2+}": [0.004, 0.02]},
            "note": (
                _TAG + f"**与 CD3 同条件对照**：Zn 库内只有 ν=1,2,3（IUPAC "
                f"Part 5），β₃°=0.18 ⟹ logK₃_c=0.18−6·D−0.05·I≈−1.45 ⟹ "
                f"x₃/z=0.035×6³≈7.6，分布 1:0.59:1.4:7.6 ⟹ "
                f"[ZnCl₃]⁻≈{xs_zn[2]:.3f}（~72%）、游离 Zn²⁺≈{z_zn:.4f}。"
                f"对照意义：同浓度下 Cd 顶点在 ν=4、Zn 顶点在 ν=3 —— "
                f"锚住 β 阶梯的**整体形状**。"),
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
        print(f"{c['name'][:50]:<50}  {'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<50} {'pH':>7} {'resid':>8} {'判读':>4}")
        for c in new:
            pr: dict = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            amt: dict = {}
            for e in (r.get("production") or []) + (r.get("final") or []):
                amt[e["name"]] = max(amt.get(e["name"], 0.0), e["mol"])
            ok = True
            for sp, lo in (c.get("has") or {}).items():
                ok = ok and amt.get(sp, 0.0) >= lo
            for sp, (lo, hi) in (c.get("has_range") or {}).items():
                ok = ok and lo <= amt.get(sp, 0.0) <= hi
            ph = r.get("final_pH")
            print(f"{c['name'][:50]:<50} "
                  f"{(ph if ph is not None else float('nan')):>7.3f} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>4}")
            for sp in list((c.get("has") or {})) + list((c.get("has_range") or {})):
                print(f"      {sp:<22} {amt.get(sp, 0.0):.5g}")
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
