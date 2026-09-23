# -*- coding: utf-8 -*-
"""第 263 轮 · 两性缓冲对锚（**答案由引擎无关的闭式解给出**）。

## 动因

第 263 轮查出 `estimate_pH` 分支 4 的「两性中点」捷径
（`pH = (pKa1+pKa2)/2`）在**两性物种与其共轭伙伴同时在账**时前提不成立
（那是缓冲体系，不是纯两性盐）。修法已落地（`speciation.estimate_state`
的两性支路改为先问精确质子条件）。本脚本为**该改动**补锚。

## 答案来源（不依赖引擎任何 pH 通路）

二元酸 H₂A，总量 `C`（mol）、外加钠 `n_Na`（mol）、体积 `V`：

    D = h² + Ka1·h + Ka1·Ka2
    [H₂A] = C·h²/D ; [HA⁻] = C·Ka1·h/D ; [A²⁻] = C·Ka1·Ka2/D
    电荷：n_Na/V + h = [HA⁻] + 2[A²⁻] + oh          （无外加游离强酸）

`f(pH)` 严格单调递减 ⟹ 唯一根。判据函数 = `tools/amph_exact.py::exact_pH`
（同一份实现，含正/负对照 `--calib`）。

## 覆盖设计

* **比例维**：把缓冲比从 0.3 扫到 0.7（1:1 已有 TC3，故取非 1:1 档）；
* **体系维**：硫系（H₂S，pKa 7.0/14.0，两性物种 HS⁻）与
  碳系（CO₂，pKa 6.4/10.3，两性物种 HCO₃⁻）——两者 pKa2−pKa1 差得远
  （7.0 vs 3.9），能把"中点式有多偏"的梯度拉开；
* **温度维**：各取一个非常温档（dH 不同，pKa1/pKa2 的相对移动不同）。

## 用法

    python tools/add_amph_buffer_cases.py            # 只打印（dry-run）
    python tools/add_amph_buffer_cases.py --write    # 写进 tests.json
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

from chemkit.data import load_tables             # noqa: E402
from chemkit.speciation import pKw_of            # noqa: E402
from tools.amph_exact import exact_pH, _pka      # noqa: E402

HALF = 0.06          # 断言半径（与 242–262 轮同口径）


def build(T) -> list[dict]:
    """返回待入库用例（`ph` 以**闭式精确解**为中心 ±HALF）。"""
    out = []
    specs = [
        # (前缀, 二元酸物种, 两性物种名, C, n_Na, T_K, 标签)
        ("AM1", "H_2S", "HS^-", 0.010, 0.003, 298.15, "硫系 3:10"),
        ("AM2", "H_2S", "HS^-", 0.010, 0.007, 298.15, "硫系 7:10"),
        ("AM3", "CO_2", "HCO_3^-", 0.010, 0.004, 298.15, "碳系 4:10"),
        ("AM4", "CO_2", "HCO_3^-", 0.010, 0.006, 298.15, "碳系 6:10"),
        ("AM5", "H_2S", "HS^-", 0.010, 0.004, 323.15, "硫系 4:10 @323.15K"),
        ("AM6", "CO_2", "HCO_3^-", 0.010, 0.006, 323.15, "碳系 6:10 @323.15K"),
        ("AM7", "H_2S", "HS^-", 0.002, 0.001, 298.15, "硫系 1:1 稀档"),
        ("AM8", "CO_2", "HCO_3^-", 0.002, 0.001, 298.15, "碳系 1:1 稀档"),
    ]
    for pre, acid, amph, C, nNa, T_K, tag in specs:
        p1 = _pka(T, acid, T_K, 0)
        p2 = _pka(T, amph, T_K, 0)
        pKw = pKw_of(T_K)
        ex = exact_pH(p1, p2, C, n_Na=nNa, pKw=pKw)
        mid = 0.5 * (p1 + p2)
        out.append({
            "name": f"{pre} {acid} {C:g}M + NaOH {nNa:g}M（两性缓冲对）",
            "subs": [[acid, C], ["NaOH", nNa]],
            "cond": {"V_L": 1.0, "T_K": T_K},
            "ph": [round(ex - HALF, 2), round(ex + HALF, 2)],
            "note": (
                f"第 263 轮扩充 · **两性缓冲对锚**（{tag}）。"
                f"答案由**引擎无关的闭式解**给出：H₂A 总量 C、外加钠 n_Na 时 "
                f"`D = h²+Ka1·h+Ka1·Ka2`、`[HA⁻] = C·Ka1·h/D`、"
                f"`[A²⁻] = C·Ka1·Ka2/D`，电荷条件 `n_Na/V + h = [HA⁻] + 2[A²⁻] + oh` "
                f"（`f(pH)` 严格单调 ⟹ 唯一根，二分）。"
                f"本档 pKa1={p1:.4f}、pKa2={p2:.4f}（T={T_K} K，含 van't Hoff），"
                f"**闭式解 {ex:.4f}**；"
                f"作为对照，两性中点式 `(pKa1+pKa2)/2 = {mid:.4f}` "
                f"（**偏 {mid - ex:+.4f}**）——中点式的推导前提是"
                f"「该两性物种是唯一质子条件物种」，而本档账本里 "
                f"{acid} 与 {amph} 同时在账 ⟹ 它是缓冲对的一员，前提不成立。"
                f"第 263 轮修法正是把这条支路改为先问精确质子条件。"
                f"（判据与正/负对照见 `tools/amph_exact.py --calib`。）"),
        })
    return out


def main(argv: list[str]) -> int:
    T = load_tables()
    new = build(T)
    path = os.path.join(ROOT, "chemkit", "data", "tests.json")
    with io.open(path, encoding="utf-8") as fh:
        db = json.load(fh)
    have = {c["name"] for c in db}
    print(f"{'name':<46} {'ph 区间':>16}  推导")
    for c in new:
        mark = "已存在" if c["name"] in have else "新增"
        print(f"{c['name'][:46]:<46} {str(c['ph']):>16}  {mark}")
    if "--check" in argv:
        # 只**测量**引擎与独立闭式解的差，不用它调断言（标准以闭式解为准）。
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<46} {'闭式解':>8} {'引擎':>8} {'Δ':>8} "
              f"{'resid':>7} {'判读':>6}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c["cond"], T, _probe=pr)
            ex = 0.5 * (c["ph"][0] + c["ph"][1])
            got = r.get("final_pH")
            d = (got - ex) if got is not None else None
            ok = "✓" if (d is not None and abs(d) <= HALF) else "✗"
            print(f"{c['name'][:46]:<46} {ex:>8.4f} {str(got):>8} "
                  f"{(f'{d:+.4f}' if d is not None else '—'):>8} "
                  f"{abs(_live(pr.get('active') or [])):>7.4f} {ok:>6}")
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
