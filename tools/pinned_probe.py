# -*- coding: utf-8 -*-
"""第 233 轮 · B 组单例可行性验证：`pinned` 能给 H45 一个自洽解吗？

第 232 轮定下：B 组 17 例 / 156 残差需 `pinned` 路线（大工程），
**启动前先做一次便宜的单例验证**——挑 `H45`（25.4，全库最大残差）。

判据：`H45` 终态账本里有 `Al(OH)_3 0.00046`（固相在场）与
`[Al(OH)_4]^- 0.8745`（铝酸根主导）。引擎给 pH 12.185，而账本电荷自洽解
10.534。第 213 轮已知 `charge_pH` **平账本**时把"总量当溶解量"⟹ 在有固相
时偏；但 `pinned` 参数正是为此而生：把 `[Al] = 10^((y·(pKw−pH) − pKsp)/x)`
作为电荷项**随 pH 连续变化**，于是固相从"分支选择"变成"方程的一项"。

本脚本直接构造 `pinned` 元组调用 `charge_pH`，看：
  ① 能否解出（不返回 None）；
  ② 解在哪个 pH；
  ③ 该 pH 下账本是否更自洽（与引擎现值 12.185 比）。

`pinned` 元素格式（`acidbase.charge_pH` docstring）：
  `(z, pKsp, x, y)` 或 `(z, pKsp, x, y, ladder)` 或 `(z, pKsp, x, y, ladder, frozen)`
  `ladder` = `((z_i, logβ_i, nu_i, 物种名), …)`（羟合梯，从属于钉住的 [M]）
  `frozen` = `((z_i, coef_i), …)`（非 OH⁻ 配合物从属项）

用法： python tools/pinned_probe.py
"""
import io
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

import chemkit.engine as eng                                # noqa: E402
import chemkit.acidbase as ab                               # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith("H45")][0]
probe = {}
eng.judge([{"name": a, "mol": b} for a, b in c["subs"]],
          c.get("cond") or {"V_L": V}, T, _probe=probe)
led = dict(probe.get("ledger") or {})
print(f"用例 {c['name']}")
print(f"  引擎 pH = {probe.get('pH_solver')}   max|S| = {probe.get('max_abs_S')}")
print(f"  账本 = { {k: round(v,6) for k,v in led.items() if k != 'H_2O'} }")

print("\n=== ① 平账本 charge_pH（已知在有固相时偏）===")
print(f"  {ab.charge_pH(dict(led), V, T, T_K)}")

# 取 Al(OH)3 的 pKsp / x / y
e_al = next((e for e in T.ksp if e["solid"] == "Al(OH)_3"), None)
print(f"\n=== ② Al(OH)_3 的 ksp 条目 ===\n  {e_al}")
x_a, y_a = eng._ksp_xy(e_al)
pk = spec._pksp(e_al, T_K)
print(f"  x={x_a} y={y_a} pKsp(T)={pk}")

# 羟合梯：Al 的 OH⁻ 配合物（从属于钉住 [Al]）
ladder = []
for b in T.beta:
    if b.get("center") == "Al^{3+}" and b.get("ligand") == "OH^-":
        ladder.append((spec.charge_of(b["complex"]), b["logb"], b["nu"],
                       b["complex"]))
print(f"\n=== ③ Al 羟合梯（{len(ladder)} 级）===")
for t in ladder:
    print(f"  z={t[0]:+d} logb={t[1]} nu={t[2]} {t[3]}")

print("\n=== ④ 逐种 pinned 组合试解 ===")
led2 = dict(led)
# 钉住前提：删掉游离 Al³⁺ 与羟合梯成员（电荷由钉住项+梯项表达）
removed = []
for sp in [t[3] for t in ladder] + ["Al^{3+}"]:
    if led2.pop(sp, None) is not None:
        removed.append(sp)
print(f"  已从账本删除（钉住取代）: {removed}")
print(f"  剩余账本 = { {k: round(v,6) for k,v in led2.items() if k != 'H_2O'} }")

cands = [
    ("裸 (z,pKsp,x,y)", (3, pk, x_a, y_a)),
    ("+梯", (3, pk, x_a, y_a, tuple(ladder))),
    ("+梯+空 frozen", (3, pk, x_a, y_a, tuple(ladder), ())),
]
for tag, pin in cands:
    try:
        r = ab.charge_pH(dict(led2), V, T, T_K, pinned=(pin,))
    except Exception as exc:                                # noqa: BLE001
        r = f"异常 {type(exc).__name__}: {exc}"
    print(f"  {tag:18s} -> {r}")

print("\n=== 判读 ===")
print("  · 若带梯的 pinned 解出一个**中性附近**（5~11）的 pH，且与")
print("    '引擎 12.185 / 平账本 10.534' 都不同 ⟹ pinned 路线可行，")
print("    值得投入实现（它把固相写进方程，是唯一自洽的做法）。")
print("  · 若返回 None/异常 ⟹ 需先补 pinned 的前置（如 frozen 项），")
print("    或该态不满足 pinned 的结构前提（固相量 0.00046 极小）。")
