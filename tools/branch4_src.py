# -*- coding: utf-8 -*-
"""第 275 轮 · **分支 4 选支来源普查**：呈现 pH 是按哪几个物种、哪一种角色定出来的，
以及**账本里哪些显著物种根本没有角色**。

## 为什么需要（§1.12 的下一步由"调阈值"改为"读账本"）

第 274 轮拿到 §1.12 的硬证据：`H43 AlCl3+NaOH 1:3.5` 的**呈现 pH = 3.26**，
而把它的**账本**代回电荷平衡手算得 **pH 12.03**（差 8.8 个单位）。第 274 轮
还查明该例的**账本本身是合理的**（0.845 mol 铝酸根 + 0.0326 Al(OH)₃ + 0.0313
Al³⁺）⟹ 错的是**选支**，不是走步。

`estimate_pH` 分支 4 的 `h_c`/`o_c` 由**角色表**（`rolemap`）扫描**在账物种**
累加而来，而 `rolemap` 只由 `acids_map ∪ bases_map ∪ hyd_map ∪ amph_pH` 构成
（speciation.py L843-848），后三者都源自 `T.pka` / `T.beta` 注册。**推论**：
账本里的碱侧配合物 `[Al(OH)₄]⁻` 若没有注册成某个 pKa 对的**碱侧**，
它就**没有角色** ⟹ `o_c` 看不见它 ⟹ 只剩痕量 `Al³⁺` 撑起 `h_c` ⟹ 选酸侧。

本工具用引擎自己的 `PH_SRC` 审计（生产路径恒 None）**把这件事量死**：
① 列出实际参与 `h_c`/`o_c` 的 (物种, 角色, 值)；② 列出账本里 ≥ 阈值的显著
物种中**没有角色**的那些（那就是选支盲区）。

用法：
    python tools/branch4_src.py H43
    python tools/branch4_src.py H43 --min=1e-4
"""

from __future__ import annotations

import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
import chemkit.speciation as spec                           # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    thresh = 1e-6
    for a in argv:
        if a.startswith("--min="):
            thresh = float(a[6:])
    pre = args[0]
    T = load_tables()
    hit = [c for c in load_cases(None) if c["name"].split()[0] == pre]
    if not hit:
        print(f"未找到用例 {pre}")
        return 1
    c = hit[0]
    cond = c.get("cond") or {}
    V = float(cond.get("V_L", 1.0))

    spec.PH_SRC = []
    pr: dict = {}
    try:
        r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T,
                      _probe=pr)
    finally:
        src = spec.PH_SRC
        spec.PH_SRC = None
    print(f"=== {c['name']} ===  呈现 pH = {r.get('final_pH')}  "
          f"（solver {pr.get('pH_solver')}）")

    # ① 实际参与选支的来源
    agg: dict = defaultdict(lambda: [0, 0.0])
    picks = []
    for e in src:
        if e[0] == "__branch4__":
            picks.append(e)
            continue
        g = agg[(e[0], e[1])]
        g[0] += 1
        g[1] = max(g[1], e[2])
    print("\n① 参与 h_c/o_c 的 (物种, 角色)：")
    print(f"   {'物种':<16}{'角色':<14}{'次数':>7}{'最大值(M)':>14}")
    for (sp, tag), (n, v) in sorted(agg.items(), key=lambda kv: -kv[1][1]):
        print(f"   {sp:<16}{tag:<14}{n:>7}{v:>14.6g}")
    if picks:
        print(f"   最终选支：{picks[-1][1]}  值 {picks[-1][2]:.6g}")

    # ② 账本里的显著物种，哪些**没有角色**（选支盲区）
    # ⚠️ 必须读**引擎真正的 `rolemap`**（从 `T._est_static[T_K]` 取），
    # 而不是 `T.pka_acid`/`T.pka_base` —— 第 276 轮起含氧酸根等角色的来源
    # 还包括 `_first_k`（β 一级水解）、`hyd_map`（Ksp 水解）与
    # **含氧酸根碱侧角色**（由 β/Ksp 推出），只看 pKa 表会漏报。
    _sc = (getattr(T, "_est_static", None) or {}).get(298.15)
    _rolemap = _sc[6] if _sc else {}
    led = dict(pr.get("ledger") or {})
    if not led:
        for e in (r.get("final") or []):
            led[e["name"]] = e["mol"]
    print(f"\n② 账本显著物种（≥ {thresh:g} M）的角色状态"
          f"（rolemap 共 {len(_rolemap)} 项）：")
    print(f"   {'物种':<18}{'M':>12}  {'in rolemap':<11}角色明细")
    blind = []
    for sp, m in sorted(led.items(), key=lambda kv: -kv[1]):
        if sp == "H_2O" or m / V < thresh:
            continue
        role = _rolemap.get(sp)
        tags = []
        if role:
            Ka, Kb, conj_p, Kh_qc, amph_v = role
            if Ka is not None:
                tags.append("Ka")
            if Kb is not None:
                tags.append(f"Kb={Kb:.4g}")
            if Kh_qc is not None:
                tags.append("Kh")
            if amph_v is not None:
                tags.append("amph")
        tag = "  ".join(tags) if tags else "**无角色（选支盲区）**"
        if role is None:
            blind.append((sp, m))
        print(f"   {sp:<18}{m / V:>12.5g}  {str(bool(role)):<11}{tag}")
    if blind:
        print(f"\n   ⟹ {len(blind)} 个显著物种没有酸碱角色："
              f"{', '.join(sp for sp, _ in blind)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
