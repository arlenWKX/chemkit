# -*- coding: utf-8 -*-
"""第 277 轮 · **角色表注入自检**：把 `estimate_state` 静态块里"含氧酸根碱侧
角色"的计算**逐行重放**一遍，看它到底在哪一步停住。

背景：第 276 轮把该角色写进了 `bases_map`（代码在 `speciation.py` 约 L880-901），
但 `tools/branch4_src.py H45` 实测 `[Al(OH)_4]^-` **不在 `rolemap`** 里
⟹ `o_c` 看不见账本里 0.5 M 的铝酸根 ⟹ 分支 4 给 pH 7.0。
本工具重放同一段逻辑并打印每一步的中间量，定位是哪一步 `continue` 掉了。

用法：python tools/role_inject_check.py
"""

from __future__ import annotations

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

from chemkit.core import charge_of                        # noqa: E402
from chemkit.data import load_tables                      # noqa: E402
from chemkit.speciation import _pksp                      # noqa: E402

T_K = 298.15


def main() -> int:
    T = load_tables()
    _ksp_oh = {}
    for _e in T.ksp:
        _c0, _a0 = _e["pair"]
        if _a0 == "OH^-" and _c0 not in _ksp_oh:
            _ksp_oh[_c0] = _e
    print(f"_ksp_oh 命中 {len(_ksp_oh)} 个 OH 型阳离子")
    print(f"  Al^{{3+}} 在内？ {('Al^{3+}' in _ksp_oh)}"
          f"   {_ksp_oh.get('Al^{3+}', {}).get('solid')}")
    _lad = {}
    n_pos = n_lad = 0
    for b in T.beta:
        if b["ligand"] != "OH^-" or b.get("m", 1) != 1:
            continue
        cx = b.get("complex")
        if not cx or cx in T.solids:
            continue
        if charge_of(cx) >= 0:
            n_pos += 1
            continue
        _lad[(b["center"], b.get("nu", 1))] = (cx, float(b["logb"]))
        n_lad += 1
    print(f"_lad 收 {n_lad} 个含氧酸根（跳掉 {n_pos} 个阳离子/中性梯级）")
    for (c, k), (cx, lb) in sorted(_lad.items()):
        e = _ksp_oh.get(c)
        z = charge_of(c)
        n = k - z
        kb = None if (e is None or n <= 0) else 10.0 ** ((_pksp(e, T_K) - lb) / n)
        print(f"  {cx:<16} 中心={c:<10} k={k} z={z} n={n} "
              f"Ksp={'有' if e else '**无**'}  Kb={kb}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
