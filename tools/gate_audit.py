"""闸门相位清单（§7 X-35）：把动力学闸门的条件按"相脆弱/相稳定"分类。

动机（第 50 轮）：EU01 的硝酸根闸门被"走步**瞬时**酸度"点燃，而 couples.json
注记把正解指向该闸门。要谈"相位一致性"，先得有清单：哪些闸门读的是**瞬时态**
（`pH_max/pH_min`——随走步摆动）、哪些读的是**相稳定量**（`c_min/c_max` 配
`c_basis:"shadow"`，即强酸/配体的总量记账，与走步相位无关）。

化学纪律：**动力学闸门是"速率"论断，应在相稳定量上求值**；瞬时量只适合
"形态"闸门（如 Mn²⁺ 游离形态仅酸性存在——那是热力学形态问题，不是速率）。
但 `pH_max` 不能从硝酸根闸门简单删掉：硝酸**盐**（Eu(NO₃)₃ 等）在账本里同样
贡献 shadow 强酸量，删掉 pH 条件会让中性介质里的硝酸盐被解锁还原 ✗。
⟹ 正解是给 pH 条件本身换一个**相稳定**的写法（如"平衡质子预算/总强酸"），
而不是删条件。本工具只做清单，不改行为。

用法：python tools/gate_audit.py [--all]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, ".")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="列出全部闸门条目")
    args = ap.parse_args()
    with open(os.path.join("chemkit/data/couples.json"), encoding="utf-8") as f:
        rows = json.load(f)
    couples = rows["couples"] if isinstance(rows, dict) else rows
    frag, stable, form = [], [], []
    for c in couples:
        g = c.get("gate")
        if not g:
            continue
        name = f"{c['ox']}/{c['red']}"
        ph_win = ("pH_max" in g) or ("pH_min" in g)
        c_basis = g.get("c_basis", "账本瞬时浓度")
        stable_c = ("c_min" in g) or ("c_max" in g)
        row = (name, {k: v for k, v in g.items()}, c_basis)
        if ph_win and stable_c:
            frag.append(row)          # 两种相位混用：本轮病灶形态
        elif ph_win:
            (form if "形态" in (c.get("note") or "") else frag).append(row)
        elif stable_c:
            stable.append(row)
    for tag, bag in (("混用（pH 窗 + 浓度条件）", frag), ("相稳定（仅浓度）", stable),
                     ("形态闸门（仅 pH 窗）", form)):
        print(f"\n== {tag}：{len(bag)} 条 ==")
        for name, g, cb in (bag if args.all else bag[:12]):
            print(f"  {name:<22} {g}   c_basis={cb}")
    print(f"\n共 {len(frag) + len(stable) + len(form)} 条带闸门电对")
    return 0


if __name__ == "__main__":
    sys.exit(main())
