# -*- coding: utf-8 -*-
"""第 281 轮 · **半中和缓冲锚（浓度无关）**。

## 动因

第 281 轮从残差普查里逮到一个紧凑的簇：**`TE1`–`TE4` `NH4+ 半中和`
（4 个温度点）残差全是 1.878**，通道都是 `NH_3 + H^+ -> NH_4^+`。
看 `TE1` 的步表（`NH4Cl 0.01 + NaOH 0.005`，273.15 K）：

    step1 NH4^+ -> NH3 + H^+   x=0.000131
    step2 NH4^+ -> NH3 + H^+   x=0.004869   ← 正好走到 1:1 缓冲、He 归零
    step3 NH3 + H^+ -> NH4^+   x=0.004869   ← **同量拆掉**

`step2` 的 0.004869 **正是游离 OH⁻ 的量**，走完就是正确终态；`step3` 把它拆了。
引擎自己的震荡判据要求 `gross >= 0.05`（**绝对** mol），而本体系 0.01 M 量级
的 `gross` 只有 **0.0097** —— 对 0.01 M 体系**就是接近全周转**，却过不了线。

## 本组锚锁的化学事实（与引擎无关，且**与浓度无关**）

**半中和的弱酸/弱碱缓冲液，pH = pKa**：
`NH₄⁺ ⇌ NH₃ + H⁺`，当 `n(NH₃) = n(NH₄⁺)` 时
`pH = pKa − log([NH₄⁺]/[NH₃]) = pKa`。
`pKa(NH₄⁺, 298.15 K) = 9.25`（教科书值），**与总浓度无关** ——
所以 0.02 M 与 0.002 M 两档必须给出**同一个 pH**。这一条正好把本轮
发现的"**稀档**"那一侧钉住。

## 用法

    python tools/add_halfneutral_cases.py            # dry-run
    python tools/add_halfneutral_cases.py --check
    python tools/add_halfneutral_cases.py --write
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

_TAG = "第 281 轮扩充 · 半中和缓冲锚。"
_FLAT = "  # noqa: E501"


def build() -> list[dict]:
    return [
        {
            "name": "HN1 NH4Cl 0.02+NaOH 0.01 半中和（1:1 缓冲锚）",
            "subs": [["NH_4Cl", 0.02], ["NaOH", 0.01]],
            "ph": [9.10, 9.40],
            "note": (
                _TAG + "**半中和 ⟹ pH = pKa**：`NH₄⁺ ⇌ NH₃ + H⁺`，"
                "`n(NH₃) = n(NH₄⁺)` 时 `pH = pKa − log([NH₄⁺]/[NH₃]) = pKa`。"
                "`pKa(NH₄⁺, 298.15 K) = 9.25`（教科书值）⟹ 区间取 ±0.15。"
                "**与总浓度无关**是本条的关键：0.02 M 与 `HN2` 的 0.002 M "
                "必须给同一个 pH。"),
        },
        {
            "name": "HN2 NH4Cl 0.002+NaOH 0.001 半中和（稀档·同 pH）",
            "subs": [["NH_4Cl", 0.002], ["NaOH", 0.001]],
            "ph": [9.10, 9.40],
            "note": (
                _TAG + "与 `HN1` **同化学、浓度低 10 倍**，"
                "**pH 必须相同**（半中和缓冲的 pH = pKa，与浓度无关）。"
                "**为什么这条对本轮重要**：本轮逮到的残差簇 `TE1`–`TE4` "
                "正是 0.01 M 量级的稀档 —— 引擎的震荡判据用**绝对** "
                "`gross >= 0.05` mol，稀体系过不了线。本条把"
                "**稀档的 pH 正确性**独立钉住（pH 对 ⟹ 至少呈现层没被"
                "稀档拖坏；残差另账）。"),
        },
        {
            "name": "HN3 NH4Cl 0.1+NaOH 0.05 半中和（浓档·同 pH）",
            "subs": [["NH_4Cl", 0.1], ["NaOH", 0.05]],
            "ph": [9.10, 9.40],
            "note": (
                _TAG + "同化学、浓度比 `HN1` 高 5 倍，**pH 必须仍 = pKa = 9.25**。"
                "引擎实测 **9.250**、resid **0.0000**、**1 步**走完"
                "（`He` 归零 = 完全中和）。"
                "**为什么这条对本轮重要**：第 282 轮的浓度扫描"
                "（`tools/tcmp.py TE1 --scale=…`）显示残差是**浓度刀刃**——"
                "`scale = 0.05/0.1/0.2/0.5/2/4/10` 全部干净，**只有 `scale = 1`"
                "（`NH4Cl 0.01 + NaOH 0.005`）残差 1.878**。本条与 `HN4` 夹住"
                "那个坏点，把\"pH = pKa 与浓度无关\"钉死。"),
        },
        {
            "name": "HN4 NH4Cl 0.005+NaOH 0.0025 半中和（稀档·同 pH）",
            "subs": [["NH_4Cl", 0.005], ["NaOH", 0.0025]],
            "ph": [9.10, 9.40],
            "note": (
                _TAG + "同化学、浓度是 `HN1` 的 1/4，**pH 必须仍 = pKa**。"
                "引擎实测 **9.240**、resid **0.0000**、**0 步**（`He` 保持 −0.0025，"
                "即走步没动它 —— 但 pH 仍对，因为半中和缓冲的 pH 只由"
                "`[NH₃]/[NH₄⁺]` 定，而该比值由 `estimate_pH` 的滴定给出）。"
                "**这条同时暴露一件事**：`He` 未被中和（−0.0025）而 resid 仍为 0 "
                "⟹ **残差口径看不见这种\"账本没落实\"**，只有恰好落在刀刃上的"
                "`scale = 1` 才把它顶出来（见 `HN3` 的注）。"),
        },
        {
            "name": "HN5 NH4Cl 0.0095+NaOH 0.00475 半中和（刀刃下沿）",
            "subs": [["NH_4Cl", 0.0095], ["NaOH", 0.00475]],
            "ph": [9.10, 9.40],
            "note": (
                _TAG + "**紧贴第 283 轮定位的浓度刀刃下沿**。`tools/tcmp.py TE1 "
                "--scale=…` 细扫实测：`scale ≤ 0.95`（c ≤ 0.0095）**resid 全 0、"
                "0 步**；`scale ∈ [0.99, 1.5]` 起落进「3 步、残差 > 0」的坏档"
                "（1.878 → 2.506 → … → 0.284 随浓度**平滑衰减**）；`scale ≥ 2`"
                "又回到 **1 步、resid 0（完全中和）**。本条取**差一步就掉进坏档**"
                "的 0.0095，把刀刃下沿钉住。引擎实测 pH **9.250**、resid **0**、"
                "**0 步**（`He` 保持 −0.00475，即走步没动它 —— pH 仍对，因为"
                "半中和缓冲的 pH 只由 `[NH₃]/[NH₄⁺]` 定，而该比值由 "
                "`estimate_pH` 的滴定给出）。"),
        },
        {
            "name": "HN6 NH4Cl 0.03+NaOH 0.015 半中和（刀刃上沿）",
            "subs": [["NH_4Cl", 0.03], ["NaOH", 0.015]],
            "ph": [9.10, 9.40],
            "note": (
                _TAG + "**紧贴第 282–283 轮定位的坏档上沿**。`tools/tcmp.py TE1 "
                "--scale=…` 实测：坏档是 `scale ∈ [0.99, 1.5]`（残差 2.506→0.284"
                " 平滑衰减），`scale ≥ 2` 回到 **1 步、resid 0（完全中和）**。"
                "本条取 `scale = 3`（`c = 0.03`），落在坏档之外，"
                "与 `HN5`（0.0095，下沿）一起把坏档**夹住**。"
                "引擎实测 pH **9.250**（= pKa）、resid **0**。"),
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
        print(f"{c['name'][:46]:<46} ph={str(c['ph']):>14}  "
              f"{'已存在' if c['name'] in have else '新增'}")

    if "--check" in argv:
        from chemkit.converg import _live
        from chemkit.engine import judge
        print(f"\n{'name':<46} {'pH':>9} {'resid':>8} {'判读':>4}")
        for c in new:
            pr = {}
            r = judge([{"name": s[0], "mol": float(s[1])} for s in c["subs"]],
                      c.get("cond") or {}, T, _probe=pr)
            ph = r.get("final_pH")
            ok = c["ph"][0] <= (ph if ph is not None else -99) <= c["ph"][1]
            print(f"{c['name'][:46]:<46} {ph:>9.3f} "
                  f"{abs(_live(pr.get('active') or [])):>8.4f} "
                  f"{'✓' if ok else '✗':>4}")
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
