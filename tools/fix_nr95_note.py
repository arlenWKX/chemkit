# -*- coding: utf-8 -*-
"""第 243 轮 · 定向修复 `NR95 AgCl+CaCO3` 的 note（张冠李戴）。

**背景**：`tools/fix_eq_standards.py` 的旧实现把四类修正合进一个 `FIX` dict，
`NR95` 同时属于"整数倍放大"与"报告口径化石"两类 ⟹ 后写的覆盖前写的；
而 `REASON` 是按 name 索引的单个字符串，脚本用
`REASON.get(name) or <第一类通用文案>` 取值时，`REASON` 里那条**是存在的**，
只是内容被写成了上一条（`P10`）的口径……最终落盘的是**第一类通用文案**，
于是 note 里出现"如 99:52 族 Fe 52≠21+18+12=51"——**与 AgCl/CaCO₃ 毫无关系**。

`eq` 值本身是对的（3:2:2:1 主通道守恒式），只有 note 文案错。
本脚本只改这一条的 note，不动 `eq`。

用法： python tools/fix_nr95_note.py [--write]
"""
import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

PATH = os.path.join(ROOT, "chemkit", "data", "tests.json")
NAME = "NR95 AgCl+CaCO3"

NOTE = (
    "【v0.5.0 标准修正·报告口径化石 + 整数倍放大，**跨两类**】"
    "原式 `497CaCO_3 + 343H_2O + 56AgCl -> …` 是引擎旧路径吃 `round(x,6)` "
    "后的净差被**量化**（两条迹量通道都被抬成主通道）再放大 56 倍印出来的；"
    "第一步先按「整数倍放大」约到 9:1 式，第二步按「报告口径化石」去掉 AgCl 通道"
    "（1.41e-5，占 CaCO₃ 通道 11%，落在主通道呈现的接受线之外），"
    "得 CaCO₃ 的溶解-水解守恒式 **3:2:2:1**（已过 eqcheck）。"
    "审计链：497:343:56 → 9:6:1 → 3:2:2:1。"
    "⚠️ 旧 note 里写的「99:52 族 Fe 52≠21+18+12=51」是**张冠李戴**"
    "（那属 N30/D32 Fe-SCN 族），系 fix_eq_standards 的 dict 覆盖缺陷所致；"
    "该工具已改为有序 list 并加了 FINAL 目标表（第 243 轮）。"
)


def main():
    write = "--write" in sys.argv
    rows = json.load(io.open(PATH, encoding="utf-8"))
    hit = [c for c in rows if c["name"] == NAME]
    if not hit:
        print(f"!! 未找到 {NAME}")
        return 2
    c = hit[0]
    print(f"当前 note 片段: …{c.get('note', '')[-120:]}")
    print(f"\n新 note:\n{NOTE}")
    if not write:
        print("\n（--write 写回）")
        return 0
    c["note"] = NOTE
    with io.open(PATH, "w", encoding="utf-8", newline="\n") as f:
        json.dump(rows, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"\n[写入] {PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
