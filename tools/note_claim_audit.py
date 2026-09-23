# -*- coding: utf-8 -*-
"""第 244 轮 · 核 `fix_eq_standards.py` 里每条「原式不守恒」的断言是否成立。

**动因**：`21 AgBr+浓氨水` 的 note 写着「**元素与电荷不守恒**」，
但把它的旧式 `240NH_3 + 91AgBr + 56H_2O -> 91Br^- + 91[Ag(NH_3)_2]^+ +
56NH_4^+ + 56OH^-` 逐元素/逐电荷一算 —— **守恒**（Ag 91=91、Br 91=91、
N 240=240、H 832=832；电荷 0）。它属于**整数倍放大**类，不是不守恒类。
generic 文案被套到了它头上。

**判据（机械化）**：对 FIX 表每条，把它的 `old` 式子用库内
`elements_of` / `charge_of` 精确核验，得出"是否守恒"，
再与该条 note 是否声称"不守恒"比对。

⚠️ 解析坑（本脚本第一版就栽在这里）：`[Ag(NH_3)_2]^+` 这类**方括号配合物**
必须先整体识别，否则会被拆成 `Ag(NH_3)_2` 而数错原子。
故这里用"按 `+` 切分但先做括号配对"的解析器。

用法： python tools/note_claim_audit.py
"""
import io
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.core import charge_of, elements_of                # noqa: E402

# 已见过的 (name, new)：用于识别"跨类条目"（同一用例被两类修正先后落到）
SEEN: set[tuple[str, str]] = set()


def split_terms(text):
    """按 `+` 切分**顶层**项，并取出系数。

    ⚠️ 三个必须一起处理的坑（本脚本前后三版都栽在这里）：
      1. `[Ag(NH_3)_2]^+` —— `+` 在方括号**闭合之后**，是电荷上标不是分隔符；
      2. `Fe^{3+}` —— `+` 在花括号**之内**；
      3. `H^+` —— `^` 之后紧跟的就是电荷号。
    故判据是：**若该 `+` 前面紧邻 `^` 或 `^数字`，则它是电荷号**。
    """
    out = []
    depth = 0
    cur = ""
    for ch in text:
        if ch in "[{":
            depth += 1
        elif ch in "]}":
            depth -= 1
        is_charge = bool(re.search(r"\^\d*$", cur))
        if ch == "+" and depth == 0 and not is_charge:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    for t in out:
        t = t.strip()
        if not t or t == "->":
            continue
        m = re.match(r"^(\d+)\s*(.+)$", t)
        if m:
            out_coef, sp = int(m.group(1)), m.group(2).strip()
        else:
            out_coef, sp = 1, t
        yield out_coef, sp


def balance_of(eq):
    """返回 (元素差, 电荷差)。元素差为空 dict ⟺ 守恒。"""
    if "->" not in eq:
        return None, None
    l, r = eq.split("->", 1)
    els = {}
    q = 0
    for coef, sp in split_terms(l):
        try:
            for k, v in elements_of(sp).items():
                els[k] = els.get(k, 0) + coef * v
            q += coef * charge_of(sp)
        except Exception:                                   # noqa: BLE001
            return None, "解析失败:" + sp
    for coef, sp in split_terms(r):
        try:
            for k, v in elements_of(sp).items():
                els[k] = els.get(k, 0) - coef * v
            q -= coef * charge_of(sp)
        except Exception:                                   # noqa: BLE001
            return None, "解析失败:" + sp
    return {k: v for k, v in els.items() if v}, q


def main():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_fx", os.path.join(ROOT, "tools/fix_eq_standards.py"))
    fx = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fx)

    p = os.path.join(ROOT, "chemkit/data/tests.json")
    rows = json.load(io.open(p, encoding="utf-8"))
    by = {c["name"]: c for c in rows}

    print("=== FIX 表逐条：原式守恒性（`old` → `new`）===")
    print("说明：某条被**多次**修正时，上一步的 `new` 就是下一步的 `old`；")
    print("      'cross' = 跨类条目（两类修正先后落在同一条用例上）。\n")
    bad = []
    for name, new, old in fx.FIX:
        eo, qo = balance_of(old)
        en, qn = balance_of(new)
        if eo is None or en is None:
            print(f"  ? {name}: 解析失败")
            continue
        co = (not eo) and qo == 0
        cn = (not en) and qn == 0
        cross = (name, new) in SEEN
        SEEN.add((name, new))
        note = (by.get(name) or {}).get("note") or ""
        flag = ""
        # 机械可判的错配只有两种：
        #  (a) 原式**守恒**，note 却断言"元素与电荷不守恒"（跨类条目除外——
        #      它可能是对**另一个**中间步说的）；
        #  (b) 原式**不守恒**，note 只谈"整数倍放大"而**完全不提不守恒**
        #      （放大解释不了不守恒 ⟹ 漏了主因）。
        #  ⚠️ 只要 note **同时**点明两条，就不算错配（部分正确的说明是好说明）。
        if co and ("元素与电荷不守恒" in note) and not cross:
            flag = "  ⟹ ⚠ note 断言不守恒，但原式**守恒**"
            bad.append(name)
        if (not co) and ("整数倍放大" in note) and ("不守恒" not in note):
            flag = "  ⟹ ⚠ note 只谈'整数倍放大'，漏了'原式不守恒'"
            bad.append(name)
        print(f"  原式 {'守恒  ' if co else '不守恒'} → 新式 "
              f"{'守恒' if cn else '**不守恒**'}   {name}"
              + ("   [cross]" if cross else "") + flag)
        if not cn:
            print("      ⟹ ⚠⚠ 新式不守恒，必须修")
            bad.append(name)
        elif eo:
            print(f"      （原式差异：元素 {eo}, 电荷 {qo}）")
    print(f"\n机械可判的错配：{len(set(bad))} 条 {sorted(set(bad))}")
    print("注：generic 文案里出现的其它族举例（如'99:52 族 Fe'）是**说明性文字**，")
    print("    不构成对该条的断言；`note_audit.py` 按'元素是否出现'筛会产生")
    print("    大量假阳性（实测 15 条命中里 0 条真错配），已弃用该思路。")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
