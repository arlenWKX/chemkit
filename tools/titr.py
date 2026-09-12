"""滴定记账审计：`_buffer_titration` 入/出对照（族总量守恒 + 虚拟账本差异）。

背景（architecture §7 L）：`_buffer_titration` 的多级堆循环每弹出一个条目就
做两件事——(a) 质量记账 (b) Henderson 定 pH——而 (a) 用**绝对写**
`ledger2[base] = m - take`，m 只是该堆条目入堆时的**局部量**。多级模式下
产物按部分量重新入堆，再次弹出时绝对写会把账本里同名物种的其余量整块抹掉
（FeCl3+Na2CO3：碳 3.000000 → 2.975767，pH 3.09 虚低）。

本工具给出**可断言的量化**：每次调用前后各质子化族的成员总量差。
守恒的实现在该列上必须恒为 0（族总量是质子化梯的基本不变量）。

用法：
    python tools/titr.py B19 E42 K02          # 指定用例前缀
    python tools/titr.py --all                # 全量（只打汇总）
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

from chemkit.acidbase import build_families                      # noqa: E402
from chemkit.data import load_tables                             # noqa: E402
from chemkit import speciation                                   # noqa: E402
from chemkit.engine import judge                                 # noqa: E402
from chemkit.testsuit import load_cases                          # noqa: E402

_ORIG = speciation._buffer_titration
_CALLS: list = []


def _wrap(ledger, H_excess, V, T, pKw, multilevel=False, T_K=298.15,
          cache=None, touch=None):
    before = dict(ledger)
    out = _ORIG(ledger, H_excess, V, T, pKw, multilevel, T_K, cache, touch)
    # 快照返回值：虚拟账本被调用方（judge 走步）继续原地改写，
    # 直接留引用会在报告期看到"来自未来的差异"（假破例）
    _CALLS.append((before, H_excess, multilevel, (out[0], out[1], dict(out[2]))))
    return out


speciation._buffer_titration = _wrap


def _fam_totals(led: dict, fams: dict) -> dict:
    tot: dict = {}
    for s, m in led.items():
        fi = fams.get(s)
        if fi is None:
            continue
        tot[fi[0]] = tot.get(fi[0], 0.0) + m
    return tot


def _losses(before: dict, after: dict, fams: dict) -> list:
    a = _fam_totals(before, fams)
    b = _fam_totals(after, fams)
    out = []
    for k in set(a) | set(b):
        d = b.get(k, 0.0) - a.get(k, 0.0)
        if abs(d) > 1e-9:
            out.append((k, a.get(k, 0.0), b.get(k, 0.0), d))
    return out


def show(prefixes: list[str]) -> None:
    T = load_tables()
    fams = build_families(T)
    grand_bad = grand_calls = 0
    for c in load_cases(None):
        if prefixes and not c["name"].startswith(tuple(prefixes)):
            continue
        _CALLS.clear()
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        judge(subs, cond, T)
        bad = []
        for i, (before, He, ml, out) in enumerate(_CALLS):
            tit, He_res, led2 = out
            ls = _losses(before, led2, fams)
            if ls:
                bad.append((i, He, ml, ls, tit))
        grand_calls += len(_CALLS)
        grand_bad += len(bad)
        if not prefixes:
            continue
        print(f"\n{'='*78}\n{c['name']}  调用 {len(_CALLS)} 次；族总量破例 {len(bad)} 次")
        for i, He, ml, ls, tit in bad[:6]:
            print(f"  #{i} He={He:+.6f} multilevel={int(ml)} pH={tit} "
                  f"族损失 {len(ls)}：")
            for k, a, b, d in sorted(ls, key=lambda t: -abs(t[3]))[:5]:
                print(f"      {k:22s} {a:12.6f} → {b:12.6f}  Δ={d:+.6g}")
        if _CALLS:
            before, He, ml, (tit, He_res, led2) = _CALLS[-1]
            print(f"  末次调用：He={He:+.6f} multilevel={int(ml)} → pH={tit} "
                  f"He_res={He_res:+.3g}")
            for s in sorted(set(before) | set(led2)):
                a, b = before.get(s, 0.0), led2.get(s, 0.0)
                if abs(a - b) > 1e-12:
                    fi = fams.get(s)
                    tag = f" [族 {fi[0]}]" if fi else ""
                    print(f"      {s:24s} {a:12.6f} → {b:12.6f}{tag}")
    print(f"\n合计：调用 {grand_calls} 次；族总量破例 {grand_bad} 次")


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("-")]
    show([] if "--all" in sys.argv else (a or ["B19"]))
