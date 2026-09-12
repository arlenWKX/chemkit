"""数据库机械审计：守恒、电子数复算、重复与冲突。

覆盖 chemkit/data 的全部数值表，检查**不依赖文献**的自洽性错误——
这类错误一旦存在就是硬错误（反应写不出来 / 电子数错 / 同一事实两个值）：

  pka.json    acid -> base + n·H⁺：非 H 元素必须相等、n 必须等于 H 数差、
              电荷差必须等于 n
  ksp.json    solid + w·H2O -> x·cat + y·an：按 H/O 配平求 w，非 H/O 元素
              必须相等；x/y 与 _ksp_xy 一致
  beta.json   center + nu·ligand -> complex：元素与电荷**精确**相等
  couples.json ox + a·H⁺ + b·H2O + n·e⁻ -> red：非 H/O 元素必须相等；
              **n 由电荷平衡独立复算**（与表内 n 对照）；b 由 O 平衡给出
  thermo.json ΔHf：键可被 elements_of 解析
  全表       同一键（acid,base）/solid/(center,ligand,nu)/(ox,red) 不得重复；
              跨表引用的物种必须可解析

用法：python tools/data_audit.py
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, ".")

from chemkit.core import (charge_of, elements_of,  # noqa: E402
                          FormulaError, parse_species)
from chemkit.candidates import _ksp_xy            # noqa: E402

DATA = "chemkit/data"
WATER = "H_2O"


def load(name: str):
    with open(os.path.join(DATA, name), encoding="utf-8") as f:
        return json.load(f)


def _el(s: str) -> dict:
    try:
        return dict(elements_of(s))
    except (FormulaError, Exception):     # noqa: BLE001
        return {}


def _chg(s: str) -> int:
    try:
        return charge_of(s)
    except Exception:                     # noqa: BLE001
        return 0


def _nonho(e: dict) -> dict:
    return {k: v for k, v in e.items() if k not in ("H", "O")}


def rep(kind: str, msg: str, bag: list) -> None:
    bag.append(f"[{kind}] {msg}")


def _solve_nw(acid: str, base: str, nmax: int = 4):
    """`acid + w·H2O -> base + n·H⁺` 的元素空间求解。

    返回 (n, w) 列表（所有可行的 n）。不能假定 n 就是表内值——那正是要
    检验的对象；也不能忽略隐含水（CO₂ + H₂O -> HCO₃⁻ + H⁺），否则
    全部含氧酸都会误报（首版即此）。
    """
    ea, eb = _el(acid), _el(base)
    els = set(ea) | set(eb)
    out = []
    for n in range(1, nmax + 1):
        w = eb.get("O", 0) - ea.get("O", 0)
        ok = True
        for el in els:
            lhs = ea.get(el, 0) + w * (2 if el == "H" else 1 if el == "O" else 0)
            rhs = eb.get(el, 0) + n * (1 if el == "H" else 0)
            if lhs != rhs:
                ok = False
                break
        if ok and _chg(acid) - _chg(base) - n != 0:
            ok = False
        if ok:
            out.append((n, w))
    return out


def audit_pka(rows, bad, dup) -> int:
    seen = defaultdict(list)
    for e in rows:
        a, b, n = e["acid"], e["base"], e.get("n", 1)
        seen[(a, b)].append(e["pka"])
        if not _el(a) or not _el(b):
            rep("解析失败", f"pka {a} / {b}", bad)
            continue
        sols = _solve_nw(a, b)
        if not sols:
            rep("无解", f"{a} (pKa {e['pka']}) 与 {b} 配不出守恒式", bad)
        elif n not in [s[0] for s in sols]:
            rep("n 值", f"{a} -> {b} + {n}H⁺：元素空间只支持 n="
                        f"{[s[0] for s in sols]}（w={[s[1] for s in sols]}）", bad)
    for k, v in seen.items():
        if len(v) > 1:
            rep("重复", f"pKa 对 {k} 出现 {len(v)} 次：{v}", dup)
    return len(rows)


def audit_ksp(rows, bad, dup) -> int:
    seen = defaultdict(list)
    for e in rows:
        solid, (cat, an), pk = e["solid"], e["pair"], e["pKsp"]
        seen[solid].append(pk)
        es, ec, ea = _el(solid), _el(cat), _el(an)
        if not es or not ec or not ea:
            rep("解析失败", f"ksp {solid}/{cat}/{an}", bad)
            continue
        x, y = _ksp_xy(e)
        # w 由 O 平衡定：O(s) + w = x·O(cat) + y·O(an)
        w = ec.get("O", 0) * x + ea.get("O", 0) * y - es.get("O", 0)
        lhs = dict(es)
        lhs["H"] = lhs.get("H", 0) + 2 * w
        lhs["O"] = lhs.get("O", 0) + w
        rhs = {}
        for sp, k in ((cat, x), (an, y)):
            for el, v in _el(sp).items():
                rhs[el] = rhs.get(el, 0) + v * k
        if _nonho(lhs) != _nonho(rhs):
            rep("元素", f"{solid} + {w}H2O -> {x}{cat} + {y}{an}："
                       f"非 H/O 元素不等 {_nonho(lhs)} vs {_nonho(rhs)}", bad)
        elif lhs.get("H", 0) != rhs.get("H", 0):
            rep("H计数", f"{solid} 溶解式 H 不平（w={w}）", bad)
        dq = x * _chg(cat) + y * _chg(an)
        if dq:
            rep("电荷", f"{solid} -> {x}{cat} + {y}{an}：电荷 {dq}", bad)
        if pk <= 0:
            rep("值域", f"{solid} pKsp={pk} ≤ 0", bad)
    for k, v in seen.items():
        if len(v) > 1:
            rep("重复", f"固相 {k} 出现 {len(v)} 次：{v}", dup)
    return len(rows)


def audit_beta(rows, bad, dup) -> int:
    seen = defaultdict(list)
    for e in rows:
        c, cen, lig, nu = e["complex"], e["center"], e["ligand"], e["nu"]
        seen[(cen, lig, nu)].append(c)
        ec, en, el_ = _el(c), _el(cen), _el(lig)
        if not ec or not en or not el_:
            rep("解析失败", f"beta {c}", bad)
            continue
        rhs = {}
        for sp, k in ((cen, 1), (lig, nu)):
            for el, v in _el(sp).items():
                rhs[el] = rhs.get(el, 0) + v * k
        if ec != rhs:
            rep("元素", f"{cen} + {nu}{lig} -> {c}：{rhs} vs {ec}", bad)
        dq = _chg(c) - _chg(cen) - nu * _chg(lig)
        if dq:
            rep("电荷", f"{cen} + {nu}{lig} -> {c}：电荷差 {dq}", bad)
        if e["logb"] <= 0:
            rep("值域", f"{c} logβ={e['logb']} ≤ 0", bad)
    for k, v in seen.items():
        if len(v) > 1:
            rep("重复", f"配位键 {k} 出现 {len(v)} 次：{v}", dup)
    return len(rows)


def _solve_half(ox: str, red: str, kmax: int = 6):
    """`k_ox·ox + a·H⁺ + b·H2O -> k_red·red` 的元素空间求解。

    电对表只存 (ox, red) 两个物种，**不带系数**：`Cl_2/Cl^-` 实为
    `Cl2 + 2e -> 2Cl⁻`。首版假定两侧系数都是 1，把一大半电对误报成
    "元素不守恒"。返回 (k_ox, k_red, a, b, n) 列表（n 由电荷平衡给）。
    """
    eo, er = _el(ox), _el(red)
    els = set(eo) | set(er)
    out = []
    for ko in range(1, kmax + 1):
        for kr in range(1, kmax + 1):
            # O 平衡定 b，H 平衡定 a（a 可负 = 右侧出 H⁺）
            b = kr * er.get("O", 0) - ko * eo.get("O", 0)
            a = (kr * er.get("H", 0) - ko * eo.get("H", 0)) - 2 * b
            ok = True
            for el in els:
                lhs = (ko * eo.get(el, 0)
                       + a * (1 if el == "H" else 0)
                       + b * (2 if el == "H" else 1 if el == "O" else 0))
                rhs = kr * er.get(el, 0)
                if lhs != rhs:
                    ok = False
                    break
            if not ok:
                continue
            n = ko * _chg(ox) + a - kr * _chg(red)
            if n > 0:
                out.append((ko, kr, a, b, n))
        if out:
            break
    return out


def audit_couples(rows, bad, dup, info) -> int:
    seen = set()
    mismatch = 0
    released = 0
    for e in rows:
        ox, red, n = e["ox"], e["red"], e["n"]
        if (ox, red) in seen:
            rep("重复", f"电对 ({ox}, {red}) 重复", dup)
        seen.add((ox, red))
        if not _el(ox) or not _el(red):
            rep("解析失败", f"couple {ox}/{red}", bad)
            continue
        # 配体/阴离子释放型：ox 比 red 多出的元素在 pair 里无处安放，但它们
        # 是**作为自由离子释放**的（AgCl + e -> Ag + Cl⁻、[Ag(NH3)2]⁺ + e ->
        # Ag + 2NH3、[PtCl6]²⁻ + e -> Pt + 6Cl⁻）——pair 模式只写 ox/red，
        # 释放物是隐含的。这类不是数据错误，单列信息项。
        extra = {el: v for el, v in _el(ox).items() if v > _el(red).get(el, 0)}
        extra = {el: v - _el(red).get(el, 0) for el, v in extra.items()}
        extra.pop("H", None)
        extra.pop("O", None)
        sols = _solve_half(ox, red)
        if not sols and extra:
            released += 1
            info.append(f"[释放型] {ox} -> {red}：释放 {extra}（pair 不含释放物）")
            continue
        if not sols:
            rep("无解", f"{ox} -> {red}：元素空间配不出半反应", bad)
            continue
        ns = {s[4] for s in sols}
        if n not in ns:
            mismatch += 1
            rep("电子数", f"{ox} -> {red}：电荷平衡给 n={sorted(ns)}，表内 n={n}"
                          f"（配平 {sols[0][:2]}，n 按式量或最小整数比的等价表达）",
                bad)
    print(f"   （释放型电对 {released} 条——pair 模式只写 ox/red，配体/阴离子"
          f"释放是隐含的，单列信息项）")
    if mismatch:
        print(f"   （电子数不符 {mismatch} 条——多为 n 按式量 vs 按最小整数比"
              f"的等价表达，需逐条判读）")
    return len(rows)


def main() -> None:
    bad: list[str] = []
    dup: list[str] = []
    info: list[str] = []
    counts = {}
    counts["pka"] = audit_pka(load("pka.json"), bad, dup)
    counts["ksp"] = audit_ksp(load("ksp.json"), bad, dup)
    counts["beta"] = audit_beta(load("beta.json"), bad, dup)
    counts["couples"] = audit_couples(load("couples.json"), bad, dup, info)
    th = load("thermo.json")
    counts["thermo"] = len(th)
    import re as _re
    for k in th:
        if _re.search(r"\((aq|s|l)\)$", k):
            rep("死键", f"thermo 键 {k}：代码只查 (g) 后缀，"
                        f"(aq)/(s)/(l) 从未被读取", bad)
        elif not _el(k):
            rep("解析失败", f"thermo 键 {k}", bad)

    print("== 审计规模 ==")
    for k, v in counts.items():
        print(f"   {k:9s} {v}")
    print(f"\n== 硬错误 {len(bad)} 条 ==")
    for m in bad[:60]:
        print("  " + m)
    if len(bad) > 60:
        print(f"  … 余 {len(bad) - 60} 条")
    print(f"\n== 重复/冲突 {len(dup)} 条 ==")
    for m in dup[:40]:
        print("  " + m)
    if info:
        print(f"\n== 信息项 {len(info)} 条（非错误，按设计）==")
        for m in info[:8]:
            print("  " + m)
        if len(info) > 8:
            print(f"  … 余 {len(info) - 8} 条")


if __name__ == "__main__":
    main()
