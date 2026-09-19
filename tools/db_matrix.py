# -*- coding: utf-8 -*-
"""数据库**跨表完备性**审计（`tools/data_audit.py` 管的是表内自洽，本工具管
"某个物种在这张表里有、在那张表里没有"）。

## 为什么需要（第 144 轮实测两例）
* `HSO_4^-` 在 `pka.json` 与 `thermo.json` **两表都没有** ⟹ 硫酸根在整个模型里
  被当成"全离解旁观离子"，而 pKa₂(HSO₄⁻) = 1.99 意味 pH 3 时约 9%、pH 1 时约 91%
  以 HSO₄⁻ 存在（补这一条实测**87 例**失败：既需要"酸碱质子化不算反应事件"的
  呈现语义，也需要约 50 例标准按分析形态重裁）。
* `Co(OH)_3` 不在 `thermo.json` 也不在 `ksp.json` ⟹ 表里没有 Co(OH)₃/Co(OH)₂
  电对，第 143 轮给 O₂/H₂O 补的 `closed_except_red: ["Co(OH)_2"]` 成了**死数据**
  （析出的 Co(OH)₂ 在引擎里没有任何氧化通道）。

两例都是"静默缺口"——不报错、不违反守恒，只是化学少了一块。本工具把它们
**机械地列出来**，并按"测试用例触达数"排序，使数据库扩充有优先级依据。

用法：
    python tools/db_matrix.py                # 全部规则 + 摘要
    python tools/db_matrix.py --rule ksp_cat_no_beta --top 30
    python tools/db_matrix.py --json .tmp_dbmatrix.json
"""
import argparse
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
sys.path.insert(0, REPO)

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.core import elements_of                     # noqa: E402
from chemkit.data import load_tables                        # noqa: E402


def _is_element(sp: str) -> bool:
    """单质（只含一种元素）——`ΔHf ≡ 0`、无 Ksp 可言，不算缺口。"""
    try:
        return len(elements_of(sp)) <= 1
    except Exception:                                       # noqa: BLE001
        return False


# 强酸阴离子：其共轭酸 pKa ≤ 0，引擎按"强酸直读"建模（`_acids_map` 明确
# 跳过 pKa ≤ 0 的物种），**本来就不该在 pKa 表里**。与 `_is_element` 同一
# 原则：缺的是"引擎不建模的东西"就不是缺口（D10：先读设计意图再记缺口）。
# 之前把 Cl⁻/Br⁻/I⁻/NO₃⁻/ClO₄⁻ 报成"配体缺 pKa"，把 161 例的假缺口顶到榜首，
# 掩盖了真缺口（第 152 轮）。
_STRONG_ANION = frozenset({
    "Cl^-", "Br^-", "I^-", "NO_3^-", "ClO_4^-", "ClO_3^-", "BrO_3^-",
    "IO_3^-", "MnO_4^-", "HSO_4^-", "SO_4^{2-}", "ClO^-", "ClO_2^-",
    "N_3^-", "SCN^-", "CN^-",
})


def _read(name: str):
    p = os.path.join(REPO, "chemkit", "data", name)
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def _case_reach() -> dict:
    """物种 -> 触达它的用例数（用例 subs/eq/断言里出现即算一次）。

    粗粒度但足够做优先级：只按物种名字符串在用例 JSON 文本里的出现计数。
    """
    cases = _read("tests.json")
    reach: dict[str, int] = {}
    for c in cases:
        blob = json.dumps(c, ensure_ascii=False)
        for m in set(re.findall(r"[A-Za-z][A-Za-z0-9_^{}\-\[\]()+]*", blob)):
            reach[m] = reach.get(m, 0) + 1
    return reach


def main() -> int:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("--rule", default=None, help="只跑某条规则")
    ap.add_argument("--top", type=int, default=15, help="每条规则最多列几条")
    ap.add_argument("--json", default=None, help="把结果写成 JSON")
    args = ap.parse_args()

    T = load_tables()
    pka, thermo = _read("pka.json"), _read("thermo.json")
    beta, ksp, couples = _read("beta.json"), _read("ksp.json"), _read("couples.json")
    reach = _case_reach()

    th_keys = {k for k in thermo if not k.endswith("(g)")}
    pka_species = {e["acid"] for e in pka} | {e["base"] for e in pka}
    beta_center = {e["center"] for e in beta}
    beta_all = {e["complex"] for e in beta} | beta_center | {e["ligand"] for e in beta}
    ksp_cat = {e["pair"][0] for e in ksp}
    # **Ksp-OH 阳离子**：只有 `pair = [M, OH^-]` 才谈得上"一级水解"。
    # 第 154 轮修正：原实现把**任意** Ksp 条目的阳离子都算进来，于是
    # `H_2SiO_3 ⇌ H⁺ + SiO₃²⁻` 让 H⁺ 以 309 例高居榜首，Na⁺/K⁺ 也因
    # 碳酸盐/高氯酸盐上榜——而它们的氢氧化物是**易溶强碱**，根本没有
    # 可补的一级水解。判据必须与 `speciation.hyd_map` 同源（那里只遍历
    # `an == "OH^-"` 的条目）。
    ksp_oh_cat = {e["pair"][0] for e in ksp if e["pair"][1] == "OH^-"}
    ksp_all = {s for e in ksp for s in (e["pair"][0], e["pair"][1], e["solid"])}
    cpl_all = {s for e in couples for s in (e["ox"], e["red"])}
    solids = set(T.solids)

    def rank(items):
        return sorted(items, key=lambda s: (-reach.get(s, 0), s))

    rules: dict[str, list] = {}

    # 1) pKa 物种缺 thermo ⟹ 没有 van't Hoff（温度依赖退化）
    rules["pka_no_thermo"] = rank(
        s for s in pka_species if s not in th_keys and s != "H_2O")

    # 2) 电对成员缺 thermo ⟹ 该电对没有 dH（只能 ΔS≈0 回退）
    miss = set()
    for e in couples:
        if e["ox"] not in th_keys or e["red"] not in th_keys:
            miss.add(f'{e["ox"]}/{e["red"]}')
    rules["couple_no_thermo"] = sorted(
        miss, key=lambda s: -reach.get(s.split("/")[0], 0))

    # 3) Ksp-OH 阳离子但没有 nu=1 的 OH⁻ β ⟹ 仍走 Kh（Ksp 派生）复合近似
    explicit1 = {e["center"] for e in beta
                 if e.get("ligand") == "OH^-" and e.get("nu") == 1}
    rules["ksp_cat_no_beta1"] = rank(
        c for c in ksp_oh_cat if c not in explicit1)

    # 4) 电对/β 里的**化合物**固相缺溶解度数据（Co(OH)_3 型：有氧化通道却无
    #    沉淀平衡）。单质不算（金属/非金属单质本来就没有 Ksp）。
    rules["solid_no_ksp"] = rank(
        s for s in (cpl_all | beta_all)
        if s in solids and s not in ksp_all and not _is_element(s))

    # 5) **化合物**固相缺 thermo ⟹ 进不了 Hess 审计（单质 ΔHf ≡ 0，不算缺口）
    rules["solid_no_thermo"] = rank(
        s for s in solids if s not in th_keys and not _is_element(s))

    # 6) β 配体缺 pKa/thermo（配体自身的酸碱行为与温度依赖都没有）。
    #    强酸阴离子（Cl⁻/Br⁻/I⁻/NO₃⁻…）排除——它们的共轭酸 pKa ≤ 0，
    #    引擎走"强酸直读"，不该进 pKa 表（见 `_STRONG_ANION` 的说明）。
    rules["ligand_no_pka"] = rank(
        {e["ligand"] for e in beta} - pka_species - {"OH^-"} - _STRONG_ANION)

    # 7) **M–Cl 络合缺口**（第 156 轮）：水解常数是 I→0/无其它配体/无多聚的
    #    单核值，真实体系（1 M GaCl₃ ⟹ 3 M Cl⁻）里氯络合与多聚会显著降低自由
    #    M^{n+}、压低 H⁺ 释放 ⟹ 引擎比真实体系更酸。实测偏差随 Cl⁻ 负载排序
    #    （Ga 3 M 最大、In 次之、Rh 硝酸盐体系最小），证实缺的是**这个过程**
    #    而非常数取值。凡有 OH⁻ β₁（会走一级水解）却没有任何 Cl⁻ β 的中心离子。
    cl_centers = {e["center"] for e in beta if e.get("ligand") == "Cl^-"}
    rules["mcl_no_beta"] = rank(
        c for c in explicit1 if c not in cl_centers)

    # 8) **死数据候选**（第 165 轮，D4 的机械化）：化合物固相既不缺数据、
    #    也没被任何表引用、且**没有任何用例投料**它 ⟹ 当前无从被触及。
    #    与 `solid_no_thermo` 分开列：那 146 条里 141 条是真缺口（136 条被表
    #    引用 + 5 条被用例投料），只有 5 条落进这一类。
    used_tables = (ksp_all | beta_all | pka_species | cpl_all)
    feed = set()
    try:
        import json as _json
        _cs = _read("tests.json")
        _cs = _cs if isinstance(_cs, list) else (
            _cs.get("cases") or _cs.get("entries") or [])
        for _c in _cs:
            for _x in (_c.get("subs") or []):
                feed.add(_x[0])
    except Exception:                                       # noqa: BLE001
        pass
    rules["solid_unused"] = rank(
        s for s in solids
        if s not in used_tables and s not in feed and not _is_element(s))

    # 9) 表内物种总数与交集（信息项）
    allsp = th_keys | pka_species | beta_all | ksp_all | cpl_all
    rules["_summary"] = [f"物种并集 {len(allsp)}"]
    for nm, st in (("thermo", th_keys), ("pka", pka_species), ("beta", beta_all),
                   ("ksp", ksp_all), ("couples", cpl_all), ("solids", solids)):
        rules["_summary"].append(f"  {nm:8s} {len(st):4d}")

    only_one = [s for s in allsp
                if sum(s in st for st in (th_keys, pka_species, beta_all,
                                          ksp_all, cpl_all)) == 1]
    rules["_summary"].append(f"只出现在一张表里的物种 {len(only_one)}")

    shown = [args.rule] if args.rule else [k for k in rules if not k.startswith("_")]
    for name in shown:
        if name not in rules:
            print(f"!! 未知规则 {name}；可用：{sorted(k for k in rules if not k.startswith('_'))}")
            continue
        items = rules[name]
        print(f"\n== {name}（{len(items)} 条）==")
        for s in items[:args.top]:
            print(f"   {reach.get(s.split('/')[0], 0):4d} 例  {s}")
        if len(items) > args.top:
            print(f"   … 另有 {len(items) - args.top} 条")

    if args.rule is None:
        print("\n== 规模摘要 ==")
        for line in rules["_summary"]:
            print("  " + line)

    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump({k: v for k, v in rules.items()}, fh,
                      ensure_ascii=False, indent=1)
        print(f"\n[写出] {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
