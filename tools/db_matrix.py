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
        c for c in ksp_cat if c not in explicit1)

    # 4) 电对/β 里的**化合物**固相缺溶解度数据（Co(OH)_3 型：有氧化通道却无
    #    沉淀平衡）。单质不算（金属/非金属单质本来就没有 Ksp）。
    rules["solid_no_ksp"] = rank(
        s for s in (cpl_all | beta_all)
        if s in solids and s not in ksp_all and not _is_element(s))

    # 5) **化合物**固相缺 thermo ⟹ 进不了 Hess 审计（单质 ΔHf ≡ 0，不算缺口）
    rules["solid_no_thermo"] = rank(
        s for s in solids if s not in th_keys and not _is_element(s))

    # 6) β 配体缺 pKa/thermo（配体自身的酸碱行为与温度依赖都没有）
    rules["ligand_no_pka"] = rank(
        {e["ligand"] for e in beta} - pka_species - {"OH^-"})

    # 7) 表内物种总数与交集（信息项）
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
