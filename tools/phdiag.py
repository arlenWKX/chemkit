"""pH 机器（启发式分支机）vs 电荷平衡精确解（acidbase.charge_pH）全量对照。

目的（0.5.x 施工顺序 ③ 的量化前置）：在**引擎自己的退出账本**上比较两条
pH 求法，回答两个问题：

  1. 分歧有多大、集中在哪些用例（决定替换风险面）；
  2. 分歧例上谁是化学正确的一方（决定替换方向）。

对照口径：probe 记录的 `ledger` 已经过 `_respeciate_strong_acids` +
`_presentation_He`，`probe["pH"]` 即机器在**同一账本**上的输出——两者
输入完全相同，唯一差别是算法，故差值即纯算法差。

用法：
    python tools/phdiag.py [-n 20] [--cases 路径] [--only 前缀,前缀]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter

sys.path.insert(0, ".")

from chemkit.acidbase import charge_pH, build_families          # noqa: E402
from chemkit.core import charge_of                              # noqa: E402
from chemkit.data import load_tables                            # noqa: E402
from chemkit.engine import judge                                # noqa: E402
from chemkit.testsuit import load_cases                         # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=20)
    ap.add_argument("--cases", default=None)
    ap.add_argument("--only", default=None, help="逗号分隔的用例名前缀")
    ap.add_argument("--tol", type=float, default=0.05)
    ap.add_argument("--classify", action="store_true",
                    help="对 |ΔpH| ≥ 1 的分歧按「质子条件前提」分类")
    args = ap.parse_args()

    T = load_tables()
    fams = build_families(T)
    # 金属水解阳离子（Ksp 对里阴离子为 OH^-）：机器的 hyd_map 快捷通道
    hyd_cats = {e["pair"][0] for e in T.ksp if e["pair"][1] == "OH^-"}
    cases = load_cases(args.cases)
    if args.only:
        pre = tuple(x.strip() for x in args.only.split(",") if x.strip())
        cases = [c for c in cases if c["name"].startswith(pre)]

    rows = []
    none_n = 0
    # 储库物种集合：Ksp 的阳/阴离子对 + 气相物种——账本里出现它们意味着
    # pH 可能由 Ksp/逸度（外部储库）决定，而不是账本内部的质子条件。
    res_sp: set = set()
    for e in T.ksp:
        res_sp.update(e["pair"])
    res_sp |= set(getattr(T, "gases", ()) or ())
    for c in cases:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        probe: dict = {}
        judge(subs, cond, T, _probe=probe)
        if "ledger" not in probe:
            continue                        # OVERRIDE 直出，无账本画像
        V = float(cond.get("V_L", 1.0) or 1.0)
        T_K = float(cond.get("T_K", 298.15) or 298.15)
        led = probe["ledger"]
        He = probe["H_excess"]
        p_m = probe["pH"]
        c_H = float(cond.get("c_H", 0.0) or 0.0)
        # 正确口径：账本自身电荷平衡定 pH（H_excess 是待求量的另一半，
        # 传进去即双重记账）；c_H 是无阴离子的强酸条件量，须补回。
        # fast=True：阻尼 Newton（与二分同根，tools/phcmp.py 全量对拍
        # 最大 |ΔpH| < 1e-12）——本工具要扫全库，省掉二分的一半时间。
        p_c = charge_pH(led, V, T, T_K, c_H=c_H, fast=True)
        if p_c is None:
            none_n += 1
        net = sum(charge_of(s) * m for s, m in led.items()
                  if not s.startswith("__") and s not in ("H^+", "OH^-"))
        fam_hit = sum(1 for s in led if s in fams)
        # 分类用画像
        solid_n = sum(1 for s, m in led.items()
                      if m > 1e-6 and s in T.solids)
        gas_n = sum(1 for s, m in led.items()
                    if m > 1e-6 and s in res_sp and s not in T.solids)
        hyd_n = sum(1 for s, m in led.items() if m > 1e-6 and s in hyd_cats)
        rows.append({
            "name": c["name"], "V": V, "T_K": T_K, "He": He, "c_H": c_H,
            "machine": p_m, "charge": p_c,
            "d": None if p_c is None else round(p_c - p_m, 3),
            "families": fam_hit, "nsp": len(led),
            "cons": abs(net + He + c_H),
            "solid": solid_n, "gas": gas_n, "hyd": hyd_n,
            "note": (c.get("note") or "")[:70],
        })

    got = [r for r in rows if r["d"] is not None]
    big = sorted(got, key=lambda r: -abs(r["d"]))
    print(f"用例 {len(rows)}；charge_pH 无括号 {none_n}；可对照 {len(got)}")
    bad = [r for r in rows if r["cons"] > 1e-6]
    print(f"账本守恒违规 |Σz·n+He|>1e-6 : {len(bad)} 例"
          f"（最大 {max((r['cons'] for r in rows), default=0):.3g}）")
    hist = Counter()
    for r in got:
        a = abs(r["d"])
        k = ("<0.01" if a < 0.01 else "<0.05" if a < 0.05 else "<0.2" if a < 0.2
             else "<0.5" if a < 0.5 else "<1" if a < 1 else ">=1")
        hist[k] += 1
    for k in ("<0.01", "<0.05", "<0.2", "<0.5", "<1", ">=1"):
        print(f"  |ΔpH| {k:6s} : {hist[k]:5d}")
    print(f"  分歧 >{args.tol} : "
          f"{sum(1 for r in got if abs(r['d']) > args.tol)}")
    print(f"\n== 分歧最大 {args.n} 例（Δ = charge − machine）==")
    for r in big[:args.n]:
        print(f"  Δ{r['d']:+8.3f}  机器 {r['machine']:6.2f} → 电荷 {r['charge']:6.2f} "
              f" He={r['He']:+10.6f} 族{r['families']}/物种{r['nsp']:3d}  {r['name'][:52]}")
    if args.classify:
        _classify(got, args.n)
    print(f"\n== 一致例（|Δ|<0.01）样例 ==")
    for r in [x for x in got if abs(x["d"]) < 0.01][:args.n]:
        print(f"  pH {r['machine']:6.2f}  He={r['He']:+10.6f}  {r['name'][:56]}")
    return 0


# 分歧分类：判定"精确质子条件在哪些用例上根本无资格发言"。
# 记账（§7 N-5/N-6）：质子条件的成立前提是**账本是闭系、电中性、
# 且 pH 相关物种全部落在 pKa 族里**。三条前提各自对应一类分歧：
BUCKETS = (
    ("B1 账本不守恒", "|Σz·n+He+c_H| > 1e-6：投料本身不电中性（裸离子）"
                      "或记账户型不一致——质子条件无对象"),
    ("B2 无族可分配", "账本里没有 pKa 族成员：质子条件没有可再分配的形态，"
                      "pH 只能由水解/两性/缓冲启发式或外部储库给"),
    ("B3 族且无储库", "有 pKa 族、无固相/水解阳离子：**质子条件唯一"
                      "有资格的档**——这里的分歧必须逐条判化学对错"),
    ("B4 族+储库·直读档", "有族也有储库，且 |He| ≥ 1e-3：机器应走自由强酸/"
                          "强碱直读（pH 只由 He 定）"),
    ("B5 族+储库·缓冲档", "有族也有储库，He 小：机器走水解/两性/缓冲启发式，"
                          "储库（Ksp/逸度）不在账本的质子条件里"),
)

# 储库判据**只算固相与水解阳离子**，不算"气体"：
# 溶解态 CO₂/H₂S/SO₂/NH₃/H₂Se 既是 T.gases 成员、也是 pKa 族成员，它们
# **在账本里**（量由走步与逸出步决定），是质子条件的一等公民。首版分类把
# "账本里有气体"当外部储库，误把 175 例碳酸/硫化物缓冲体系排除在 B3 之外
# （本轮自查发现——分类器自己就是一处"近似掩盖问题"）。真正在账本之外的
# 自由度只有两个：固相持有质量、金属阳离子的水解由 Ksp 定（hyd_map）。
RESERVOIR_SP = "T.solids ∪ Ksp-OH 阳离子"


def _bucket(r: dict) -> str:
    if r["cons"] > 1e-6:
        return "B1 账本不守恒"
    if r["families"] == 0:
        return "B2 无族可分配"
    if not (r["solid"] or r["hyd"]):
        return "B3 族且无储库"
    if abs(r["He"]) >= 1e-3:
        return "B4 族+储库·直读档"
    return "B5 族+储库·缓冲档"


def _classify(got: list, n: int) -> None:
    """按"质子条件的前提"给分歧分档，并给两档阈值的计数表。

    这张表就是"精确 pH 闭环能不能入联立求解器"的判据来源：**只有当
    B3（族且无储库）在分歧里有实质占比时，替换 pH 定义才有收益**——
    否则分歧全部落在质子条件**无资格**的档（非电中性账本 / 外部储库）。
    """
    soft = [r for r in got if abs(r["d"]) > 0.05]
    hard = [r for r in got if abs(r["d"]) >= 1.0]
    print(f"\n== 分歧分类（>0.05：{len(soft)} 例；≥1：{len(hard)} 例）==")
    print(f"  {'档':22s} {'>0.05':>7s} {'≥1':>7s}   含义")
    for name, desc in BUCKETS:
        ns = sum(1 for r in soft if _bucket(r) == name)
        nh = sum(1 for r in hard if _bucket(r) == name)
        print(f"  {name:22s} {ns:7d} {nh:7d}   {desc}")
    for name, _desc in BUCKETS:
        sub = [r for r in hard if _bucket(r) == name]
        if not sub:
            continue
        print(f"\n  --- [{name}] ≥1 档 {len(sub)} 例 ---")
        for r in sorted(sub, key=lambda r: -abs(r["d"]))[:n]:
            print(f"     Δ{r['d']:+8.3f} 机器 {r['machine']:6.2f} → 电荷 "
                  f"{r['charge']:6.2f}  He={r['He']:+9.5f} 族{r['families']} "
                  f"固{r['solid']} 水解{r['hyd']} 气{r['gas']}  {r['name'][:44]}")
    print(f"\n  --- [B3 族且无储库] >0.05 档全列（唯一有资格的档）---")
    for r in sorted((r for r in soft if _bucket(r) == "B3 族且无储库"),
                    key=lambda r: -abs(r["d"])):
        print(f"     Δ{r['d']:+8.3f} 机器 {r['machine']:6.2f} → 电荷 "
              f"{r['charge']:6.2f}  He={r['He']:+9.5f} 族{r['families']} "
              f"物种{r['nsp']:3d}  {r['name'][:48]}")


if __name__ == "__main__":
    sys.exit(main())
