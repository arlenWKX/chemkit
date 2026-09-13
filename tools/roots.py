"""求根审计：保括号假位法（Illinois）能否替换主求解的二分？

**为什么问这个**（architecture §7 W-5）：二分在主求解分支固定迭代
（容差落地后 ~35 次 f 求值/步，f = pH 估计 + S_of 一遍），而全量
`S_of` 调用 844k 中约四成花在这里——这是 §7 W-5 之后剩下的最大杠杆。
Illinois 在光滑单调括号上 ~10 次求值即可到同一精度。

**唯一的语义风险是根选定**：括号内 f 非单调（多根）时，快方法可能收敛到
另一个根，而"选中哪个根"等于"实际执行多少"等于化学。断言锁化学不等于
允许换根——换根就是换化学，必须避免。故本工具在同一括号上审计两问：

  ① 括号内 f 是否单调（33 点网格的符号翻转次数）；
  ② 两法选出的根差多大（`|x_ill − x_bis| / max(1, x_max)`）。

用法：
    python tools/roots.py              # 每 500 个括号抽样一个
    python tools/roots.py -e 50        # 抽样间隔
    python tools/roots.py --cases NaClO T78   # 只跑这些用例，**每个**括号都审计
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from chemkit import engine                        # noqa: E402
from chemkit.data import load_tables              # noqa: E402
from chemkit.testsuit import load_cases           # noqa: E402


def _pct(xs: list[float], q: float) -> float:
    if not xs:
        return 0.0
    return xs[min(len(xs) - 1, int(len(xs) * q))]


def main() -> None:
    every = 500
    if "-e" in sys.argv:
        every = int(sys.argv[sys.argv.index("-e") + 1])
    want: tuple[str, ...] = ()
    if "--cases" in sys.argv:
        want = tuple(sys.argv[sys.argv.index("--cases") + 1:])
        every = 1

    T = load_tables()
    cases = load_cases(None)
    if want:
        cases = [c for c in cases if c["name"].startswith(want)]
    engine.ROOT_AUDIT = {"n": 0, "every": every, "rec": [], "case": None}
    print(f"审计：{len(cases)} 例，每 {every} 个主求解括号抽 1 个"
          f"（33 点网格 + Illinois 复算）…")
    for c in cases:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        engine.ROOT_AUDIT["case"] = c["name"]
        engine.judge(subs, c.get("cond") or {"V_L": 1.0}, T)
    rec = engine.ROOT_AUDIT["rec"]
    brackets = engine.ROOT_AUDIT["n"]
    engine.ROOT_AUDIT = None

    solid = [r for r in rec if not r["micro"] and r["f0"] > 0]
    micro = [r for r in rec if r["micro"]]
    bad = [r for r in rec if r["f0"] <= 0]
    print(f"\n主求解括号总数 {brackets}；抽样 {len(rec)}"
          f"（完整二分档 {len(solid)}、微步快通道档 {len(micro)}）"
          f"；f(0)≤0 的非法括号 {len(bad)}（f(0)={len(bad)} 例，"
          f"二分在该括号上不成立）")

    if not solid:
        return
    nb = sum(r["n_bis"] for r in solid)
    ni = sum(r["n_ill"] for r in solid)
    print(f"\n① 求值次数（完整二分档）：二分 {nb} → Illinois {ni}"
          f"  = {nb / max(ni, 1):.2f}×  "
          f"（均值 {nb / len(solid):.1f} → {ni / len(solid):.1f} 次/括号）")

    flips = sorted(r["flips"] for r in solid)
    multi = [r for r in solid if r["flips"] > 1]
    print(f"\n② 括号内多根（33 点网格符号翻转 >1）：{len(multi)}/{len(solid)}"
          f" = {len(multi) / len(solid) * 100:.2f}%"
          f"；翻转数 p50={_pct(flips, 0.5)} p90={_pct(flips, 0.9)}"
          f" max={flips[-1]}")

    d = sorted(abs(r["x_ill"] - r["x_bis"]) / max(1.0, r["x_max"]) for r in solid)
    n_loose = [r for r in solid
               if abs(r["x_ill"] - r["x_bis"]) / max(1.0, r["x_max"]) > 1e-9]
    n_hi = [r for r in solid
            if abs(r["x_ill"] - r["x_bis"]) / max(1.0, r["x_max"]) > 1e-6]
    print(f"\n③ 两法根差 / max(1,x_max)：p50={_pct(d, 0.5):.2e}"
          f" p90={_pct(d, 0.9):.2e} p99={_pct(d, 0.99):.2e} max={d[-1]:.2e}")
    print(f"   > 1e-9（远超走步判据所需精度）：{len(n_loose)} 例"
          f"；> 1e-6（会改变呈现量级）：{len(n_hi)} 例")

    if n_hi:
        print("\n   换根的括号（>1e-6）：")
        for r in sorted(n_hi, key=lambda r: -abs(r["x_ill"] - r["x_bis"])
                        / max(1.0, r["x_max"]))[:15]:
            print(f"     Δ/max={abs(r['x_ill'] - r['x_bis']) / max(1.0, r['x_max']):.2e}"
                  f" 翻转={r['flips']} x_max={r['x_max']:.3g}"
                  f" 二分={r['x_bis']:.6g} Illinois={r['x_ill']:.6g}"
                  f"  {str(r['case'])[:26]} | {str(r['eq'])[:52]}")
        # 形状诊断：把最差几例的 f 在 [0, 1.3·两法较大根] 上细扫一遍——
        # 看它是"两次穿越"（真多根）还是"平台/台阶"（pH 机器分支造成的伪结构）；
        # 同时给出 pH 的斜率符号行（对齐同一网格）：若 f 的口袋正好落在 pH
        # 的"折点"上，就证实口袋来自 pH 机器而非化学。
        print("\n   形状诊断（f 细扫，+ 表示 f>0；下一行是 pH 斜率：+升 -降 .平）：")
        for r in sorted(n_hi, key=lambda r: -abs(r["x_ill"] - r["x_bis"])
                        / max(1.0, r["x_max"]))[:5]:
            hi = 1.3 * max(r["x_bis"], r["x_ill"], 1e-30)
            n = 61
            row = "".join("+" if r["f"](hi * k / (n - 1)) > 0 else "."
                          for k in range(n))
            print(f"     x∈[0,{hi:.4g}]  {row}")
            tr = r.get("trace")
            if tr:
                phs = [t[1] for t in tr if t[1] is not None]
                slope = "".join(
                    "+" if b - a > 0.02 else ("-" if a - b > 0.02 else ".")
                    for a, b in zip(phs, phs[1:]))
                print(f"        pH 斜率     {slope}")
                jump = max(((abs(b - a), i, a, b)
                            for i, (a, b) in enumerate(zip(phs, phs[1:]))),
                           default=(0, 0, 0, 0))
                print(f"        pH {phs[0]:.2f}→{phs[-1]:.2f}；最大单步跳变 "
                      f"{jump[0]:.2f} @ 第 {jump[1]}/60 点（{jump[2]:.2f}→{jump[3]:.2f}）")
                # 分支序列：跳变处换的是哪条返回路径（§7 X 的归因）
                tags = [t[3] if len(t) > 3 else "?" for t in tr]
                seq: list[str] = []
                for k, tg in enumerate(tags):
                    if k == 0 or tg != tags[k - 1]:
                        seq.append(f"{tg}@{k}")
                print(f"        分支序列    {' → '.join(seq[:12])}")
                if jump[0] > 0.2:
                    j = jump[1]
                    print(f"        跳变归因    第 {j} 点 {tags[j]} → 第 {j + 1} 点 "
                          f"{tags[j + 1] if j + 1 < len(tags) else '?'}")
            print(f"       二分根 {r['x_bis']:.6g} @ {r['x_bis'] / hi * (n - 1):.1f}/60"
                  f"；Illinois 根 {r['x_ill']:.6g} @ {r['x_ill'] / hi * (n - 1):.1f}/60"
                  f"  [{str(r['case'])[:22]}]")

    ok = [r for r in solid if r["flips"] <= 1]
    if ok:
        d_ok = sorted(abs(r["x_ill"] - r["x_bis"]) / max(1.0, r["x_max"])
                      for r in ok)
        n_ok_bad = sum(1 for x in d_ok if x > 1e-6)
        nb_ok = sum(r["n_bis"] for r in ok)
        ni_ok = sum(r["n_ill"] for r in ok)
        print(f"\n④ 只看网格单调的括号（可安全用快方法的那一档）："
              f"{len(ok)}/{len(solid)} = {len(ok) / len(solid) * 100:.1f}%；"
              f"其中根差 >1e-6 的 {n_ok_bad} 例；"
              f"求值 {nb_ok} → {ni_ok} = {nb_ok / max(ni_ok, 1):.2f}×")

    # ⑤ pH(x) 连续性普查：口袋只是"pH 不连续"的一个后果，先把病本身量出来
    census(rec, label="pH 连续性与换根普查")


def census(rec: list, label: str = "pH 连续性普查") -> dict:
    """pH(x) 连续性 + 换根普查（`tools/phcont.py` 复用；返回统计字典）。

    只统计**完整二分档**（微步快通道档的括号本就注定丢弃）。
    """
    solid = [r for r in rec if not r["micro"] and r["f0"] > 0]
    out: dict = {"n": len(solid)}
    traced = [r for r in solid if r.get("trace")]
    jumps = []
    for r in traced:
        phs = [t[1] for t in r["trace"] if t[1] is not None]
        if len(phs) < 2:
            continue
        jumps.append((max(abs(b - a) for a, b in zip(phs, phs[1:])), r))
    if not jumps:
        return out
    jumps.sort(key=lambda t: -t[0])
    vals = sorted(j for j, _ in jumps)
    n02 = sum(1 for j in vals if j > 0.2)
    n10 = sum(1 for j in vals if j > 1.0)
    n30 = sum(1 for j in vals if j > 3.0)
    div = [r for r in solid
           if abs(r["x_ill"] - r["x_bis"]) / max(1.0, r["x_max"]) > 1e-6]
    out.update(n_ph=len(vals), jump_p50=_pct(vals, 0.5),
               jump_p90=_pct(vals, 0.9), jump_max=vals[-1],
               n_gt02=n02, n_gt10=n10, n_gt30=n30, n_div=len(div))
    print(f"\n⑤ {label}（{len(vals)} 个括号，61 点网格内相邻点的最大 pH 变化）："
          f"p50={_pct(vals, 0.5):.3f} p90={_pct(vals, 0.9):.3f} max={vals[-1]:.2f}")
    print(f"   >0.2：{n02} 例（{n02 / len(vals) * 100:.1f}%）；"
          f">1.0：{n10} 例（{n10 / len(vals) * 100:.1f}%）；"
          f">3.0：{n30} 例（{n30 / len(vals) * 100:.1f}%）")
    print(f"   换根（两法根差 >1e-6·max(1,x_max)）：{len(div)}/{len(solid)}"
          f" = {len(div) / max(len(solid), 1) * 100:.1f}%")
    # 电荷自洽性：精确质子条件（charge_pH）的**前提**。不一致 ⟹ 账本内部
    # 电荷平衡不是真实约束（§7 F/N 的 L08 反例），统一方程不能盲用。
    res = sorted(abs(t[4]) for r in traced for t in r["trace"] if len(t) > 4)
    if res:
        n6 = sum(1 for v in res if v > 1e-6)
        n3 = sum(1 for v in res if v > 1e-3)
        out.update(resid_p50=_pct(res, 0.5), resid_max=res[-1],
                   n_res_gt6=n6, n_res_gt3=n3)
        print(f"   探头态电荷自洽 |Σz·n + He|：p50={_pct(res, 0.5):.2e} "
              f"p90={_pct(res, 0.9):.2e} max={res[-1]:.2e}；"
              f">1e-6 的 {n6}/{len(res)} = {n6 / len(res) * 100:.1f}%，"
              f">1e-3 的 {n3}（{n3 / len(res) * 100:.1f}%）")
    print("   跳变最大的 12 例：")
    for j, r in jumps[:12]:
        print(f"     ΔpH={j:5.2f}  x_max={r['x_max']:.3g}"
              f" 翻转={r['flips']}  根差={abs(r['x_ill'] - r['x_bis']) / max(1.0, r['x_max']):.1e}"
              f"  {str(r['case'])[:24]} | {str(r['eq'])[:44]}")
    return out


if __name__ == "__main__":
    main()
