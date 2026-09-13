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
        # 看它是"两次穿越"（真多根）还是"平台/台阶"（pH 机器分支造成的伪结构）
        print("\n   形状诊断（f 细扫，+ 表示 f>0）：")
        for r in sorted(n_hi, key=lambda r: -abs(r["x_ill"] - r["x_bis"])
                        / max(1.0, r["x_max"]))[:5]:
            hi = 1.3 * max(r["x_bis"], r["x_ill"], 1e-30)
            n = 61
            row = "".join("+" if r["f"](hi * k / (n - 1)) > 0 else "." for k in range(n))
            print(f"     x∈[0,{hi:.4g}]  {row}")
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


if __name__ == "__main__":
    main()
