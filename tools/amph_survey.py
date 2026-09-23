# -*- coding: utf-8 -*-
"""第 263 轮 · 分支 4 标签普查（**每一次** `estimate_pH` 调用都计数）。

目的：量出"两性中点"这条捷径在**热路径**（solve_extent 二分探针）上的
真实触发频率与代价面 —— 决定"把它换成精确解"的爆破半径。

做法：`speciation._tag` 是模块级函数、`estimate_state` 按全局名调用它，
故把 `_sp._tag` 换成计数器即可**零侵入**统计（不改 chemkit/ 一行）。

输出 `logs/amph_survey.json`：
  * `tags`      —— 全局标签直方图（全部探针）
  * `per_case`  —— 逐例：末态 ledger 上的 est/chg、触发过的标签集合
  * `amph_partner` —— "两性物种与其共轭伙伴同时在账"的例（第 263 轮的核心族）

⚠️ 串行运行（~200 s）。运行期间**冻结 chemkit/**。
"""
from __future__ import annotations

import collections
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

from chemkit import speciation as _sp            # noqa: E402
from chemkit.data import load_tables             # noqa: E402
from chemkit.speciation import (                 # noqa: E402
    charge_pH, estimate_state, pKw_of)


TAGS = collections.Counter()


def _counting_tag(why: str) -> None:
    TAGS[why] += 1


def main() -> int:
    _sp._tag = _counting_tag
    T = load_tables()
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    from chemkit.engine import judge

    # 两性物种集合与其"共轭酸/共轭碱"伙伴（与 estimate_state 同源口径）
    amph = set(T.pka_acid) & set(T.pka_base)
    conj_acid = {}       # amph -> 其共轭酸（来自 pka_base 侧）
    for b, entries in T.pka_base.items():
        if b not in amph:
            continue
        e1 = [e for e in entries if e["n"] == 1] or entries
        best = max(e1, key=lambda e: e["pka"])
        # pka_base[base] 的条目里 acid 就是它的共轭酸
        conj_acid[b] = best.get("acid")
    conj_base = {}       # amph -> 其共轭碱
    for a, entries in T.pka_acid.items():
        if a not in amph:
            continue
        e1 = min(entries, key=lambda e: e["pka"])
        conj_base[a] = e1.get("base")

    per_case = []
    before = sum(TAGS.values())
    for i, c in enumerate(cases, 1):
        subs = [{"name": s[0], "mol": float(s[1])} for s in c["subs"]]
        cond = c.get("cond") or {}
        pr = {}
        try:
            r = judge(subs, cond, T, _probe=pr)
        except Exception as exc:
            per_case.append({"name": c["name"], "err": repr(exc)})
            continue
        led = {e["name"]: e["mol"] for e in (r.get("final") or [])}
        V = float(cond.get("V_L", 1.0))
        T_K = float(cond.get("T_K", 298.15))
        He = float(pr.get("H_excess", 0.0) or 0.0)
        _sp.PH_TAGS = []
        try:
            ph, _vl, _hr = estimate_state(led, He, V, T, T_K)
            tags = list(_sp.PH_TAGS)
        finally:
            _sp.PH_TAGS = None
        try:
            ex = charge_pH(led, V, T, T_K)
        except Exception:
            ex = None
        # 共轭伙伴是否同时在账（元素量显著）
        part = []
        for sp, m in led.items():
            if m <= 1e-9 or sp not in amph:
                continue
            for other, kind in ((conj_acid.get(sp), "共轭酸"),
                                (conj_base.get(sp), "共轭碱")):
                if other and led.get(other, 0.0) > 1e-9:
                    part.append(f"{sp}+{other}({kind})")
        per_case.append({
            "name": c["name"], "suite_pH": r.get("final_pH"),
            "pH_est": round(ph, 4),
            "pH_chg": (round(ex, 4) if ex is not None else None),
            "d": (round(ph - ex, 4) if ex is not None else None),
            "resid": abs(pr.get("resid_live") or 0.0),
            "tags": tags, "amph_partner": part,
            "n_species": len(led),
        })
        if i % 100 == 0:
            print(f"  {i}/{len(cases)} …", flush=True)

    out = {
        "tags": dict(TAGS.most_common()),
        "n_calls": sum(TAGS.values()),
        "per_case": per_case,
    }
    with io.open(os.path.join(ROOT, "logs", "amph_survey.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(f"\n总标签计数 {sum(TAGS.values())}（本进程首例前基线 {before}）")
    for k, v in TAGS.most_common(20):
        print(f"   {k:<28} {v}")
    tagd = collections.Counter()
    for r in per_case:
        for t in r.get("tags") or []:
            tagd[t] += 1
    print("\n末态分支归属（逐例）：")
    for k, v in tagd.most_common():
        print(f"   {k:<28} {v}")
    bad = [r for r in per_case if r.get("d") is not None and abs(r["d"]) > 0.1]
    print(f"\n末态 |Δ(est, chg)| > 0.1 的例：{len(bad)}")
    for r in sorted(bad, key=lambda r: -abs(r["d"]))[:25]:
        print(f"   {str(r['name'])[:34]:<34} est={r['pH_est']:>7} "
              f"chg={r['pH_chg']:>7} Δ={r['d']:>7} resid={r['resid']:>6.3f} "
              f"{'/'.join(r['tags'])}")
    print(f"\n两性物种 + 共轭伙伴同时在账的例："
          f"{sum(1 for r in per_case if r.get('amph_partner'))}")
    print("\n-> logs/amph_survey.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
