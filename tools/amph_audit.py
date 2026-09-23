# -*- coding: utf-8 -*-
"""第 263 轮 · `estimate_pH` 分支归属普查 + 两性中点分支的合法性检验。

问题（第 263 轮起点）：3 条我自己加的 H₂S 半中和锚
（TC1/TJ3/TC3）实测 `pH 7.35` 而 `pH_solver 11.084` ——
**呈现 pH 对、走步内 pH 机器错 3.7 个单位**，且带 3.5 的残差。

本工具做三件事：

1. **分支归属**：对给定 ledger 直接调 `estimate_state`，打开 `PH_TAGS` 与
   `PH_SRC`，打印"赢下 pH 的是哪一支、哪些来源参与竞争"。
2. **精确对照**：同一 ledger 调 `charge_pH`（账本电荷平衡的精确解，
   第 260–262 轮已验证 ±0.0005）作为独立判据。
3. **合法性判据**：两性中点 `(pKa1+pKa2)/2` 的推导前提是
   **该两性物种是唯一质子条件物种**（由 `[H₂A] = [A²⁻]` 得 `h² = Ka1·Ka2`）。
   若账本里**已经存在该两性物种的共轭酸**（`conj[amph]`）且量与它可比，
   前提被破坏 ⟹ 中点式不适用。

用法：
  python tools/amph_audit.py --demo            # 固定几个已知体系
  python tools/amph_audit.py --suite           # 全库普查（走 judge 取末态 ledger）
"""
from __future__ import annotations

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
from chemkit.candidates import H_ION, WATER      # noqa: E402


def probe(T, ledger: dict, He: float, V: float, T_K: float,
          label: str = "") -> dict:
    """返回该 (ledger, He) 上的分支归属 + 三路 pH。"""
    pKw = pKw_of(T_K)
    _sp.PH_TAGS = []
    _sp.PH_SRC = []
    try:
        ph, vled, He_res = estimate_state(ledger, He, V, T, T_K)
        tags = list(_sp.PH_TAGS)
        src = list(_sp.PH_SRC)
    finally:
        _sp.PH_TAGS = None
        _sp.PH_SRC = None
    try:
        ex = charge_pH(ledger, V, T, T_K)
    except Exception as exc:                     # pragma: no cover
        ex = None
        tags.append(f"charge_pH 抛异常: {exc!r}")
    out = {
        "label": label, "pH_est": round(ph, 4),
        "pH_chg": (round(ex, 4) if ex is not None else None),
        "tags": tags, "He_res": He_res, "pKw": round(pKw, 4),
        "ledger": {k: round(v, 8) for k, v in sorted(ledger.items())
                   if v and k != WATER},
        "src": [(s, t, round(v, 8)) for s, t, v in src],
    }
    if ex is not None:
        out["d_est_chg"] = round(ph - ex, 4)
    return out


def _show(o: dict) -> None:
    print(f"--- {o['label']} ---")
    print(f"    tags={o['tags']}   pH_est={o['pH_est']}  "
          f"pH_chg={o['pH_chg']}  Δ={o.get('d_est_chg')}")
    led = ", ".join(f"{k}={v:g}" for k, v in o["ledger"].items())
    print(f"    ledger: {led}")
    top = sorted(o["src"], key=lambda t: -t[2])[:6]
    print("    src(top6): " + "  ".join(f"{s}[{t}]={v:g}" for s, t, v in top))


def demo(T) -> int:
    """正/负对照：判据必须能分开"纯两性盐"与"共轭缓冲对"。"""
    cases = [
        # (标签, ledger, He, V, T_K, 独立判据说明)
        ("TC1 H₂S 半中和 273.15K（H₂S/HS⁻ 1:1 缓冲）",
         {"H_2S": 0.005, "HS^-": 0.005, "Na^+": 0.005}, 0.0, 1.0, 273.15, ""),
        ("TC3 H₂S 半中和 298.15K",
         {"H_2S": 0.005, "HS^-": 0.005, "Na^+": 0.005}, 0.0, 1.0, 298.15, ""),
        ("纯 NaHS 0.01 M（真两性盐，中点式应适用）",
         {"HS^-": 0.01, "Na^+": 0.01}, 0.0, 1.0, 298.15, ""),
        ("AB03 H₂S+少量NaOH",
         {"H_2S": 0.009, "HS^-": 0.001, "Na^+": 0.001}, 0.0, 1.0, 298.15, ""),
        ("纯 NaHCO₃ 0.01 M（教科书两性，pH≈8.3）",
         {"HCO_3^-": 0.01, "Na^+": 0.01}, 0.0, 1.0, 298.15, ""),
        ("CO₂/NaHCO₃ 1:1 缓冲（共用 CO₂/HCO₃⁻）",
         {"CO_2": 0.005, "HCO_3^-": 0.005, "Na^+": 0.005},
         0.0, 1.0, 298.15, ""),
        ("纯 NaH₂PO₄ 0.01 M（两性，pKa 2.1/7.2）",
         {"H_2PO_4^-": 0.01, "Na^+": 0.01}, 0.0, 1.0, 298.15, ""),
    ]
    for label, led, He, V, T_K, _ in cases:
        _show(probe(T, led, He, V, T_K, label))
        print()
    return 0


def suite(T) -> int:
    """全库普查：对每例末态 ledger 问"走步 pH 机器 vs 精确解 vs 呈现 pH"。"""
    with io.open(os.path.join(ROOT, "chemkit", "data", "tests.json"),
                 encoding="utf-8") as fh:
        cases = json.load(fh)
    from chemkit.engine import judge
    rows = []
    for c in cases:
        subs = [{"name": s[0], "mol": float(s[1])} for s in c["subs"]]
        cond = c.get("cond") or {}
        pr = {}
        try:
            r = judge(subs, cond, T, _probe=pr)
        except Exception as exc:
            rows.append({"name": c["name"], "err": repr(exc)})
            continue
        led = {}
        for e in (r.get("final") or []):
            led[e["name"]] = e["mol"]
        V = float(cond.get("V_L", 1.0))
        T_K = float(cond.get("T_K", 298.15))
        He = float(pr.get("H_excess", 0.0) or 0.0)
        o = probe(T, led, He, V, T_K, c["name"])
        o["suite_pH"] = r.get("final_pH")
        o["resid"] = abs(pr.get("resid_live") or 0.0)
        rows.append(o)
    with io.open(os.path.join(ROOT, "logs", "amph_audit.json"), "w",
                 encoding="utf-8") as fh:
        json.dump(rows, fh, ensure_ascii=False, indent=1)
    amphi = [r for r in rows if "tags" in r and "两性中点" in r["tags"]]
    print(f"全库 {len(rows)} 例；分支 4 走**两性中点**的末态 {len(amphi)} 例")
    print()
    print(f"{'name':<34} {'pH_est':>7} {'pH_chg':>7} {'Δ':>7} {'resid':>7} "
          f"{'suite':>6}")
    for r in sorted(amphi, key=lambda r: -abs(r.get("d_est_chg") or 0.0)):
        print(f"{str(r['label'])[:34]:<34} {r.get('pH_est'):>7} "
              f"{str(r.get('pH_chg')):>7} {str(r.get('d_est_chg')):>7} "
              f"{r.get('resid', 0):>7.3f} {str(r.get('suite_pH')):>6}")
    print("\n-> logs/amph_audit.json")
    return 0


def main(argv: list[str]) -> int:
    T = load_tables()
    if "--suite" in argv:
        return suite(T)
    return demo(T)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
