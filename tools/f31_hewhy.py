# -*- coding: utf-8 -*-
"""第 203 轮 · F31 取值定影：走步过程中 H_excess 与各 pH 估计器的实际取值。

核心矛盾（已实测）：
  · normalize 起始态：ledger={Cl-3, Na+3, Ga3+1}, He=-3.0, estimate_pH=14.477 ✓ 化学正确
  · judge  终态：He=0.0, pH=1.62，ledger={Ga3+0.2534, [Ga(OH)4]-0.7431, [Ga(OH)]2+0.0034}
  · 但同一账本的 Σz·n = -2.976 ⟹ 若 He 真是"净质子账本"，应有 n(H+)=2.976
    ⟹ pH≈0.5，而不是 1.62。**(He=0.0, pH=1.62) 二者不可能同时成立。**

本脚本勾住 4 个模块级函数，记录每次调用的入参/返回值，定位：
  · He 是**何时、被哪个反应**从 -3.0 变成 0.0 的；
  · pH 估计器的输入 He 与输出 pH 是否自洽（0 净酸不该给 pH 1.62）；
  · `_presentation_He`（呈现口径"幻影碱归零"闸）是否把 He 改写。

用法： python tools/f31_hewhy.py [用例前缀，默认 F31]
"""
import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith(PRE + " ")][0]

LOG = []


def _hook(name):
    fn = getattr(eng, name, None)
    if fn is None:
        print(f"  [缺] engine.{name}")
        return
    def w(*a, **k):
        try:
            res = fn(*a, **k)
        except Exception as exc:                            # noqa: BLE001
            LOG.append({"fn": name, "err": f"{type(exc).__name__}: {exc}"})
            raise
        rec = {"fn": name}
        if name in ("closed_pH", "estimate_pH", "presentation_pH",
                    "_presentation_He"):
            rec["He_in"] = (a[1] if len(a) > 1 else k.get("H_excess"))
            led = a[0] if a else k.get("ledger")
            if isinstance(led, dict):
                rec["led_sig"] = {s: round(m, 6) for s, m in led.items()
                                  if m > 1e-9 and s != "H_2O"}
        rec["out"] = (list(res) if isinstance(res, tuple)
                      else (round(res, 6) if isinstance(res, float) else res))
        LOG.append(rec)
        return res
    setattr(eng, name, w)


for nm in ("closed_pH", "estimate_pH", "presentation_pH",
           "_presentation_He"):
    _hook(nm)

probe = {}
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
r = eng.judge(subs, c.get("cond") or {"V_L": 1.0}, T, _probe=probe)

print(f"用例 {c['name']}  调用记录 {len(LOG)} 条")
print("\n=== 各函数：输入 He -> 输出 ===")
seen = {}
for rec in LOG:
    fn = rec["fn"]
    seen.setdefault(fn, []).append(rec)
for fn, rs in seen.items():
    print(f"\n--- {fn}  ({len(rs)} 次) ---")
    for rec in rs[:6]:
        if "err" in rec:
            print(f"   ERR {rec['err']}")
            continue
        he = rec.get("He_in")
        print(f"   He_in={he if not isinstance(he, float) else round(he, 6)}"
              f" -> {rec['out']}")
        if fn == "_presentation_He" and rec.get("led_sig"):
            print(f"      led={rec['led_sig']}")
    if len(rs) > 6:
        last = rs[-1]
        print(f"   ... 末次 He_in={last.get('He_in')} -> {last['out']}")

print("\n=== 关键：He 的轨迹（每步执行后的 He）===")
print("  走步 steps（judge 结果）:")
for i, st in enumerate((r.get("steps") or []), 1):
    eq = st.get("equation") if isinstance(st, dict) else str(st)
    ex = st.get("extent") if isinstance(st, dict) else None
    print(f"   {i:2d}. ext={ex} {str(eq)[:80]}")
print(f"\n  终态 probe He = {probe.get('H_excess')}  pH = {probe.get('pH')}")

print("\n=== 手工核算：由终态账本反推 He（He ≡ Σz·n，排除 H+/OH-）===")
import chemkit.speciation as spec                           # noqa: E402
led = probe.get("ledger") or {}
s = 0.0
for sp, m in led.items():
    z = spec._charge_cached(sp)
    if z and sp not in ("H^+", "OH^-"):
        s += z * m
        print(f"   {sp:22s} {m:11.5g} z={z:+d} -> {z * m:+.6g}")
print(f"   Σz·n(除 H+/OH-) = {s:+.6f} mol")
print(f"   ⟹ 若 He 定义成立，n(H+) = {-s:.6f} mol ⟹ pH = {-s and -__import__('math').log10(abs(-s)):.3f}")
print(f"   引擎自报 He = {probe.get('H_excess')}，pH = {probe.get('pH')}")

out = os.path.join(ROOT, "logs", f"f31hewhy-{PRE}.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"case": c["name"], "log": LOG, "probe": probe,
               "he_from_ledger": s, "steps": r.get("steps")},
              f, ensure_ascii=False, indent=1, default=str)
print(f"\n已写 {out}")
