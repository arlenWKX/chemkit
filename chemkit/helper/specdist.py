# -*- coding: utf-8 -*-
"""**形态分布呈现探针**：把「账本里有什么」与「呈现层报了什么」并列打出来。

## 为什么需要它

第 292 轮查明"高阶氯络合物族"的失败形态是：**账本正确、呈现层把弱池折叠掉** ——
`ZN6 ZnCl2 0.1+NaCl 6` 实测账本净差

    0.24583Cl^- + 0.089491Zn^{2+}  ->  0.0721183[ZnCl_3]^- + 0.0121036[ZnCl_2]
                                       + 0.0052677[ZnCl]^+ + 1.39e-06[Zn(OH)]^+

（`[ZnCl₃]⁻` 占 Zn 的 72%，与库值独立不动点推导一致），
但同一例的**呈现层**给 `degree=0 / changed=False / reacted=False`，
`has` 断言读到的 `[ZnCl_3]^-` 产量 = **0**。

⟹ 判"是引擎错还是标准错"之前，必须先分清**三层**各自说了什么：

1. **账本**（`probe["ledger"]`，`judge(_probe=…)`）：化学真值所在；
2. **净差**（`账本净差`，引擎自己渲染的 initial→final 之差）；
3. **呈现**（`Reaction(r).consumption / production / final`，`testsuit` 的
   `has`/`has_not`/`has_range` 断言**只读这一层**）。

本工具把三层逐物种并列，差异处标 `← 只在账本` / `← 只在呈现`。

用法：

    python chemkit/helper/specdist.py CD3 ZN6 [更多前缀…]
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.candidates import WATER            # noqa: E402
from chemkit.data import load_tables            # noqa: E402
from chemkit.engine import judge                # noqa: E402
from chemkit.system import Reaction             # noqa: E402
from chemkit.testsuit import load_cases         # noqa: E402


def _amt(seq) -> dict:
    out: dict = {}
    for e in (seq or []):
        out[e["name"]] = e.get("mol", 0.0)
    return out


def main(keys: list[str]) -> int:
    T = load_tables()
    cases = load_cases(None)
    picks = [c for c in cases if c["name"].startswith(tuple(keys))]
    if not picks:
        print("未匹配到用例")
        return 2
    for c in picks:
        pr: dict = {}
        r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {}, T, _probe=pr)
        led = {k: v for k, v in (pr.get("ledger") or {}).items()
               if v > 1e-9 and k != WATER and not k.startswith("__")}
        init = {e["name"]: e.get("mol", 0.0) for e in (r.get("initial") or [])}
        rx = Reaction(r)
        prod = dict(rx.production or {}) if isinstance(rx.production, dict) \
            else _amt(rx.production)
        print("=" * 78)
        print(f"{c['name']}")
        print(f"  degree={r.get('degree')} changed={r.get('changed')} "
              f"reacted={r.get('reacted')} pH={r.get('final_pH')}")
        keys_all = sorted(set(led) | set(init) | set(prod),
                          key=lambda k: -max(led.get(k, 0.0), init.get(k, 0.0),
                                             abs(prod.get(k, 0.0))))
        print(f"  {'物种':<22}{'账本(终)':>12}{'初态':>12}{'呈现Δ':>12}   备注")
        for k in keys_all:
            a, b, d = led.get(k, 0.0), init.get(k, 0.0), prod.get(k, 0.0)
            tag = ""
            if abs(a - b) > 1e-6 and abs(d) < 1e-9:
                tag = "← 只在账本（呈现层丢了）"
            elif abs(d) > 1e-9 and abs(a - b) < 1e-9:
                tag = "← 只在呈现"
            print(f"  {k:<22}{a:>12.5f}{b:>12.5f}{d:>12.5f}   {tag}")
        if c.get("has") or c.get("has_range") or c.get("has_not"):
            print(f"  断言: has={c.get('has')} range={c.get('has_range')} "
                  f"not={c.get('has_not')}")
        print("  ── 呈现层原始清单 ──")
        for k in ("initial", "consumption", "production", "final"):
            seq = r.get(k)
            if seq is None:
                print(f"    {k:<12} (缺)")
                continue
            items = [(e["name"], e.get("mol", 0.0)) for e in seq]
            print(f"    {k:<12} n={len(items):<3} "
                  + ", ".join(f"{n}={m:.5g}" for n, m in items[:8]))
        print(f"    keys: {sorted(r.keys())}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["ZN6"]))
