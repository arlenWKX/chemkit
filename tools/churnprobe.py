# -*- coding: utf-8 -*-
"""**空转/张力族的走步诊断探针**：把一个用例的"走步为什么会空转"所需的**全部**
计数器一次打出来 —— 免得为了看某一个计数器反复跑套件。

## 为什么需要它

第 287–288 轮查 `EU01` 的空转（2815 步）时，我是在 `dev.py case`（步表按 extent 取
top-N）、`CHEM_TRACE`（逐键次数）、`fscan`（求根）、`extent.py`（终态候选）之间来回跑，
每换一个角度就重跑一次 —— 违反 "跑一次读多次"。而且有几个关键计数器
（`joint_tries` / `joint_ok` / 失败签名黑名单 / `windows` 的（净移, 周转）对）
**只在 `_probe` 里**，`dev.py case` 不打。

## 打出来的东西（与引擎内部口径同名，不做二次解释）

* `iters` / `exit` / `degree` / 终态 pH / `resid_live`
* **联立梯子**：`joint_tries` / `joint_ok` ⟹ "第一级（联立）到底有没有被问过"
* **冻结梯子**：`freeze_events` / `frozen_at` / `freeze_revive` / `freeze_expired`
* **微步**：`micro_steps`；**循环**：`cycle_jumps` / `cycle_freeze` / `cycle_dx`
* **周转窗**：`windows` 里的 `(净移, 周转)` 对 —— 判"高周转 + 零净移"的**形状**
  是否成立（若形状成立而检测器没开火，那就是**门槛**问题，不是判据问题）

用法：

    python tools/churnprobe.py EU01 [更多前缀…]
"""
from __future__ import annotations

import sys

sys.path.insert(0, ".")

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.data import load_tables          # noqa: E402
from chemkit.engine import judge              # noqa: E402
from chemkit.testsuit import load_cases       # noqa: E402


def main(keys: list[str]) -> int:
    T = load_tables()
    cases = load_cases(None)
    picks = [c for c in cases if c["name"].startswith(tuple(keys))]
    if not picks:
        print("未匹配到用例")
        return 2
    for c in picks:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        pr: dict = {}
        r = judge(subs, c.get("cond") or {}, T, _probe=pr)
        print("=" * 76)
        print(f"{c['name']}")
        print(f"  pH={r.get('final_pH')}  degree={r.get('degree')}  "
              f"changed={r.get('changed')}  steps={len(r.get('steps') or [])}")
        for k in ("joint_tries", "joint_ok", "freeze_events", "freeze_revive",
                  "freeze_expired", "micro_steps", "cycle_jumps", "cycle_freeze",
                  "cycle_dx", "iters", "active"):
            if k in pr:
                v = pr[k]
                if k == "active":
                    v = f"{len(v)} 条活跃键"
                print(f"  {k:<15} {v}")
        fa = pr.get("frozen_at")
        if fa:
            print(f"  frozen_at       {len(fa)} 个键")
            for kk, vv in list(fa.items())[:6]:
                print(f"      {str(kk)[:56]}  <- {str(vv)[:46]}")
        w = pr.get("windows") or []
        if w:
            print(f"  windows（净移, 周转）共 {len(w)} 窗；末 5 窗：")
            for dr, to in w[-5:]:
                rr = (dr / to) if to else 0.0
                print(f"      net={dr:<12.6g} turnover={to:<12.6g} "
                      f"net/turnover={rr:.4f}  {'← 形状成立（高周转零净移）' if to > 0 and rr <= 0.02 else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["EU01"]))
