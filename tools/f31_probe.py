# -*- coding: utf-8 -*-
"""第 203 轮 · F31 走步层诊断 v2（D15 形核帧缺陷）：用 `judge(_probe=)` 官方钩子。

第 203 轮教训：`_joint_fire` / `_joint_collect` 是 `engine.run` 内部的**嵌套闭包**
（engine.py L1122/L1155），不是模块级函数 ⟹ 旧版 `setattr(eng, name, w)` 勾不住，
打印"[缺] 不存在"是勾子写错，不是引擎缺机制。正确入口是 `judge(..., _probe={})`
（engine.py L2140 `_probe_exit`），它导出 joint_tries/joint_ok/cycle_jumps/
freeze_events/frozen_at/micro_steps + 逐候选 S。

用法： python tools/dev.py run tools/f31_probe.py [用例前缀，默认 F31]
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
cases = {c["name"]: c for c in load_cases(None)}
hit = [v for n, v in cases.items() if n.startswith(PRE + " ")]
if not hit:
    print("无用例匹配", PRE)
    sys.exit(1)
c = hit[0]
print(f"用例: {c['name']}")
print(f"  subs = {c['subs']}")
print(f"  cond = {c.get('cond')}")
for k in ("changed", "reacted", "degree", "has", "has_not", "has_range",
          "ann", "override"):
    if k in c:
        print(f"  期望 {k} = {c[k]}")

probe = {}
subs = [{"name": n, "mol": m} for n, m in c["subs"]]
r = eng.judge(subs, c.get("cond") or {"V_L": 1.0}, T, _probe=probe)
print(f"\n结果: changed={r['changed']} reacted={r['reacted']} "
      f"degree={r['degree']} pH={r.get('final_pH')}")

print("\n=== ③ 联立/循环机制是否被问过（决定'归中接力'为何不触发）===")
for k in ("exit", "iters", "hist", "steps_n", "cycle_jumps", "joint_tries",
          "joint_ok", "freeze_events", "freeze_expired", "freeze_revive",
          "micro_steps", "cycle_freeze", "cycle_dx"):
    if k in probe:
        print(f"  {k:15s} = {probe[k]}")
print("  -- 判读: cycle_jumps=0 且 joint_tries=0 ⟹ 触发条件(it>=64 且近24步"
      "全微步<0.02)从未成立，机制**没被问过**；>0 才是'问了但失败'。")

print("\n=== ① 最终账本（这是终态化学事实）===")
led = probe.get("ledger") or {}
for s_, m in sorted(led.items(), key=lambda kv: -kv[1]):
    if m > 1e-6:
        print(f"  {s_:22s} {m:.6g}")

print("\n=== ① 终态候选残留（按 |S| 降序，只看两侧在场）===")
act = [a for a in (probe.get("active") or []) if a.get("two_sided")]
act.sort(key=lambda a: -abs(a["S"]))
print(f"  两侧在场候选 {len(act)} 条；max_abs_S={probe.get('max_abs_S')}")
for a in act[:12]:
    print(f"  S={a['S']:+9.3f} ext={a['ext_max']:<11.3g} "
          f"froz={int(a['frozen'])} dis={int(a['dis_fwd'])}/{int(a['dis_rev'])} "
          f"{a['kind']:9s} {a['eq'][:62]}")
    if a.get("dis_why"):
        print(f"        why={str(a['dis_why'])[:88]}")

print("\n=== ① 单侧在场（产物/反应物缺席 ⟹ 不判残差，但要看有没有镓酸根/固）===")
one = [a for a in (probe.get("active") or []) if not a.get("two_sided")]
for a in one:
    if any(t in a["eq"] for t in ("Ga", "H3Ga", "Ga(OH)", "GaO")):
        print(f"  S={a['S']:+9.3f} ext={a['ext_max']:<11.3g} "
              f"dis={int(a['dis_fwd'])}/{int(a['dis_rev'])} {a['eq'][:70]}")

print(f"\n  frozen_n={probe.get('frozen_n')} disabled_n={probe.get('disabled_n')}")
print(f"  frozen_at = {json.dumps(probe.get('frozen_at') or {}, ensure_ascii=False)[:400]}")

out = os.path.join(ROOT, "logs", f"f31probe-{PRE}.json")
with io.open(out, "w", encoding="utf-8") as f:
    json.dump({"case": c["name"], "spec": c, "probe": probe,
               "result": {k: v for k, v in r.items()
                          if k in ("changed", "reacted", "degree",
                                   "final_pH", "annotations", "steps")}},
              f, ensure_ascii=False, indent=1, default=str)
print(f"\n已写 {out}")
