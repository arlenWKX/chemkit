# -*- coding: utf-8 -*-
"""第 203 轮 · F31 走步层诊断（D15 形核帧缺陷）：回答三问。

用法： python tools/dev.py run tools/f31_diag.py [用例前缀，默认 F31]

三问：
  ① 沉淀候选在各轮 pick 中是否存在、S 多少、是否被禁用/限幅；
  ② derived 步的 `pin_mode` 冻结值与 f(x) 的 pH 轨迹（ROOT_AUDIT 细扫）；
  ③ 联合步机制是否触发过（`_joint_fire` / `_joint_collect`）。
"""
import io
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

import chemkit.engine as eng                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases, run_case           # noqa: E402

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
T = load_tables()
cases = {c["name"]: c for c in load_cases(None)}
c = [v for n, v in cases.items() if n.startswith(PRE + " ")][0]
print(f"用例: {c['name']}  subs={c['subs']}")

# ---- 走步轨迹：抓事件 ----
EV = []


def _wrap(name, extra=None):
    fn = getattr(eng, name, None)
    if fn is None:
        print(f"  [缺] engine.{name} 不存在")
        return
    def w(*a, **k):
        r = fn(*a, **k)
        EV.append((name, a, k, r))
        return r
    setattr(eng, name, w)
    print(f"  [勾住] engine.{name}")


_orig = {n: getattr(eng, n, None) for n in
         ("_joint_fire", "_joint_collect")}
for n in ("_joint_fire", "_joint_collect"):
    _wrap(n)

pr = {}
try:
    ok = run_case(c, T, verbose=False)
except Exception as exc:                                    # noqa: BLE001
    print("  运行异常:", type(exc).__name__, exc)
    ok = False
print(f"断言: {'PASS' if ok else 'FAIL'}")

print(f"\n=== ③ 联合步事件（{len(EV)} 次）===")
for name, a, k, r in EV[:10]:
    print(f"  {name} -> {str(r)[:70]}")

# ---- 履历：用 dev 工具的能力补足 ----
print("\n=== ① 沉淀候选在场性 / ② f(x) 细扫 ===")
print("  用 ROOT_AUDIT 抓 f(x) 轨迹（61 点）")
eng.ROOT_AUDIT = {"case": c["name"], "n": 0, "every": 1, "rec": []}
pr2 = {}
try:
    run_case(c, T, verbose=False)
except Exception as exc:                                    # noqa: BLE001
    print("  二次运行异常:", type(exc).__name__, exc)
rec = eng.ROOT_AUDIT.get("rec") or []
print(f"  记录到 {len(rec)} 个括号")
for r in rec[:3]:
    tr = r.get("trace") or []
    if tr:
        print(f"  --- 括号 x_max={r.get('x_max')} 方案={r.get('dir')} "
              f"trace {len(tr)} 点 ---")
        for row in tr[:8]:
            print(f"     {row}")
        break
