# -*- coding: utf-8 -*-
"""第 238 轮 · 直接调 `joint_solve`：B 组停滞态上它到底返回什么？

`tools/joint_visible.py` 已否证"单侧通道不可见"这一线索：
B 组终态上联立**看得见 14 条两侧在场的强驱动方程**（H45 最大 |S|=25.4），
远超 `JOINT_MIN_M=2` 与驱动阈值 0.02。⟹ **它看得到方程，却解不动。**

本脚本在终态账本上**直接调 `joint_solve`**（与 `_joint_fire` 同参数），
打印返回的 `(status, x, resid)`，以确定失败形态：
  · `fail`   —— Newton 未收敛 / 残差大
  · `boundary` —— 不动点在物理域外（数据张力型）
  · `ok`     —— 竟然能解（那说明接线问题仍在别处）

用法： python tools/joint_call.py [用例前缀...]
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

import chemkit.engine as eng                                # noqa: E402
import chemkit.joint as joint                               # noqa: E402
from chemkit.data import load_tables                        # noqa: E402
from chemkit.templates import enumerate_candidates          # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
PREFIXES = sys.argv[1:] or ["H45", "F31", "H43", "E52", "E41"]
cases = {c["name"]: c for c in load_cases(None)}

for pre in PREFIXES:
    hit = [v for n, v in cases.items() if n.startswith(pre)]
    if not hit:
        continue
    c = hit[0]
    probe = {}
    eng.judge([{"name": a, "mol": b} for a, b in c["subs"]],
              c.get("cond") or {"V_L": V}, T, _probe=probe)
    led = dict(probe.get("ledger") or {})
    pH = probe.get("pH_solver")
    He = probe.get("H_excess") or 0.0
    cands = enumerate_candidates(led, He, pH or 7.0, V, T_K, T, True)
    # 复刻 _joint_collect(None) 的过滤（两侧在场 + 按 |S| 截前 JOINT_MAX_M）
    seen = {}
    for cd in cands:
        ps = cd.pres_specs
        if not all(led.get(s, 0.0) > eng.X_MIN for s in ps[0] + ps[1]):
            continue
        if cd.meta.get("slow") or cd.meta.get("deferred"):
            continue
        nk = (cd.netkey_fwd if cd.netkey_fwd <= cd.netkey_rev
              else cd.netkey_rev)
        try:
            S = eng.S_of(cd, led, V, pH or 7.0, T_K, T, frozenset(),
                         eng.P_EXT_KPA, True, {})
        except Exception:                                   # noqa: BLE001
            continue
        if nk not in seen or abs(S) > abs(seen[nk][1]):
            seen[nk] = (cd, S)
    act = sorted(((cd, 1 if S > 0 else -1, S) for cd, S in seen.values()),
                 key=lambda e: -abs(e[2]))[:joint.JOINT_MAX_M]
    print(f"\n{'=' * 74}")
    print(f"{c['name']}   pH={pH}  m(actives)={len(act)}")
    if len(act) < joint.JOINT_MIN_M:
        print("  actives 不足 JOINT_MIN_M，联立不会启动")
        continue
    print("  actives:")
    for cd, d, S in act:
        print(f"    d={d:+d} |S|={abs(S):7.3f} {cd.kind:9s} {cd.r} -> {cd.pr}")
    try:
        st, x, res = joint.joint_solve(
            led, He, [(cd, d) for cd, d, _S in act], V, T_K, T,
            frozenset(), eng.P_EXT_KPA, True, eng.S_of,
            eng.H_ION, eng.WATER)
        print(f"  ⟹ joint_solve: status=**{st}**  resid={res:.4f}")
        print(f"     x = {[round(v, 5) for v in x]}")
    except Exception as exc:                                # noqa: BLE001
        print(f"  ⟹ joint_solve 异常 {type(exc).__name__}: {exc}")

print("\n=== 判读 ===")
print("  · status=fail/boundary ⟹ **联立看得到方程但解不动**：")
print("    这是第四次否证的最后一环，应停止在该轴投入。")
print("  · status=ok 且 resid 小 ⟹ 接线问题仍在别处（如触发时机/落实路径）。")
