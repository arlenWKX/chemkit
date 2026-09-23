# -*- coding: utf-8 -*-
"""第 203 轮 · `estimate_pH` 是否**依赖隐藏状态**（缓存/预热）？

矛盾点：同一 (ledger, He) 组合，诊断工具与新进程给出 pH=12.3803，
而引擎走步过程中在同一账本上得到 1.3109 / 1.6197（相差 10.8 单位）。

若 `estimate_pH` 是**纯函数**，两者必须相等 ⟹ 矛盾只能用"引擎走步时的
账本/He 与探针导出的不同"解释（例如 respeciate 改了 He）。
若 `estimate_pH` **依赖 cache/touch 等隐藏状态**，则走步序会改变 pH 估计
⟹ 这是比 F31 更根本的缺陷（走步结果不可复现）。

本脚本三种调用方式对比同一 (ledger, He)：
  ① 冷调用（cache=None）
  ② 预热 cache（先算几个邻居再算目标）
  ③ 走步实跑后（judge 之后再算同一账本）
并逐 He 打印，看是否存在 cache 依赖。

用法： python tools/phcache.py
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
from chemkit.data import load_tables                        # noqa: E402
from chemkit.testsuit import load_cases                     # noqa: E402

T = load_tables()
V, T_K = 1.0, 298.15
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith("F31 ")][0]
subs = [{"name": n, "mol": m} for n, m in c["subs"]]

probe = {}
r = eng.judge(subs, c.get("cond") or {"V_L": V}, T, _probe=probe)
LED = dict(probe.get("ledger") or {})
LED["H_2O"] = 55.6
print("终态账本 =", {k: round(v, 6) for k, v in LED.items() if k != "H_2O"})
print("引擎 probe pH =", probe.get("pH"), " He =", probe.get("H_excess"))

print("\n=== ① 冷调用 vs ② 预热 cache vs ③ 不同 He ===")
print(f"  {'He':>12} {'冷(cache=None)':>16} {'预热 cache':>14} {'差':>8}")
for he in (-3.0, -1.0, -0.1, -0.05, -0.03, -0.024007, -0.02, -0.01, 0.0,
           0.001, 0.1):
    cold = eng.estimate_pH(LED, he, V, T, T_K)
    cache = {}
    for h2 in (-3.0, -1.0, -0.5, -0.1, -0.02, 0.0, 0.5, 1.0):
        try:
            eng.estimate_pH(LED, h2, V, T, T_K, cache=cache)
        except Exception:                                   # noqa: BLE001
            pass
    warm = eng.estimate_pH(LED, he, V, T, T_K, cache=cache)
    d = "" if abs(warm - cold) < 1e-9 else f"{warm - cold:+.4f} ⚠"
    print(f"  {he:+12.6g} {cold:16.6f} {warm:14.6f} {d:>8}")

print("\n=== ③ 用 judge 内部同一路径复算（He 取 respeciate 后的值）===")
_respeciate_strong_acids = eng._respeciate_strong_acids      # noqa: E402
for he0 in (0.0, -0.024007):
    for tag, hev in (("原样", he0),
                     ("respeciate", _respeciate_strong_acids(
                         dict(LED), he0, V, T))):
        try:
            ph = eng.estimate_pH(LED, hev, V, T, T_K)
            print(f"  He={he0:+.6f} {tag:12s} -> He'={hev:+.6f}  pH={ph:.4f}")
        except Exception as exc:                            # noqa: BLE001
            print(f"  He={he0:+.6f} {tag:12s} 异常 {type(exc).__name__}")

print("\n=== 判读 ===")
print("  冷/预热两列若一致 ⟹ estimate_pH 是纯函数，矛盾只能在'He 被改写'；")
print("  若不一致 ⟹ pH 估计依赖隐藏 cache，走步结果不可复现（更严重）。")
print("  注意 `touch`/`no_base`/`no_acid`/`pin_mode` 参数：它们是**接缝开关**，")
print("  若走步与探针传了不同的组合，同一账本会得到不同 pH（D14 类）。")
