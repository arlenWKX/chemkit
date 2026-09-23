# -*- coding: utf-8 -*-
"""第 223 轮 · Pitzer 数据勘察：结构、覆盖面、与 SIT 的接入点对照。

`chemkit/data/pitzer.json` 已有 426 条参数但引擎零消费（D4 死数据）。
本脚本回答实现前必须知道的四件事：
  ① 条目结构（字段名、取值形态、是否带温度依赖）；
  ② 各 `kind` 的条数与覆盖的离子对；
  ③ **与库内物种命名的匹配率**（Pitzer 用 PHREEQC 命名，账本用 chemkit 命名）；
  ④ 现 SIT 层的接入点（`speciation.sit_logK` / `logK_T`）以便对称接线。

用法： python tools/pitzer_survey.py
"""
import io
import json
import os
import sys
import collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.data import load_tables                        # noqa: E402

T = load_tables()
pj = json.load(io.open(os.path.join(ROOT, "chemkit", "data",
                                    "pitzer.json"), encoding="utf-8"))
print("=== ① 顶层 ===")
for k, v in pj.items():
    if k == "entries":
        print(f"  {k}: {len(v)} 条")
    else:
        print(f"  {k}: {str(v)[:150]}")

ents = pj["entries"]
print("\n=== ② 条目结构 ===")
print(f"  样例: {json.dumps(ents[0], ensure_ascii=False)[:300]}")
fields = collections.Counter()
for e in ents:
    for k in e:
        fields[k] += 1
print(f"  字段出现次数: {dict(fields)}")

print("\n=== ③ kind 分布 + 值形态 ===")
kinds = collections.Counter(e.get("kind") for e in ents)
for k, v in kinds.most_common():
    ex = next(e for e in ents if e.get("kind") == k)
    print(f"  {str(k):10s} {v:4d}  例: {json.dumps(ex, ensure_ascii=False)[:180]}")

print("\n=== ④ 与库内物种命名的匹配率 ===")
# chemkit 端可能出现的物种名
cands = set()
for b in T.beta:
    cands.add(b["center"])
    cands.add(b["complex"])
cands |= set(T.pka_acid) | set(T.pka_base)
cands |= {e["pair"][0] for e in T.ksp} | {e["pair"][1] for e in T.ksp}
for e in T.thermo:
    if isinstance(e, dict) and e.get("species"):
        cands.add(e["species"])


def names_of(e):
    """Pitzer 条目的离子名（字段是 `ions`，不是 `ion`）。"""
    v = e.get("ions")
    if isinstance(v, str):
        return [v]
    if isinstance(v, (list, tuple)):
        return [x for x in v if isinstance(x, str)]
    return []


hit = miss = 0
missed = collections.Counter()
for e in ents:
    for nm in names_of(e):
        if nm in cands:
            hit += 1
        else:
            miss += 1
            missed[nm] += 1
print(f"  名称命中库内物种: {hit}   未命中: {miss}")
print(f"  未命中最多者（前 15）:")
for nm, c in missed.most_common(15):
    print(f"    {nm:26s} ×{c}")

print("\n=== ⑤ 现 SIT 接入点（对称接线参考）===")
from chemkit import speciation as spec                      # noqa: E402
import inspect                                              # noqa: E402
for fn in ("sit_logK", "ionic_strength", "ionic_strength_reaction"):
    f = getattr(spec, fn, None)
    if f is None:
        print(f"  {fn}: 缺")
        continue
    try:
        print(f"  {fn}{inspect.signature(f)}")
    except Exception:                                       # noqa: BLE001
        print(f"  {fn}: ?")
from chemkit import candidates as cand                      # noqa: E402
f = getattr(cand, "logK_T", None)
if f is not None:
    print(f"  candidates.logK_T{inspect.signature(f)}")
print(f"  speciation.SIT_ALL = {getattr(spec, 'SIT_ALL', '?')}")
print(f"  speciation._SIT_A   = {getattr(spec, '_SIT_A', '?')}")

print("\n=== 判读 ===")
print("  · 若未命中率高 ⟹ 需要一个 PHREEQC→chemkit 的名称映射层（参考")
print("    tools/xcheck.py 的既有映射）。")
print("  · 若条目多为 b0/b1 而缺 Cφ ⟹ 只能做简化式，需注明近似的边界。")
