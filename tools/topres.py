"""残差排行（含口径分类明细）。用法：python tools/topres.py [-n 20] [dump.json]"""

from __future__ import annotations

import json
import sys

path = "converg-baseline.json"
n = 20
args = sys.argv[1:]
if "-n" in args:
    n = int(args[args.index("-n") + 1])
    args = args[:args.index("-n")] + args[args.index("-n") + 2:]
if args:
    path = args[0]

doc = json.load(open(path, encoding="utf-8"))
cases = doc["cases"]
print(f"{path}: n={len(cases)}")
row = ("{r:8.3f} {e:9s} it={i:4d} frz={f:7.2f} slow={s:7.2f} "
       "blk={b:7.2f} trc={t:7.2f}  {nm}")
print("\n== live 残差排行（真实收敛债务）==")
for x in sorted(cases, key=lambda c: -c.get("resid_live", 0))[:n]:
    print(row.format(r=x.get("resid_live", 0), e=x.get("exit", "?"),
                     i=x.get("iters", 0), f=x.get("resid_frozen", 0),
                     s=x.get("resid_slow", 0), b=x.get("resid_blocked", 0),
                     t=x.get("resid_trace", 0), nm=x["name"][:44]))
print("\n== 被排除类的最大者（看得见但不计入）==")
for k in ("resid_frozen", "resid_slow", "resid_blocked", "resid_trace"):
    top = max(cases, key=lambda c: c.get(k, 0))
    print(f"  {k:14s} {top.get(k, 0):8.3f}  {top['name'][:50]}")
