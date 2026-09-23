# -*- coding: utf-8 -*-
"""第 264 轮 · 拾取修法的影响面与回归例逐条核对（只读档）。"""
import io
import json
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

A = {r["name"]: r for r in json.load(io.open("logs/suite-r263-1301.json",
                                             encoding="utf-8"))}
B = {r["name"]: r for r in json.load(io.open("logs/suite-latest.json",
                                             encoding="utf-8"))}
drift = [k for k in A
         if abs((A[k].get("resid_live") or 0) - (B[k].get("resid_live") or 0))
         > 1e-6 or A[k].get("ok") != B[k].get("ok")]
print(f"受影响例数（残差或通过性变化）：{len(drift)} / {len(A)}")
flip_up = [k for k in drift if not A[k]["ok"] and B[k]["ok"]]
flip_dn = [k for k in drift if A[k]["ok"] and not B[k]["ok"]]
print(f"  通过性 FAIL->PASS {len(flip_up)}；PASS->FAIL {len(flip_dn)}")
print(f"  仅残差变化 {len(drift) - len(flip_up) - len(flip_dn)}")

print("\n-- 残差上升最多的 12 例 --")
rows = sorted(drift,
              key=lambda k: -((B[k].get("resid_live") or 0)
                              - (A[k].get("resid_live") or 0)))
for k in rows[:12]:
    d = (B[k].get("resid_live") or 0) - (A[k].get("resid_live") or 0)
    if d <= 0:
        break
    print(f"  {k[:38]:<38} ok {A[k]['ok']!s:<5}->{B[k]['ok']!s:<5} "
          f"resid {A[k].get('resid_live'):>7} -> {B[k].get('resid_live'):>7} "
          f"({d:+.4f})  pH {A[k].get('pH')} -> {B[k].get('pH')}  "
          f"kind {A[k].get('resid_src_kind')}->{B[k].get('resid_src_kind')}")

print("\n-- 残差下降最多的 12 例 --")
for k in rows[::-1][:12]:
    d = (B[k].get("resid_live") or 0) - (A[k].get("resid_live") or 0)
    if d >= 0:
        break
    print(f"  {k[:38]:<38} ok {A[k]['ok']!s:<5}->{B[k]['ok']!s:<5} "
          f"resid {A[k].get('resid_live'):>7} -> {B[k].get('resid_live'):>7} "
          f"({d:+.4f})  pH {A[k].get('pH')} -> {B[k].get('pH')}")
