"""`estimate_pH`（走步内 pH 机器）与呈现 pH 的分歧普查（第 263 轮）。

**只读档**。判据：`|pH_solver − pH| > 0.1` ⟹ 走步用的 pH 与呈现给用户的 pH
不是同一个数 —— 走步的 S 评估、通道推进全建立在 `pH_solver` 上，
故分歧大的例必然是残差来源候选。

用法：
    python tools/phdiv.py                 # 读 logs/suite-latest.json
    python tools/phdiv.py logs/suite-r262-1293.json
"""
from __future__ import annotations

import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def main(argv: list[str]) -> int:
    path = argv[1] if len(argv) > 1 else "logs/suite-latest.json"
    with io.open(path, encoding="utf-8") as fh:
        import json
        rows = json.load(fh)
    bad = []
    miss = 0
    for r in rows:
        ph, sol = r.get("pH"), r.get("pH_solver")
        if ph is None or sol is None:
            miss += 1
            continue
        try:
            d = abs(float(sol) - float(ph))
        except (TypeError, ValueError):
            miss += 1
            continue
        if d > 0.1:
            bad.append((d, r))
    bad.sort(key=lambda t: -t[0])
    print(f"档 {path}   例数 {len(rows)}   缺字段 {miss}   "
          f"|pH_solver − pH| > 0.1 的例数 {len(bad)}")
    print()
    print(f"{'#':>3} {'name':<36} {'ΔpH':>7} {'pH':>7} {'solver':>8} "
          f"{'resid':>8} {'src_kind':<12} exit")
    for i, (d, r) in enumerate(bad, 1):
        print(f"{i:>3} {str(r.get('name'))[:36]:<36} {d:>7.3f} "
              f"{str(r.get('pH')):>7} {str(r.get('pH_solver')):>8} "
              f"{abs(r.get('resid_live') or 0.0):>8.3f} "
              f"{str(r.get('resid_src_kind'))[:12]:<12} {str(r.get('exit'))[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
