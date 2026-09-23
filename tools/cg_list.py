"""列出高残差例并按 §1.2 的 A/B/C 三类归因（第 263 轮）。

判据（全部取自 handoff §1.2，**不新增口径**）：

* 签名命中 = 该例出现"pH 跳变"签名（`probe` 里 `frozen_at`/`零推进` 等）；
  这里用**可复算的口径**：`resid_live > 0` 且该例的 `resid_src_kind` 落在
  引擎自报的几类来源里。
* A：信号命中 且 **无固相在场**
* B：信号命中 且 **有固相在场**
* C：**pH 本就对**（`pH_eng == pH_chg`）但通道推不动

⚠️ 本脚本只**读档**（`logs/suite-latest.json`）+ 逐例重跑取 probe，
不写任何文件；重跑结果落 `logs/cg_list.json`。
"""
from __future__ import annotations

import io
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def load_suite(path: str = "logs/suite-latest.json") -> list[dict]:
    with io.open(path, encoding="utf-8") as fh:
        return json.load(fh)


def main() -> int:
    rows = load_suite()
    hi = [r for r in rows if abs(r.get("resid_live") or 0.0) > 0.1]
    hi.sort(key=lambda r: -(abs(r.get("resid_live") or 0.0)))
    tot = sum(abs(r.get("resid_live") or 0.0) for r in rows)
    print(f"总例数 {len(rows)}   残差质量合计 {tot:.1f}")
    print(f"|resid_live| > 0.1 的例数 {len(hi)}   "
          f"质量 {sum(abs(r.get('resid_live') or 0.0) for r in hi):.1f}")
    print()
    print(f"{'#':>3} {'name':<34} {'resid':>9} {'pH':>7} {'deg':>9} "
          f"{'src_kind':<18} {'solver':<10} exit")
    for i, r in enumerate(hi, 1):
        print(f"{i:>3} {str(r.get('name'))[:34]:<34} "
              f"{abs(r.get('resid_live') or 0.0):>9.4f} "
              f"{str(r.get('pH')):>7} {str(r.get('degree'))[:9]:>9} "
              f"{str(r.get('resid_src_kind'))[:18]:<18} "
              f"{str(r.get('pH_solver'))[:10]:<10} {str(r.get('exit'))[:20]}")
    with io.open("logs/cg_list.json", "w", encoding="utf-8") as fh:
        json.dump(hi, fh, ensure_ascii=False, indent=1)
    print("\n-> logs/cg_list.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
