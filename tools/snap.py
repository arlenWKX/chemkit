"""逐例快照（pH / He / 残差 / 判定 / 摘要）→ JSON，供跨版本对照。

用法：python tools/snap.py out.json
"""
from __future__ import annotations

import json
import sys

sys.path.insert(0, ".")

from chemkit.converg import _result_digest                # noqa: E402
from chemkit.data import load_tables                      # noqa: E402
from chemkit.engine import judge                          # noqa: E402
from chemkit.testsuit import load_cases                   # noqa: E402


def main(out: str) -> None:
    T = load_tables()
    recs = {}
    for c in load_cases(None):
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        probe: dict = {}
        r = judge(subs, cond, T, _probe=probe)
        recs[c["name"]] = {
            "pH": r.get("final_pH"), "deg": r["degree"],
            "changed": r["changed"], "digest": _result_digest(r),
            "He": probe.get("H_excess"), "exit": probe.get("exit"),
            "iters": probe.get("iters"), "resid": probe.get("max_abs_S"),
            "steps": len(r.get("steps", [])),
            "net": r.get("net_equation"),
            "led": {k: v for k, v in (probe.get("ledger") or {}).items()
                    if not k.startswith("__")},
        }
    with open(out, "w", encoding="utf-8") as f:
        json.dump(recs, f, ensure_ascii=False)
    print(f"{len(recs)} 例 → {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "_snap.json")
