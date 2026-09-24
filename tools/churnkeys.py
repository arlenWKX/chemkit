# -*- coding: utf-8 -*-
"""EU01 空转对的**事实清单**（第 289 轮修正用）：

  ① 步数按净键计数（谁在空转、各走了多少量）；
  ② 终态活跃候选里 resid_live 那条的全部判据字段
     （two_sided / slow / blocked / ext_max / frozen / S）——
     回答"`_freeze` 到底会不会拒冻它、拒冻的话卡在哪一条"。

用法：python tools/churnkeys.py EU01
"""
from __future__ import annotations

import sys
from collections import Counter

sys.path.insert(0, ".")

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from chemkit.data import load_tables          # noqa: E402
from chemkit.engine import judge              # noqa: E402
from chemkit.testsuit import load_cases       # noqa: E402


def main(keys: list[str]) -> int:
    T = load_tables()
    cases = load_cases(None)
    picks = [c for c in cases if c["name"].startswith(tuple(keys))]
    for c in picks:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        pr: dict = {}
        r = judge(subs, c.get("cond") or {}, T, _probe=pr)
        steps = r.get("steps") or []
        cnt = Counter(s["equation"] for s in steps)
        ext = Counter()
        for s in steps:
            ext[s["equation"]] += s["extent"]
        print("=" * 76)
        print(f"{c['name']}  steps={len(steps)}")
        rs = pr.get("resid_src")
        if rs:
            print(f"  resid_src: S={rs.get('S'):+.3f} two={rs.get('two_sided')} "
                  f"slow={rs.get('slow')} blocked={rs.get('blocked')} "
                  f"frozen={rs.get('frozen')} ext_max={rs.get('ext_max'):.3g} "
                  f"dis={rs.get('disabled')} why={rs.get('dis_why')} "
                  f"{str(rs.get('eq'))[:60]}")
        print("  -- 步数 top6（净量 / 毛量）--")
        for eq, n in cnt.most_common(6):
            print(f"  {n:>6}  net={0.0:+.3g}  gross={ext[eq]:.4g}  {eq[:60]}")
        print("  -- 终态活跃候选（|S| 降序 top8）--")
        act = pr.get("active") or []
        act = sorted(act, key=lambda a: -abs(a.get("S", 0.0)))
        for a in act[:8]:
            print(f"  S={a.get('S'):+.3f} two={a.get('two_sided')} "
                  f"slow={a.get('slow')} blocked={a.get('blocked')} "
                  f"frozen={a.get('frozen')} ext_max={a.get('ext_max'):.3g} "
                  f"{str(a.get('eq'))[:56]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["EU01"]))
