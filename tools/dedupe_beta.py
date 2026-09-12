"""合并 beta.json 中 (center, ligand, nu, complex, logb) 完全相同的重复条目。

只合并**键与数值都相同**的条目（差异仅在 note 文字），保留先出现者的
note 并追加后者的独有说明；数值不同的真冲突一律不合并、原样报出。

用法：python tools/dedupe_beta.py [--write]
"""
from __future__ import annotations

import json
import sys
from collections import Counter

PATH = "chemkit/data/beta.json"


def main(write: bool) -> None:
    rows = json.load(open(PATH, encoding="utf-8"))
    key = lambda e: (e["center"], e["ligand"], e["nu"], e["complex"], e["logb"])  # noqa: E731
    cnt = Counter(key(e) for e in rows)
    dupes = {k for k, v in cnt.items() if v > 1}
    conflicts = {k for k in dupes
                 if len({(e["complex"], e["logb"]) for e in rows
                         if (e["center"], e["ligand"], e["nu"]) == k[:3]}) > 1}
    print(f"总条目 {len(rows)}；完全重复键 {len(dupes)}；"
          f"其中数值真冲突 {len(conflicts)}")
    for k in sorted(conflicts):
        print("   真冲突（不合并）:", k)

    out, seen = [], {}
    for e in rows:
        k = key(e)
        if k in dupes and k not in conflicts:
            if k in seen:
                prev = out[seen[k]]
                n2 = e.get("note", "")
                if n2 and n2 not in (prev.get("note") or ""):
                    prev["note"] = ((prev.get("note") or "") + "；" + n2).strip("；")
                continue
            seen[k] = len(out)
        out.append(e)
    print(f"合并后 {len(out)} 条（删去 {len(rows) - len(out)} 条纯冗余）")
    if write:
        with open(PATH, "w", encoding="utf-8") as f:
            json.dump(out, f, ensure_ascii=False, indent=1)
            f.write("\n")
        print(f"已写回 {PATH}")


if __name__ == "__main__":
    main("--write" in sys.argv)
