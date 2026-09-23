# -*- coding: utf-8 -*-
"""第 243 轮 · 直接导出任意投料的走步画像（不需要它是套件用例）。

回答"引擎为什么没析出 X"这一类问题：把每一步的
候选集合 / 每个候选的 S / 退出原因 / 迭代数 全打出来。

用法：
  python tools/walkprobe.py LaCl3 0.01
  python tools/walkprobe.py LaCl3 0.01 --json logs/walk_la.json
"""
import io
import json
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


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) < 2:
        print("用法: python tools/walkprobe.py <投料名> <mol> [V_L]")
        return 2
    name, mol = args[0], float(args[1])
    V = float(args[2]) if len(args) > 2 else 1.0
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()
    pr = {}
    r = judge([{"name": name, "mol": mol}], {"V_L": V}, T, _probe=pr)

    print(f"=== {name} {mol} mol / {V} L ===")
    print(f"pH={r.get('final_pH')}  degree={r.get('degree')}  "
          f"changed={r.get('changed')}")
    print(f"exit={pr.get('exit')}  iters={pr.get('iters')}  "
          f"max|S|={pr.get('max_abs_S')}  pH_solver={pr.get('pH_solver')}")
    print("\n-- steps --")
    for s in r.get("steps") or []:
        print(f"   {s['equation']}   x={s['extent']:.6g}")
    print("\n-- production --")
    for p in r.get("production") or []:
        print(f"   {p['name']:24s} {p['mol']:.6g}")
    print("\n-- final (top 12) --")
    fin = sorted((e for e in (r.get("final") or [])),
                 key=lambda e: -e["mol"])[:12]
    for e in fin:
        print(f"   {e['name']:24s} {e['mol']:.6g}")

    act = pr.get("active") or []
    print(f"\n-- active 候选（{len(act)}）--")
    for a in act:
        if isinstance(a, dict):
            print("   " + json.dumps(a, ensure_ascii=False)[:200])
        else:
            print("   " + str(a)[:200])

    out = "--json" in sys.argv
    if out:
        i = sys.argv.index("--json")
        p = sys.argv[i + 1] if i + 1 < len(sys.argv) else "logs/walk.json"
        with io.open(os.path.join(ROOT, p), "w", encoding="utf-8") as f:
            json.dump({"probe": pr, "result": {
                "pH": r.get("final_pH"), "degree": r.get("degree"),
                "steps": r.get("steps"), "production": r.get("production"),
                "final": r.get("final")}}, f, ensure_ascii=False, indent=1,
                default=str)
        print(f"\n[写入] {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
