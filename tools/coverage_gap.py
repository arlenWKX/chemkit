# -*- coding: utf-8 -*-
"""第 242 轮 · 覆盖缺口普查：现库用例 vs 库内可用数据。

**只读**。回答两个问题：
  ① 哪些元素/体系**已有** hydrolysis 常数（β₁）但**零用例**覆盖？
  ② 哪些体系的用例数远低于其数据丰度（薄覆盖）？

用法：
  python tools/coverage_gap.py            # 全部清单
  python tools/coverage_gap.py beta       # 只看 beta 中心覆盖率
  python tools/coverage_gap.py thin       # 只看薄覆盖族
"""
import io
import json
import os
import re
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def _load(p):
    with io.open(os.path.join(ROOT, p), encoding="utf-8") as f:
        return json.load(f)


def norm(name):
    """把物种名归一成"元素符号"便于统计：Ga^{3+}/[Ga(OH)_4]^- -> Ga"""
    m = re.match(r"^\[?([A-Z][a-z]?)", name or "")
    return m.group(1) if m else None


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    cases = _load("chemkit/data/tests.json")
    beta = _load("chemkit/data/beta.json")
    ksp = _load("chemkit/data/ksp.json") if os.path.exists(
        os.path.join(ROOT, "chemkit/data/ksp.json")) else []

    # 用例覆盖：把 subs（投料）与 name 里的元素都算上
    case_elems = defaultdict(int)
    for c in cases:
        els = set()
        for n, _m in c.get("subs") or []:
            e = norm(n)
            if e:
                els.add(e)
        for e in re.findall(r"\b([A-Z][a-z]?)[a-z0-9]*", c.get("name") or ""):
            if len(e) >= 1 and e[0].isupper():
                els.add(e)
        for e in els:
            case_elems[e] += 1

    # beta 中心（仅 OH⁻ 水解 = "某金属会水解"的证据）
    hyd = defaultdict(list)
    for e in beta:
        if e.get("ligand") == "OH^-":
            hyd[norm(e["center"])].append(e.get("nu"))

    # Ksp 阳离子
    ksp_cat = defaultdict(int)
    for e in ksp:
        pair = e.get("pair") or []
        if pair:
            k = norm(pair[0])
            if k:
                ksp_cat[k] += 1

    if mode in ("all", "beta"):
        print("== ① 有 OH⁻ 水解常数但用例覆盖极少的体系 ==")
        print("   %-6s %-10s %-10s %s" % ("el", "nu(OH-)", "用例数", "Ksp数"))
        rows = []
        for el in sorted(hyd, key=lambda x: case_elems.get(x, 0)):
            rows.append((case_elems.get(el, 0), el, sorted(hyd[el]),
                         ksp_cat.get(el, 0)))
        for n, el, ns, nk in rows:
            flag = "  <<< 零覆盖" if n == 0 else ("  <- 薄" if n <= 2 else "")
            print("   %-6s %-10s %-10d %-6d%s"
                  % (el, ns, n, nk, flag))

    if mode in ("all", "thin"):
        print("\n== ② 元素覆盖总表（取用例数最少的 30 个）==")
        for el, n in sorted(case_elems.items(), key=lambda kv: kv[1])[:30]:
            print("   %-6s %d 例" % (el, n))
        print("\n   用例总数 %d；涉及元素 %d 个"
              % (len(cases), len(case_elems)))

    if mode in ("all", "ksp"):
        print("\n== ③ Ksp 阳离子 vs 用例覆盖 ==")
        for el, nk in sorted(ksp_cat.items()):
            n = case_elems.get(el, 0)
            if n <= 3:
                print("   %-6s Ksp=%-3d 用例=%d" % (el, nk, n))


if __name__ == "__main__":
    main()
