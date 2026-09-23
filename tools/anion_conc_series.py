# -*- coding: utf-8 -*-
"""第 246 轮 · 阴离子不变性的**浓度序列**（0.001 / 0.01 / 0.1 M）。

第 245 轮只测了 0.01 M 一档，且用 `final_pH`（round 2）做差 —— 读数精度不足。
本轮两处升级：
  ① 取数改用 `probe["pH"]`（3 位小数），分辨率 0.001；
  ② 加浓度维度，检验"Δ=0"是不是只在某一档偶然成立。

判据仍然来自数据库自身：`beta.json` 内该金属若无 M–Y 配合物条目，
则该阴离子不参与任何含 H⁺/OH⁻ 的平衡 ⟹ 同阳离子、同计量的两盐
**应当给出同一个 pH**。

用法： python tools/anion_conc_series.py [--write]
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

BETA = json.load(io.open(os.path.join(ROOT, "chemkit/data/beta.json"),
                         encoding="utf-8"))
CS = (0.001, 0.01, 0.1)

# center -> (基准盐 NO3-, {阴离子键: 对照盐式})
SETS = {
    "Mg^{2+}": ("Mg(NO3)2", {"Br^-": "MgBr2", "I^-": "MgI2",
                             "ClO_4^-": "Mg(ClO4)2", "SO_4^{2-}": "MgSO4"}),
    "Ca^{2+}": ("Ca(NO3)2", {"Br^-": "CaBr2", "I^-": "CaI2",
                             "ClO_4^-": "Ca(ClO4)2"}),
    "Ni^{2+}": ("Ni(NO3)2", {"Br^-": "NiBr2", "I^-": "NiI2",
                             "ClO_4^-": "Ni(ClO4)2", "SO_4^{2-}": "NiSO4"}),
    "Mn^{2+}": ("Mn(NO3)2", {"Br^-": "MnBr2", "I^-": "MnI2",
                             "ClO_4^-": "Mn(ClO4)2", "SO_4^{2-}": "MnSO4"}),
    "La^{3+}": ("La(NO3)3", {"Br^-": "LaBr3", "I^-": "LaI3"}),
    "Gd^{3+}": ("Gd(NO3)3", {"Br^-": "GdBr3", "I^-": "GdI3"}),
    "Lu^{3+}": ("Lu(NO3)3", {"Br^-": "LuBr3", "I^-": "LuI3"}),
    "Al^{3+}": ("Al(NO3)3", {"Br^-": "AlBr3", "I^-": "AlI3"}),
}


def main():
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()

    def ph3(salt, c):
        try:
            pr = {}
            r = judge([{"name": salt, "mol": c}], {"V_L": 1.0}, T, _probe=pr)
            if not (r.get("final") or []):
                return None
            return pr.get("pH")
        except Exception:                                   # noqa: BLE001
            return None

    rows = []
    print("=== 阴离子不变性 × 浓度（probe pH，3 位小数）===")
    for center, (ref, pairs) in sorted(SETS.items()):
        n_ent = {}
        for e in BETA:
            if e.get("center") == center:
                n_ent[e["ligand"]] = n_ent.get(e["ligand"], 0) + 1
        for c in CS:
            b = ph3(ref, c)
            if b is None:
                continue
            for key, salt in sorted(pairs.items()):
                o = ph3(salt, c)
                if o is None:
                    continue
                d = o - b
                rows.append((center, salt, c, n_ent.get(key, 0), b, o, d))
    print("  %-9s %-12s %-8s %-5s %-9s %-9s %s"
          % ("center", "对照盐", "C/M", "M–Y", "基准", "对照", "ΔpH"))
    bad = 0
    for center, salt, c, n, b, o, d in rows:
        flag = ""
        if n == 0 and abs(d) > 0.001:
            flag = "  ⚠⚠ 应不变却变了"
            bad += 1
        elif n:
            flag = "  （有配体条目，可变）"
        print("  %-9s %-12s %-8g %-5d %-9.3f %-9.3f %+.3f%s"
              % (center, salt, c, n, b, o, d, flag))
    tot0 = sum(1 for r in rows if r[3] == 0)
    print(f"\n库内无 M–Y 的配对共 {tot0} 个（跨 {len(CS)} 个浓度）；"
          f"违反不变性 {bad} 个")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "anion_conc_series.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"center": r[0], "salt": r[1], "C": r[2],
                        "n_complex": r[3], "ph_ref": r[4], "ph_other": r[5],
                        "dph": r[6]} for r in rows], f,
                      ensure_ascii=False, indent=1)
        print(f"[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
