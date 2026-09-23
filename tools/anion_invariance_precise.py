# -*- coding: utf-8 -*-
"""第 246 轮 · 重测阴离子不变性（**用未取整的 raw pH**）。

**为什么必须重测**（第 246 轮自查发现）：
第 245 轮的矩阵拿 `r["final_pH"]` 做差，而 `final_pH = round(..., 2)`
（engine.py L2642）——**0.01 的格子会把真差 ≤0.01 全抹成 0.0000**。
所以当时"30 个配对全 Δ=+0.0000（4 位小数）"是**读数精度不足**造成的，
不是测出来的结论。教训见 lessons「分层取数：别用被取整的下游量做差分」。

本脚本改用 `judge(..., _probe={})` 的 `probe["pH"]`（`round(...,3)`）
做差分，分辨率 0.001；并同时打印 raw 的 `estimate_pH` 复算值供核对。

用法： python tools/anion_invariance_precise.py [--write]
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


def complexes_of(center):
    out = {}
    for e in BETA:
        if e.get("center") == center:
            out.setdefault(e["ligand"], []).append((e["nu"], e["logb"]))
    return out


# (center, 基准盐 NO3-, [(阴离子键, 对照盐式, 投料 mol)])；统一 0.01 M 金属
PAIRS = {
    "Mg^{2+}": ("Mg(NO3)2", [("Br^-", "MgBr2", 0.01),
                             ("I^-", "MgI2", 0.01),
                             ("ClO_4^-", "Mg(ClO4)2", 0.01),
                             ("SO_4^{2-}", "MgSO4", 0.01)]),
    "Ca^{2+}": ("Ca(NO3)2", [("Br^-", "CaBr2", 0.01),
                             ("I^-", "CaI2", 0.01),
                             ("ClO_4^-", "Ca(ClO4)2", 0.01)]),
    "Ni^{2+}": ("Ni(NO3)2", [("Br^-", "NiBr2", 0.01),
                             ("I^-", "NiI2", 0.01),
                             ("ClO_4^-", "Ni(ClO4)2", 0.01),
                             ("SO_4^{2-}", "NiSO4", 0.01)]),
    "Mn^{2+}": ("Mn(NO3)2", [("Br^-", "MnBr2", 0.01),
                             ("I^-", "MnI2", 0.01),
                             ("ClO_4^-", "Mn(ClO4)2", 0.01),
                             ("SO_4^{2-}", "MnSO4", 0.01)]),
    "Zn^{2+}": ("Zn(NO3)2", [("ClO_4^-", "Zn(ClO4)2", 0.01),
                             ("SO_4^{2-}", "ZnSO4", 0.01)]),
    "Cd^{2+}": ("Cd(NO3)2", [("ClO_4^-", "Cd(ClO4)2", 0.01),
                             ("SO_4^{2-}", "CdSO4", 0.01)]),
    "Cu^{2+}": ("Cu(NO3)2", [("SO_4^{2-}", "CuSO4", 0.01)]),
    "La^{3+}": ("La(NO3)3", [("Br^-", "LaBr3", 0.01),
                             ("I^-", "LaI3", 0.01)]),
    "Ce^{3+}": ("Ce(NO3)3", [("Br^-", "CeBr3", 0.01),
                             ("I^-", "CeI3", 0.01)]),
    "Nd^{3+}": ("Nd(NO3)3", [("Br^-", "NdBr3", 0.01),
                             ("I^-", "NdI3", 0.01)]),
    "Gd^{3+}": ("Gd(NO3)3", [("Br^-", "GdBr3", 0.01),
                             ("I^-", "GdI3", 0.01)]),
    "Lu^{3+}": ("Lu(NO3)3", [("Br^-", "LuBr3", 0.01),
                             ("I^-", "LuI3", 0.01)]),
    "Y^{3+}": ("Y(NO3)3", [("Br^-", "YBr3", 0.01),
                           ("I^-", "YI3", 0.01)]),
    "Al^{3+}": ("Al(NO3)3", [("Br^-", "AlBr3", 0.01),
                             ("I^-", "AlI3", 0.01)]),
}


def main():
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()

    def ph3(salt, mol):
        """用 probe['pH']（round 3）取数；返回 (pH3, raw, err)。"""
        try:
            pr = {}
            r = judge([{"name": salt, "mol": mol}], {"V_L": 1.0}, T,
                      _probe=pr)
            fin = r.get("final") or []
            if not fin:
                return None, None, "无 final"
            return pr.get("pH"), pr.get("pH"), None
        except Exception as exc:                            # noqa: BLE001
            return None, None, f"{type(exc).__name__}"

    print("=== 阴离子不变性（**probe['pH']，3 位小数**；基准 NO₃⁻）===")
    print("  %-9s %-13s %-6s %-10s %-10s %-9s %s"
          % ("center", "对照盐", "M–Y", "基准 pH", "对照 pH", "ΔpH", "判读"))
    rows, bad = [], []
    for center, (ref, pairs) in sorted(PAIRS.items()):
        cx = complexes_of(center)
        b, _raw, err = ph3(ref, 0.01)
        if b is None:
            print(f"  !! {center} 基准 {ref} 失败: {err}")
            continue
        for key, salt, mol in pairs:
            n = len(cx.get(key, []))
            o, _r2, err = ph3(salt, mol)
            if o is None:
                print(f"  {center:<9} {salt:<13} {n:<6} —  {err}")
                continue
            d = o - b
            if n == 0:
                v = "一致" if abs(d) <= 0.001 else "**⚠ 应不变却变了**"
                if abs(d) > 0.001:
                    bad.append((center, salt, d))
            else:
                v = "可变（库内有该配体条目）"
            rows.append((center, salt, n, b, o, d))
            print("  %-9s %-13s %-6d %-10.3f %-10.3f %+.3f    %s"
                  % (center, salt, n, b, o, d, v))
    print(f"\n库内无 M–Y 配体、但 ΔpH > 0.001 的：{len(bad)} 条")
    for c, s, d in bad:
        print(f"   {c} {s}  Δ={d:+.3f}")
    n_inv = sum(1 for r in rows if r[2] == 0)
    n_bad = len(bad)
    print(f"\n汇总：库内无 M–Y 的配对 {n_inv} 个，其中违反不变性 {n_bad} 个")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "anion_invariance_precise.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"center": r[0], "salt": r[1], "n_complex": r[2],
                        "ph_ref": r[3], "ph_other": r[4], "dph": r[5]}
                       for r in rows], f, ensure_ascii=False, indent=1)
        print(f"[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
