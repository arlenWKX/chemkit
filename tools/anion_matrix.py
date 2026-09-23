# -*- coding: utf-8 -*-
"""第 245 轮 · 阴离子配对矩阵（**显式盐名，不做字符串替换**）。

第一版用 `cl_salt.replace("Cl", suf)` 造盐名，造出 `AlNO33` / `CdClO42`
这种**不存在的式子**，引擎只好给出"无反应 pH≈7"的假结果，
于是矩阵里出现"库内无该配体却变了 +3.4"的大片假警报。
**教训：造化学式不能用字符串替换**（见 lessons「造式别用字符串替换」）。

原理：若库内**没有**该 金属–阴离子 的配合物条目，则该阴离子在模型里
完全不参与 ⟹ `MXₙ` 与 `MYₙ` 的 pH 应当相同；有则该变。

用法： python tools/anion_matrix.py [--write]
"""
import io
import json
import os
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

BETA = json.load(io.open(os.path.join(ROOT, "chemkit/data/beta.json"),
                         encoding="utf-8"))


def complexes_of(center):
    out = defaultdict(list)
    for e in BETA:
        if e.get("center") == center:
            out[e["ligand"]].append((e["nu"], e["logb"]))
    return out


# center -> {阴离子键: (盐式, 投料 mol)}；碱盐统一 0.01 M 金属
PAIRS = {
    "Mg^{2+}": {"Br^-": ("MgBr2", 0.01), "I^-": ("MgI2", 0.01),
                "NO_3^-": ("Mg(NO3)2", 0.01), "ClO_4^-": ("Mg(ClO4)2", 0.01),
                "SO_4^{2-}": ("MgSO4", 0.01)},
    "Ca^{2+}": {"Br^-": ("CaBr2", 0.01), "I^-": ("CaI2", 0.01),
                "NO_3^-": ("Ca(NO3)2", 0.01), "ClO_4^-": ("Ca(ClO4)2", 0.01)},
    "Ni^{2+}": {"Br^-": ("NiBr2", 0.01), "I^-": ("NiI2", 0.01),
                "NO_3^-": ("Ni(NO3)2", 0.01), "ClO_4^-": ("Ni(ClO4)2", 0.01),
                "SO_4^{2-}": ("NiSO4", 0.01)},
    "Zn^{2+}": {"Br^-": ("ZnBr2", 0.01), "I^-": ("ZnI2", 0.01),
                "NO_3^-": ("Zn(NO3)2", 0.01), "ClO_4^-": ("Zn(ClO4)2", 0.01),
                "SO_4^{2-}": ("ZnSO4", 0.01)},
    "Mn^{2+}": {"Br^-": ("MnBr2", 0.01), "I^-": ("MnI2", 0.01),
                "NO_3^-": ("Mn(NO3)2", 0.01), "ClO_4^-": ("Mn(ClO4)2", 0.01),
                "SO_4^{2-}": ("MnSO4", 0.01)},
    "Cd^{2+}": {"Br^-": ("CdBr2", 0.01), "I^-": ("CdI2", 0.01),
                "NO_3^-": ("Cd(NO3)2", 0.01), "ClO_4^-": ("Cd(ClO4)2", 0.01),
                "SO_4^{2-}": ("CdSO4", 0.01)},
    "Cu^{2+}": {"Br^-": ("CuBr2", 0.01), "NO_3^-": ("Cu(NO3)2", 0.01),
                "SO_4^{2-}": ("CuSO4", 0.01)},
    "La^{3+}": {"Br^-": ("LaBr3", 0.01), "I^-": ("LaI3", 0.01),
                "NO_3^-": ("La(NO3)3", 0.01)},
    "Ce^{3+}": {"Br^-": ("CeBr3", 0.01), "I^-": ("CeI3", 0.01),
                "NO_3^-": ("Ce(NO3)3", 0.01)},
    "Nd^{3+}": {"Br^-": ("NdBr3", 0.01), "I^-": ("NdI3", 0.01),
                "NO_3^-": ("Nd(NO3)3", 0.01)},
    "Gd^{3+}": {"Br^-": ("GdBr3", 0.01), "I^-": ("GdI3", 0.01),
                "NO_3^-": ("Gd(NO3)3", 0.01)},
    "Lu^{3+}": {"Br^-": ("LuBr3", 0.01), "I^-": ("LuI3", 0.01),
                "NO_3^-": ("Lu(NO3)3", 0.01)},
    "Y^{3+}": {"Br^-": ("YBr3", 0.01), "I^-": ("YI3", 0.01),
               "NO_3^-": ("Y(NO3)3", 0.01)},
    "Al^{3+}": {"Br^-": ("AlBr3", 0.01), "I^-": ("AlI3", 0.01),
                "NO_3^-": ("Al(NO3)3", 0.01)},
}


def main():
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()

    def ph(salt, mol):
        try:
            r = judge([{"name": salt, "mol": mol}], {"V_L": 1.0}, T)
            fin = {e["name"]: e["mol"] for e in (r.get("final") or [])}
            if not fin:
                return None, "无 final（式子可能不存在）"
            return r.get("final_pH"), None
        except Exception as exc:                            # noqa: BLE001
            return None, f"{type(exc).__name__}"

    # 基准必须用**对所有金属都不配位**的阴离子（NO₃⁻）。
    # ⚠️ 第一版拿 Cl⁻ 当基准，于是"Cl⁻ 会配位"这件事被算进了 Δ，
    #    使 Cd（有 Cd–Br/Cd–I）看起来"应改变却变了"。见 lessons。
    REF = ("NO_3^-", "NO3")
    rows = []
    for center, pairs in sorted(PAIRS.items()):
        cx = complexes_of(center)
        ref_salt = REF_SALT.get(center)
        if ref_salt is None:
            continue
        base, err = ph(ref_salt, 0.01)
        if base is None:
            print(f"  !! {center} 基准盐 {ref_salt} 失败: {err}")
            continue
        n_ref = len(cx.get(REF[0], []))
        for key, (salt, mol) in sorted(pairs.items()):
            if key == REF[0]:
                continue
            n_ent = len(cx.get(key, []))
            other, err = ph(salt, mol)
            if other is None:
                rows.append((center, salt, n_ent, base, None, err, n_ref))
                continue
            rows.append((center, salt, n_ent, base, other, None, n_ref))

    print("=== 阴离子配对矩阵（0.01 M 金属；基准 = NO₃⁻，库内无 M–NO₃ 配体）===")
    print("  %-9s %-13s %-6s %-9s %-9s %-9s %s"
          % ("center", "对照盐", "M–Y", "NO3-基准", "对照 pH", "ΔpH", "判读"))
    bad = []
    for center, salt, n_ent, base, other, err, n_ref in rows:
        if err:
            print(f"  {center:<9} {salt:<13} {n_ent:<6} —  {err}")
            continue
        d = other - base
        if n_ent == 0:
            verdict = "一致" if abs(d) <= 0.05 else "**⚠ 应不变却变了**"
            if abs(d) > 0.05:
                bad.append((center, salt, d))
        else:
            verdict = "可变（库内有该配体条目）"
        print("  %-9s %-13s %-6d %-9.4f %-9.4f %+.4f  %s"
              % (center, salt, n_ent, base, other, d, verdict))
    print(f"\n库内无 M–Y 配体、但 pH 变化 >0.05 的：{len(bad)} 条")
    for c, s, d in bad:
        print(f"   {c} {s}  Δ={d:+.4f}")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "anion_matrix.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"center": r[0], "salt": r[1], "n_complex": r[2],
                        "ph_ref": r[3], "ph_other": r[4], "err": r[5]}
                       for r in rows], f, ensure_ascii=False, indent=1)
        print(f"[写入] {os.path.relpath(p, ROOT)}")
    return 0


REF_SALT = {
    "Mg^{2+}": "Mg(NO3)2", "Ca^{2+}": "Ca(NO3)2", "Ni^{2+}": "Ni(NO3)2",
    "Zn^{2+}": "Zn(NO3)2", "Mn^{2+}": "Mn(NO3)2", "Cd^{2+}": "Cd(NO3)2",
    "Cu^{2+}": "Cu(NO3)2", "La^{3+}": "La(NO3)3", "Ce^{3+}": "Ce(NO3)3",
    "Nd^{3+}": "Nd(NO3)3", "Gd^{3+}": "Gd(NO3)3", "Lu^{3+}": "Lu(NO3)3",
    "Y^{3+}": "Y(NO3)3", "Al^{3+}": "Al(NO3)3",
}

if __name__ == "__main__":
    sys.exit(main())
