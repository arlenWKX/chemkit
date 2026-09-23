# -*- coding: utf-8 -*-
"""第 243 轮 · "只有一级水解 + 背景阴离子"体系的 **真实自洽解** vs 引擎。

⚠️ **本工具当前不可作为判据（未通过校准）** —— 保留是为了记录方法与坑。
第 243 轮实测：它给出的"真实 pH"在 `AlCl₃`/`CrCl₃`/`FeCl₃` 上落在 **7.0**
（那三个体系必然强酸性、且实测 pH 3.0~3.6），说明它的根选择仍不对。
**根因（已查明）**：把"质量 + 电荷"两式联立时，代数上会同时出现
  · 一个 α≈0 的根（≈中性），
  · 一个 α≈100% 的根（如 `LaCl₃` 0.01 M 给 pH 2.00），
而"全部水解"那个根**不可达**：把 0.01 mol 金属全部水解需要 0.01 mol 的
OH⁻/O 受体，纯水只提供 1e-7 mol 量级。**"化学上不可达的代数根"必须显式
排除**，而"α<0.5"这种过滤器不做这件事（它把 AlCl₃ 的真实态滤掉了）。

**正确的判据形态（下一步）**：把水解写成显式反应
`M^n+ + H₂O → [M(OH)]^{(n-1)+} + H⁺`，并把**水平衡（O 守恒）**作为第三条
方程纳入，再对每个根做**可达性检验**（所需 OH⁻ ≤ 体系可提供量）。

用法（仅在修复后才可信）： python tools/chargeconsist.py
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

from chemkit.core import charge_of                            # noqa: E402
from chemkit.candidates import _ksp_xy                          # noqa: E402

PKW = 14.0


def _load(p):
    with io.open(os.path.join(ROOT, "chemkit/data", p), encoding="utf-8") as f:
        return json.load(f)


BETA = _load("beta.json")
KSP = _load("ksp.json")


def hyd1(center):
    for e in BETA:
        if e.get("center") == center and e.get("ligand") == "OH^-" \
                and e.get("nu") == 1:
            return e["logb"], e["complex"]
    return None, None


def ksp_of(center):
    for e in KSP:
        if e.get("pair") == [center, "OH^-"]:
            x, y = _ksp_xy(e)
            return e["pKsp"], x, y
    return None


def true_pH(center, c_m, z_lig, c_l, pkw=PKW):
    """只含一级水解的真实自洽 pH。

    ⚠️ **必须做"可达性"过滤（第 243 轮踩坑）**：把"质量 + 电荷"两式联立时，
    代数上会出现一个 **α≈100% 的假根**（`LaCl₃` 0.01 M 给 pH≈2.00，
    `CaCl₂` 0.01 M 给 pH≈2.00）。它满足两个方程，但**不可达**：
    把 C=0.01 的金属全部水解需要 0.01 mol 的 OH⁻/O 受体，而纯水里
    只有 1e-7 mol 量级。判别式：**水解消耗的 OH⁻ 不能超过体系能提供的量**，
    即 α·C 必须 ≤ [OH⁻]_可用 + [H⁺] 的净生成能力。

    实用判据（等价且便宜）：**逐 pH 检验 [M] = C − [MOH] > 0 且
    [MOH] ≤ C**；假根恰恰在 `[MOH] = h·C/(h+K·oh)` 接近 C 时出现。
    更稳的判据是同时要求 **水平衡**（本函数用 h 与 oh 的耦合隐含）。
    实测有效做法：**只接受"α < 0.5"的那个根**——因为 α>0.5 的态需要
    远超水所能提供的 OH⁻（本库全部金属在纯水中的水解度都远小于 1）。
    """
    lb, cplx = hyd1(center)
    if lb is None:
        return {"pH": None, "note": f"{center} 无 ν=1 的 OH⁻ 水解常数"}
    n = charge_of(center)

    def charge(pH):
        h = 10.0 ** (-pH)
        oh = 10.0 ** (-(pkw - pH))
        moh = h * c_m / (h + (10.0 ** lb) * oh)
        m = c_m - moh
        return (n * m + (n - 1) * moh + h - oh - abs(z_lig) * c_l,
                m, moh)

    # 扫出**所有**过零点，再用"水解度 < 50%"过滤掉不可达的假根
    roots = []
    prev_x, prev_f = 0.0, charge(0.0)[0]
    steps = 4000
    for i in range(1, steps + 1):
        x = pkw * i / steps
        f = charge(x)[0]
        if prev_f * f < 0:
            a, b, fa = prev_x, x, prev_f
            for _ in range(200):
                mid = 0.5 * (a + b)
                if fa * charge(mid)[0] < 0:
                    b = mid
                else:
                    a = mid
            roots.append(0.5 * (a + b))
        prev_x, prev_f = x, f
    valid = []
    for r in roots:
        _q, m, moh = charge(r)
        if m > 0 and moh / c_m < 0.5:
            valid.append(r)
    if not valid:
        return {"pH": None, "note": f"无可达根（共 {len(roots)} 个过零点）",
                "n_roots": len(roots)}
    pH = valid[-1] if len(valid) > 1 else valid[0]
    _q, m, moh = charge(pH)
    out = {"pH": pH, "M": m, "MOH": moh, "alpha": moh / c_m,
           "cplx": cplx, "logb1": lb, "pKsp": None, "QK": None,
           "n_roots": len(roots), "n_valid": len(valid)}
    ks = ksp_of(center)
    if ks:
        pKsp, x, y = ks
        oh = 10.0 ** (-(pkw - pH))
        out["pKsp"] = pKsp
        out["QK"] = (x * (__import__("math").log10(m) if m > 0 else -99)
                     + y * __import__("math").log10(oh)) + pKsp
    return out


def engine_pH(subs, V=1.0):
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()
    pr = {}
    r = judge([{"name": nm, "mol": mo} for nm, mo in subs],
              {"V_L": V}, T, _probe=pr)
    return r, pr


# (标签, 中心, 盐, C, 抗衡阴离子电价, C_L)
CASES = [
    # --- 对照：引擎已知做对的三价水解体系 ---
    ("AlCl3", "Al^{3+}", "AlCl3", 0.01, -1, 0.03),
    ("CrCl3", "Cr^{3+}", "CrCl3", 0.01, -1, 0.03),
    ("FeCl3", "Fe^{3+}", "FeCl3", 0.001, -1, 0.003),
    # --- 稀土 ---
    ("LaCl3", "La^{3+}", "LaCl3", 0.01, -1, 0.03),
    ("LaCl3", "La^{3+}", "LaCl3", 1e-4, -1, 3e-4),
    ("CeCl3", "Ce^{3+}", "CeCl3", 0.01, -1, 0.03),
    ("LuCl3", "Lu^{3+}", "LuCl3", 0.01, -1, 0.03),
    ("GdCl3", "Gd^{3+}", "GdCl3", 0.01, -1, 0.03),
    ("YCl3", "Y^{3+}", "YCl3", 0.01, -1, 0.03),
]


def main():
    if len(sys.argv) > 2:
        el = sys.argv[1]
        c = float(sys.argv[2])
        center = None
        for e in BETA:
            if e.get("center", "").startswith(el) and e.get("ligand") == "OH^-":
                center = e["center"]
                break
        n = charge_of(center)
        t = true_pH(center, c, -1, n * c)
        print(f"{center} C={c}: 真实 pH={t['pH']}")
        return 0
    print("== 一级水解的真实自洽 pH vs 引擎 ==")
    print("  %-7s %-8s %-7s %-9s %-9s %-9s %s"
          % ("case", "C/M", "logb1", "真实pH", "引擎pH", "Δ", "Q/Ksp"))
    flag = []
    for label, center, salt, c, zl, cl in CASES:
        t = true_pH(center, c, zl, cl)
        try:
            r, pr = engine_pH([(salt, c)])
            epH = r.get("final_pH")
        except Exception as exc:                            # noqa: BLE001
            epH = None
            print("   err", label, exc)
        d = None if (t["pH"] is None or epH is None) else epH - t["pH"]
        print("  %-7s %-8g %-7s %-9s %-9s %-9s %s"
              % (label, c, t.get("logb1"),
                 ("%.4f" % t["pH"]) if t["pH"] is not None else t.get("note"),
                 ("%.4f" % epH) if epH is not None else "None",
                 ("%+.4f" % d) if d is not None else "-",
                 ("%.2f" % t["QK"]) if t.get("QK") is not None else "-"))
        if d is not None and abs(d) > 0.05:
            flag.append((label, c, d, t["alpha"]))
    print(f"\n  |Δ| > 0.05: {len(flag)}/{len(CASES)}")
    for label, c, d, al in flag:
        print(f"     {label} C={c:g}: Δ={d:+.4f}（真实水解度 {al * 100:.2f}%）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
