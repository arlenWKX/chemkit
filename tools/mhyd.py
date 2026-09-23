# -*- coding: utf-8 -*-
"""第 242 轮 · 单金属水解 **独立解**（不依赖引擎）+ 与引擎对拍。

对"M^n+ 盐 + 水"体系，pH 可**独立算出**：未知量 ([M], [L], pH)。
给定 pH 时 [M] 由金属质量平衡显式给出、[L] 由配体质量平衡给出
⟹ 外层二分 pH、内层求游离配体不动点。

这个独立解只用**库内常数**、与引擎代码无关 ⟹ 可当**判据**。
但它有**适用域**，四点必须一起满足（第 242 轮四次修正换来的）：

1. **配体质量平衡不可省**：把投料浓度当游离 [L] 会把 `ScCl₃` 误判成
   差 4.46 个 pH 单位（`[ScCl]^{2+}` logβ=0.55 看着弱，但必须计入）。
2. **先定相态、再定根**：顺序不能颠倒。液相/固相是两条不同曲线，
   不分相态就在 [0,14] 二分会得到混合假根（`LaCl₃` 液相分支根本不与零轴相交）。
3. **电荷要用相对量**：绝对电荷两项各 ~1e-2、真值 ~1e-5，相减只剩 3 位有效数字
   ⟹ 扫根漏掉过零点、表现为"无解"。故用 `Σz[i] / Σ|z|[i]` 并在 pOH 上参数化。
4. **固相分支只钉游离 [M]**：溶解总量 = [M]·(1+Σβ[OH]^ν)，其余进固相；
   写错则电荷恒负、无根。

用法：
  python tools/mhyd.py --calib    # 与引擎对拍（配位强/弱 × 饱和/未饱和）

**校准现状（第 242 轮实测）**：17 个体系里 **9 个稀土 + Mg** 与引擎落在
≤0.03 pH 单位内。**未通过校准的体系，其 Δ 不可当缺陷证据。**

已知限制：
  · **固相分支不完整**：pH_sat 之上、需要"析出量"参与电荷平衡的体系
    （`CaCl₂` 0.01 M 即此类，当前报 −5.37）**给出错误值**；
    **液相（未饱和）分支才是可信的**——稀土那批用例正是靠它。
  · 只支持"单金属中心 + 单配体"，不支持多配体 / 多金属 / 氧化还原。
"""
import io
import json
import math
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


def hyd_of(center):
    return sorted((e["nu"], e["logb"], e["complex"]) for e in BETA
                  if e.get("center") == center and e.get("ligand") == "OH^-")


def lig_of(center, lig):
    return sorted((e["nu"], e["logb"], e["complex"]) for e in BETA
                  if e.get("center") == center and e.get("ligand") == lig)


def ksp_of(center):
    """该金属中心 + OH⁻ 的氢氧化物 Ksp（若有）：返回 (pKsp, x, y)。"""
    for e in KSP:
        if e.get("pair") == [center, "OH^-"]:
            x, y = _ksp_xy(e)
            return e["pKsp"], x, y
    return None


def independent(center, c_m, lig=None, c_l=0.0, pkw=PKW):
    """独立解。返回 dict：pH / 物种 / 相态 / pH_sat。

    算法（第 242 轮第三次修正才定下来，**顺序不能颠倒**）：

    ① **先求饱和边界 pH_sat**：在"无固相"分支上找 Q(pH) = Ksp 的根。
       La(OH)₃ 实测 pH_sat = 7.43。
    ② **在 pH_sat 处算电荷不平衡量**，用它选分支 —— 这一步是关键：
       · 不平量为正 ⟹ 未饱和区就够电中性 ⟹ 根在**液相分支**（≤ pH_sat）；
       · 不平量为负 ⟹ 液相怎么调都补不上负电荷 ⟹ 金属**必须析出**，
         根在**固相分支**（此时游离 [M] 被 Ksp 钉住）。
    ③ 在选定分支上求电荷方程的第一个根。

    ⚠️ 为什么不能对 [0,14] 直接二分（前两版的错）：
      液相分支的电荷函数在 LaCl₃ 上**单调**（pH 3 时 −0.029 → pH 10 时
      −0.0398），根本不在 pH 5 附近变号；直接全区间二分/粗扫会落到
      **pH≈12 的假根**上（那里靠 [OH⁻] 自身平衡凑数），或干脆漏掉真根。
      **"先定相态、再定根"是正确顺序。**
    """
    zc = charge_of(center)
    zl = charge_of(lig) if lig else 0
    ks = ksp_of(center)

    def spec(pH, solid):
        oh = 10.0 ** (-(pkw - pH))
        h = 10.0 ** (-pH)
        hyd = hyd_of(center)
        ligs = lig_of(center, lig) if lig else []
        # 游离配体不动点
        if ligs and c_l > 0:
            lo, hi = 0.0, c_l
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                den = 1.0
                for nu, lb, _nm in hyd:
                    den += (10.0 ** lb) * oh ** nu
                for nu, lb, _nm in ligs:
                    den += (10.0 ** lb) * mid ** nu
                m_t = c_m / den if not (solid and ks) else None
                if m_t is None:
                    pKsp, x, y = ks
                    m_t = 10.0 ** (-pKsp / x) / oh ** (y / x)
                bound = sum(nu * m_t * (10.0 ** lb) * mid ** nu
                            for nu, lb, _nm in ligs)
                if mid + bound > c_l:
                    hi = mid
                else:
                    lo = mid
                if hi - lo < 1e-25:
                    break
            free_l = 0.5 * (lo + hi)
        else:
            free_l = 0.0
        if solid and ks:
            # ⚠️ 固相分支：被 Ksp 钉住的是**游离** [M]；溶解总量
            # = [M]·(1 + Σβ_OH[OH]^ν + Σβ_L[L]^ν)，其余进固相。
            # （此前把"溶解总量"也写成游离值 ⟹ 电荷恒负 ⟹ 找不到根。）
            pKsp, x, y = ks
            m = 10.0 ** (-pKsp / x) / oh ** (y / x)
        else:
            den = 1.0
            for nu, lb, _nm in hyd:
                den += (10.0 ** lb) * oh ** nu
            for nu, lb, _nm in ligs:
                den += (10.0 ** lb) * free_l ** nu
            m = c_m / den
        sp = {center: m}
        for nu, lb, nm in hyd:
            sp[nm] = m * (10.0 ** lb) * oh ** nu
        for nu, lb, nm in ligs:
            sp[nm] = m * (10.0 ** lb) * free_l ** nu
        return m, free_l, sp, h, oh

    def charge(pH, solid):
        _m, _fl, sp, h, oh = spec(pH, solid)
        q = sum(charge_of(k) * v for k, v in sp.items()) + h - oh
        return q - (zc * c_m + zl * c_l)

    def qlog(pH):
        """Q/Ksp 的常用对数（**液相分支**，用该 pH 下的游离配体）。"""
        if not ks:
            return None
        pKsp, x, y = ks
        _m, _fl, _sp, _h, oh = spec(pH, False)
        den = 1.0
        for nu, lb, _nm in hyd_of(center):
            den += (10.0 ** lb) * oh ** nu
        for nu, lb, _nm in (lig_of(center, lig) if lig else []):
            den += (10.0 ** lb) * _fl ** nu
        m = c_m / den
        return (x * math.log10(m) + y * math.log10(oh)) + pKsp

    def root(f, a, b, n=4000):
        """在 [a,b] 上找**第一个**过零点并二分细化。"""
        prev_x = a
        prev = f(a)
        for i in range(1, n + 1):
            x = a + (b - a) * i / n
            cur = f(x)
            if prev == 0.0:
                return prev_x
            if prev * cur < 0:
                lo, hi = prev_x, x
                flo = prev
                for _ in range(200):
                    mid = 0.5 * (lo + hi)
                    fm = f(mid)
                    if flo * fm < 0:
                        hi = mid
                    else:
                        lo, flo = mid, fm
                return 0.5 * (lo + hi)
            prev_x, prev = x, cur
        return None

    def charge_rel(pOH):
        """**相对电荷不平衡** Σz·[i] / Σ|z|·[i]（**在 pOH 上参数化**）。

        ⚠️ 为什么要相对量 + pOH 参数化（第 242 轮第四次修正）：
        绝对电荷 = (阳离子电荷 − 阴离子电荷)，两项各约 10⁻² 而真值约 10⁻⁵
        ⟹ 相减只剩 3 位有效数字，扫根时会**漏掉过零点**（表现为"无解"）。
        改参数化为 pOH 后，`10**-pOH` 不会像 `10**-(pKw-pOH)` 那样把
        pH 的舍入放大到 [OH⁻] 上，再取相对量即恢复满精度。
        """
        pH = pkw - pOH
        _m, _fl, sp, h, oh = spec(pH, False)
        pos = 0.0
        neg = 0.0
        for k, v in sp.items():
            z = charge_of(k)
            if z > 0:
                pos += z * v
            else:
                neg += -z * v
        pos += h
        neg += oh
        bg = zc * c_m + zl * c_l
        if bg > 0:
            neg += bg
        else:
            pos += -bg
        tot = pos + neg
        return (pos - neg) / tot if tot > 0 else 0.0

    def qlog_p(pOH):
        pH = pkw - pOH
        _m, _fl, _sp, _h, oh = spec(pH, False)
        if not ks:
            return None
        pKsp, x, y = ks
        den = 1.0
        for nu, lb, _nm in hyd_of(center):
            den += (10.0 ** lb) * oh ** nu
        for nu, lb, _nm in (lig_of(center, lig) if lig else []):
            den += (10.0 ** lb) * _fl ** nu
        m = c_m / den
        return (x * math.log10(m) + y * math.log10(oh)) + pKsp

    def charge_rel_solid(pOH):
        pH = pkw - pOH
        _m, _fl, sp, h, oh = spec(pH, True)
        pos = neg = 0.0
        for k, v in sp.items():
            z = charge_of(k)
            if z > 0:
                pos += z * v
            else:
                neg += -z * v
        pos += h
        neg += oh + zc * c_m + zl * c_l
        tot = pos + neg
        return (pos - neg) / tot if tot > 0 else 0.0

    out = {"pH": None, "sp": {}, "phase": None, "pH_sat": None,
           "solid": (f"x={ks[1]} y={ks[2]} pKsp={ks[0]}" if ks else None),
           "dbg": {}}
    # ① 饱和边界
    pH_sat = root(qlog, 0.0, 14.0) if ks else None
    out["pH_sat"] = pH_sat
    # ② 在边界处判分支
    c_sat = charge(pH_sat, False) if pH_sat is not None else None
    solid = (c_sat is not None and c_sat < 0)
    out["dbg"] = {"charge_at_sat": c_sat, "solid": solid}
    out["dbg"]["chg_true_8"] = charge(8.0, True)
    out["dbg"]["chg_true_10"] = charge(10.0, True)
    out["dbg"]["chg_false_8"] = charge(8.0, False)
    # ③ 在选定分支求根（液相限制在 ≤pH_sat 内找，避开假根区）
    if solid:
        pH = root(lambda x: charge(x, True), 0.0, 14.0)
        if pH is not None:
            _m, _fl, sp, _h, _o = spec(pH, True)
            out.update({"pH": pH, "sp": sp, "phase": "固相饱和"})
            return out
    lo_b = (pH_sat - 1e-9) if pH_sat is not None else 14.0
    pH = root(lambda x: charge(x, False), 0.0, max(lo_b, 1e-9))
    if pH is not None:
        _m, _fl, sp, _h, _o = spec(pH, False)
        out.update({"pH": pH, "sp": sp, "phase": "无固相"})
        return out
    out["phase"] = "无解"
    return out


def engine_pH(subs, V=1.0):
    from chemkit.data import load_tables
    from chemkit.engine import judge
    T = load_tables()
    pr = {}
    r = judge([{"name": n, "mol": m} for n, m in subs],
              {"V_L": V}, T, _probe=pr)
    return r.get("final_pH"), pr


CASES = [
    ("LaCl3", "La^{3+}", 0.01, "Cl^-", 0.03),
    ("CeCl3", "Ce^{3+}", 0.01, "Cl^-", 0.03),
    ("NdCl3", "Nd^{3+}", 0.01, "Cl^-", 0.03),
    ("GdCl3", "Gd^{3+}", 0.01, "Cl^-", 0.03),
    ("HoCl3", "Ho^{3+}", 0.01, "Cl^-", 0.03),
    ("TmCl3", "Tm^{3+}", 0.01, "Cl^-", 0.03),
    ("LuCl3", "Lu^{3+}", 0.01, "Cl^-", 0.03),
    ("EuCl3", "Eu^{3+}", 0.01, "Cl^-", 0.03),
    ("YCl3", "Y^{3+}", 0.01, "Cl^-", 0.03),
    ("ScCl3", "Sc^{3+}", 0.01, "Cl^-", 0.03),
    ("FeCl3", "Fe^{3+}", 0.001, "Cl^-", 0.003),
    ("AlCl3", "Al^{3+}", 0.001, "Cl^-", 0.003),
    ("MgCl2", "Mg^{2+}", 0.01, "Cl^-", 0.02),
    ("CaCl2", "Ca^{2+}", 0.01, "Cl^-", 0.02),
    ("NiCl2", "Ni^{2+}", 0.01, "Cl^-", 0.02),
    ("CoCl2", "Co^{2+}", 0.01, "Cl^-", 0.02),
    ("MnCl2", "Mn^{2+}", 0.01, "Cl^-", 0.02),
]


def run(group, title):
    print(f"== {title} ==")
    print("  %-8s %8s %9s %9s %9s %9s %s"
          % ("case", "C/M", "独立pH", "pH_sat", "相态", "引擎pH", "Δ"))
    bad = []
    for label, center, c, lig, cl in group:
        ind = independent(center, c, lig, cl)
        try:
            epH, _pr = engine_pH([(label, c)])
        except Exception as exc:                            # noqa: BLE001
            epH = None
            print("   engine error", label, type(exc).__name__, exc)
        ipH = ind["pH"]
        d = None if (ipH is None or epH is None) else epH - ipH
        print("  %-8s %8g %9s %9s %9s %9s %s"
              % (label, c,
                 ("%.4f" % ipH) if ipH is not None else "None",
                 ("%.3f" % ind["pH_sat"]) if ind["pH_sat"] is not None else "-",
                 ind["phase"] or "-",
                 ("%.4f" % epH) if epH is not None else "None",
                 ("%+.4f" % d) if d is not None else "-"))
        if ipH is not None and d is not None and abs(d) > 0.05:
            bad.append((label, d, ind["phase"]))
    print(f"  |Δ| > 0.05: {len(bad)}/{len(group)}")
    for label, d, ph in bad:
        print(f"     {label} ({ph}): Δ={d:+.4f}")
    print()
    return bad


def main():
    if "--calib" in sys.argv or not sys.argv[1:]:
        run(CASES, "校准组：独立解 vs 引擎（配位强/弱 × 饱和/未饱和）")
    else:
        print("用法: python tools/mhyd.py --calib | --lanth")


if __name__ == "__main__":
    main()
