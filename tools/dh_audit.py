# -*- coding: utf-8 -*-
"""第 247 轮 · **温度外推可信度**审计（van't Hoff 的 dH 是否物理可信）。

**动因（第 247 轮实测）**：一元弱酸半中和时 pH = pKa(T)，而
`pKa(T) = pKa(298) − dH/(R·ln10)·(1/298.15 − 1/T)`。
对大多数体系这条链路**逐位正确**（HCN/NH₄⁺/H₂S/HClO 实测 Δ ≤ 0.003），
但对少数体系它会给出**物理上不可能**的 pKa：

    HIO_4  pKa298=1.6  dH=+399.0 ⟹ pKa(373) = −12.4（比高氯酸还"强" 2 个数量级）
    HVO_4^{2-} dH=+243.0        ⟹ pKa(273) = 16.4

根因**不在 van't Hoff 公式**，而在输入 ΔfH：
  · `HIO_4` 的 ΔfH = −550 kJ/mol 是**固态偏高碘酸**的值
    （NIST：HIO₄(cr) ≈ −536），而它被当成**水相物种**参与酸碱平衡；
  · `HClO_4` 的 ΔfH = −10 kJ/mol，文献 HClO₄(aq) ≈ −40 ⟹ dH 被抬高约 3 倍。

**判据（本脚本）**：一元质子解离的 ΔH 物理量级
  参照：H₂O +55.8 / HF +12 / HSO₄⁻ +21 / 醋酸 −0.4（|ΔH| 一般 < 60 kJ/mol）。
  另加**后果判据**：|ΔpKa(273→373)| 若 > 3，则该条目在液态水温度域内
  会把 pKa 推出常见酸区 ⟹ **温度外推不可信**。

用法： python tools/dh_audit.py [--write]
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

R_LN10 = 8.31446261815324e-3 * math.log(10.0)
DH_SANE = 60.0          # kJ/mol，一元质子解离的经验上限
DPKA_MAX = 3.0          # 液态水温度域内可接受的 pKa 漂移


def main():
    from chemkit.data import load_tables
    T = load_tables()
    rows = []
    for e in T.pka:
        dH = e.get("dH")
        if dH is None or e.get("n") != 1:
            continue
        dpka = (dH / R_LN10 * (1 / 298.15 - 1 / 373.15)
                - dH / R_LN10 * (1 / 298.15 - 1 / 273.15))
        rows.append((abs(dH), abs(dpka), e["acid"], e["base"], e["pka"], dH))
    rows.sort(reverse=True)

    print(f"一元弱酸 pKa 条目（带 dH）：{len(rows)} 条")
    print(f"判据：|dH| > {DH_SANE} kJ/mol 或 |ΔpKa(273→373)| > {DPKA_MAX}\n")
    print("  %-16s %-8s %-9s %-11s %s"
          % ("acid", "pKa298", "dH", "ΔpKa", "判读"))
    bad = []
    for adh, adp, acid, base, pk, dH in rows[:20]:
        flag = ""
        if adh > DH_SANE:
            flag += " ⚠|dH|过大"
        if adp > DPKA_MAX:
            flag += " ⚠pKa漂移过大"
        if flag:
            bad.append((acid, pk, dH, adp))
        print("  %-16s %-8s %+9.1f %-11.2f%s" % (acid, pk, dH, adp, flag))
    print(f"\n合计需复核：{len(bad)} 条")
    for acid, pk, dH, adp in bad:
        print(f"   {acid:<16s} pKa298={pk:<7} dH={dH:+.1f}  ΔpKa={adp:.2f}")
    print("\n处置建议：**不自行改 ΔfH**（须逐条查文献并核相态），"
          "先记账；在这些条目上**不要写温度断言**。")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "dh_audit.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"acid": a, "pKa298": pk, "dH": dH, "dpKa": dp}
                       for a, pk, dH, dp in bad], f,
                      ensure_ascii=False, indent=1)
        print(f"[写入] {os.path.relpath(p, ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
