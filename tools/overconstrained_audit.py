# -*- coding: utf-8 -*-
"""第 257 轮 · 普查「过度约束态」族：**账本内同族物种的比值是否与库内常数一致**。

**判据（第 256 轮定下）**：一个末态若满足"质量守恒 + 电荷守恒"，但**账本里同族
两个物种的比值**与库内 `Ka`（及该 pH）给出的比值不符，则该态**不满足解离平衡** ——
即"被多条约束同时要求，其中必有一条不成立"。

**算法**：
  1. 逐例 `judge(..., _probe={})`，取 `probe["ledger"]` 与 `probe["pH"]`；
  2. 用 `build_families(T)` 把账本物种归族；
  3. 对**同一族内同时出现 ≥2 个物种**的例，按库内 `pKa` 算"应有比值"
     `w_i/w_0 = Π 10^(sgn_k·(pKa_k − pH))`，与**账本实际比值**比；
  4. 相对偏差 > 1% 即记为"不一致"。

**为什么要普查**：第 243/254/256 轮分别在 `LaCl₃`/`H₂O₂`/`HIO₄` 上撞到同一结构，
但**从未量过全库有多少例**。本脚本补这一刀。

用法： python tools/overconstrained_audit.py [--write]
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

TOL_REL = 0.01          # 比值相对偏差阈值


def main():
    from chemkit.data import load_tables
    from chemkit.engine import judge
    from chemkit.acidbase import build_families, eff_pka
    import chemkit.testsuit as ts

    T = load_tables()
    fams = build_families(T)
    cases = ts.load_cases(None)

    rows = []
    nerr = 0
    for c in cases:
        cond = c.get("cond") or {"V_L": 1.0}
        tk = cond.get("T_K") or 298.15
        pr = {}
        try:
            judge([{"name": n, "mol": m} for n, m in c["subs"]], cond, T,
                  _probe=pr)
        except Exception:                                   # noqa: BLE001
            nerr += 1
            continue
        led = pr.get("ledger")
        pH = pr.get("pH")
        if not led or pH is None:
            continue
        # 归族：族 -> {物种: 量}
        groups = {}
        for s, m in led.items():
            if m <= 0.0 or s.startswith("__") or s in ("H_2O", "H^+", "OH^-"):
                continue
            info = fams.get(s)
            if info is None:
                continue
            groups.setdefault(info[0], {})[s] = m
        for fid, members in groups.items():
            if len(members) < 2:
                continue
            # 取该族的 order/pka/sgn
            info = fams[next(iter(members))]
            order, pka, dh, sgn = info[1], info[3], info[4], info[5]
            ek = eff_pka(pka, dh, tk)
            # 应有权重
            logs = {order[0]: 0.0}
            acc = 0.0
            for k in range(len(ek)):
                acc += (ek[k] - pH) * (1 if sgn is None else sgn[k])
                logs[order[k + 1]] = acc
            # 账本里出现的成员，两两比值
            present = [sp for sp in order if sp in members]
            if len(present) < 2:
                continue
            bad = 0.0
            for i in range(len(present)):
                for j in range(i + 1, len(present)):
                    sp_i, sp_j = present[i], present[j]
                    if members[sp_j] <= 0:
                        continue
                    ratio_led = members[sp_i] / members[sp_j]
                    ratio_th = 10.0 ** (logs[sp_i] - logs[sp_j])
                    rel = abs(ratio_led - ratio_th) / max(ratio_th, 1e-300)
                    bad = max(bad, rel)
            rows.append((c["name"], tk, fid, present,
                         {sp: members[sp] for sp in present}, bad,
                         bool(c.get("ok", True))))

    print(f"用例 {len(cases)}；有 ledger 的检查完成（异常 {nerr}）")
    print(f"**同族内出现 ≥2 物种**的例：{len(rows)}")
    viol = [r for r in rows if r[5] > TOL_REL]
    print(f"其中**比值与库内常数不一致**（相对偏差 > {TOL_REL:.0%}）的："
          f"**{len(viol)}**\n")
    print("  %-46s %-8s %-16s %-9s %s"
          % ("用例", "T/K", "族", "最大相对偏差", "通过?"))
    for nm, tk, fid, present, mem, bad, ok in sorted(viol,
                                                     key=lambda r: -r[5])[:25]:
        print("  %-46s %-8.2f %-16s %-9.3g %s"
              % (nm[:46], tk, str(fid)[:16], bad, "PASS" if ok else "FAIL"))
    if len(viol) > 25:
        print(f"  … 另 {len(viol) - 25} 例")
    npass = sum(1 for r in viol if r[6])
    print(f"\n  不一致例中：PASS {npass}、FAIL {len(viol) - npass}")
    print(f"  不一致例占全部'同族多样物种'例的 "
          f"{100.0 * len(viol) / max(1, len(rows)):.1f}%")
    if "--write" in sys.argv:
        p = os.path.join(ROOT, "logs", "overconstrained_audit.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump([{"name": r[0], "T": r[1], "family": str(r[2]),
                        "species": r[3], "amounts": r[4],
                        "max_rel_dev": r[5], "ok": r[6]} for r in viol],
                      f, ensure_ascii=False, indent=1)
        print(f"\n[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
