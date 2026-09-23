# -*- coding: utf-8 -*-
"""第 248 轮 · 库内 ΔfH vs llnl.dat 权威值（找"同族里对不上的那一个"）。

**动因**：第 247 轮查出 8 条 pKa 的派生 dH 物理上不可信（`HIO_4` dH=+399）。
追到输入 ΔfH：`thermo.json` 的 `IO_4^-` = −151.0、`ClO_4^-` = −129.33，
与 llnl 的**计算焓**逐位一致（−36.2 kcal = −151.4；−30.91 kcal = −129.3），
**但 `HIO_4`/`HClO_4` 这些中性酸在 llnl 里根本没有对应物种**。

⟹ 假设：** anions 是校准过的，中性酸是可疑的**（`HIO_4` = −550 是
固态偏高碘酸的值）。本脚本把这个假设变成普查：
把 llnl 每个反应块的 `# Enthalpy of formation: X kcal/mol` 连同它上一行的
反应式一起解出来，做成"llnl 权威 ΔfH 表"，再与 `thermo.json` 比对。

判据：**同一张表的其它条目都对得上、只有某几条差得远** ⟹ 那几条要复核
（"同族一致性"判据，不依赖我对外部文献的记忆）。

用法： python tools/thermo_xcheck.py [--write]
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

KCAL = 4.184          # kcal -> kJ


def parse_llnl(path):
    """返回 {反应式: ΔfH(kJ/mol)}（llnl 注释里的 Calculated enthalpy of formation）。

    ⚠️ 分段必须按 **PHREEQC 的段落约定：段名行不缩进**。
    第一版用 `re.split(r"\\n(?=\\S)")`，于是注释行 `#\\tEnthalpy of formation:`
    （`#` 是非空白）也被当成新段开头，段头全变成注释行 ⟹ **匹配 0 条**。
    """
    t = io.open(path, encoding="utf-8", errors="replace").read()
    out = {}
    cur = None
    for line in t.splitlines():
        s = line.strip()
        if not s:
            continue                           # 空行不改变当前段
        # 段名行 = 不缩进且不是注释（注释行也可能顶格，不能当段名）
        if not line[0].isspace() and not s.startswith("#"):
            cur = s
            continue
        m = re.search(r"Enthalpy of formation:\s*([-\d.]+)\s*kcal/mol", line)
        if m and cur:
            # ⚠️ llnl 对"算不出来的"条目也会写一行 `Enthalpy of formation:`
            # 但值为 0 或空 ⟹ 必须滤掉，否则会造出 −1453 kJ/mol 这种假差异。
            if "Not possible" in line:
                continue
            try:
                v = float(m.group(1)) * KCAL
            except ValueError:
                continue
            if abs(v) < 1e-9:
                continue
            out[cur] = v
    return out


def main():
    from chemkit.data import load_tables
    T = load_tables()
    lp = os.path.join(ROOT, "thirdparty/phreeqc-3.8.6-17100/database/llnl.dat")
    if not os.path.exists(lp):
        print(f"!! 找不到 {lp}")
        return 2
    ll = parse_llnl(lp)
    print(f"llnl 解出带 ΔfH 的反应块：{len(ll)} 个")
    print(f"库内 thermo 条目：{len(T.thermo)}\n")

    # 库内 ΔfH 单位应为 kJ/mol；建立"物种名 -> ΔfH"
    ours = {k: v for k, v in T.thermo.items() if isinstance(v, (int, float))}
    print(f"thermo 里是数值的条目：{len(ours)}")

    # 直接比对：把 llnl 的 ΔfH 归给"出现在反应式里的那个'非常规试剂'物种"
    # 这里只做**可判定的窄比对**：单元素阴离子/中性物种，且库内有同名条目。
    norm = {}
    for head, dh in ll.items():
        parts = [p.strip() for p in head.split("=")]
        if len(parts) != 2:
            continue
        lhs, rhs = parts
        # 去掉 O2 / H2O / H+ / e- 这些"构成项"
        skip = {"O2", "H2O", "H+", "e-", "H2", "O2(g)", "CO2", "H2O(l)"}
        cand = [p for p in re.split(r"\s*\+\s*", rhs) if p.strip() not in skip]
        if len(cand) == 1:
            norm.setdefault(cand[0].strip(), dh)

    # 库内键名 -> llnl 键名 的宽松映射（只查可唯一对应的）
    def key_candidates(k):
        yield k
        yield k.replace("^", "")
        yield k.replace("_", "")
    matched, diffs, only_ours = 0, [], []
    for k, v in sorted(ours.items()):
        hit = None
        for kk in key_candidates(k):
            if kk in norm:
                hit = norm[kk]
                break
        if hit is None:
            only_ours.append(k)
            continue
        matched += 1
        d = v - hit
        if abs(d) > 20.0:                 # >20 kJ/mol 即值得复核
            diffs.append((abs(d), k, v, hit, d))

    print(f"\n可与 llnl 对上名的：{matched} 条")
    diffs.sort(reverse=True)
    print(f"其中差 > 20 kJ/mol 的：{len(diffs)} 条\n")
    print("  %-16s %-12s %-12s %-10s" % ("物种", "库内", "llnl", "差"))
    for _a, k, v, h, d in diffs[:25]:
        print("  %-16s %-12.1f %-12.1f %+.1f" % (k, v, h, d))
    if "--write" in sys.argv:
        import json
        p = os.path.join(ROOT, "logs", "thermo_xcheck.json")
        with io.open(p, "w", encoding="utf-8") as f:
            json.dump({"matched": matched,
                       "diffs": [{"species": k, "ours": v, "llnl": h,
                                  "d": d} for _a, k, v, h, d in diffs]},
                      f, ensure_ascii=False, indent=1)
        print(f"\n[写入] {os.path.relpath(p, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
