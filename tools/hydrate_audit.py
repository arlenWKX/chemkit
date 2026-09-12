"""水合物/无水盐歧义审计 + 溶解度边界敏感性扫描。

两件事：
  A. **水合物歧义**：扫描 ksp.json，找出固相稳定析出形态是**水合物**、
     而条目给的是无水式（或未声明晶型）的记录。判据：note 提及水合/结晶
     水/晶型，或落在常见"文献 Ksp 按水合物报道"的清单里。
     引擎账本只有无水式 + 水（溶剂活度 1），所以**不需要新物种**——
     需要的是选对**相**的 Ksp 并在数据里写明晶型。

  B. **边界敏感性**：对每个用例，列出析出量很小的固相。析出量小 ⟹ 该固相
     在平衡面上"刚好成立"，±0.5 个 pKsp 就能让它消失 ⟹ 用水合物还是无水
     盐的 Ksp 会直接翻转结论（V16 即此：pQ=7.207 vs pKsp=7.2，析出
     0.0556 mol；换成二水盐的 6.5 则不析出）。

用法：python tools/hydrate_audit.py [-n 30]
"""
from __future__ import annotations

import json
import re
import sys

sys.path.insert(0, ".")

from chemkit.data import load_tables          # noqa: E402
from chemkit.engine import judge              # noqa: E402
from chemkit.testsuit import load_cases       # noqa: E402

# 常见"文献 Ksp 多按水合物报道"的难溶盐（无水式入库时最易踩坑）
KNOWN_HYDRATE = {
    "CaSO_3": "CaSO3·2H2O（pKsp 6.5）vs 无水 CaSO3（7.2）",
    "CaSO_4": "CaSO4·2H2O 石膏（4.6）vs 硬石膏（4.4–5.0）",
    "CaC_2O_4": "CaC2O4·H2O（8.6）vs 无水（8.0）",
    "CaCO_3": "方解石/文石无水（8.3–8.5）vs CaCO3·6H2O 六水（ikaite）",
    "MgNH_4PO_4": "MgNH4PO4·6H2O 鸟粪石",
    "Fe(OH)_3": "Fe(OH)3·xH2O 水合氧化铁（无定形 vs 针铁矿）",
    "Al(OH)_3": "Al(OH)3·xH2O 无定形 vs 三水铝石/勃姆石",
    "Zn(OH)_2": "无定形 vs ε-Zn(OH)2",
    "Mg(OH)_2": "brucite（无水）",
    "BaSO_4": "重晶石（无水）",
    "Na_2CO_3": "Na2CO3·10H2O（<32℃）/·H2O（>32℃）vs 无水",
    "NaHCO_3": "无水（注意：极易溶，pKsp<0 表形态）",
    "Na_2SO_4": "Na2SO4·10H2O 芒硝 vs 无水（32.4℃ 转变）",
    "Na_2S_2O_3": "Na2S2O3·5H2O 海波",
    "Li_2CO_3": "无水（590℃ 以上才生成，室温无水）",
    "MnCO_3": "无水菱锰矿",
    "CuSO_4": "CuSO4·5H2O 胆矾",
    "MgCO_3": "MgCO3·3H2O vs 菱镁矿",
    "NiCO_3": "NiCO3·xH2O",
    "CoCO_3": "CoCO3·xH2O",
}
HYD_WORDS = ("水合", "结晶水", "晶型", "水化物", "·H2O", "水合物", "hydrate")


def main(n: int) -> None:
    T = load_tables()
    rows = json.load(open("chemkit/data/ksp.json", encoding="utf-8"))

    print("== A. 水合物/晶型歧义（ksp.json）==")
    flag = []
    for e in rows:
        s, note = e["solid"], e.get("note") or ""
        hit = [w for w in HYD_WORDS if w in note]
        if hit or s in KNOWN_HYDRATE:
            flag.append((s, e["pKsp"], hit, KNOWN_HYDRATE.get(s, ""), note[:70]))
    print(f"   命中 {len(flag)} / {len(rows)} 条：")
    for s, pk, hit, known, note in flag:
        print(f"   · {s:18s} pKsp={pk:<7} {'注中提及'+str(hit) if hit else '注未提，但属常见水合物'}")
        if known:
            print(f"       已知形态差异：{known}")
        if note:
            print(f"       note: {note}")
    print("\n   其中**文献 Ksp 按水合物报道、而条目未声明晶型**的高风险项：")
    for s, pk, hit, known, _ in flag:
        if s in KNOWN_HYDRATE and not any(w in ("".join(x[4] for x in flag if x[0] == s))
                                          for w in ("水合", "晶型")):
            print(f"   · {s}  pKsp={pk}   {known}")

    print("\n== B. 边界敏感性（析出量小的固相）==")
    out = []
    for c in load_cases(None):
        subs = [{"name": x, "mol": m} for x, m in c["subs"]]
        r = judge(subs, c.get("cond") or {"V_L": 1.0}, T)
        solids = [(e["name"], e["mol"]) for e in r.get("final", [])
                  if e["name"] in T.solids and e["mol"] > 0]
        feed = sum(m for _x, m in c["subs"])
        for s, m in solids:
            if m < 0.2 and m < 0.2 * max(feed, 1e-9):
                out.append((m / max(feed, 1e-9), c["name"], s, m, feed,
                            s in KNOWN_HYDRATE))
    out.sort()
    print(f"   析出量 < 20% 投料的固相实例 {len(out)} 条，最小 {n} 条：")
    for frac, name, s, m, feed, risky in out[:n]:
        print(f"   {frac*100:6.2f}%  {s:16s} {m:.5f}/{feed:g} mol  "
              f"{'⚠水合物歧义' if risky else '          '}  {name[:44]}")


if __name__ == "__main__":
    nn = 30
    a = sys.argv[1:]
    if "-n" in a:
        nn = int(a[a.index("-n") + 1])
    main(nn)
