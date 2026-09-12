"""账本电荷守恒审计。

**原理**：候选反应逐个都是配平的（原子 + 电荷守恒），且引擎把游离质子
统一记为带符号标量 `H_excess`（正 = 游离强酸 mol）。因此恒等式

        Σ_s z_s·n_s(final) + H_excess(final) = 0

必须在任意经候选步推进的状态上成立。唯一可能破坏它的是"不经候选反应"
的记账环节：normalize 的酸碱中和、分子态强酸重组、以及投料条件
c_H/c_OH/pH。

**第一版审计（已作废）的错误**：直接判 `Σz·n + He ≠ 0` 为违规，却没
把 `cond["c_H"]`（初始强酸浓度 mol/L）计入预期游离质子——于是
E51(c_H=9)、N09(c_H=14)、E12(c_H=6)、E02(c_H=4) 这类**带酸条件**的
用例全被误报成"凭空生电"。把"工具不看的条件"当成"引擎的错"，是审计
工具自身的设计缺陷。

用法：python tools/charge_audit.py [--top N]
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chemkit.data import load_tables                    # noqa: E402
from chemkit.core import charge_of                      # noqa: E402
from chemkit.engine import judge                        # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _q(led: dict) -> float:
    return sum(charge_of(s) * m for s, m in led.items()
               if not s.startswith("__"))


def main(argv):
    top = 30
    if "--top" in argv:
        top = int(argv[argv.index("--top") + 1])
    T = load_tables()
    judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, T)
    cases = json.load(open(os.path.join(HERE, "chemkit/data/tests.json"),
                           encoding="utf-8"))
    rows = []
    for c in cases:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        cond = c.get("cond") or {"V_L": 1.0}
        r = judge(subs, cond, T)
        led = {e["name"]: e["mol"] for e in r["final"]}
        led.pop("H^+", None)
        led.pop("OH^-", None)
        he = r["H_excess"]
        dev = _q(led) + he          # 名义恒等式：Σz·n(final) = −He
        rows.append((abs(dev), c["name"], _q(led), he, dev,
                     bool(cond.get("c_H") or cond.get("c_OH")
                          or cond.get("pH") is not None)))
    rows.sort(reverse=True)
    bad = [x for x in rows if x[0] > 1e-3]
    noacid = [x for x in bad if not x[5]]
    print(f"总例 {len(rows)}；|Σz·n + He| > 1e-3: {len(bad)} "
          f"（{100.0 * len(bad) / len(rows):.1f}%）"
          f"；其中无酸/碱条件者 {len(noacid)}")
    for x in rows[:top]:
        print(f"  偏差={x[0]:9.4g}  Σz·n={x[2]:+10.4g}  He={x[3]:+10.4g}"
              f"  {'[带酸碱条件]' if x[5] else ''}  {x[1][:40]}")


if __name__ == "__main__":
    main(sys.argv[1:])
