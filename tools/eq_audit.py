"""测试库 `eq`/`eq_has` 方程式的守恒审查（含大系数舍入容差）。

判定分两级：
  · **精确**：元素差向量与净电荷全为零；
  · **舍入容差**：|差值| ≤ max|系数|·2e-3（测试比对签名按"最大系数归一
    + round(…,3)"生成，大系数混合通道式（银氨/硫氰配位族）真实比为
    无理数，取整后残差随最大系数线性放大——这类是**呈现层的诚实近似**，
    不是化学错误）。
不落入任何一级者，即为测试标准自身的化学错误。

用法：python tools/eq_audit.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chemkit.core import elements_of, charge_of        # noqa: E402
from chemkit.system import _parse_equation             # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(eq: str):
    """返回 (level, detail)：level ∈ {"exact", "round", "bad"}。"""
    try:
        r, p = _parse_equation(eq)
    except Exception as exc:                            # noqa: BLE001
        return "bad", f"解析失败 {exc}"
    if not r or not p:
        return "bad", "单侧为空"
    el = {}
    q = 0
    mx = 0.0
    for d, sg in ((r, -1), (p, 1)):
        for sp, nu in d.items():
            mx = max(mx, abs(nu))
            q += sg * nu * charge_of(sp)
            for e, c in elements_of(sp).items():
                el[e] = el.get(e, 0) + sg * nu * c
    dev_el = max((abs(v) for v in el.values()), default=0.0)
    dev = max(dev_el, abs(q))
    if dev == 0.0:
        return "exact", ""
    tol = mx * 2e-3
    if dev <= tol:
        return "round", f"残差 {dev:g} ≤ {tol:g}（最大系数 {mx:g} 的舍入界）"
    return "bad", f"元素差 {[ (k, v) for k, v in el.items() if v ]} 净电荷 {q:+g}"


def main():
    cases = json.load(open(os.path.join(HERE, "chemkit/data/tests.json"),
                           encoding="utf-8"))
    n_eq = n_round = 0
    bad = []
    for c in cases:
        for tag in ("eq", "eq_has"):
            v = c.get(tag)
            if v is None:
                continue
            items = v if isinstance(v, list) else [v]
            for eq in items:
                if not isinstance(eq, str):
                    continue
                n_eq += 1
                lv, det = check(eq)
                if lv == "round":
                    n_round += 1
                elif lv == "bad":
                    bad.append((c["name"], tag, eq, det))
    print(f"期望方程式共 {n_eq} 条：精确守恒 {n_eq - n_round - len(bad)}、"
          f"舍入容差内 {n_round}、**不守恒 {len(bad)}**")
    for nm, tag, eq, det in bad:
        print(f"  [{tag}] {nm[:40]:42s} {eq}")
        print(f"        → {det}")


if __name__ == "__main__":
    main()
