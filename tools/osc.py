"""走步过程审计：反复逼近（锯齿/乒乓）的量化画像。

**为什么需要这个指标**（用户判据："反复逼近才达到最终产物不是好兆头"）：
walk 是逐候选贪心（Gauss-Seidel 式坐标松弛），每一步把**单个**平衡解到
Q=K（其余物种当冻结）。耦合平衡因此靠"锯齿式"几何收敛——同一对净反应
来回改写、幅度按 ρ<1 递减。终态对，但：

  · 步数 = 收敛轮数 × O(1/ln ρ)，量级直接决定性能；
  · 停在收敛判据上时残差 ~ρ^n，正是"系数差 1e-3"的来源（N15 的
    `1.999CO_2`）；判据一放松就落到错态；
  · ρ→1 的张力体系（Ag32/J06 型）能爬到上千步仍不收敛。

指标：
  n_steps      走步落实的步数
  n_key        不同净反应的个数（步数 / n_key = 平均重复度）
  max_rep      单个净反应的最大重复次数
  alt          锯齿对数：相邻两步互为逆反应（同净键反向）的次数
  rep_ge3      重复 ≥3 次的净反应个数（真正的"反复逼近"签名）

用法：
    python tools/osc.py                # 全量：汇总 + 最差 25 例
    python tools/osc.py N15 H89 P11    # 指定用例前缀
    python tools/osc.py -n 40          # 最差 N 例
"""
from __future__ import annotations

import sys
from collections import Counter

sys.path.insert(0, ".")

from chemkit.data import load_tables                # noqa: E402
from chemkit.engine import judge                     # noqa: E402
from chemkit.testsuit import load_cases              # noqa: E402


def _key(eq: str) -> tuple:
    """净键：箭头两侧的物种集合（方向归一，逆反应同键）。"""
    if " -> " not in eq:
        return (eq,)
    a, b = eq.split(" -> ", 1)
    left = frozenset(x.strip() for x in a.split(" + "))
    right = frozenset(x.strip() for x in b.split(" + "))
    return (left, right) if repr(left) <= repr(right) else (right, left)


def profile(r: dict) -> dict:
    steps = [s for s in r.get("steps", []) if s.get("kind") != "neutralize"]
    keys = [_key(s["equation"]) for s in steps]
    cnt = Counter(keys)
    # 锯齿：相邻两步**同净键、反方向**（equation 文本两侧互换）——这正是
    # "反复逼近"的签名：同一条平衡被来回推，幅度按 ρ<1 递减。
    zig = sum(1 for i in range(1, len(steps))
              if keys[i] == keys[i - 1]
              and steps[i]["equation"] != steps[i - 1]["equation"])
    return {
        "n_steps": len(steps),
        "n_key": len(cnt),
        "max_rep": max(cnt.values(), default=0),
        "rep_ge3": sum(1 for v in cnt.values() if v >= 3),
        "zig": zig,
    }


def main() -> int:
    argv = sys.argv[1:]
    top_n = 25
    if "-n" in argv:
        i = argv.index("-n")
        top_n = int(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    pre = tuple(x for x in argv if not x.startswith("-"))
    T = load_tables()
    rows = []
    hist = Counter()
    for c in load_cases(None):
        if pre and not c["name"].startswith(pre):
            continue
        r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
        p = profile(r)
        p["name"] = c["name"]
        rows.append(p)
        e = p["n_steps"]
        hist["1" if e <= 1 else "2-3" if e <= 3 else "4-6" if e <= 6
             else "7-10" if e <= 10 else "11-20" if e <= 20
             else "21-50" if e <= 50 else "51+"] += 1
    print(f"用例 {len(rows)}")
    print("步数分布： " + "  ".join(f"{k}:{hist[k]}" for k in
                                 ("1", "2-3", "4-6", "7-10", "11-20", "21-50",
                                  "51+") if hist[k]))
    z = [r for r in rows if r["zig"] or r["rep_ge3"]]
    print(f"锯齿/重复≥3 的体系：{len(z)} 例"
          f"（占 {100.0 * len(z) / max(len(rows), 1):.1f}%）")
    print(f"总步数 {sum(r['n_steps'] for r in rows)}；"
          f"最多步 {max(r['n_steps'] for r in rows)}")
    print(f"\n== 最差 {top_n} 例（按 锯齿数+重复度 排序）==")
    rows.sort(key=lambda r: -(r["zig"] * 3 + r["rep_ge3"] * 2 + r["n_steps"]))
    for r in rows[:top_n]:
        print(f"  步{r['n_steps']:4d} 净键{r['n_key']:3d} 最大重复{r['max_rep']:3d} "
              f"重复≥3:{r['rep_ge3']:3d} 锯齿{r['zig']:3d}  {r['name'][:56]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
