"""断言脆弱性审计：对**物理上无意义的数值扰动**敏感的断言 = 锁轨迹的断言。

**背景**（architecture §7 V）：全量差分实测显示，把 `solve_extent` 的二分
提前 1e-11 收敛（状态差 ~1e-11·x_max，比任何判据阈值小 5 个数量级）会让
**42/1173 例**改变走步路径——这些用例处在判据边界（分岔点）。若断言锁的是
"某一轮数值轨迹恰好走出来的那条路"（具体的分步路线、整数配比、贴着阈值的
产量），那么在**物理等价**的实现变更下就会翻红，从而把算法改进挡在门外。

**判据（本工具的审计口径）**：把每个用例的投料量按 `(1±ε)` 相对扰动
（ε = 1e-11 / 1e-9，**远低于任何有化学意义的差别**：1e-9 相对 = 1 mol 投料
差 1 nmol），重跑断言。**仍然翻红的断言 = 对无意义扰动敏感 = 锁轨迹**，
应当改写成锁化学的形式（守恒 + 产量区间 / 净方程 / pH 区间），或者如实
记录"该用例确实处在化学边界上"。

用法：
    python tools/fragility.py               # ε=1e-11 与 1e-9，两种符号
    python tools/fragility.py -e 1e-9       # 只用给定 ε
"""
from __future__ import annotations

import copy
import sys

sys.path.insert(0, ".")

from chemkit.data import load_tables                      # noqa: E402
from chemkit import testsuit as ts                        # noqa: E402


def _kind(err: str) -> str:
    """错误串 → 断言类别（用于聚合"哪一类断言在翻红"）。"""
    if "净方程不符" in err:
        return "eq"
    if "多步方程缺" in err:
        return "eq_has"
    if "不在 [" in err:
        return "has_range"
    if "不应生成" in err:
        return "has_not"
    if "产量" in err:
        return "has"
    if err.startswith("pH") or " pH" in err or "不在" in err:
        return "ph"
    if "changed" in err or "reacted" in err or "degree" in err:
        return "判定层"
    if "缺标注" in err or "override" in err:
        return "标注层"
    return "其他"


def _run(cases: list[dict], eps: float) -> dict:
    """按 (1±ε) 扰动投料重跑；返回 {用例名: [错误串]}。"""
    out: dict[str, list[str]] = {}
    for c0 in cases:
        c = copy.deepcopy(c0)
        c["subs"] = [[n, m * (1.0 + eps)] for n, m in c["subs"]]
        n0 = len(ts.RESULTS)
        ts.run_case(c, T, verbose=False)
        if len(ts.RESULTS) > n0:
            if not ts.RESULTS[-1]["ok"]:
                out[c0["name"]] = list(ts.RESULTS[-1]["errors"])
            del ts.RESULTS[n0:]
    return out


T = None

if __name__ == "__main__":
    argv = sys.argv[1:]
    eps_list = [1e-11, 1e-9]
    if "-e" in argv:
        eps_list = [float(argv[argv.index("-e") + 1])]
    T = load_tables()
    ts.FAILS.clear()
    ts.PASS_N = 0
    cases = ts.load_cases(None)
    ts.judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, T)   # 预热
    base_fails = set()
    for c in cases:
        n0 = len(ts.RESULTS)
        ts.run_case(c, T, verbose=False)
        if len(ts.RESULTS) > n0:
            if not ts.RESULTS[-1]["ok"]:
                base_fails.add(c["name"])
            del ts.RESULTS[n0:]
    print(f"基线：{len(cases)} 例，红 {len(base_fails)} 例 {sorted(base_fails)}")
    agg: dict[str, set] = {}
    for eps in eps_list:
        for sign in (+1.0, -1.0):
            e = eps * sign
            ts.PASS_N = 0
            res = _run(cases, e)
            for nm, errs in res.items():
                if nm in base_fails:
                    continue
                agg.setdefault(nm, set()).update(
                    f"{_kind(x)}({'+' if sign > 0 else '-'}ε={eps:g})"
                    for x in errs)
            print(f"  ε={e:+.0e}: 新翻红 {len(res) - len(set(res) & base_fails)} 例")
    print(f"\n=== 对物理无意义扰动敏感的用例：{len(agg)} 例 ===")
    for nm in sorted(agg, key=lambda k: -len(agg[k])):
        print(f"  {nm[:46]:48s} {sorted(agg[nm])}")
