"""候选数据一致性审计：同一净变换的**正反写法** logK 必须互为相反数。

**为什么需要**（architecture §7 X-12 的副产物）：T34/MX01 一类"反复乒乓"
（`[Fe(SCN)₃] → Fe³⁺+3SCN⁻` 与 `Fe³⁺+3SCN⁻ → [Fe(SCN)₃]` 交替执行）在
走步里表现为张力例。乒乓若来自**数据**（同一对正反写法由两条不同数据源
生成、logK 不相消），那是数据缺陷；若来自**模型**（pH 悬崖让同一写法的根
在两侧来回），那是算法缺陷。本工具把前者一次性排掉。

判据（精确、无容差）：
  ① 同一净键的正反两条候选：`logK(正) + logK(反) == 0`；
  ② 同一净键出现多条候选（重复生成）：logK 必须一致；
  ③ 派生候选（kind=derived）的 logK 必须等于其基候选的线性组合
     （Hess 一致性；`tools/hess_audit.py` 已覆盖，这里只做交叉引用）。

用法：
    python tools/cand_audit.py             # 全量用例的初始态
    python tools/cand_audit.py --cases T34 MX01
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from chemkit.candidates import WATER, H_ION, logK_T          # noqa: E402
from chemkit.data import load_tables                          # noqa: E402
from chemkit.engine import _fmt                              # noqa: E402
from chemkit.normalize import normalize                       # noqa: E402
from chemkit.speciation import estimate_pH                    # noqa: E402
from chemkit.templates import enumerate_candidates            # noqa: E402
from chemkit.testsuit import load_cases                       # noqa: E402


def _eq(c) -> str:
    r = " + ".join(f"{_fmt(nu)}{s}" for s, nu in c.r.items() if s != WATER)
    p = " + ".join(f"{_fmt(nu)}{s}" for s, nu in c.pr.items() if s != WATER)
    return f"{r} -> {p}"


def main() -> None:
    want: tuple[str, ...] = ()
    if "--cases" in sys.argv:
        want = tuple(sys.argv[sys.argv.index("--cases") + 1:])
    T = load_tables()
    T_K = 298.15
    cases = load_cases(None)
    if want:
        cases = [c for c in cases if c["name"].startswith(want)]

    n_cand = 0
    bad_sym: list = []      # 正反不相消
    dup: list = []          # 同键重复且 logK 不一致
    seen_states = 0
    for case in cases:
        subs = [{"name": n, "mol": m} for n, m in case["subs"]]
        cond = case.get("cond") or {"V_L": 1.0}
        try:
            led, He, _esc, _unk = normalize(subs, cond, T)
        except Exception:
            continue
        V = float(cond.get("V_L", 1.0) or 1.0)
        T_Kc = float(cond.get("T_K", 298.15) or 298.15)
        seen_states += 1
        pH = estimate_pH(led, He, V, T, T_Kc)
        bykey: dict = {}
        for c in enumerate_candidates(led, He, pH, V, T_Kc, T, True):
            n_cand += 1
            lk = logK_T(c, T_Kc)
            k = c.netkey_fwd
            prev = bykey.get(k)
            if prev is None:
                bykey[k] = (lk, c)
            elif prev[0] != lk:
                dup.append((case["name"], _eq(prev[1]), prev[0], _eq(c), lk))
        # 正反检查：反向键存在时，两者 logK 必须相消
        for k, (lk, c) in bykey.items():
            rev = bykey.get((k[1], k[0]))
            if rev is not None and lk + rev[0] != 0.0:
                bad_sym.append((case["name"], _eq(c), lk, _eq(rev[1]), rev[0]))

    print(f"扫描 {seen_states} 个初始态、{n_cand} 条候选")
    print(f"\n① 正反写法 logK 不相消：{len(bad_sym)} 例")
    for nm, e1, k1, e2, k2 in bad_sym[:20]:
        print(f"   [{nm[:22]}] {e1[:56]} logK={k1:+.6f}")
        print(f"       反向 {e2[:56]} logK={k2:+.6f}  和={k1 + k2:+.3g}")
    print(f"\n② 同净键重复且 logK 不一致：{len(dup)} 例")
    for nm, e1, k1, e2, k2 in dup[:20]:
        print(f"   [{nm[:22]}] {e1[:56]} logK={k1:+.6f}")
        print(f"       重复 {e2[:56]} logK={k2:+.6f}  差={k1 - k2:+.3g}")


if __name__ == "__main__":
    main()
