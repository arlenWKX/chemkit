"""chemkit.testsuit：统一测试模块（默认不随包加载，__init__ 不导出）。

运行：
    python -m chemkit.testsuit                # 默认 ./data/tests.json + 环闭合检查
    python -m chemkit.testsuit path/to.json   # 指定用例库

内容：
  1) 用例运行器：加载 tests.json（默认包内 ./data/tests.json），逐例调 judge
     并校验 reacted/degree/ann/override/has/has_not/has_range/has_any/ph，
     以及离子方程式断言 eq（净离子方程式，规范化比对）与
     eq_has（多步离子方程式中必须出现的步骤式）。
  2) 温度域/气压校验：273.15–373.15K 域外、非正气压必须报错不静默。
  3) 热力学环闭合检查：任何数据不许单点存在，必须能在 Hess 网络里自洽。

用例字段（断言维度）：
  subs     [["化学式", mol], ...]            投料
  cond     {"V_L":..., "T_K"/"T_C":..., "p_kpa":...}  条件（可选）
  changed  bool                              体系是否发生显著净变化（宽口径，
                                            含纯溶解/电离/水解等形态变化）
  reacted  bool                              狭义化学反应（不含纯溶解/电离/水解）
  degree   int/str/[int|str, ...]            程度（可给可接受列表；int=2/1/0，
                                            str=complete/incomplete/hardly/none）
  has      {化学式: 最少 mol}                产物/终态下限
  has_not  {化学式: 上限 mol}                不应生成（默认上限外的量即失败）
  has_range {化学式: [lo, hi]}
  has_initial {化学式: 最少 mol}             初态下限（post-normalize）
  has_in_initial_not {化学式}                初态不应含
  ph       [lo, hi]
  ann      ["slow"/"blocked"...]
  override OVERRIDE id
  eq       净离子方程式（如 "Zn + 2H^+ -> H_2 + Zn^{2+}"）；
           两侧内部物种顺序无关，系数必须一致。null 表示要求无净方程
  eq_has   [离子方程式, ...]                 每个都必须出现在 equations 中
  note     备注（不校验）
"""
from __future__ import annotations
import json
import os
import sys
import time

from .data import load_tables, Tables, _half_balance
from .engine import judge
from .equations import Equation
from .system import Reaction, _parse_equation

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
DEFAULT_CASES = os.path.join(DATA_DIR, "tests.json")

PASS_N = 0
FAILS: list[str] = []
TIMES: list[tuple[float, str]] = []   # 逐例计时（秒，用例名）
# 逐例结构化结果（name/ok/ms/errors/...）——供 `--out` 落盘，
# 便于完整保留一轮测试结果并随时读取比对（不必重跑）。
RESULTS: list[dict] = []
# 辅助检查电池（T 区间 / 酸碱锚点 / API 边界）的**逐块**结果：每项
# {name, ok, passed, fails:[文本]}。用例走 RESULTS，电池走这里——
# 两者都进摘要与结构化产物（否则"摘要含电池、产物只有用例"必然打架）。
CHECKS: list[dict] = []


def load_cases(path: str | None = None) -> list[dict]:
    """加载测试用例库（默认包内 ./data/tests.json）。"""
    with open(path or DEFAULT_CASES, encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------- 方程式比对

def _eq_signature(eq: str) -> tuple | None:
    """方程式的规范化签名：({反应物: 系数}, {产物: 系数})，两侧内部顺序无关。

    系数按比例归一（以最大系数为 1），容差 1e-3 比较——断言关注化学计量
    关系而非书写顺序。无法解析返回 None。

    **H⁺/OH⁻/H₂O 规范形（v0.5.0）**：净方程写成 H⁺ 还是 OH⁻ 是**呈现规范**
    而非化学——两者只差水的自电离关系 `H⁺ + OH⁻ → H₂O`：

        H₂S + OH⁻ → HS⁻ + H₂O   ≡   H₂S → HS⁻ + H⁺        （Q05）
        CO₂ + OH⁻ → HCO₃⁻       ≡   CO₂ + H₂O → HCO₃⁻ + H⁺ （H39）

    用例库两种写法都有（B12/H78 期望解离式、H39 期望中和式、E46/N24 期望
    水解式），引擎的呈现政策（终态 pH > 7 用 OH⁻）不可能同时满足——把它当
    正确性判据就是把格式选择当成化学事实。故比较前**先消去 OH⁻**：
    `OH⁻(消耗 c) → H₂O(消耗 c) + H⁺(生成 c)`（反向同理），再跨侧抵消。
    规范形对元素/电荷守恒与氧化还原配平**完全保真**（只挪动 H/O 记账位置），
    真正的错解（少一个原子、电荷不平）照样不匹配。
    """
    r, p = _parse_equation(eq)
    if not r or not p:
        return None
    # 消去 OH⁻（把水自电离的规范自由度固定下来）
    for side, other in ((r, p), (p, r)):
        c = side.pop("OH^-", 0.0)
        if c:
            side["H_2O"] = side.get("H_2O", 0.0) + c
            other["H^+"] = other.get("H^+", 0.0) + c
    # 跨侧抵消（规范形可能把物种挪到另一侧）
    for sp in set(r) & set(p):
        x = min(r[sp], p[sp])
        r[sp] -= x
        p[sp] -= x
    r = {s: v for s, v in r.items() if v > 1e-9}
    p = {s: v for s, v in p.items() if v > 1e-9}
    if not r or not p:
        # 规范形为空 = 该方程**就是**水的自电离关系本身（`H⁺ + OH⁻ → H₂O`
        # 及等价写法）——纯中和体系的标准式正是它。返回空签名（可比较），
        # 与"解析失败"（None）区分开：前者是合法方程，后者不是。
        return ({}, {}) if not r and not p else None
    mx = max(list(r.values()) + list(p.values()))
    if mx <= 0:
        return None
    return ({s: round(nu / mx, 3) for s, nu in r.items()},
            {s: round(nu / mx, 3) for s, nu in p.items()})


def _eq_match(expected: str, actual: str | None) -> bool:
    """方程式匹配：方向无关（正反应与逆反应视为同一）。

    引擎按记账顺序产出步骤方程，方向可能与教材期望相反（如
    'CH_3COOH -> CH_3COO^- + H^+' vs 'Ac^- + H^+ -> HAc'）。化学上等价，
    断言应放过。"""
    if actual is None:
        return False
    es = _eq_signature(expected)
    as_ = _eq_signature(actual)
    if es is None or as_ is None:
        return False
    return es == as_ or (es[0], es[1]) == (as_[1], as_[0])


def _eq_vec(eq: str) -> dict | None:
    """方程式 → 净向量（消耗为正、生成为负；`Fraction` 精确系数）。

    用于 `eq_has` 的**张成判定**（见 `_step_spanned`）。"""
    from fractions import Fraction
    r, p = _parse_equation(eq)
    if not r or not p:
        return None
    v: dict = {}
    for s, nu in r.items():
        v[s] = v.get(s, Fraction(0)) + Fraction(nu).limit_denominator(10 ** 6)
    for s, nu in p.items():
        v[s] = v.get(s, Fraction(0)) - Fraction(nu).limit_denominator(10 ** 6)
    return {s: x for s, x in v.items() if x != 0}


def _step_spanned(need: dict, steps: list[dict]) -> bool:
    """`need` 能否由实际步骤的**非负线性组合**得到（化学等价分解的判据）。

    **为什么需要它**（architecture §7 V/W）：`eq_has` 原本要求被断言的步骤
    **逐字出现**在分步叙述里，于是断言锁住了"某一轮数值轨迹恰好走出来的
    那条路线"。实测：把投料做 ±1e-11 相对扰动（1 mol 差 1 nmol，物理上无
    意义），T13/P11 的 `eq_has` 就翻红；把二分提前 1e-11 收敛，N19/Z31 同样
    翻红。而路线完全可以**化学等价地分解**——例如

        H⁺ + [Al(OH)₄]⁻ → Al(OH)₃ + H₂O
        ≡  (4H⁺ + [Al(OH)₄]⁻ → 4H₂O + Al³⁺) + (3OH⁻ + Al³⁺ → Al(OH)₃)
           + 3(H⁺ + OH⁻ → H₂O)

    两者描述同一净变换。判据应当是"叙述**覆盖**了这个变换"，而不是"逐字
    出现这一条"。实现：把各步与需求都写成净向量（精确有理数），解
    Σ cᵢ·stepᵢ = need；有解且 cᵢ ≥ 0 即通过（负系数 = 要倒着走某一步，
    那不是一条合法叙述）。纯 stdlib 精确算术，无容差。
    """
    cols = [_eq_vec(a) for a in steps]
    cols = [c for c in cols if c]
    if not cols:
        return False
    # **规范方向**（§7 O-5）：H⁺ + OH⁻ → H₂O 是呈现规范而非化学，故叙述可以
    # 自由地"加上/减去"这一条——提供 (H⁺ + OH⁻ − H₂O) 与 −(·) 两个方向 ⟹
    # `4H⁺ + [Al(OH)₄]⁻ → 4H₂O + Al³⁺` 与 `3OH⁻ + Al³⁺ → Al(OH)₃` 的组合
    # 即可覆盖 `H⁺ + [Al(OH)₄]⁻ → Al(OH)₃ + H₂O`（N19 型路线）。
    from fractions import Fraction
    gauge = {"H^+": Fraction(1), "OH^-": Fraction(1), "H_2O": Fraction(-1)}
    n_step = len(cols)
    cols.append(gauge)
    cols.append({k: -v for k, v in gauge.items()})
    species = sorted(set(need) | {s for c in cols for s in c})
    rows = [[c.get(s, 0) for c in cols] for s in species]
    rhs = [need.get(s, 0) for s in species]
    n = len(cols)
    # 增广矩阵上做精确高斯消元（列主元）
    aug = [rows[i] + [rhs[i]] for i in range(len(species))]
    piv_col: list[int] = []
    r = 0
    for cidx in range(n):
        piv = None
        for i in range(r, len(aug)):
            if aug[i][cidx] != 0:
                piv = i
                break
        if piv is None:
            continue
        aug[r], aug[piv] = aug[piv], aug[r]
        pv = aug[r][cidx]
        aug[r] = [x / pv for x in aug[r]]
        for i in range(len(aug)):
            if i != r and aug[i][cidx] != 0:
                f = aug[i][cidx]
                aug[i] = [a - f * b for a, b in zip(aug[i], aug[r])]
        piv_col.append(cidx)
        r += 1
        if r == len(aug):
            break
    # 一致性：全零行而右端非零 ⟹ 无解
    for i in range(r, len(aug)):
        if all(x == 0 for x in aug[i][:n]) and aug[i][n] != 0:
            return False
    sol = [0] * n
    for i, cidx in enumerate(piv_col):
        sol[cidx] = aug[i][n]
    # 非负约束只施加于**真实步骤**：规范方向（H⁺+OH⁻⇌H₂O）是自由的记账
    # 自由度，允许负系数——否则 `A + B − 3·gauge` 型分解会被误判为不合法。
    if any(x < 0 for x in sol[:n_step]):
        return False          # 需要倒着走某一步 ⟹ 不是合法叙述
    # 代回校验（数值上由精确算术保证，防御性再算一遍）
    for s in species:
        got = sum(c.get(s, 0) * k for c, k in zip(cols, sol))
        if got != need.get(s, 0):
            return False
    return True


def check_equations(c: dict, rxn: Reaction, errs: list[str]) -> None:
    """离子方程式断言：eq（净方程，可为 null 要求无净方程）与
    eq_has（多步方程式必须**覆盖**的变换，见 `_step_spanned`）。"""
    if "eq" in c:
        exp = c["eq"]
        if exp is None:
            if rxn.net_equation is not None:
                errs.append(f"不应有净方程，实际 {rxn.net_equation}")
        elif not _eq_match(exp, rxn.net_equation):
            errs.append(f"净方程不符：实际 {rxn.net_equation} 期望 {exp}")
    for exp in c.get("eq_has", []):
        if any(_eq_match(exp, a) for a in rxn.equations):
            continue                      # 逐字出现（快路径）
        need = _eq_vec(exp)
        if need is not None and _step_spanned(need, rxn.equations):
            continue                      # 被叙述的步骤组合覆盖（化学等价分解）
        errs.append(f"多步方程缺 {exp}（实际 {rxn.equations}）")


def run_case(c: dict, T: Tables, verbose: bool = True) -> bool:
    """运行单条用例并校验全部断言维度；返回是否通过。"""
    global PASS_N
    name = c["name"]
    subs = [{"name": n, "mol": m} for n, m in c["subs"]]
    t0 = time.time()
    # 探针（收敛画像）顺带取回：`judge` 本来只多填一个 dict，成本可忽略；
    # 存到模块级供并行 worker 回传（第 214 轮），避免为拿残差再跑一遍 judge。
    _pr: dict = {}
    r = judge(subs, c.get("cond") or {"V_L": 1.0}, T, _probe=_pr)
    _LAST_PROBE.clear()
    _LAST_PROBE.update(_pr)
    TIMES.append((time.time() - t0, name))
    errs = []
    # 量值校验统一用 最终态∪产物 的合并视图
    # v0.4.0 弱池口径：断言池成员（Zn²⁺/Ni²⁺/Cu²⁺ 及弱卤配形态）时按池
    # 总量聚合——O02 has Zn²⁺ 0.95 是池口径（溶解态总量），物种级 free
    # Zn²⁺ 0.74 属形态细节。has_not/has_range 维持物种级（上限/区间断言
    # 锁的是具体形态量）
    pool_map = r.get("pool_map") or {}
    amt = {}
    amt_sp = {}   # 物种级视图（has_not/has_range 用：上限/区间锁具体形态量）
    pool_tot: dict[str, float] = {}
    for e in r["production"] + r["final"]:
        sp = e["name"]
        amt_sp[sp] = max(amt_sp.get(sp, 0.0), e["mol"])
        a = pool_map.get(sp)
        if a is not None:
            pool_tot[a] = pool_tot.get(a, 0.0) + e["mol"]
        else:
            amt[sp] = max(amt.get(sp, 0.0), e["mol"])
    for a, v in pool_tot.items():
        amt[a] = max(amt.get(a, 0.0), v)
    # 初态视图（post-normalize，强电解质电离 + 气体处理 + 中和已记账）
    init_map: dict[str, float] = {e["name"]: e["mol"]
                                  for e in r.get("initial", [])}
    if "changed" in c and r["changed"] != c["changed"]:
        errs.append(f"changed={r['changed']} 期望{c['changed']}")
    if "reacted" in c and r["reacted"] != c["reacted"]:
        errs.append(f"reacted={r['reacted']} 期望{c['reacted']}")
    if "degree" in c:
        exp_raw = c["degree"]
        exp_list = exp_raw if isinstance(exp_raw, list) else [exp_raw]
        # 同时支持 str（complete/...）与 int（2/1/0）；引擎输出为 int
        _DSTR_TO_INT = {"complete": 2, "incomplete": 1,
                        "hardly": 0, "none": 0}
        exp_ints = set()
        for v in exp_list:
            if isinstance(v, int):
                exp_ints.add(v)
            else:
                exp_ints.add(_DSTR_TO_INT.get(v, -1))
        if r["degree"] not in exp_ints:
            errs.append(f"degree={r['degree']} 期望{exp_list}")
    for a in c.get("ann", []):
        if a not in r["annotations"]:
            errs.append(f"缺标注 {a}")
    if "override" in c and r.get("override") != c["override"]:
        errs.append(f"override={r.get('override')} 期望{c['override']}")
    for sp, lo in c.get("has", {}).items():
        m = amt.get(sp, 0.0)
        if m < lo:
            errs.append(f"{sp} 产量 {m:.4g} < {lo}")
    for sp, hi in c.get("has_not", {}).items():
        m = amt_sp.get(sp, 0.0)
        if m > hi:
            errs.append(f"{sp} 不应生成 {m:.4g} > {hi}")
    for sp, (lo, hi) in c.get("has_range", {}).items():
        m = amt_sp.get(sp, 0.0)
        if not (lo <= m <= hi):
            errs.append(f"{sp}={m:.4g} 不在 [{lo},{hi}]")
    for sp, lo in c.get("has_any", {}).items():
        m = amt.get(sp, 0.0)
        if m < lo:
            errs.append(f"{sp}(any) 产量 {m:.4g} < {lo}")
    for sp, lo in c.get("has_initial", {}).items():
        m = init_map.get(sp, 0.0)
        if m < lo:
            errs.append(f"初态 {sp}={m:.4g} < {lo}")
    for sp in c.get("has_in_initial_not", {}):
        m = init_map.get(sp, 0.0)
        if m > 1e-6:
            errs.append(f"初态不应含 {sp}（实际 {m:.4g}）")
    if "ph" in c and r["final_pH"] is not None:
        if not (c["ph"][0] <= r["final_pH"] <= c["ph"][1]):
            errs.append(f"ph={r['final_pH']} 不在 {c['ph']}")
    # 离子方程式断言（仅当用例声明了 eq/eq_has 时才构建 Reaction）
    if "eq" in c or "eq_has" in c:
        check_equations(c, Reaction(r), errs)
    if errs:
        FAILS.append(name)
        RESULTS.append({"index": c.get("__gidx") or (len(RESULTS) + 1),
                        "name": name, "ok": False,
                        "ms": round((time.time() - t0) * 1000, 2), "errors": errs,
                        "note": c.get("note") or "",
                        "pH": r.get("final_pH"), "degree": r.get("degree"),
                        "changed": r.get("changed"),
                        "annotations": list(r.get("annotations") or []),
                        "net_equation": (
                            None if Reaction(r).net_equation is None
                            else Reaction(r).net_equation.plain())
                        if ("eq" in c or "eq_has" in c) else None})
        if verbose:
            print(f"[FAIL] {name}  -- {'; '.join(errs)}"
                  + (f"  ({c['note']})" if c.get("note") else ""))
            print("   steps:", [(s["equation"], s["extent"]) for s in r["steps"][:6]])
            print("   production:", [(p["name"], p["mol"]) for p in r["production"]],
                  "| pH:", r["final_pH"], "| degree:", r["degree"],
                  "| ann:", r["annotations"])
        return False
    PASS_N += 1
    RESULTS.append({"index": c.get("__gidx") or (len(RESULTS) + 1),
                    "name": name, "ok": True,
                    "ms": round((time.time() - t0) * 1000, 2), "errors": [],
                    "note": c.get("note") or "",
                    "pH": r.get("final_pH"), "degree": r.get("degree"),
                    "changed": r.get("changed"),
                    "annotations": list(r.get("annotations") or []),
                    "net_equation": None})
    if verbose:
        print(f"[PASS] {name}" + (f"  ({c['note']})" if c.get("note") else ""))
    return True


def case_T_range(T: Tables) -> int:
    """温度域与气压校验：
       - 273.15–373.15K 之外（非常压液态水）必须报错不静默
       - p_kpa <= 0 必须报错（负压/零压无物理意义）
    返回通过数。"""
    global PASS_N
    for tk, ok in [(250.0, False), (400.0, False), (273.15, True), (373.15, True)]:
        try:
            judge([], {"T_K": tk}, T)
            if ok:
                PASS_N += 1
            else:
                FAILS.append(f"T_K={tk} 未报错")
        except ValueError:
            if ok:
                FAILS.append(f"T_K={tk} 被误拒")
            else:
                PASS_N += 1
    # p_kpa 边界：0 / 负数必须报错；正常气压应通过
    for pk, ok in [(0.0, False), (-10.0, False), (50.0, True), (101.3, True), (500.0, True)]:
        try:
            judge([{"name": "NaCl", "mol": 0.01}], {"p_kpa": pk}, T)
            if ok:
                PASS_N += 1
            else:
                FAILS.append(f"p_kpa={pk} 未报错")
        except ValueError:
            if ok:
                FAILS.append(f"p_kpa={pk} 被误拒")
            else:
                PASS_N += 1
    return PASS_N


# ============================================================ 酸碱地基检查
# 教科书对照（溶液 pH 的解析值或文献值；浓度与投料组成均写全）。
ACIDBASE_ANCHORS: list[tuple] = [
    ("HCl 1M",       {"Cl^-": 1.0},                        0.00),
    ("NaOH 1M",      {"Na^+": 1.0},                        14.00),
    ("NaAc 1M",      {"Na^+": 1.0, "CH_3COO^-": 1.0},       9.38),
    ("NH4Cl 0.1M",   {"NH_4^+": 0.1, "Cl^-": 0.1},          5.12),
    ("NH3 0.1M",     {"NH_3": 0.1},                        11.12),
    ("NaHCO3 0.1M",  {"Na^+": 0.1, "HCO_3^-": 0.1},         8.34),
    ("Na2CO3 0.1M",  {"Na^+": 0.2, "CO_3^{2-}": 0.1},      11.62),
    ("H3BO3 0.1M",   {"H_3BO_3": 0.1},                      5.12),
    ("NaHSO3 0.1M",  {"Na^+": 0.1, "HSO_3^-": 0.1},         4.50),
    ("NaH2PO4 0.1M", {"Na^+": 0.1, "H_2PO_4^-": 0.1},       4.65),
    ("Na2HPO4 0.1M", {"Na^+": 0.2, "HPO_4^{2-}": 0.1},      9.75),
    ("TlCl 饱和",     {"Tl^+": 0.0135, "Cl^-": 0.0135},      7.00),
]


def acidbase_check(T: Tables) -> int:
    """精确酸碱地基（`chemkit.acidbase`）的常驻不变量检查。返回通过数。

    两道检查，都是"数据自身的化学约束"，不依赖任何引擎输出：

    1. **族分布必须复现条目自己的 Ka**（本轮抓到的缺陷类）：pKa 条目的定义
       是 `A ⇌ B + H⁺` ⟹ `[B]/[A] = Ka/[H⁺] = 10^(pH − pKa)`。族内逐级
       权重若在**标记酸恰是少质子一侧**的条目（H₃BO₃/[B(OH)₄]⁻、Tl⁺/TlOH、
       CO₂/HCO₃⁻…）上取错方向，整个族的分布会倒过来——实测硼酸 pH 7 下
       `[B(OH)₄]⁻/[H₃BO₃] = 174`（真值 0.0058）、TlCl 溶液被算成 pH 1.87。
       该检查对全部 n==1 条目逐条比对，是全库级的（不是抽样）。
    2. **教科书锚点**：12 个溶液的 pH 与解析/文献值对照（±0.15）。

    这两条都落在"数据 + 精确地基"层，与判定引擎的走步无关，改数据或改
    `build_families`/`dist_charge` 都会被立刻抓住。
    """
    global PASS_N
    from .acidbase import build_families, charge_pH, dist_charge, eff_pka
    from .core import charge_of
    fams = build_families(T)
    T_K = 298.15
    pH_ref = 7.0
    n_ok = 0
    for e in T.pka:
        if e.get("n", 1) != 1:
            continue
        a, b, pka = e.get("acid"), e.get("base"), e.get("pka")
        if a is None or b is None or a not in fams or b not in fams:
            continue
        fa = fams[a]
        if fa[0] != fams[b][0] or a not in fa[1] or b not in fa[1]:
            continue
        order, pkaT, dh, sgn = fa[1], fa[3], fa[4], fa[5]
        ek = eff_pka(pkaT, dh, T_K)
        logw = [0.0]
        acc = 0.0
        for k in range(len(ek)):
            acc += (ek[k] - pH_ref) * sgn[k]
            logw.append(acc)
        mx = max(logw)
        ws = [10.0 ** (v - mx) for v in logw]
        # 注：此处原有 `sw = sum(ws)`——**死变量**，从未被消费（第 203 轮清理）。
        # 比值用 ws 的两项直接相除，`mx` 平移已在比值中相消 ⟹ 归一化本来就不需要。
        got = ws[order.index(b)] / ws[order.index(a)]
        want = 10.0 ** (pH_ref - pka)
        if abs(got / want - 1.0) > 0.02:
            FAILS.append(f"族分布不复现 Ka：{a} / {b} pKa={pka} "
                         f"族比 {got:.4g} 应为 {want:.4g}")
        else:
            n_ok += 1
    for name, led, want in ACIDBASE_ANCHORS:
        got = charge_pH(led, 1.0, T, T_K)
        if got is None or abs(got - want) > 0.15:
            FAILS.append(f"酸碱锚点 {name}：解 {got} 期望 {want}")
        else:
            n_ok += 1
    _ = dist_charge, charge_of
    PASS_N += n_ok
    return n_ok


# ============================================================ 环闭合检查
PKW = 14.0
K = 0.05916


def consistency(T: Tables, verbose: bool = True) -> list[str]:
    """热力学环闭合检查（数据纪律：任何数据不许单点存在）。返回失败列表。"""
    fails = []

    def check(name, lhs, rhs, tol=0.05):
        if abs(lhs - rhs) > tol:
            fails.append(f"[FAIL] {name}: {lhs:.3f} vs {rhs:.3f} "
                         f"(差 {abs(lhs - rhs):.3f} > {tol})")
        elif verbose:
            print(f"[ok] {name}: {lhs:.3f} ≈ {rhs:.3f}")

    # 1) h 形电对与 oh 文献参考值的换算闭合：E_h − E_oh = k·(h/n)·pKw
    #    （h 不入库，运行时配平重算）
    for c in T.couples:
        if "oh_ref" in c:
            bal = _half_balance(c["ox"], c["red"])
            expect = K * ((bal[1] if bal else 0) / c["n"]) * PKW
            check(f"电对 h/oh 换算 {c['ox']}/{c['red']}",
                  c["E0"] - c["oh_ref"], expect, 0.02)

    # 2) 多元酸分级 pKa 之和 = 合并条目
    for e_sum in T.pka:
        if e_sum["n"] >= 2:
            chain = [e for e in T.pka if e["n"] == 1 and
                     (e["acid"] == e_sum["acid"] or e["base"] == e_sum["base"])]
            if len(chain) >= e_sum["n"]:
                s = sum(sorted([e["pka"] for e in chain])[:e_sum["n"]])
                check(f"多元酸 pKa 和 {e_sum['acid']}→{e_sum['base']}",
                      s, e_sum["pka"], 0.3)

    # 3) 配位电对 E0 与 logβ 自洽：E(complex/M) = E(M^n+/M) − k·logβ/n
    for c in T.couples:
        if c["ox"] in T.beta_by_complex:
            b = T.beta_by_complex[c["ox"]]
            ref = next((x for x in T.couples
                        if x["ox"] == b["center"] and x["red"] == c["red"]), None)
            if ref:
                check(f"配位电对 {c['ox']}/{c['red']} ↔ logβ",
                      c["E0"], ref["E0"] - K * b["logb"] / c["n"], 0.05)

    # 4) Fe(OH)3/Fe(OH)2 与 Fe3+/Fe2+ 及两侧 Ksp 自洽
    fe = next(c for c in T.couples
              if c["ox"] == "Fe(OH)_3" and c["red"] == "Fe(OH)_2")
    fe_ion = next(c for c in T.couples
                  if c["ox"] == "Fe^{3+}" and c["red"] == "Fe^{2+}")
    k3 = T.ksp_by_solid["Fe(OH)_3"]["pKsp"]
    k2 = T.ksp_by_solid["Fe(OH)_2"]["pKsp"]
    check("Fe(OH)3/Fe(OH)2 oh 形 ↔ Ksp 网络",
          fe["oh_ref"], fe_ion["E0"] - K * (k3 - k2), 0.05)

    # 5) 氢氧化物 Ksp 与"阳离子水解 pKa"等价（单一数据源声明，打印备查）
    if verbose:
        for e in T.ksp:
            if e["pair"][1] == "OH^-":
                from .engine import _ksp_xy
                x, y = _ksp_xy(e)
                print(f"[info] {e['solid']}: 水解 logK = pKsp − n·pKw = "
                      f"{e['pKsp'] - y * PKW:.1f} "
                      f"（{x}{e['pair'][0]} + {y}H2O ⇌ {e['solid']} + {y}H+）")

    # 6) 氨合配离子溶解度梯度（Ksp⊗β 派生 logK：AgCl 微负、AgBr 更负、AgI 极负）
    if verbose:
        for solid in ("AgCl", "AgBr", "AgI"):
            e = T.ksp_by_solid[solid]
            b = T.beta_by_complex["[Ag(NH_3)_2]^+"]
            print(f"[info] {solid} + 2NH3 ⇌ 配离子 + 卤离子: "
                  f"logK = {b['logb'] - e['pKsp']:.2f}")

    # 7) 配合物电荷/原子一致性（beta schema 单一配体纪律）：
    #    complex 电荷 = center + nu×ligand；complex 原子组成 = center + nu×ligand。
    #    混配体配合物（如 [Co(NH3)5Cl]2+ 含非配体元素 Cl）不符合
    #    center+nu*ligand 生成路径，一律拒绝（引擎配位步将守恒被破坏）
    from .core import charge_of, elements_of
    for e in T.beta:
        q_ok = (charge_of(e["complex"])
                == charge_of(e["center"]) + e["nu"] * charge_of(e["ligand"]))
        a_c = elements_of(e["complex"])
        a_ref = dict(elements_of(e["center"]))
        for el, cnt in elements_of(e["ligand"]).items():
            a_ref[el] = a_ref.get(el, 0) + e["nu"] * cnt
        if not q_ok:
            fails.append(f"[FAIL] beta 电荷不守恒 {e['complex']}: "
                         f"{charge_of(e['complex'])} ≠ {charge_of(e['center'])} "
                         f"+ {e['nu']}×{charge_of(e['ligand'])}")
        elif a_c != a_ref:
            extra = set(a_c) ^ set(a_ref)
            fails.append(f"[FAIL] beta 原子不守恒 {e['complex']} "
                         f"(混配体/非单一配体条目: 差异元素 {sorted(extra)}; "
                         f"{a_c} vs {a_ref})")
        elif verbose:
            print(f"[ok] beta 守恒 {e['complex']} "
                  f"(q={charge_of(e['complex'])})")

    # 8) 同名配合物双注册检查（不同 center 的重名条目互相覆盖，产生
    #    complex/decomplex 幻影振荡，如 [SbCl6]- 曾同时注册 Sb3+/Sb5+ 中心）
    seen: dict = {}
    for e in T.beta:
        if e["complex"] in seen and seen[e["complex"]] != e["center"]:
            fails.append(f"[FAIL] beta 同名双中心 {e['complex']}: "
                         f"{seen[e['complex']]} vs {e['center']}")
        seen[e["complex"]] = e["center"]

    return fails


def dH_coverage(T: Tables) -> None:
    """dH 运行时派生覆盖率报告（原 build_dH.py 职责并入）。"""
    print("\n---- dH 覆盖率 ----")
    for name, tbl in (("couples", T.couples), ("pka", T.pka),
                      ("ksp", T.ksp), ("beta", T.beta)):
        hit = sum(1 for e in tbl if e.get("dH") is not None)
        print(f"  {name}: {hit}/{len(tbl)}")


def case_api() -> None:
    """System/Reaction/Engine 高层 API 自检：对象语义、累计投料再平衡、
    raw 过程保留、net_equation 净离子方程式、H+/OH-/H2O 显式出现、外界气压调节，
    Engine 对象三层级平级方法（v0.3.6），以及核心字段：initial/degree/
    consumption/production/reacted/changed。"""
    import chemkit
    global PASS_N
    sys = chemkit.System(V=1.0)
    sys.add("NaOH", 0.1)
    r2 = sys.add("HCl", 0.1)
    ok1 = (isinstance(r2, chemkit.Reaction) and r2.changed
           and r2.degree == 2 and r2.pH is not None
           and 6.0 <= r2.pH <= 8.0 and len(sys.history) == 2
           and sys.feeds.get("NaOH") == 0.1)
    r3 = chemkit.react({"Ca(OH)_2": 1.0, "CO_2": 0.5}, V=1.0)
    ok2 = (r3.production.get("CaCO_3", 0.0) >= 0.45
           and isinstance(r3.raw.get("steps"), list))
    ok3 = chemkit.System({"NaCl": 0.1}).result is not None
    # 净离子方程式：Zn + H2SO4 → Zn + 2H+ -> H2 + Zn2+
    r4 = chemkit.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0)
    # v0.5.0：净方程式**先结构化后渲染**——`net_equation` 是 Equation
    # （left/right 系数 + reversible 布尔），字符串只是视图：
    # `str()`/`.tex()` = TeX（进 Markdown $$ 块），`.plain()` = 无标记 fallback。
    _n4 = r4.net_equation
    ok4 = (isinstance(_n4, Equation)
           and _n4.left == {"H^+": 2, "Zn": 1}
           and _n4.right == {"H_2": 1, "Zn^{2+}": 1}
           and _n4.reversible is False
           and _n4.tex() == ("2\\mathrm{H}^{+} + \\mathrm{Zn} \\rightarrow "
                            "\\mathrm{H_{2}} + \\mathrm{Zn}^{2+}")
           and _n4.plain() == "2H+ + Zn -> H2 + Zn2+")
    # H+/OH-/H2O 显式出现在 consumption/production：
    # NaOH+HCl → consumption 含 H+ 和 OH-，production 含 H2O
    r5 = chemkit.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0)
    ok5 = ("H^+" in r5.consumption and "OH^-" in r5.consumption
           and "H_2O" in r5.production)
    # 外界气压调节：低压下气体更易逸出（不应报错）
    r6 = chemkit.react({"Na_2CO_3": 0.1, "HCl": 0.2}, V=1.0, p=50.0)
    ok6 = r6.changed and "CO_2" in r6.production
    # degree / consumption / production / net_equation
    ok7 = (r4.degree == 2  # Zn + H2SO4 完全反应
           and r4.consumption is not None
           and r4.production is not None
           and r4.net_equation is not None)
    # initial（初态，post-normalize）
    # SO3 + H2O → H+ + HSO4-；初态应含 HSO4- 与 H+（强电解质已电离 + 气体已反应）
    r8 = chemkit.react({"SO_3": 0.1}, V=1.0)
    ok8 = (r8.initial.get("HSO_4^-", 0.0) >= 0.09
           and r8.initial.get("H^+", 0.0) >= 0.09
           and "SO_3" not in r8.initial)
    # reacted（狭义化学反应）
    # NaCl 溶解：changed=False（NaCl 已在 normalize 阶段完全电离，账本无 NaCl 残留），
    # reacted=False（纯溶解不算化学反应）
    r9 = chemkit.react({"NaCl": 0.1}, V=1.0)
    ok9 = (r9.changed is False and r9.reacted is False
           and r9.degree == 0)
    # Zn + H2SO4：changed=True, reacted=True（redox）
    ok10 = (r4.reacted is True and r4.degree == 2)
    # 绝热耦合（isothermal=False）：中和热教材值 55.84 kJ/mol 量热式验证
    rt = chemkit.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0, isothermal=False)
    ok11 = (rt.heat_kJ is not None and 5.0 <= rt.heat_kJ <= 6.2
            and rt.dT_K is not None and 0.8 <= rt.dT_K <= 1.8
            and 298.9 <= rt.T_final_K <= 299.9
            and rt.thermal.get("converged") is True
            and len(rt.thermal.get("trace", [])) >= 1
            and rt.degree == 2)
    # 默认恒温：单遍求解、不计热效应（判定引擎主用途）
    ok12 = (r5.thermal.get("mode") == "isothermal" and r5.heat_kJ is None
            and r5.dT_K is None)
    # 绝热耦合收敛（大热效应多轮）+ 化学结果仍在自洽终温下正确
    rz = chemkit.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0, isothermal=False)
    # 【第 147 轮】此检查一度翻红：`thermo.analyze` 的"缺 ΔHf"判别只看
    # ≥TRACE_MOL 的物种，而求和遍历全体 ⟹ 痕迹量且 ΔHf 未知的新物种
    # （[Zn(OH)]⁺）触发 `float * None`，整条绝热耦合降级成 heat_kJ=None。
    # 真因是引擎缺陷（判别与求和集合不一致），已在 thermo.py 根治 ⟹ 本检查
    # 恢复严格形式（热值必须落在教材区间）。
    ok13 = (rz.heat_kJ is not None and 140.0 <= rz.heat_kJ <= 165.0
            and rz.T_final_K is not None and 325.0 <= rz.T_final_K <= 345.0
            and rz.thermal.get("converged") is True
            and rz.degree == 2 and "H_2" in rz.production)
    # v0.3.5 热层新覆盖：thermo.json 统一后新增固体的溶解焓路径
    # BaF2 沉淀放热（溶解焓 +4.2 kJ/mol → 沉淀放热）
    rbf = chemkit.react({"BaCl_2": 0.1, "NaF": 0.2}, V=0.1, isothermal=False)
    ok14 = (rbf.heat_kJ is not None and 0.25 <= rbf.heat_kJ <= 0.55
            and rbf.dT_K is not None and 0.5 <= rbf.dT_K <= 1.5
            and rbf.thermal.get("converged") is True
            and rbf.production.get("BaF_2", 0.0) >= 0.09)
    # ZnCO3 沉淀吸热（溶解焓 -18.2 → 沉淀吸热，dT 为负——罕见但真实的
    # 吸热沉淀路径，热层符号正确性检验）
    rzc = chemkit.react({"ZnCl_2": 0.1, "Na_2CO_3": 0.1}, V=0.1, isothermal=False)
    ok15 = (rzc.heat_kJ is not None and -2.2 <= rzc.heat_kJ <= -1.4
            and rzc.dT_K is not None and -5.5 <= rzc.dT_K <= -3.5
            and rzc.thermal.get("converged") is True)
    # CaCO3+HCl 气液两相拆分：气相 CO2（gas）按 ΔHf(g) 计价、残留溶解态
    # 按水溶值（比 v0.3.4 全气相口径更准 0.07 kJ 量级）
    rcc = chemkit.react({"CaCO_3": 0.05, "HCl": 0.1}, V=0.1, isothermal=False)
    ok16 = (rcc.heat_kJ is not None and 0.6 <= rcc.heat_kJ <= 1.0
            and rcc.gas.get("CO_2", 0.0) >= 0.04
            and rcc.thermal.get("converged") is True
            and rcc.dT_K is not None and 1.5 <= rcc.dT_K <= 2.5)
    # v0.3.6 Engine 对象：三个层级是平级方法，与函数式入口同语义
    eng = chemkit.Engine()
    re1 = eng.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0)
    ok17 = (isinstance(re1, chemkit.Reaction) and re1.degree == 2
            and re1.net_equation == r5.net_equation
            and re1.thermal.get("mode") == "isothermal")
    s2 = eng.system(V=1.0)
    s2.add("NaOH", 0.1)
    re2 = s2.add("HCl", 0.1)
    ok18 = (s2._tables is eng.tables and len(s2.history) == 2
            and re2.degree == 2 and 6.0 <= re2.pH <= 8.0)
    re3 = eng.react({"NaOH": 0.1, "HCl": 0.1}, V=1.0, isothermal=False)
    ok19 = (re3.heat_kJ is not None and 5.0 <= re3.heat_kJ <= 6.2
            and re3.thermal.get("converged") is True)
    # v0.5.2：用户侧删去与 react 重叠的 judge()，raw 账本改从 Reaction.raw 取
    # ——通道测的是"引擎 raw dict 契约仍可达"，入口换成唯一入口。
    raw4 = eng.react({"Zn": 1.0, "H_2SO_4": 1.0}, V=1.0).raw
    ok20 = (isinstance(raw4, dict) and raw4["reacted"] is True
            and raw4["degree"] == 2 and raw4["override"] is None
            and {"steps", "H_excess", "cond"} <= set(raw4)
            and raw4["cond"]["V_L"] == 1.0)
    # v0.5.2：judge 的独有能力（无阴离子的质子条件量 c_H/c_OH/pH）不能在
    # "删重复入口"时被顺手削掉——它必须从存活入口 react/System 可达。
    rc_h = chemkit.react({"Zn": 0.1}, V=1.0, c_H=0.5)
    rs_h = chemkit.System(V=1.0, c_H=0.5)
    rs_h.add("Zn", 0.1)
    ok21 = (rc_h.pH == rs_h.result.pH and rc_h.degree == rs_h.result.degree
            and rc_h.raw["H_excess"] > 0.0)
    ok22 = not hasattr(chemkit, "judge") and not hasattr(eng, "judge")
    for tag, ok in (("API System 累计投料再平衡", ok1),
                    ("API react 一步式+raw 过程", ok2),
                    ("API System 建立即反应", ok3),
                    ("API Reaction.net_equation 净离子方程式", ok4),
                    ("API H+/OH-/H2O 显式出现", ok5),
                    ("API 外界气压 p 调节", ok6),
                    ("API degree/consumption/production/net_equation", ok7),
                    ("API initial 初态含 HSO4-（SO3+H2O 反应）", ok8),
                    ("API reacted 纯溶解=False", ok9),
                    ("API reacted redox=True", ok10),
                    ("API 绝热耦合中和热/dT/trace", ok11),
                    ("API 默认恒温 isothermal", ok12),
                    ("API 绝热耦合大热效应收敛", ok13),
                    ("API 热层新覆盖 BaF2 沉淀放热", ok14),
                    ("API 热层新覆盖 ZnCO3 沉淀吸热", ok15),
                    ("API 热层气相拆分 CaCO3+HCl", ok16),
                    ("API Engine.react 与函数式同语义", ok17),
                    ("API Engine.system 共享表与缓存", ok18),
                    ("API Engine 绝热耦合", ok19),
                    ("API Engine.react().raw 账本直通", ok20),
                    ("API c_H/c_OH/pH 质子条件量可达（react/System）", ok21),
                    ("API 用户侧无 judge 重复入口", ok22)):
        if ok:
            PASS_N += 1
        else:
            FAILS.append(tag)


# 性能分界线（ms）——README「性能」节的 5 档固定口径，改动需同步文档
PERF_EDGES = (10.0, 20.0, 50.0, 100.0, 500.0)


def _perf_buckets(ms: list[float]) -> dict:
    """5 档性能分布（口径与 README 性能表一致）+ 分位数。"""
    n = max(len(ms), 1)
    srt = sorted(ms)
    return {
        "edges_ms": list(PERF_EDGES),
        "counts": {f">{e:g}ms": sum(1 for v in ms if v > e) for e in PERF_EDGES},
        "pct": {f">{e:g}ms": round(sum(1 for v in ms if v > e) / n * 100, 1)
                for e in PERF_EDGES},
        "mean_ms": round(sum(ms) / n, 2),
        "p50_ms": round(srt[n // 2], 2),
        "p90_ms": round(srt[int(n * 0.9)], 2),
        "max_ms": round(srt[-1], 2),
        "total_ms": round(sum(ms), 1),
    }


def write_report(path: str, extra: dict | None = None) -> str:
    """把本轮测试结果**结构化落盘**（JSON），供保留与离线读取比对。

    结构：
      meta     版本/时间/规模/总墙钟
      summary  pass/fail + 5 档性能分布（含 >100ms 档的逐例编号）
      cases    逐例：index/name/ok/ms/errors/note/pH/degree/changed/
               annotations/net_equation
      checks   附加检查（T 区间/API/环闭合）结果
    """
    ms = [t * 1000.0 for t, _ in TIMES]
    perf = _perf_buckets(ms)
    over100 = [{"index": i + 1, "name": r["name"], "ms": r["ms"]}
               for i, r in enumerate(RESULTS) if r["ms"] > PERF_EDGES[3]]
    over100.sort(key=lambda x: -x["ms"])
    # 质量口径（第 214 轮）：并行路径也把探针带来的收敛残差写进留档，
    # 于是"通过性 + 耗时 + 残差 + 归因"一份 JSON 齐备 ⟹ 不必再跑来换信息。
    _enrich_results_with_probe(RESULTS, _PROBES)
    doc = {
        "meta": {
            "at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "version": _version(),
            "cases": len(RESULTS),
            "elapsed_s": round(sum(t for t, _ in TIMES), 2),
        },
        "summary": {
            "pass": PASS_N, "fail": len(FAILS),
            "total": PASS_N + len(FAILS),
            # 分开报（同一份留档里两种口径必须一致，见本文件顶部说明）
            "cases_pass": sum(1 for r in RESULTS if r["ok"]),
            "cases_fail": sum(1 for r in RESULTS if not r["ok"]),
            "checks_pass": sum(c["passed"] for c in CHECKS),
            "checks_fail": sum(len(c["fails"]) for c in CHECKS),
            "failed_names": list(FAILS),
            "perf": perf,
            "over_100ms": over100,
        },
        "cases": RESULTS,
        "checks": extra or {},
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
        f.write("\n")
    return path


def _version() -> str:
    try:
        from . import __version__
        return __version__
    except Exception:            # noqa: BLE001
        return "?"


def _perf_report() -> None:
    """性能分档统计（5 档）+ >100ms 档的用例编号。

    口径与 README 性能表一致（>10/20/50/100/500 ms），便于跨版本对照。
    **>100 ms 档给出编号**：这一档是个位数到几十例的量级，逐例可追踪
    （回归时点名即可定位）；>10/20/50 档只给计数与占比（量大，逐例列出
    无信息量）。编号 = 用例在全量清单中的序号（1 起）+ 用例名。
    """
    if not TIMES:
        return
    ms = [t * 1000.0 for t, _ in TIMES]
    n = len(ms)
    print("\n---- 性能分界（5 档，单例墙钟）----")
    print(f"  {'分界线':>9s} {'用例数':>7s} {'占比':>8s}")
    for e in PERF_EDGES:
        k = sum(1 for v in ms if v > e)
        print(f"  {'> ' + format(e, 'g') + ' ms':>9s} {k:7d} {k / n * 100:7.1f}%")
    srt = sorted(ms)
    print(f"  均值 {sum(ms) / n:.1f} ms；中位 {srt[n // 2]:.1f} ms；"
          f"P90 {srt[int(n * 0.9)]:.1f} ms；最值 {srt[-1]:.1f} ms")
    over = [(v, i + 1, nm) for i, ((_t, nm), v) in enumerate(zip(TIMES, ms))
            if v > PERF_EDGES[3]]
    if over:
        over.sort(reverse=True)
        print(f"  >100 ms 档（{len(over)} 例，编号＝清单序号）：")
        for v, idx, nm in over:
            print(f"    #{idx:<5d} {v:9.1f} ms  {nm}")


def _par_init():
    """并行 worker 启动：自己载表 + **预热**。

    `spawn` 下 worker 不继承父进程状态 ⟹ 必须自己载表（实测 182 ms）。
    **预热不可省**（第 206 轮实测）：否则每个 worker 的首个用例要独自承担
    全部惰性缓存构建，使套件最靠前的简单用例（`1 HCl+NaOH` 等）显示成
    12~24 s —— 看着像算法病态，其实同一用例预热后只要 9~270 ms。
    `main()` 本就有预热，并行路径必须一致，否则"最慢用例"榜被冷启动污染。
    """
    global _PAR_T
    _T = load_tables()
    judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, _T)
    _PAR_T = _T


_PAR_T = None
# 并行路径收集的探针（收敛画像），供 write_report 落盘 `resid_live` 等
# 质量口径字段——否则并行留档只有通过性/耗时，残差还得另跑 converg（第 214 轮）。
_PROBES: dict = {}
# `run_case` 每次调用会把本次 judge 的探针存到这里（`_probe` 传入），
# 并行 worker 取走后再回传父进程。**避免了并行路径重复调用 judge**
# （第 214 轮踩到：为拿 probe 又调一次 judge，墙钟直接翻倍到 108s）。
_LAST_PROBE: dict = {}


def _par_run(case: dict):
    """并行 worker：跑单个用例，回传 (name, ok, ms, errors, result_dict, probe)。

    与 `run_case` 共用同一份判定逻辑（**不复制**）：清空 worker 内的全局
    残余，调用 `run_case` 取本次追加的失败项与 RESULTS 条目，再把全局复原。

    **第 214 轮：同时回传 `probe`**（`judge(..., _probe=...)` 的收敛画像）。
    此前并行路径丢弃 probe ⟹ `resid_live`（质量口径）拿不到，只能另跑
    `converg` 串行全量。现在一份留档同时含：通过性、耗时、**与收敛残差**，
    满足"跑一次、读多次"。probe 只含 float/list/dict，跨进程可 pickle。
    """
    global PASS_N
    name = case["name"]
    n_f0, n_r0, n_p0 = len(FAILS), len(RESULTS), PASS_N
    t0 = time.time()
    try:
        ok = run_case(case, _PAR_T, verbose=False)
        errs = list(FAILS[n_f0:])
    except Exception as exc:                                # noqa: BLE001
        ok = False
        FAILS.append(name)
        RESULTS.append({"index": 0, "name": name, "ok": False,
                        "ms": round((time.time() - t0) * 1000, 2),
                        "errors": [f"运行异常 {type(exc).__name__}: {exc}"],
                        "note": case.get("note") or "", "pH": None,
                        "degree": None, "changed": None,
                        "annotations": [], "net_equation": None})
    ms = round((time.time() - t0) * 1000, 2)
    errs = list(FAILS[n_f0:])
    rec = dict(RESULTS[n_r0]) if len(RESULTS) > n_r0 else None
    if rec is not None:
        rec["ms"] = ms
    del FAILS[n_f0:]
    del RESULTS[n_r0:]
    PASS_N = n_p0
    return (name, bool(ok), ms, errs, rec, _LAST_PROBE)


def run_cases_parallel(cases: list[dict], jobs: int,
                       progress: bool = True) -> tuple[int, float]:
    """并行跑用例，把结果**按原顺序**填入 RESULTS/TIMES/FAILS/PASS_N。

    返回 (通过数, 墙钟秒)。用例之间互相独立（各跑自己的 judge），
    故并行安全；已实测**逐例结果（通过性 + 错误文本）与串行完全一致**。

    `progress=True` 时**边完成边打进度**（每例一行 `[ok]/[FAIL]`，含累计
    与预计剩余时间）—— 第 296 轮加：长跑时不再让人空等。
    进度只写 **stderr**，保证 stdout 的摘要格式不被污染。
    """
    global PASS_N
    from concurrent.futures import ProcessPoolExecutor
    if jobs <= 1:
        raise ValueError("jobs<=1 应走串行分支")
    # 长尾负载均衡：**先跑慢例**（静态排序表；没有表就按原顺序）。
    # 实测最慢单例 ~20 s，若落到最后启动，墙钟直接被它拖住。
    order = list(range(len(cases)))
    slow = _SLOW_FIRST
    if slow:
        rank = {n: i for i, n in enumerate(slow)}
        order.sort(key=lambda i: rank.get(cases[i]["name"], len(rank)))
    t0 = time.perf_counter()
    got: dict[int, tuple] = {}
    n_done = n_ok_so_far = 0
    total = len(cases)
    _next_pct = 10               # **每 10% 报一次**（第 296 轮：不再让人空等）
    with ProcessPoolExecutor(max_workers=jobs,
                             initializer=_par_init) as ex:
        futs = {ex.submit(_par_run, cases[i]): i for i in order}
        from concurrent.futures import as_completed
        for fu in as_completed(futs):
            got[futs[fu]] = fu.result()
            n_done += 1
            _nm, _ok, _ms, _errs, _rec, _pr = got[futs[fu]]
            n_ok_so_far += 1 if _ok else 0
            pct = n_done * 100 // total
            # **每跨过 10% 打一行进度**；失败例另外单打一行 ✗（不带计数，
            # 免得看起来像多报了几次进度）。
            if pct >= _next_pct:
                el = time.perf_counter() - t0
                eta = (el / n_done) * (total - n_done) if n_done else 0.0
                print(f"  [{pct:3d}%] {n_done:>4}/{total}  通过 {n_ok_so_far}"
                      f"  {el:6.1f}s  剩~{eta:4.0f}s", file=sys.stderr, flush=True)
                _next_pct = (pct // 10 + 1) * 10
            if not _ok:
                print(f"    ✗ {_nm[:56]} -- {'; '.join(_errs or [])[:78]}",
                      file=sys.stderr, flush=True)
    wall = time.perf_counter() - t0
    PASS_N = 0
    global _PROBES
    _PROBES = {}
    for i, c in enumerate(cases):
        _nm, ok, ms, errs, rec, probe = got[i]
        if probe:
            _PROBES[c["name"]] = probe
        if rec is None:
            rec = {"index": i + 1, "name": c["name"], "ok": bool(ok),
                   "ms": ms, "errors": errs, "note": c.get("note") or "",
                   "pH": None, "degree": None, "changed": None,
                   "annotations": [], "net_equation": None}
        rec["index"] = i + 1
        # **分片要用全局序号**（第 296 轮）：分片跑时 `i` 是子集内的局部序号，
        # 各片都从 1 开始 ⟹ 合并时大量 index 重复、无法归位。用例对象里带
        # `__gidx`（由 dev.py 分片时写入）时以它为准。
        if c.get("__gidx") is not None:
            rec["index"] = c["__gidx"]
        RESULTS.append(rec)
        TIMES.append((ms / 1000.0, c["name"]))
        if ok:
            PASS_N += 1
        else:
            FAILS.append(c["name"])
    return PASS_N, wall


# 已知慢例（仅用于并行**调度**，不影响判定）。数值取自
# `logs/suite-parallel-latest.json` 的逐例留档；表过期只会让负载均衡变差，
# 不会错判。
#
# ⚠️ 第 241 轮实测：**同一用例的单例耗时会随 worker 数大幅波动**
# （`Y05` 在 4 workers 下 14.7 s、6 workers 18.9 s、8 workers 20.2 s）——
# 因为超订物理核后单例内部就在跟别的 worker 抢 CPU。故任何形如
# "Y05 耗时 20 s 是算法病态"的结论，**必须先说明当时的 worker 数**；
# 本表只当作**排序提示**，不当作性能结论。
_SLOW_FIRST = (
    "Y05 BaCl2+Na2SO3 白沉", "Z07 MnS+醋酸 溶解",
    "N19 Na[Al(OH)4]+CO2过量", "X10 Cu+Hg(NO3)2 置换",
    "M01 AlCl3+少量NaOH", "Amp14 AlCl3+少量NaOH",
    "Z04 Ag+稀硝酸 放NO", "MX01 Fe3+ + SCN- 显色平衡",
    "Cu31 CuSO4+NaCl+Cu 还原", "DR6 MgCl2+Na2CO3",
    "H42 AlCl3+NaOH 1:2.5", "Co41 CoCl2+NaHCO3",
)


def _enrich_results_with_probe(results: list, probes: dict) -> None:
    """把并行收集的探针（收敛画像）并进逐例结果，**原地**改 `results`。

    产出字段（质量口径，与 `converg._live` 同一函数，口径只此一处）：
      `resid_live`（值得解且 walk 会解的两侧平衡最大 |S|）、`resid_max`、
      `iters`、`exit`、`pH_solver`、`resid_src_eq/kind`（残差归因）。
    这样一份 `logs/suite-latest.json` 同时含通过性、耗时、残差与归因 ⟹
    **不必为看残差另跑 converg 串行全量**（第 214 轮动因）。
    """
    if not probes:
        return
    try:
        from .converg import _live
    except Exception:                                       # noqa: BLE001
        return
    for r in results:
        p = probes.get(r["name"])
        if not p:
            continue
        try:
            act = p.get("active") or []
            src = p.get("resid_src") or {}
            r["resid_live"] = _live(act)
            r["resid_max"] = p.get("max_abs_S")
            r["iters"] = p.get("iters")
            r["exit"] = p.get("exit")
            r["pH_solver"] = p.get("pH_solver")
            r["resid_src_eq"] = src.get("eq")
            r["resid_src_kind"] = src.get("kind")
        except Exception:                                   # noqa: BLE001
            continue


def main(cases_path: str | None = None, out_path: str | None = None,
         jobs: int = 1, cases: list[dict] | None = None) -> int:
    """跑套件。

    `cases` 直接给用例列表（**分片用**）：走 `cases_path` 时用例要先序列化到
    临时 JSON，而 `json.dump` 会丢掉引擎侧加的 `__gidx` 之类非 JSON 原生键的
    语义（读回来是普通键，且分片各片局部编号）⟹ 合并时 index 重复。直接传对象
    可保住全局序号（第 296 轮）。
    """
    T = load_tables()
    judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, T)   # 预热：模板/缓存冷启动不计入首例
    TIMES.clear()
    RESULTS.clear()
    FAILS.clear()
    _cases = cases if cases is not None else load_cases(cases_path)
    if jobs and jobs > 1:
        global _PAR_T
        _PAR_T = T                       # 串行分支复用
        n_ok, wall = run_cases_parallel(_cases, jobs)
        print(f"[并行 {jobs} workers] 墙钟 {wall:.1f}s")
    else:
        # 串行也报进度（**每 10% 一档**，失败例立即报）——第 296 轮：分片在
        # 外层并行时，片内走的就是这条路径，没有进度就完全看不到动静。
        _t0 = time.perf_counter()
        _total = len(_cases)
        _ok_n = 0
        _next = 10
        for _i, c in enumerate(_cases, 1):
            run_case(c, T)
            if RESULTS and RESULTS[-1].get("ok"):
                _ok_n += 1
            _pct = _i * 100 // _total if _total else 100
            if _pct >= _next:
                _el = time.perf_counter() - _t0
                _eta = (_el / _i) * (_total - _i) if _i else 0.0
                print(f"  [{_pct:3d}%] {_i:>4}/{_total}  通过 {_ok_n}"
                      f"  {_el:6.1f}s  剩~{_eta:4.0f}s", file=sys.stderr, flush=True)
                _next = (_pct // 10 + 1) * 10
            if RESULTS and not RESULTS[-1].get("ok"):
                print(f"    ✗ {c['name'][:56]} -- "
                      f"{'; '.join(RESULTS[-1].get('errors') or [])[:78]}",
                      file=sys.stderr, flush=True)
    # ===== 质量口径富化（第 214 轮）：把探针带来的收敛残差并进逐例结果 =====
    # 放在 `main` 里（而不是只放 `write_report`）——`dev.py suite` 调用的是
    # `main(tmp, None)`，`out_path=None` 时 `write_report` 根本不执行，
    # 而 `dev.py` 直接 dump `RESULTS` ⟹ 只写 write_report 会**静默丢失**。
    _enrich_results_with_probe(RESULTS, _PROBES)
    CHECKS.clear()
    for _nm, _fn in (("T_range", lambda: case_T_range(T)),
                     ("acidbase", lambda: acidbase_check(T)),
                     ("api", lambda: case_api())):
        _f0, _p0 = len(FAILS), PASS_N
        _fn()
        _new = list(FAILS[_f0:])
        CHECKS.append({"name": _nm, "ok": not _new,
                       "passed": PASS_N - _p0 - len(_new), "fails": _new})
    n_case = len(RESULTS)
    n_case_ok = sum(1 for r in RESULTS if r["ok"])
    n_chk = sum(c["passed"] + len(c["fails"]) for c in CHECKS)
    n_chk_ok = sum(c["passed"] for c in CHECKS)
    total = PASS_N + len(FAILS)
    print(f"\n===== 用例 {n_case_ok}/{n_case} · 辅助检查 {n_chk_ok}/{n_chk}"
          f" · 合计 {PASS_N}/{total} PASS =====")
    if FAILS:
        print("失败:", FAILS)

    slow = sorted((t, n) for t, n in TIMES if t > 0.5)
    if slow:
        print(f"\n---- 慢用例（>{0.5}s，共 {len(slow)} 例）----")
        for t, n in reversed(slow):
            print(f"  {t:6.2f}s {n}")
    print(f"总耗时 {sum(t for t, _ in TIMES):.1f}s")
    _perf_report()
    dH_coverage(T)

    print("\n---- 环闭合检查 ----")
    cfail = consistency(T)
    if cfail:
        print("\n".join(cfail))
    else:
        print("全部环闭合检查通过。")
    if out_path:
        p = write_report(out_path, extra={
            "batteries": CHECKS, "consistency_fail": cfail,
        })
        print(f"\n结构化结果已写入 {p}"
              f"（cases/summary/checks；summary.perf 为 5 档性能分布）")
    return 0 if (not FAILS and not cfail) else 1


if __name__ == "__main__":
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:                                # pragma: no cover
            pass
    args = sys.argv[1:]
    if "-h" in args or "--help" in args:
        # v0.5.2：此前 `--help` 会被当成用例库路径去 open() ⟹ 抛栈退出；
        # 文档里写着这个入口，就得让它打印用法而不是崩。
        print("用法: python -m chemkit.testsuit [用例库.json] [--out 结果.json]\n"
              "  省略用例库 = 包内 data/tests.json + 温度域/酸碱地基/高层 API\n"
              "  自检 + 环闭合检查；--out 另存结构化结果（cases/summary/checks）。")
        sys.exit(0)
    out = None
    if "--out" in args:
        i = args.index("--out")
        if i + 1 >= len(args):
            print("!! --out 需要一个文件名参数")
            sys.exit(2)
        out = args[i + 1]
        args = args[:i] + args[i + 2:]
    sys.exit(main(args[0] if args else None, out))
