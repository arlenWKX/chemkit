"""**不变量审计**：抓"引擎在什么情况下于物理上不可能对"的三类错误。

第 304 轮建立。动机来自第 303 轮的全面审查（`_audit/SUMMARY.md`）：
套件 1283/1378 这个数字里混着三类**结构性假阳性**，而当时的检查一条都抓不到。

三条不变量（都是**可以在不看真值的前提下**判定对错的）：

  ① **尺度不变性**（`scale`）—— 同一浓度、同一配比的溶液，
     改体积不能改变任何**强度量**（pH、各物种浓度、degree）。
     体积是热力学的无关变量。反例：`GaCl3+4NaOH` 在 V=0.1…10 之间
     pH 漂 6.1 个单位（`acidbase.dist_charge` 量纲错）。

  ② **进程确定性**（`det`）—— 同一份代码、同一份数据、同一输入，
     在不同 `PYTHONHASHSEED` 下必须给出**逐位相同**的结果。
     反例：`normalize._split_salt` 遍历未排序的 `T.anions`，
     `H_3PO_4` 在 seed=0 拆成 `H+ + H2PO4-`、seed=1 拆成 `3H+ + PO4^3-`。

  ③ **输出守恒**（`cons`）—— 引擎**自己给出**的净方程式 / 分步方程式
     必须元素与电荷守恒。注意这与 `eqcheck` 不同：`eqcheck` 校验的是
     `tests.json` 里的**期望式**，本工具校验的是**引擎的实际输出**。
     反例：`D43 Cl2+石灰乳` 输出 `Ca(OH)2 + Cl2 -> ClO- + Ca2+`，Cl 差 1。

  ④ **绝热放热恒定**（`heat`）—— 固定投料摩尔数、改体积 ⟹ 绝热放热
     `heat_kJ` 必须恒定（广延量）。**锁封闭体系**（`gas_escape=False`）：
     开放体系的逸出焓差是真实物理（CO₂ 脱气 20.29 kJ/mol），逸出量随
     泡点（浓度）变，总放热随之而变是正确的；封闭口径把该项归零，
     剩下的漂移就是账本缺陷。
     反例（第 306 轮实测，已修）：`CaCO3 0.05 + HCl 0.1` 封闭体系
     1.196→1.672 kJ（V=0.1→2.0），根因是缺 ΔHf 的 `[CaCl]⁺` 被静默丢弃。

用法：
    python -m chemkit.helper.invariants              # 四条全跑（用 tests.json）
    python -m chemkit.helper.invariants scale        # 只跑尺度
    python -m chemkit.helper.invariants det --n 40   # 只跑确定性，抽 40 例
    python -m chemkit.helper.invariants cons         # 只跑守恒
    python -m chemkit.helper.invariants heat         # 只跑放热恒定
    python -m chemkit.helper.invariants cons --jobs 8 # 并行（默认 min(4, cpu)）

退出码：0 = 全部通过；1 = 有违规（供 CI / 发版门槛使用）。

**为什么这几条值得单独立工具**：它们检查的是"引擎能不能自相矛盾"，
而不是"引擎与某个外部真值差多少"。后者需要权威数据、会随数据更新而变；
前者是逻辑同一性，永不过期，而且一旦违反，无论真值是什么都说明有 bug。
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:                                    # pragma: no cover
    pass

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from chemkit.core import charge_of, elements_of       # noqa: E402

# ---------------------------------------------------------------- 公共

_TOL_PH = 0.02       # pH 分辨率：引擎 pH 保留 3 位小数，0.02 足够宽松
_TOL_REL = 0.05      # 浓度相对容差（5%）——吸收数值路径的合理差异
# 痕量门的**两层**设计（第 304 轮实测教训）：
#   · `_TRACE`（绝对）= 1e-6 M：低于此浓度不参与比较。
#   · `_TRACE_REL`（相对）= 1e-3：相对主量低于 0.1% 的物种也不比较。
_TRACE = 1e-6
_TRACE_REL = 1e-3


def load_cases() -> list[dict]:
    """读包内 tests.json（与 testsuit 同一份真值源）。"""
    p = os.path.join(_ROOT, "chemkit", "data", "tests.json")
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def _subs_of(case: dict) -> dict:
    raw = case.get("subs") or []
    return dict(raw) if isinstance(raw, dict) else {k: v for k, v in raw}


def _react(case: dict, V: float | None = None, T=None,
           gas_escape: bool | None = None):
    """按用例跑一次；`V` 给定时覆盖用例自带的体积。

    `gas_escape=False` 用于封闭体系口径（④ 放热恒定检查：隔离逸出焓差
    这一真实物理效应，只留账本缺陷）。
    """
    import chemkit
    subs = _subs_of(case)
    cond = case.get("cond") or {}
    v = V if V is not None else cond.get("V_L", 1.0)
    tk = T if T is not None else cond.get("T_K", 298.15)
    kw = {"V": v, "T": tk}
    if gas_escape is not None:
        kw["gas_escape"] = gas_escape
    return chemkit.react(subs, **kw)


def _conc(final: dict, V: float) -> dict:
    """呈现层 `final`（mol）→ 浓度（mol/L），滤掉痕量。"""
    return {s: m / V for s, m in (final or {}).items() if m / V > _TRACE}


def _jobs() -> int:
    """并行 worker 数：默认 min(4, cpu)；`--jobs N` 或 `CHEMKIT_JOBS` 覆盖。

    与 testsuit 同一范式（`ProcessPoolExecutor`）。scale/cons 的每条用例
    互相独立（各跑自己的 react），并行安全；已实测逐例结果与串行一致。
    """
    import argparse
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--jobs", type=int, default=0)
    try:
        ns, _ = p.parse_known_args()
        if ns.jobs > 0:
            return ns.jobs
    except SystemExit:                           # pragma: no cover
        pass
    env = os.environ.get("CHEMKIT_JOBS")
    if env:
        try:
            if int(env) > 0:
                return int(env)
        except ValueError:
            pass
    return max(1, min(4, os.cpu_count() or 1))


# ------------------------------------------------- ① 尺度不变性

def check_scale(cases: list[dict], factors=(0.1, 1.0, 10.0),
                limit: int | None = None) -> list[dict]:
    """同一**浓度**体系，体积放大/缩小 ⟹ 强度量必须不变。

    做法：把用例的投料量与体积**同比缩放**（保持浓度不变），比较
    pH、degree、changed 与各物种**浓度**。总摩尔数（广延量）理应随
    体积线性变化，故不比较 mol，只比较 mol/L。

    注意：气相也要比——泡点是**浓度**界（`c_sat = H·p`），同浓度就该同泡点，
    与总体积无关。这正是它能抓住 `_equil_gas_phase` 里
    `cap`（mol）与浓度混用那类错误的原因。
    """
    bad: list[dict] = []
    use = cases[:limit] if limit else cases
    jobs = _jobs()
    if jobs > 1 and len(use) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for row in ex.map(_scale_one, use, [factors] * len(use)):
                bad.extend(row)
        return bad
    for c in use:
        bad.extend(_scale_one(c, factors))
    return bad


def _scale_one(case: dict, factors=(0.1, 1.0, 10.0)) -> list[dict]:
    """单例尺度检查（worker 入口；`factors` 必须显式传入——模块级函数
    看不到 `check_scale` 的局部变量，Windows spawn 下闭包也不可用）。"""
    bad: list[dict] = []
    subs0 = _subs_of(case)
    if not subs0:
        return bad
    cond = case.get("cond") or {}
    V0 = cond.get("V_L", 1.0)
    T0 = cond.get("T_K", 298.15)
    base = None
    for f in factors:
        subs = {s: m * f for s, m in subs0.items()}
        try:
            r = _react({"subs": list(subs.items()),
                        "cond": {"V_L": V0 * f, "T_K": T0}})
        except Exception as e:                   # 引擎崩了也是缺陷
            bad.append({"case": case["name"], "factor": f,
                        "kind": "异常", "detail": f"{type(e).__name__}: {e}"})
            continue
        cur = {"pH": r.pH, "degree": r.degree, "changed": r.changed,
               "conc": _conc(r.final, V0 * f)}
        if base is None:
            base = (cur, f)
            continue
        b, bf = base
        if (cur["pH"] is None) != (b["pH"] is None):
            bad.append({"case": case["name"], "factor": f, "kind": "pH 有无",
                        "detail": f"f={bf} → {b['pH']}, f={f} → {cur['pH']}"})
        elif cur["pH"] is not None and abs(cur["pH"] - b["pH"]) > _TOL_PH:
            bad.append({"case": case["name"], "factor": f, "kind": "pH 漂移",
                        "detail": f"f={bf} → {b['pH']:.3f}, "
                                  f"f={f} → {cur['pH']:.3f}"})
        if cur["degree"] != b["degree"] or cur["changed"] != b["changed"]:
            bad.append({"case": case["name"], "factor": f, "kind": "判定翻转",
                        "detail": f"degree {b['degree']}→{cur['degree']}, "
                                  f"changed {b['changed']}→{cur['changed']}"})
        sp = set(cur["conc"]) | set(b["conc"])
        _scale = max([abs(b["conc"].get(s, 0.0)) for s in sp]
                     + [abs(cur["conc"].get(s, 0.0)) for s in sp] + [1e-12])
        diff = []
        for s in sorted(sp):
            x, y = b["conc"].get(s, 0.0), cur["conc"].get(s, 0.0)
            if max(abs(x), abs(y)) < _TRACE_REL * _scale:
                continue
            if abs(x - y) > _TOL_REL * max(abs(x), abs(y), 1e-12):
                diff.append(f"{s}: {x:.4g}→{y:.4g}")
        if diff:
            bad.append({"case": case["name"], "factor": f, "kind": "浓度漂移",
                        "detail": "; ".join(diff[:5])})
    return bad


# ------------------------------------------------- ② 进程确定性

_DET_PROBE = r"""
import json, sys
sys.path.insert(0, {root!r})
import chemkit
from chemkit.converg import _result_digest
out = []
for c in json.loads({payload!r}):
    try:
        r = chemkit.react(dict(c["subs"]), V=c.get("V", 1.0), T=c.get("T", 298.15))
        out.append({{"name": c["name"], "digest": _result_digest(r.raw),
                     "pH": r.pH, "degree": r.degree}})
    except Exception as e:
        out.append({{"name": c["name"], "err": f"{{type(e).__name__}}: {{e}}"}})
print(json.dumps(out))
"""


def check_det(cases: list[dict], seeds=("0", "1", "2"), limit: int | None = None,
              digest_only: bool = False) -> list[dict]:
    """同一输入在多个 `PYTHONHASHSEED` 下必须给出同一个结果。

    **必须分进程跑**：`PYTHONHASHSEED` 在解释器启动时冻结，同进程改无效。
    """
    use = cases[:limit] if limit else cases
    payload = []
    for c in use:
        cond = c.get("cond") or {}
        subs = _subs_of(c)
        payload.append({"name": c["name"],
                        "subs": [list(kv) for kv in subs.items()],
                        "V": cond.get("V_L", 1.0), "T": cond.get("T_K", 298.15)})
    src = _DET_PROBE.format(root=_ROOT, payload=json.dumps(payload))
    runs: dict = {}
    for sd in seeds:
        env = dict(os.environ)
        env["PYTHONHASHSEED"] = sd
        p = subprocess.run([sys.executable, "-X", "utf8", "-c", src],
                           capture_output=True, text=True, encoding="utf-8",
                           env=env, cwd=_ROOT)
        if p.returncode != 0:
            return [{"case": "<启动失败>", "kind": "子进程错误",
                     "detail": (p.stderr or "")[-400:]}]
        runs[sd] = json.loads(p.stdout)
    bad: list[dict] = []
    base_seed = seeds[0]
    for i, row in enumerate(runs[base_seed]):
        for sd in seeds[1:]:
            other = runs[sd][i]
            name = row.get("name", "?")
            if row.get("err") or other.get("err"):
                continue
            if digest_only:
                if row["digest"] != other["digest"]:
                    bad.append({"case": name, "kind": "digest 不同",
                                "detail": f"{base_seed}:{row['digest']} vs "
                                          f"{sd}:{other['digest']}"})
                continue
            if row["digest"] != other["digest"]:
                d = []
                if abs((row.get("pH") or 0) - (other.get("pH") or 0)) > 1e-9:
                    d.append(f"pH {row.get('pH')} vs {other.get('pH')}")
                if row.get("degree") != other.get("degree"):
                    d.append(f"degree {row.get('degree')} vs {other.get('degree')}")
                bad.append({"case": name,
                            "kind": "数值不同" if d else "仅摘要顺序不同",
                            "detail": "; ".join(d) or
                                      f"{row['digest'][:12]} vs {other['digest'][:12]}"})
    return bad


# ------------------------------------------------- ③ 输出守恒

def _side_atoms(side: dict) -> tuple[dict, float]:
    atoms: dict = {}
    q = 0.0
    for sp, coef in side.items():
        for el, n in elements_of(sp).items():
            atoms[el] = atoms.get(el, 0.0) + coef * n
        q += coef * charge_of(sp)
    return atoms, q


def _conservation(ne) -> list[str]:
    """净方程（Equation 结构）的元素与电荷守恒检查。返回违规描述。"""
    la, lq = _side_atoms(ne.left)
    ra, rq = _side_atoms(ne.right)
    out = []
    for el in sorted(set(la) | set(ra)):
        a, b = la.get(el, 0.0), ra.get(el, 0.0)
        # 相对容差：大系数允许 0.1% 的浮点残差，但绝不容忍"缺一个原子"。
        if abs(a - b) > max(1e-6, 1e-3 * max(abs(a), abs(b))):
            out.append(f"{el} 左{a:.4g} 右{b:.4g}（Δ{a - b:+.4g}）")
    if abs(lq - rq) > max(1e-6, 1e-3 * max(abs(lq), abs(rq))):
        out.append(f"电荷 左{lq:+.4g} 右{rq:+.4g}")
    return out


def check_cons(cases: list[dict], limit: int | None = None) -> list[dict]:
    """引擎**输出**的净方程与分步方程必须元素 + 电荷守恒。

    与 `helper.eqcheck` 的分工：eqcheck 验 `tests.json` 的**期望式**，
    本函数验**引擎实际产出的式子**。
    """
    use = cases[:limit] if limit else cases
    bad: list[dict] = []
    jobs = _jobs()
    if jobs > 1 and len(use) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for row in ex.map(_cons_one, use):
                bad.extend(row)
        return bad
    for c in use:
        bad.extend(_cons_one(c))
    return bad


def _cons_one(case: dict) -> list[dict]:
    """单例守恒检查（worker 入口：进程池需要可 pickle 的模块级函数）。"""
    bad: list[dict] = []
    try:
        r = _react(case)
    except Exception as e:                           # noqa: BLE001
        return [{"case": case["name"], "kind": "异常",
                 "detail": f"{type(e).__name__}: {e}"}]
    ne = r.net_equation
    if ne is not None:
        v = _conservation(ne)
        if v:
            bad.append({"case": case["name"], "kind": "净方程不守恒",
                        "detail": "; ".join(v)})
    for i, eq in enumerate(getattr(r, "equations", None) or [], 1):
        v = _conservation(eq)
        if v:
            bad.append({"case": case["name"],
                        "kind": f"分步方程#{i} 不守恒",
                        "detail": "; ".join(v)})
            break          # 每例只报第一条，避免同一病灶刷屏
    return bad


# ------------------------------------------------- ④ 绝热放热恒定

def _heat_one(case: dict, factors=(0.1, 1.0, 10.0)) -> list[dict]:
    """单例放热恒定检查（worker 入口）。

    **判据**：固定投料摩尔数，改变体积 ⟹ 绝热放热 `heat_kJ` 必须恒定
    （广延量：同摩尔数的同一反应，放热与体积无关）。

    **必须锁封闭体系**（`gas_escape=False`）：开放体系的逸出焓差
    （CO₂ 脱气 20.29 kJ/mol）是**真实物理**——逸出量随泡点（浓度）变，
    于是总放热随体积变是**正确的**，不是缺陷。封闭体系把逸出项归零，
    剩下的漂移就是账本缺陷（T2 主因：缺 ΔHf 的配离子被静默丢弃，
    第 306 轮已修）。

    历史（第 306 轮实测）：`CaCO₃ 0.05 + HCl 0.1` 封闭体系修前
    1.196→1.672 kJ（V=0.1→2.0，漂 0.48）；修后应恒定。
    """
    bad: list[dict] = []
    subs0 = _subs_of(case)
    if not subs0:
        return bad
    cond = case.get("cond") or {}
    V0 = cond.get("V_L", 1.0)
    T0 = cond.get("T_K", 298.15)
    base = None
    for f in factors:
        try:
            r = _react({"subs": list(subs0.items()),
                        "cond": {"V_L": V0 * f, "T_K": T0}},
                       gas_escape=False)
        except Exception as e:                   # 引擎崩了也是缺陷
            bad.append({"case": case["name"], "factor": f,
                        "kind": "异常", "detail": f"{type(e).__name__}: {e}"})
            continue
        th = r.thermal or {}
        q = th.get("heat_kJ")
        if q is None:
            # 无热值（H/O 不闭合/无变换）：跨体积必须一致地无
            if base is not None and base[0] is not None:
                bad.append({"case": case["name"], "factor": f,
                            "kind": "热值有无",
                            "detail": f"f={base[1]} → {base[0]}, "
                                      f"f={f} → None"})
            if base is None:
                base = (None, f)
            continue
        if base is None:
            base = (q, f)
            continue
        b, bf = base
        if b is None:
            bad.append({"case": case["name"], "factor": f, "kind": "热值有无",
                        "detail": f"f={bf} → None, f={f} → {q}"})
            continue
        # 相对容差 2%：ΔHf 表值 3 位有效数字 + 配位分配残差
        if abs(q - b) > 0.02 * max(abs(q), abs(b), 1e-9):
            bad.append({"case": case["name"], "factor": f, "kind": "放热漂移",
                        "detail": f"f={bf} → {b:.4f}, f={f} → {q:.4f} "
                                  f"（Δ{q - b:+.4f} kJ）"})
    return bad


def check_heat(cases: list[dict], factors=(0.1, 1.0, 10.0),
               limit: int | None = None) -> list[dict]:
    """固定投料改体积 ⟹ 绝热放热恒定（封闭体系口径）。"""
    use = cases[:limit] if limit else cases
    bad: list[dict] = []
    jobs = _jobs()
    if jobs > 1 and len(use) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for row in ex.map(_heat_one, use, [factors] * len(use)):
                bad.extend(row)
        return bad
    for c in use:
        bad.extend(_heat_one(c, factors))
    return bad


# ---------------------------------------------------------------- CLI

_KINDS = {
    "scale": ("① 尺度不变性", check_scale),
    "det": ("② 进程确定性", check_det),
    "cons": ("③ 输出守恒", check_cons),
    "heat": ("④ 绝热放热恒定", check_heat),
}


def main(argv: list[str]) -> int:
    args = list(argv)
    limit = None
    if "--n" in args:
        i = args.index("--n")
        limit = int(args[i + 1])
        del args[i:i + 2]
    jobs = _jobs()
    for flag in ("--jobs",):
        if flag in args:
            i = args.index(flag)
            jobs = int(args[i + 1])
            del args[i:i + 2]
    picks = [a for a in args if a in _KINDS] or list(_KINDS)
    verbose = "--all" in args
    cases = load_cases()
    if limit:
        cases = cases[:limit]
    total_bad = 0
    import time
    t0 = time.perf_counter()
    print(f"不变量审计：{len(cases)} 例 · 检查 {', '.join(picks)}"
          f" · jobs={jobs}")
    for key in picks:
        title, fn = _KINDS[key]
        t1 = time.perf_counter()
        bad = fn(cases)
        el = time.perf_counter() - t1
        total_bad += len(bad)
        flag = "✓ 通过" if not bad else f"✗ {len(bad)} 条违规"
        print(f"\n{title}  {flag}  [{el:.1f}s]")
        show = bad if verbose else bad[:12]
        for b in show:
            print(f"  · {b['case'][:34]:36} [{b['kind']}] {b['detail'][:100]}")
        if len(bad) > len(show):
            print(f"  … 另有 {len(bad) - len(show)} 条（--all 全列）")
    print(f"\n===== 违规合计 {total_bad} 条 "
          f"（总耗时 {time.perf_counter() - t0:.1f}s）=====")
    return 0 if total_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
