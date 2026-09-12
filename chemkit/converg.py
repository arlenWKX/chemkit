"""chemkit.converg：收敛质量基准与全量差分工具（v0.3.8）。

背景（worklog Task 12 教训）：两次深水区试探都被 1247 全量差分否决
（Co32/33 轨迹重排）——动 walk 语义前必须先立"收敛质量"基准。本模块
把基准制度化：

  ① 残差指标：walk 退出点上对"两侧均在场且未被冻结"的平衡评
     S = logK − logQ，热力学终态应全部为零；max|S| 即该用例的
     欠收敛度（log 单位）。冻结平衡（frozen）是引擎的既定语义
     （宣告平衡止震），不计入质量指标；disabled 平衡计入——
     J06 型欠收敛正是"微步禁用把仍有驱动的平衡锁死"。
  ② 全量 dump：1218 用例 × {结果摘要 + 残差画像 + 计时}，
     JSON 快照作差分基线（sha256 摘要比对，bit 级语义守护）。
  ③ 差分报告：逐例对比两份 dump，分类语义位移
     （A1 changed / degree / pH / 量值 / 残差方向）。

用法：
    python -m chemkit.converg dump  [-o out.json]   # 全量 dump
    python -m chemkit.converg diff  a.json b.json   # 差分报告
    python -m chemkit.converg top   [-n 30] [json]  # 慢例/欠收敛排行
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

from .candidates import ANN_MIN_EXTENT
from .data import load_tables
from .engine import judge
from .testsuit import load_cases, DEFAULT_CASES


# ---------------------------------------------------------- 摘要与差分

def _result_digest(r: dict) -> str:
    """结果字典的语义摘要（量值与判定维度；steps 计数入摘要、轨迹不入
    ——轨迹是求解过程自由度，摘要锁的是对外呈现语义）。"""
    keys = ("changed", "reacted", "degree", "annotations", "consumption",
            "production", "final", "escaped", "ionize", "net_equation",
            "equations", "final_pH", "unknown")
    core = {k: r.get(k) for k in keys if k in r}
    blob = json.dumps(core, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def _live(active: list, only: str | None = None) -> float:
    """质量残差 = **值得解且 walk 会解**的两侧平衡的最大 |S|。

    口径原则（不是"调到好看"，而是"只统计求解器承诺要解的东西"）：

    **第一层：|S| 大 ≠ 欠收敛。** `S = logK − logQ` 是 log 尺度驱动力，
    痕量物种对可以有巨大的 S 而只能走 ~1e-6 mol——那是"无关"而非"没
    收敛"。实测：D38 Pb(NO3)2+K2CrO4 的 `resid_live` 长期报 **33.398**，
    来自 PbO₂ 2.9e-06 / Cr³⁺ 1.9e-06 mol 的一对（PbCrO₄ 已沉 0.99996）。
    引擎自身的"显著程度"判据是 `ANN_MIN_EXTENT`（slow 标注同用：
    "仅当慢反应可达显著程度才标注"），故口径取
    `ext_max ≥ ANN_MIN_EXTENT`（ext_max = 驱动方向反应物的化学计量上限）。

    **第二层：引擎既定语义不算欠收敛。** walk 的候选评估循环对下列三类
    直接 `continue`，它们**从未进入求解承诺**：
      · `frozen`  —— 引擎宣告平衡止震（极限环/实测仲裁）；
      · `slow`    —— 动力学层判定永不执行、只作标注；
      · `blocked` —— 致密膜抑制溶剂氧化通道（钝化模型）。
    实测代价（旧口径只排除 frozen）：TS04 报 **17.975**（实为
    `SO_4^{2-}/S_2O_3^{2-}` 已动力学封闭的幻影硫酸盐通道）、
    E24 Al+NaOH 报 **156.12**（膜封锁的 Al/H₂O 通道）。

    `disabled` 仍**计入**（J06 型欠收敛正是"微步禁用把仍有驱动的平衡
    锁死"），这是已知真实病灶，不能被口径优化掉。

    `only` 给类别名时改报**该类被排除者**的 |S|（透明化：不统计，但
    必须看得见）。行末 `ext_max` 恒为真，故对 trace 类用 `only="trace"`。
    """
    def _ok(a: dict) -> bool:
        return (a["two_sided"] and not a["frozen"]
                and not a.get("slow") and not a.get("blocked")
                and a.get("ext_max", float("inf")) >= ANN_MIN_EXTENT)

    def _pick(a: dict) -> bool:
        if only == "trace":
            return (a["two_sided"] and not a["frozen"]
                    and not a.get("slow") and not a.get("blocked")
                    and a.get("ext_max", float("inf")) < ANN_MIN_EXTENT)
        return bool(a["two_sided"] and a.get(only))

    f = _pick if only is not None else _ok
    return round(max((abs(a["S"]) for a in active if f(a)), default=0.0), 3)


def run_all(cases_path: str | None = None, with_probe: bool = True) -> list[dict]:
    """跑全量用例，返回逐例记录（结果摘要 + 收敛画像 + 计时）。"""
    cases = load_cases(cases_path)
    T = load_tables()
    out = []
    for c in cases:
        subs = [{"name": n, "mol": m} for n, m in c["subs"]]
        probe = {} if with_probe else None
        t0 = time.perf_counter()
        r = judge(subs, c.get("cond") or {"V_L": 1.0}, T, _probe=probe)
        dt = (time.perf_counter() - t0) * 1000.0
        rec = {
            "name": c["name"], "ms": round(dt, 1),
            "deg": r["degree"], "changed": r["changed"],
            "reacted": r.get("reacted"), "pH": r.get("final_pH"),
            "steps_n": len(r.get("steps", [])), "digest": _result_digest(r),
        }
        if probe is not None and probe:   # OVERRIDE 直出路径无探针画像
            rec["exit"] = probe["exit"]
            rec["iters"] = probe["iters"]
            rec["resid"] = round(probe["max_abs_S"], 3)
            # 质量口径：**值得解且 walk 会解**的两侧平衡最大 |S|
            # （|S| 大而可达程度 < ANN_MIN_EXTENT 的痕量方向不算欠收敛；
            #  frozen/slow/blocked 三类既定语义排除）。下列各项把每类
            # 被排除者的最大 |S| 单列——不统计，但必须看得见。
            rec["resid_live"] = _live(probe["active"])
            rec["resid_frozen"] = _live(probe["active"], "frozen")
            rec["resid_slow"] = _live(probe["active"], "slow")
            rec["resid_blocked"] = _live(probe["active"], "blocked")
            rec["resid_trace"] = _live(probe["active"], "trace")
        out.append(rec)
    return out


def dump(path_out: str = "converg-baseline.json",
         cases_path: str | None = None) -> str:
    """全量 dump：逐例记录 + 总摘要（计时分位数/残差分布/哈希指纹）。"""
    recs = run_all(cases_path)
    ms = sorted(x["ms"] for x in recs)
    res = sorted(x.get("resid_live", 0.0) for x in recs)
    n = len(recs)
    iters = sorted(x.get("iters", 0) for x in recs)
    summary = {
        "n": n,
        "ms_mean": round(sum(ms) / n, 2),
        "ms_p50": round(ms[n // 2], 2),
        "ms_p90": round(ms[int(n * 0.9)], 2),
        "ms_max": round(ms[-1], 2),
        "n_gt50": sum(1 for m in ms if m > 50),
        "n_gt100": sum(1 for m in ms if m > 100),
        "n_gt500": sum(1 for m in ms if m > 500),
        "iters_total": sum(iters),
        "iters_p99": iters[int(n * 0.99)],
        "resid_p50": round(res[n // 2], 3),
        "resid_p90": round(res[int(n * 0.9)], 3),
        "resid_max": round(res[-1], 3),
        "n_resid_gt1": sum(1 for x in res if x > 1.0),
        "n_resid_gt0p1": sum(1 for x in res if x > 0.1),
        "digest_all": hashlib.sha256(
            json.dumps([x["digest"] for x in recs]).encode()).hexdigest()[:16],
    }
    doc = {"summary": summary, "cases": recs}
    with open(path_out, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
    return path_out


def _load(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def diff(path_a: str, path_b: str, verbose: bool = True) -> dict:
    """两份 dump 的差分报告。分类：语义位移（digest/deg/pH/量）与
    纯性能变化（ms/iters 变了但 digest 不变）。"""
    A, B = _load(path_a)["cases"], _load(path_b)["cases"]
    ma = {x["name"]: x for x in A}
    stats = {"sem": [], "perf_only": [], "new_pass_like": []}
    for b in B:
        a = ma.get(b["name"])
        if a is None:
            stats["new_pass_like"].append(b["name"])
            continue
        sem = (a["digest"] != b["digest"] or a["deg"] != b["deg"]
               or a.get("changed") != b.get("changed")
               or a.get("pH") != b.get("pH"))
        if sem:
            stats["sem"].append((a, b))
        elif (a["ms"], a.get("iters")) != (b["ms"], b.get("iters")):
            stats["perf_only"].append((a, b))
    if verbose:
        print(f"语义位移 {len(stats['sem'])} 例 / 纯性能 {len(stats['perf_only'])} 例")
        for a, b in stats["sem"][:40]:
            print(f"  [sem] {a['name']}: deg {a['deg']}→{b['deg']} "
                  f"pH {a.get('pH')}→{b.get('pH')} "
                  f"resid_live {a.get('resid_live')}→{b.get('resid_live')} "
                  f"ms {a['ms']}→{b['ms']}")
        tot_a = sum(x["ms"] for x in A)
        tot_b = sum(x["ms"] for x in B)
        print(f"总时长 {tot_a:.0f}ms → {tot_b:.0f}ms "
              f"({(tot_b - tot_a) / max(tot_a, 1) * 100:+.1f}%)")
    return stats


def top(n: int = 30, path: str | None = None) -> None:
    """慢例 / 欠收敛例 / 迭代数排行。"""
    doc = _load(path or "converg-baseline.json")
    cases = doc["cases"]
    print("== 最慢 ==")
    for x in sorted(cases, key=lambda c: -c["ms"])[:n]:
        print(f"  {x['ms']:8.1f}ms iters={x.get('iters', 0):5d} "
              f"resid={x.get('resid_live', 0):6.2f} exit={x.get('exit', '?'):9s} "
              f"{x['name'][:60]}")
    print("== 欠收敛（live 残差）==")
    for x in sorted(cases, key=lambda c: -c.get("resid_live", 0))[:n]:
        print(f"  resid={x.get('resid_live', 0):6.2f} ms={x['ms']:8.1f} "
              f"exit={x.get('exit', '?'):9s} {x['name'][:60]}")


def main(argv: list[str]) -> None:
    cmd = argv[0] if argv else "dump"
    if cmd == "dump":
        args = argv[1:]
        out = "converg-baseline.json"
        if "-o" in args:
            out = args[args.index("-o") + 1]
        cases = None
        for a in args:
            if a.endswith(".json") and a != out and "-o" not in args[max(0, args.index(a) - 1):args.index(a)]:
                cases = a
        p = dump(out, cases)
        doc = _load(p)
        print(json.dumps(doc["summary"], ensure_ascii=False, indent=1))
        print(f"dump → {p}")
    elif cmd == "diff":
        diff(argv[1], argv[2])
    elif cmd == "top":
        n = 30
        if "-n" in argv:
            i = argv.index("-n")
            n = int(argv[i + 1])
            argv = argv[:i] + argv[i + 2:]
        args = [a for a in argv[1:] if not a.startswith("-")]
        top(n, args[0] if args else None)
    else:
        print(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
