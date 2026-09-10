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
            # 质量口径：非冻结两侧平衡的最大 |S|（冻结=宣告平衡，语义豁免）
            rec["resid_live"] = round(
                max((abs(a["S"]) for a in probe["active"]
                     if a["two_sided"] and not a["frozen"]), default=0.0), 3)
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
