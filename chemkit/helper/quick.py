"""快速子集测试：只跑指定用例（按名称前缀/序号），带单例计时。

用途：针对某一处改动的验证，**只跑受影响的用例**，避免每次全量 1173 例
的开销。典型用法是先把上一轮报错项喂进来：

    python tools/quick.py --fails            # 读 _fails.txt 里的上轮失败项
    python tools/quick.py V16 Z31 N30        # 按名称前缀
    python tools/quick.py --list 12 73 501   # 按清单序号（同 testsuit 的编号）

约定：失败项落盘到 `_fails.txt`（一行一个用例名），供下一轮 `--fails` 复用。
"""
from __future__ import annotations

import sys
import time

sys.path.insert(0, ".")

from chemkit import testsuit as ts        # noqa: E402
from chemkit.data import load_tables      # noqa: E402

FAILS_PATH = "_fails.txt"


def _match(cases: list[dict], keys: list[str]) -> list[dict]:
    """按**完整名**优先、否则前缀匹配（前缀会撞车：N09 同时命中
    `N09 K2Cr2O7+FeSO4酸` 与 `N09 BeCl2+过量NaOH`）。"""
    exact = [c for c in cases if c["name"] in keys]
    if exact:
        return exact
    return [c for c in cases if c["name"].startswith(tuple(keys))]


def main(argv: list[str]) -> int:
    T = load_tables()
    cases = ts.load_cases(None)
    picks: list[dict] = []
    if "--fails" in argv:
        try:
            names = [ln.strip() for ln in open(FAILS_PATH, encoding="utf-8")
                     if ln.strip()]
        except FileNotFoundError:
            print(f"没有 {FAILS_PATH}；先跑一次 tools/quick.py <名称…>")
            return 2
        picks = [c for c in cases if c["name"] in set(names)]
        print(f"上一轮失败项 {len(names)} 个，匹配到 {len(picks)} 例")
    elif "--list" in argv:
        idx = {int(x) for x in argv[argv.index("--list") + 1:] if x.isdigit()}
        picks = [c for i, c in enumerate(cases, 1) if i in idx]
    else:
        keys = [x for x in argv if not x.startswith("-")]
        if not keys:
            print(__doc__)
            return 2
        picks = _match(cases, keys)
    if not picks:
        print("未匹配到用例")
        return 2

    ts.judge_warm(T) if hasattr(ts, "judge_warm") else _warm(T)
    ts.FAILS.clear()
    ts.PASS_N = 0
    ts.TIMES.clear()
    t0 = time.perf_counter()
    for c in picks:
        ts.run_case(c, T)
    wall = time.perf_counter() - t0
    n = ts.PASS_N + len(ts.FAILS)
    print(f"\n===== 子集 {ts.PASS_N}/{n} PASS（{wall:.1f}s，"
          f"{wall / max(n, 1) * 1000:.1f} ms/例）=====")
    with open(FAILS_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(ts.FAILS))
    if ts.FAILS:
        print(f"失败项已写入 {FAILS_PATH}（下一轮用 --fails 复跑）")
    mx = sorted(ts.TIMES, reverse=True)[:5]
    if mx and mx[0][0] > 0.1:
        print("最慢 5 例：", "  ".join(f"{t*1000:.0f}ms {nm[:34]}"
                                     for t, nm in mx))
    return 1 if ts.FAILS else 0


def _warm(T) -> None:
    from chemkit.engine import judge
    judge([{"name": "NaCl", "mol": 0.1}], {"V_L": 1.0}, T)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
