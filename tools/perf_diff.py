# -*- coding: utf-8 -*-
"""**同用例对照**的性能/通过性差分：读两份 `logs/suite-*.json`，只在**共同用例**
上比较逐例 `iters` / `ms`，并列出通过性翻转。

## 为什么需要它（而不是 `dev.py perf`）

`dev.py perf` 比的是**受控基线** `.tmp_dev_converg.json`，那是**上一次刷新时**
的用例集：一旦本轮增删了用例，它的"位移"清单里就会混进 `n` 这类纯规模量，
而且 `iters_total`/`sof_total` 的增减**分不清是"用例变多"还是"引擎变慢"**。

本工具做的是**配对比较**：只取两份档案里**同名**的用例，逐例算
`Δiters`、`Δms`，因此：

* 用例数变化被自动剔除；
* 能直接定位"**是哪几例**变慢了多少"（算法病态 vs 全局变慢）。

用法：

    python tools/perf_diff.py <基线档> <对照档> [-n 8]

纪律（`agents/lessons.md` §2「跑一次读多次」）：性能数字一律从**已存的档案**里读，
不要为了"换个角度看同一批数据"重跑套件。
"""
from __future__ import annotations

import json
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def _load(p: str) -> dict:
    with open(p, encoding="utf-8") as fh:
        return {r["name"]: r for r in json.load(fh)}


def _i(r: dict) -> int:
    return int(r.get("iters") or 0)


def _ms(r: dict) -> float:
    return float(r.get("ms") or 0.0)


def main(argv: list[str]) -> int:
    top = 8
    if "-n" in argv:
        i = argv.index("-n")
        top = int(argv[i + 1])
        del argv[i:i + 2]
    if len(argv) < 2:
        print(__doc__)
        return 2
    A, B = _load(argv[0]), _load(argv[1])
    common = [n for n in A if n in B]
    only_b = [n for n in B if n not in A]
    ia, ib = sum(_i(A[n]) for n in common), sum(_i(B[n]) for n in common)
    ma, mb = sum(_ms(A[n]) for n in common), sum(_ms(B[n]) for n in common)
    print(f"共同用例 {len(common)}（另有仅 B 有 {len(only_b)}）")
    print(f"  迭代合计  {ia:>8} -> {ib:>8}   {ib / ia - 1:+.1%}" if ia else "")
    print(f"  耗时合计  {ma / 1000:>7.1f}s -> {mb / 1000:>7.1f}s  "
          f"{mb / ma - 1:+.1%}" if ma else "")
    ka = sum(1 for n in common if A[n].get("ok"))
    kb = sum(1 for n in common if B[n].get("ok"))
    print(f"  通过      {ka:>8} -> {kb:>8}   {kb - ka:+d}")

    print(f"\nΔiters 最大 {top} 例：")
    for n in sorted(common, key=lambda x: _i(B[x]) - _i(A[x]), reverse=True)[:top]:
        print(f"   {n[:42]:<44} iters {_i(A[n]):>6} -> {_i(B[n]):>6} "
              f"({_i(B[n]) - _i(A[n]):+6d})   ms {_ms(A[n]):>8.0f} -> {_ms(B[n]):>8.0f}")

    print(f"\n耗时最大 {top} 例（对照档）：")
    for n in sorted(common, key=lambda x: _ms(B[x]), reverse=True)[:top]:
        print(f"   {n[:42]:<44} {_ms(B[n]):>9.0f} ms   iters {_i(B[n]):>6}"
              f"   （基线 {_ms(A[n]):>8.0f} ms / {_i(A[n])} iters）")

    flips = [n for n in common if A[n].get("ok") != B[n].get("ok")]
    print(f"\n通过性翻转 {len(flips)} 例（翻绿 "
          f"{sum(1 for n in flips if B[n].get('ok'))} / 翻红 "
          f"{sum(1 for n in flips if not B[n].get('ok'))}）：")
    for n in flips:
        tag = "绿" if B[n].get("ok") else "红"
        err = (B[n].get("errors") or [""])[0][:60]
        print(f"   {tag}  {n[:44]:<46} {err}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
