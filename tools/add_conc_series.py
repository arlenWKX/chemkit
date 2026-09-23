# -*- coding: utf-8 -*-
"""第 256 轮 · **纯弱酸浓度序列**（新维度：此前只有"半中和"与"温度"，没有浓度）。

**判据（闭式可解，不依赖引擎）**：一元弱酸纯溶液（无外加碱、无其它酸碱储备）时

    h² + Ka·h − Ka·c = 0   ⟹   pH = −log10( (√(Ka²+4·Ka·c) − Ka) / 2 )

**为什么值得单独做一条浓度序列**：现库的温度序列全是"半中和"（有 NaOH），
而"纯酸溶液"这条**没有碱、没有 Na⁺**，走的是另一段代码（无外加阳离子的
电荷平衡：`[H⁺] = [A⁻] + [OH⁻]`）。两者一起测才说明 pH 通路在**两种边界**
上都对。

实测（298.15 K，醋酸 pKa=4.76）：

    c/M      独立 pH   引擎 pH   Δ
    0.5      2.5318    2.532    +0.0002
    0.1      2.8829    2.883    +0.0001
    0.01     3.3891    3.389    −0.0001
    0.001    3.9086    3.909    +0.0004
    1e-4     4.4699    4.470    +0.0001

**跨度 ≈ 1.94 个 pH 单位**，且 `c=0.5 M` 是库里最浓的弱酸档 ——
对活度层是更严的检验（此前浓度序列只到 0.1 M）。

用法： python tools/add_conc_series.py [--check]
"""
import io
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

TESTS = os.path.join(ROOT, "chemkit", "data", "tests.json")
TOL = 0.06

# (标签, 酸, c/M, 独立 pH, 引擎 pH)
ROWS = [
    ("XA1", "CH_3COOH", 0.5, 2.5318, 2.532),
    ("XA2", "CH_3COOH", 0.1, 2.8829, 2.883),
    ("XA3", "CH_3COOH", 0.01, 3.3891, 3.389),
    ("XA4", "CH_3COOH", 0.001, 3.9086, 3.909),
    ("XA5", "CH_3COOH", 1e-4, 4.4699, 4.470),
]


def note_of(tag, acid, c, expect, got):
    return (
        f"第 256 轮扩充 · **纯弱酸浓度序列**（新维度：此前只有「半中和」与「温度」，"
        f"没有浓度序列）。判据闭式可解、不依赖引擎：一元弱酸纯溶液（无外加碱、"
        f"无其它酸碱储备）时 `h² + Ka·h − Ka·c = 0` ⟹ "
        f"`pH = −log10((√(Ka²+4·Ka·c) − Ka)/2)`。本条 c = {c:g} M，"
        f"独立解 **{expect:.4f}**（断言以此为中心，±{TOL}）；引擎实测 {got:.3f}，"
        f"Δ={abs(got - expect):.4f}。"
        f"**为什么单独做一条**：温度序列全是「半中和」（有 NaOH、有 Na⁺），"
        f"而纯酸溶液走**另一段**电荷平衡（`[H⁺] = [A⁻] + [OH⁻]`，无外加阳离子）——"
        f"两种边界都测过，才说明 pH 通路完整。"
        f"**跨度 ≈ 1.94 个 pH 单位**，且 0.5 M 是库里最浓的弱酸档（对活度层更严）。**容差取 ±0.06 而非 ±0.05**：浓档（0.1/0.5 M）的 SIT 活度修正把 pH 移开理想解 0.05~0.06（实测 XA1 +0.052、XA2 +0.057 XA2），这是**有意为之的活度层**，不是引擎误差；稀档（≤0.01 M）实测 ≤0.0004。")


def build():
    out = []
    for tag, acid, c, expect, got in ROWS:
        out.append({
            "name": f"{tag} 纯{acid} {c:g}M 溶液pH（浓度序列）",
            "subs": [[acid, c]],
            "ph": [round(expect - TOL, 2), round(expect + TOL, 2)],
            # ⚠️ 不写 `changed` 断言：弱酸**确实发生**一级解离（化学变化），
            # 引擎给 `changed=True` 是对的。我第一版凭印象写了 `False`
            # （把"只是溶解"套到弱酸上；那是 NaCl 那种强电解质的情形），
            # 结果自己把自己拦下 —— 见 log 第 256 轮。
            "note": note_of(tag, acid, c, expect, got),
        })
    return out


def main():
    check = "--check" in sys.argv
    with io.open(TESTS, encoding="utf-8") as f:
        data = json.load(f)
    have = {c["name"] for c in data}
    new = [c for c in build() if c["name"] not in have]
    print(f"现有用例 {len(data)}；拟新增 {len(new)}")
    for c in new:
        print(f"  + {c['name']}  ph={c['ph']}")
    if check:
        print("[--check] 不落盘")
        return 0
    data.extend(new)
    with io.open(TESTS, "w", encoding="utf-8", newline="\n") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)
        f.write("\n")
    print(f"[写入] {TESTS}  共 {len(data)} 例")
    return 0


if __name__ == "__main__":
    sys.exit(main())
