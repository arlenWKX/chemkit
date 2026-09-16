"""tools/dev.py — 一轮工作的**固定工具链**（单一入口，控制台安全，输出有界）。

为什么要有它：0.5.x 的每一轮都要跑"全量套件 + 确定性指标 + 单例深探 +
锚定开关 + 守恒核验 + 补丁"，此前每次都用临时脚本/内联 `python -c` 现拼，
踩过的坑全部沉淀成本工具的**默认行为**：

  · PowerShell 不支持 heredoc（`<<'PY'`）⟹ 一律走文件，不再内联多行脚本；
  · 控制台 GBK ⟹ 本工具内置 UTF-8 重配置，且**只打印 ASCII 记号**
    （方程一律 `.plain()`：`->`/`<=>`；不接受 TeX 串进控制台）；
  · 输出淹没 ⟹ `suite` 只打摘要 + 失败明细（全套 1200+ 行 `[ok]` 丢弃）；
  · "这轮到底是哪个 pKw 约定" ⟹ **每个子命令首行强制打印约定**；
  · `converg.dump()` 会覆盖受控基线文件 ⟹ `perf` 默认写临时文件并**与
    `HEAD` 的基线逐项对比**（确定性指标 Δ + 判定行），只有 `--write` 才落盘；
  · `git stash` 做前后对比容易丢改动 ⟹ `snapshot`/`cmp` 用 JSON 快照对比；
  · 手写补丁脚本易错 ⟹ `patch` 走声明式 `PATCHES`，先全量校验再原子落盘，
    且保留原换行风格；
  · 行尾噪声（`git status` 显示 M 但 `git diff` 为空）⟹ `hygiene` 检出并
    `--fix` 一键还原。

用法（仓库根目录）：
  python tools/dev.py guide                      # 打印本轮协议（10 行）
  python tools/dev.py suite [前缀...]            # 全量/子集套件：摘要 + 失败明细
  python tools/dev.py perf                       # 收敛基线：确定性指标 vs HEAD
  python tools/dev.py case 21 T89 [--all]        # 单例深探：步表/账本/两版净方程/断言判定
  python tools/dev.py snapshot .tmp_a.json       # 全库关键输出快照（含净方程）
  python tools/dev.py cmp .tmp_a.json .tmp_b.json# 只打差异（前后对比用）
  python tools/dev.py anchor on|off|status       # pKw 锚定一键切换（幂等）
  python tools/dev.py eqcheck                    # 全库方程守恒精确核验
  python tools/dev.py hygiene [--fix]            # 行尾/临时文件卫生
  python tools/dev.py patch spec.py [--check]    # 声明式补丁（原子 + 计数断言）
"""
from __future__ import annotations

import io
import json
import os
import subprocess
import sys
from contextlib import redirect_stdout

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                    # pragma: no cover
        pass

_TMP = ".tmp_dev_"
# 墙钟字段：只报不比对（同一份代码在不同时刻能差 30%+，见 §7 性能纪律）
_MS_KEYS = ("ms_mean", "ms_p50", "ms_p90", "ms_max",
            "n_gt50", "n_gt100", "n_gt500")


def _run(cmd: list[str], cwd: str = ROOT) -> str:
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    return (r.stdout or "") + (r.stderr or "")


def convention() -> str:
    """当前 pKw 约定（每个子命令首行都打，避免"这轮跑的是哪个约定"）。"""
    from chemkit.core import pKw_of
    v = pKw_of(298.15)
    tag = "未锚定" if abs(v - 14.0) > 1e-9 else "锚定 14.0"
    return f"[约定] pKw(298.15) = {v:.6f}  ({tag})   pKw(373) = {pKw_of(373.15):.4f}"


def _banner(verb: str) -> None:
    print(f"== dev.py {verb} ==")
    print(convention())


# --------------------------------------------------------------- suite
def cmd_suite(prefixes: list[str], keep: int = 40) -> int:
    _banner("suite" + (f" {' '.join(prefixes)}" if prefixes else "（全量）"))
    from chemkit import testsuit
    tmp = None
    if prefixes:
        cases = [c for c in testsuit.load_cases(None)
                 if c["name"].startswith(tuple(prefixes))]
        if not cases:
            print("!! 没有匹配的用例")
            return 2
        tmp = _TMP + "suite.json"
        with io.open(tmp, "w", encoding="utf-8") as f:
            json.dump(cases, f, ensure_ascii=False)
        print(f"[子集] {len(cases)} 例")
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = testsuit.main(tmp, None)
    out = buf.getvalue()
    lines = out.split("\n")
    for ln in lines:                      # 摘要 + 结论行（其余明细丢弃）
        if ln.startswith("=====") or ln.startswith("失败:") \
                or ln.startswith("总耗时") or "环闭合" in ln:
            print(ln[:400])
    bad = [r for r in testsuit.RESULTS if not r.get("ok")]   # 失败明细（套件本体不打）
    for r in bad[:12]:
        print(f"  [FAIL] {r['name']} -- {'; '.join(r.get('errors') or [])[:300]}")
    if len(bad) > 12:
        print(f"  ... 另有 {len(bad) - 12} 例失败（全部见 `失败:` 清单）")
    if tmp and os.path.exists(tmp):
        os.remove(tmp)
    return rc


# --------------------------------------------------------------- perf
def cmd_perf(baseline: str = "converg-baseline.json", write: bool = False) -> int:
    _banner("perf")
    from chemkit import converg
    out = baseline if write else _TMP + "converg.json"
    buf = io.StringIO()
    with redirect_stdout(buf):
        converg.dump(out)
    new = json.load(io.open(out, encoding="utf-8"))
    head = _run(["git", "show", f"HEAD:{baseline}"])
    if not head.strip():
        print("!! 取不到 HEAD 基线，跳过对比")
        return 1
    old = json.loads(head)
    o, n = old["summary"], new["summary"]
    keys = [k for k in n if k not in _MS_KEYS]
    print(f"{'指标':<16}{'HEAD':>16}{'本轮':>16}   Δ")
    moved = []
    for k in keys:
        a, b = o.get(k), n.get(k)
        same = a == b
        if not same and k not in ("digest_all",):
            moved.append(k)
        mark = "" if same else "  ← 位移"
        if k == "digest_all":
            mark = "" if same else "  ← digest 位移（辅助信号，非门槛）"
        print(f"{k:<16}{str(a):>16}{str(b):>16}{mark}")
    print(f"[墙钟] mean={n['ms_mean']}ms p90={n['ms_p90']}ms "
          f"max={n['ms_max']}ms（机器噪声，不作门槛）")
    if not moved:
        print("[判定] 确定性指标与 HEAD 完全一致（纯呈现层改动应如此）")
    else:
        print(f"[判定] 确定性指标位移：{moved}（须逐条给出化学理由）")
    if not write:
        print(f"[提示] 未覆盖受控基线（写临时 {out}）；确需刷新基线加 --write")
    else:
        print(f"[提示] 已覆盖 {baseline}（只需刷新墙钟时才这么做）")
    return 0


# --------------------------------------------------------------- case
def cmd_case(prefixes: list[str], show_all: bool = False, limit: int = 14) -> int:
    _banner("case " + " ".join(prefixes))
    from chemkit.data import load_tables
    from chemkit.engine import judge
    from chemkit.system import Reaction
    from chemkit.testsuit import RESULTS, load_cases, run_case
    T = load_tables()
    hit = 0
    for c in load_cases(None):
        if not c["name"].startswith(tuple(prefixes)):
            continue
        hit += 1
        r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
        R = Reaction(r)
        print(f"\n-- {c['name']}  subs={c['subs']} cond={c.get('cond')}")
        print(f"   degree={r['degree']} changed={r['changed']} "
              f"reacted={r['reacted']} pH={r.get('final_pH')} "
              f"He={r.get('H_excess_initial')}->{r.get('H_excess')}")
        ne = r.get("net_exact") or {}
        if ne:
            cc = {k: v for k, v in ne["c"].items() if v > 1e-9}
            pp = {k: v for k, v in ne["p"].items() if v > 1e-9}
            fmt = lambda d: " + ".join(f"{v:.6g}{k}" for k, v in  # noqa: E731
                                       sorted(d.items(), key=lambda kv: -kv[1]))
            print(f"   账本净差: {fmt(cc)}  ->  {fmt(pp)}")
        steps = r.get("steps") or []
        order = steps if show_all else sorted(steps, key=lambda s: -s.get("extent", 0))[:limit]
        print(f"   步表（{len(steps)} 步{'，按 extent 前 %d' % limit if not show_all else ''}）:")
        for s in order:
            print(f"     [{s.get('kind','?'):9s}] ext={s.get('extent', 0.0):<10.5g} "
                  f"logK={s.get('logK')!s:<7.4} S={s.get('S')!s:<7.4} "
                  f"conv={s.get('conversion')!s:<5.3} {(s.get('equation') or '')[:70]}")
        n_eq, n_raw = R.net_equation, R.net_equation_raw
        print(f"   精编: {None if n_eq is None else n_eq.plain()}")
        print(f"   原始: {None if n_raw is None else n_raw.plain()}")
        print(f"   同对象: {n_eq is n_raw}")
        with redirect_stdout(io.StringIO()):    # 套件自身会打一行 PASS/FAIL+note
            ok = run_case(c, T)
        rec = RESULTS[-1] if RESULTS else {}
        print("   判定: " + ("PASS" if ok else "FAIL -- " +
                            "; ".join(rec.get("errors") or [])[:400]))
    if not hit:
        print("!! 没有匹配的用例")
        return 2
    return 0


# --------------------------------------------------------------- snapshot / cmp
def _snap_one(r: dict) -> dict:
    from chemkit.system import Reaction
    R = Reaction(r)
    return {
        "degree": r.get("degree"), "changed": r.get("changed"),
        "reacted": r.get("reacted"), "pH": r.get("final_pH"),
        "net": None if R.net_equation is None else R.net_equation.plain(),
        "raw": None if R.net_equation_raw is None else R.net_equation_raw.plain(),
        "steps": [f"{s.get('equation')}@{s.get('extent')}" for s in (r.get("steps") or [])],
    }


def cmd_snapshot(path: str, prefixes: list[str]) -> int:
    _banner("snapshot")
    from chemkit.data import load_tables
    from chemkit.engine import judge
    from chemkit.testsuit import load_cases
    T = load_tables()
    data = {}
    for c in load_cases(None):
        if prefixes and not c["name"].startswith(tuple(prefixes)):
            continue
        r = judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T)
        data[c["name"]] = _snap_one(r)
    with io.open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, sort_keys=True)
    print(f"[快照] {len(data)} 例 -> {path}")
    return 0


def cmd_cmp(a: str, b: str, limit: int = 60) -> int:
    _banner("cmp")
    da = json.load(io.open(a, encoding="utf-8"))
    db = json.load(io.open(b, encoding="utf-8"))
    diff = 0
    for name in sorted(set(da) | set(db)):
        x, y = da.get(name), db.get(name)
        if x == y:
            continue
        diff += 1
        if diff > limit:
            continue
        if x is None or y is None:
            print(f"\n* {name}: {'仅 A' if y is None else '仅 B'}")
            continue
        print(f"\n* {name}")
        for k in ("degree", "changed", "reacted", "pH", "net", "raw"):
            if x.get(k) != y.get(k):
                print(f"    {k}: {x.get(k)!r}\n      -> {y.get(k)!r}")
        if x.get("steps") != y.get("steps"):
            print(f"    steps: {len(x.get('steps') or [])} -> {len(y.get('steps') or [])} 条")
    print(f"\n[差异] {diff} 例（上限 {limit} 例明细）")
    return 0 if diff == 0 else 1


# --------------------------------------------------------------- anchor
_ANCHOR_PLAIN = "        v = 4471.0 / T_K - 6.09 + 0.0171 * T_K\n"
_ANCHOR_FIXED = ("        v = (4471.0 / T_K - 6.09 + 0.0171 * T_K\n"
                 "             - (4471.0 / 298.15 - 6.09 + 0.0171 * 298.15"
                 " - 14.0))\n")


def cmd_anchor(mode: str) -> int:
    p = os.path.join(ROOT, "chemkit", "core.py")
    s = io.open(p, encoding="utf-8", newline="").read()
    # 行尾风格跟随文件（避免 CRLF/LF 噪声，见 hygiene）
    nl = "\r\n" if "\r\n" in s else "\n"
    plain = _ANCHOR_PLAIN.replace("\n", nl)
    fixed = _ANCHOR_FIXED.replace("\n", nl)
    state = "on" if fixed in s else ("off" if plain in s else "?")
    print(f"[锚定] 当前 = {state}")
    if mode == "status" or state == "?":
        if state == "?":
            print("!! core.py 里的 pKw 经验式既非锚定形也非未锚定形，请手工检查")
            return 2
        _banner("anchor status")
        return 0
    want_on = mode == "on"
    if (state == "on") == want_on:
        print(f"[锚定] 已是 {mode}，无改动（幂等）")
    else:
        src, dst = (plain, fixed) if want_on else (fixed, plain)
        assert s.count(src) == 1, "锚定行匹配数 != 1"
        io.open(p, "w", encoding="utf-8", newline="").write(s.replace(src, dst))
        print(f"[锚定] 已切换 -> {mode}")
    _banner(f"anchor {mode}")
    return 0


# --------------------------------------------------------------- eqcheck
def cmd_eqcheck() -> int:
    _banner("eqcheck")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import eqcheck                                     # type: ignore
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = eqcheck.main([])
    tail = [ln for ln in buf.getvalue().split("\n") if ln.strip()][-3:]
    print("\n".join(tail))
    print("[判定] 全库 eq/eq_has 精确守恒" if rc == 0 else "[判定] 存在违规")
    return rc


# --------------------------------------------------------------- hygiene
def cmd_hygiene(fix: bool = False) -> int:
    _banner("hygiene")
    st = _run(["git", "status", "--porcelain"]).split("\n")
    mod = [ln[3:].strip() for ln in st if ln[:2] == " M"]
    unt = [ln[3:].strip() for ln in st if ln[:2] == "??"]
    eol_only = []
    for f in mod:
        if not os.path.exists(f):
            continue
        # `git diff --quiet` 干净 = 无**内容**差异（工作区只是换行风格被改写）
        if subprocess.run(["git", "diff", "--quiet", "--", f],
                          cwd=ROOT).returncode == 0:
            eol_only.append(f)
    scratch = [f for f in unt
               if f.startswith("_") and f.split(".")[-1] in ("py", "txt", "json")]
    print(f"[改动] {len(mod)} 个已跟踪文件；[未跟踪] {len(unt)} 个")
    if eol_only:
        print(f"[行尾噪声] {len(eol_only)} 个文件 git diff 为空（仅换行风格差异）：")
        for f in eol_only:
            print(f"    {f}")
        if fix:
            for f in eol_only:
                _run(["git", "checkout", "--", f])
            print("[行尾噪声] 已 git checkout 还原")
        else:
            print("    （加 --fix 一键还原，避免提交噪声）")
    else:
        print("[行尾噪声] 无")
    if scratch:
        print(f"[临时文件] {len(scratch)} 个未被 .gitignore 覆盖：")
        for f in scratch[:12]:
            print(f"    {f}")
    else:
        print("[临时文件] 无（`_*.py`/`_*.txt`/`_*.json` 建议进 .gitignore）")
    return 0


# --------------------------------------------------------------- patch
def cmd_patch(spec: str, check: bool = False) -> int:
    _banner("patch")
    if not os.path.exists(spec):
        print(f"!! 找不到 spec：{spec}")
        print("spec 格式：PATCHES = [(相对路径, 旧串, 新串[, 期望出现次数]), ...]")
        return 2
    g: dict = {}
    exec(compile(io.open(spec, encoding="utf-8").read(), spec, "exec"), g)
    patches = g.get("PATCHES")
    if not patches:
        print("!! spec 里没有 PATCHES")
        return 2
    # 先全量校验（原子性：任一条不满足则一行都不落盘）
    # 换行风格统一：spec 通常是 LF，而 Windows 检出/`write` 产出的文件可能是
    # CRLF——旧串按**归一化文本**匹配（匹配、计数、替换都在归一化文本上做，
    # 落盘时再还原原风格），免去"旧串出现 0 次"这类与内容无关的失败。
    staged: dict[str, str] = {}
    nl_of: dict[str, str] = {}
    for item in patches:
        path, old, new = item[0], item[1], item[2]
        n = item[3] if len(item) > 3 else 1
        full = os.path.join(ROOT, path)
        raw = staged.get(path)
        if raw is None:
            raw = io.open(full, encoding="utf-8", newline="").read()
            nl_of[path] = "\r\n" if "\r\n" in raw else "\n"
        crlf = nl_of[path] == "\r\n"
        cur = raw.replace("\r\n", "\n")
        old_n = old.replace("\r\n", "\n")
        new_n = new.replace("\r\n", "\n")
        got = cur.count(old_n)
        if got != n:
            print(f"!! {path}: 旧串出现 {got} 次，期望 {n} 次 -> 中止（未写入任何文件）")
            return 1
        nxt = cur.replace(old_n, new_n)
        staged[path] = nxt.replace("\n", "\r\n") if crlf else nxt
        print(f"   ok {path}: {got} 处替换（{len(old)}B -> {len(new)}B）"
              + ("  [CRLF]" if crlf else ""))
    if check:
        print("[--check] 校验通过，未写入")
        return 0
    for path, text in staged.items():
        io.open(os.path.join(ROOT, path), "w", encoding="utf-8",
                newline="").write(text)
    print(f"[写入] {len(staged)} 个文件")
    return 0


# --------------------------------------------------------------- run
def cmd_run(script: str, args: list[str]) -> int:
    """在 UTF-8 控制台环境下跑任意脚本（探针 / 审计工具 / 补丁 spec）。

    存在的理由：Windows 控制台默认 GBK，临时探针脚本里一个 `H⁺` 就能让
    整个脚本以 UnicodeEncodeError 崩掉（本工具自身已重配置 stdout，但
    被 `python xxx.py` 直接跑的脚本没有）。统一从这里跑，这类错误消失。
    """
    _banner(f"run {script}")
    if not os.path.exists(script):
        print(f"!! 找不到脚本：{script}")
        return 2
    import runpy
    argv0 = sys.argv
    sys.argv = [script] + args
    try:
        runpy.run_path(script, run_name="__main__")
    except SystemExit as e:                              # 脚本自带退出码
        return int(e.code or 0)
    finally:
        sys.argv = argv0
    return 0


# --------------------------------------------------------------- guide
_GUIDE = """\
本轮协议（固定动作，别即兴）：
  1) python tools/dev.py suite                 # 起点必须绿（记下 N/N）
  2) 只做**一处**根治性修改（优先 `edit` 工具；结构性改动写 spec 走 patch）
  3) python tools/dev.py case <受影响用例>      # 化学对账：步表/账本/两版净方程/断言
  4) python tools/dev.py suite                 # 回归；再 perf 看确定性指标
  5) python tools/dev.py eqcheck && python tools/dev.py hygiene --fix
  6) 文档（architecture §7 + changelog）→ 提交推送 → 汇报性能表
禁令（都被坑过）：
  · 不用 shell heredoc / 内联 `python -c` 拼多行与引号：写文件再跑；
  · 不用 shell 的 grep/find/type：用 read/grep/glob 工具；
  · 控制台只打 ASCII 记号（方程用 .plain()），TeX 串留给 Markdown；
  · 长输出一律截断（Select-Object -Last N / 本工具自带上限）；
  · 断言是化学事实的门槛：**先裁决化学**，再决定改引擎还是改标准；
  · 前后对比用 snapshot/cmp（不要 git stash，容易丢改动）。"""


def cmd_guide() -> int:
    _banner("guide")
    print(_GUIDE)
    return 0


def main(argv: list[str]) -> int:
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(__doc__)
        return 0
    verb, rest = argv[0], argv[1:]
    flags = {a for a in rest if a.startswith("-")}
    args = [a for a in rest if not a.startswith("-")]
    if verb == "guide":
        return cmd_guide()
    if verb == "run":
        return cmd_run(args[0], args[1:])
    if verb == "suite":
        return cmd_suite(args)
    if verb == "perf":
        base = (args[0] if args and args[0].endswith(".json")
                else "converg-baseline.json")
        return cmd_perf(base, write="--write" in flags)
    if verb == "case":
        return cmd_case(args, show_all="--all" in flags)
    if verb == "snapshot":
        return cmd_snapshot(args[0], args[1:])
    if verb == "cmp":
        return cmd_cmp(args[0], args[1])
    if verb == "anchor":
        return cmd_anchor(args[0] if args else "status")
    if verb == "eqcheck":
        return cmd_eqcheck()
    if verb == "hygiene":
        return cmd_hygiene(fix="--fix" in flags)
    if verb == "patch":
        return cmd_patch(args[0], check="--check" in flags)
    print(f"!! 未知子命令 {verb}\n")
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
