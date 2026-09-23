# -*- coding: utf-8 -*-
"""第 203 轮 · "同一含义、多个变量名"普查（使用者指出的代码混乱主形态）。

判据：把代码里所有**标识符**（变量/参数/属性）按"物理含义簇"归类，找出
同一含义在**同一文件内**出现了多少个不同名字。名字越多，"这是同一个量"
就越难看出来，改动也越容易只改一处。

方法：AST 解析（不是文本 grep——注释与文档会污染计数），只统计
真实绑定名：赋值目标、函数参数、for 目标、with 目标、comprehension 变量、
局部/全局 def 名。

用法： python tools/namescan.py [包目录，默认 chemkit]
"""
import ast
import io
import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

PKG = sys.argv[1] if len(sys.argv) > 1 else "chemkit"

# 含义簇：键 = 人话描述；值 = 匹配该含义的名字**判别式**（小写后匹配）
CLUSTERS = {
    "pH（同一个 pH 值）": (
        "ph", "ph0", "ph_x", "ph_f", "ph_solver", "ph_pres", "ph_res",
        "ph_est", "phj", "ph_before", "ph_after", "ph_init", "ph_final",
        "_ph", "_ph0", "_cls", "_cls0", "_cls_s", "_ex", "ph_s", "phi"),
    "质子过剩 He / 酸量（同一个 H_excess）": (
        "he", "he0", "he_res", "he_raw", "he_x", "he_use", "he_ok",
        "hexcess", "h_excess", "h_excess0", "he_p", "he_f", "_he",
        "_he0", "_he_ok", "_he_res", "he_solver"),
    "离子强度 I / μ": ("i", "mu", "mu_x", "i_eff", "ie", "_ie", "ionic",
                     "is_", "i_tot", "mu_calc"),
    "电荷 / z": ("z", "zc", "zi", "_z", "_zc", "charge", "chg", "_chg",
                "net", "net_chg", "_net", "q", "qtot"),
    "体积 V": ("v", "v_l", "vol", "volume", "vl"),
    "温度 T_K / T_C": ("t_k", "tk", "t_c", "tc", "temp", "t"),
    "浓度 c / 摩尔量 n / 质量 m": ("c", "conc", "n", "mol", "moles", "m",
                            "amt", "amount", "c0"),
    "物种名 sp / s": ("sp", "s", "sp_", "spec", "species", "name", "nm",
                     "s_", "key", "k"),
    "候选 c / cd / cand": ("c", "cd", "cand", "cands", "pick", "cc"),
    "总浓度/总量 tot": ("tot", "total", "sum", "s_tot", "n_tot", "ctot"),
    "logK / 常数": ("logk", "lk", "log_k", "k_eff", "logkeff", "kb"),
    "饱和指数 S": ("s_idx", "si", "sat", "s_", "sf", "s_f", "sd"),
    "方向 d / sign": ("d", "dd", "dir", "direction", "sign", "sgn"),
    "程度 x / extent / ext": ("x", "ext", "extent", "xj", "xi", "x_max",
                            "xmax", "xstar", "x_star"),
}


def bound_names(tree):
    """返回该文件的 {名字: 出现次数}（只统计真实绑定）。"""
    cnt = defaultdict(int)
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            cnt[node.id] += 1
        elif isinstance(node, ast.arg):
            cnt[node.arg] += 1
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                               ast.ClassDef)):
            cnt[node.name] += 1
        elif isinstance(node, ast.Attribute) and isinstance(node.ctx,
                                                           ast.Store):
            cnt["." + node.attr] += 1
        elif isinstance(node, ast.alias):
            cnt[(node.asname or node.name).split(".")[0]] += 1
    return cnt


files = []
for dirpath, dirnames, filenames in os.walk(os.path.join(ROOT, PKG)):
    dirnames[:] = [d for d in dirnames if d != "__pycache__"]
    for fn in sorted(filenames):
        if fn.endswith(".py"):
            files.append(os.path.join(dirpath, fn))

print(f"扫描 {len(files)} 个文件（包 {PKG}）\n")
grand = defaultdict(lambda: defaultdict(int))
for path in files:
    rel = os.path.relpath(path, ROOT)
    try:
        tree = ast.parse(io.open(path, encoding="utf-8").read(), filename=path)
    except SyntaxError as exc:
        print(f"  [语法错误] {rel}: {exc}")
        continue
    cnt = bound_names(tree)
    for cluster, pats in CLUSTERS.items():
        hits = {n: c for n, c in cnt.items() if n.lower() in pats}
        if len(hits) >= 2:
            grand[cluster][rel] = hits

print("=" * 78)
print("① 同一含义在同一文件内的**不同名字数**（≥2 即报告；越多越危险）")
print("=" * 78)
for cluster in CLUSTERS:
    per = grand.get(cluster) or {}
    if not per:
        continue
    print(f"\n### {cluster}")
    rows = sorted(per.items(), key=lambda kv: -len(kv[1]))
    for rel, hits in rows[:6]:
        names = ", ".join(f"{n}×{c}" for n, c in
                          sorted(hits.items(), key=lambda kv: -kv[1]))
        print(f"  {rel:28s} {len(hits):2d} 个名字: {names}")

print("\n" + "=" * 78)
print("② 全包范围：每个含义簇出现过多少个不同名字")
print("=" * 78)
allnames = defaultdict(set)
for cluster in CLUSTERS:
    for rel, hits in (grand.get(cluster) or {}).items():
        allnames[cluster] |= set(hits)
    print(f"  {cluster:34s} {len(allnames[cluster]):3d} 个不同名字")

print("\n" + "=" * 78)
print("③ 一个名字身兼多义（最危险：改一处以为改完了）")
print("=" * 78)
# 注意：这里**不**做 lower()——`S` 与 `s`、`T` 与 `t` 在化学代码里
# 恰恰是两种含义（饱和指数 vs 物种名；温度 vs 局部计数）。
_amb = defaultdict(list)
for cl, pats in CLUSTERS.items():
    for p in pats:
        _amb[p].append(cl)
_amb2 = defaultdict(list)
for cl, pats in CLUSTERS.items():
    for p in pats:
        _amb2[p.lower()].append((p, cl))
for name in sorted(set(list(_amb) + list(_amb2)),
                   key=lambda n: -len(set(c for _p, c in _amb2[n.lower()]))):
    clusters_of = sorted({c for _p, c in _amb2[name.lower()]})
    if len(clusters_of) < 2:
        continue
    per = []
    for path in files:
        rel = os.path.relpath(path, ROOT)
        try:
            tr = ast.parse(io.open(path, encoding="utf-8").read())
        except SyntaxError:
            continue
        cc = bound_names(tr)
        if cc.get(name):
            per.append((rel, cc[name]))
    if not per:
        continue
    print(f"  '{name}' 身兼 {len(clusters_of)} 义:")
    for cl in clusters_of:
        print(f"      - {cl}")
    print("      绑定: " + ", ".join(f"{r}×{c}" for r, c in
                                  sorted(per, key=lambda kv: -kv[1])[:5]))

print("\n=== 判读 ===")
print("  同一含义名字越多，'这是同一个量'越难看出 ⟹ 修一处漏一处。")
print("  优先治理：名字数最多、且**跨文件**出现的那几簇（跨文件 = 语义分裂）。")
print("  ③ 的身兼多义是最高危项：同一个短名在同一文件里既当索引又当物理量，")
print("  静态读代码无法分辨，改名/改动必然漏。")
