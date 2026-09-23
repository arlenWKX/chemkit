# -*- coding: utf-8 -*-
"""第 226 轮 · 量化「**同一作用域**内同一含义多名字」——真实危害面。

第 205 轮的 `tools/namescan.py` 统计的是**文件级**名字数（pH 12 个、He 9 个…），
但那是**风格总量**，不是危害面：同一文件里出现 `pH` 与 `pH_x` 很正常，
真正会让人改漏的是**同一函数内**用多个名字装同一个量。

判别式：在**同一个函数/模块作用域**内，属于同一含义簇的不同名字 ≥ 2
（且**都在该作用域被绑定**）。这类地方改一处极易漏另一处。

用法： python tools/scope_names.py [包目录，默认 chemkit]
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

# 含义簇（与 namescan.py 同源，只保留判别力强的）
CLUSTERS = {
    "pH": ("ph", "ph0", "ph_x", "ph_f", "ph_solver", "_ph", "_ph0", "ph_est",
           "_cls", "_cls0", "_cls_s", "_ex", "ph_res", "ph_buf", "ph_before",
           "ph_after", "ph_init"),
    "He(质子过剩)": ("he", "he0", "he_res", "he_raw", "he_x", "he_ok", "h_excess",
                   "h_excess0", "_he", "_he0", "_he_ok", "_he_res", "hexcess"),
    "离子强度I": ("i", "i_eff", "ie", "_ie", "ionic", "mu", "mu_x", "i_tot"),
    "电荷z": ("z", "zc", "_z", "_zc", "charge", "chg", "_chg", "net", "_net"),
    "体积V": ("v", "v_l", "vol", "volume", "vl"),
    "温度T": ("t_k", "tk", "t_c", "tc", "temp", "t"),
    "浓度/量": ("c", "conc", "n", "mol", "moles", "amt", "amount", "m", "c0"),
    "物种名": ("sp", "s", "spec", "species", "name", "nm", "key", "k"),
    "候选": ("c", "cd", "cand", "cands", "pick", "cc"),
    "程度x": ("x", "ext", "extent", "xj", "xi", "x_max", "xmax"),
}


def bound_in_scope(node):
    """收集该作用域**直接绑定**的名字（不递归进嵌套函数）。

    ⚠️ 第 226 轮修：第一版在遇到嵌套 `FunctionDef` 时 `continue`，
    但**仍会继续遍历它的 `args`/`decorator_list`**（那是 `iter_child_nodes`
    的其它子节点）⟹ 把嵌套函数的**参数名**也算进外层作用域，
    于是报出"`solve_extent` 里同时有 `_cls` 与 `_ph0`"这种**不存在**的冲突
    （`_cls` 属于它内部的 `_ph_closed` 闭包所在分支）。
    **诊断工具自身出错会伪装成代码缺陷**——这正是 lessons 里
    "诊断工具说数据缺先怀疑工具"的同类。现在遇到嵌套 def/class **整棵跳过**。
    """
    names = set()

    def walk(n, top=True):
        for ch in ast.iter_child_nodes(n):
            if isinstance(ch, (ast.FunctionDef, ast.AsyncFunctionDef,
                               ast.ClassDef)):
                if top:
                    names.add(ch.name)      # 只记名字，**不进它的子节点**
                continue
            if isinstance(ch, ast.Name) and isinstance(ch.ctx, ast.Store):
                names.add(ch.id)
            elif isinstance(ch, ast.arg):
                names.add(ch.arg)
            elif isinstance(ch, ast.alias):
                names.add((ch.asname or ch.name).split(".")[0])
            walk(ch, top=False)
    walk(node)
    return names


files = []
for dp, dn, fn in os.walk(os.path.join(ROOT, PKG)):
    dn[:] = [d for d in dn if d != "__pycache__"]
    for f in sorted(fn):
        if f.endswith(".py"):
            files.append(os.path.join(dp, f))

hits = []
for path in files:
    rel = os.path.relpath(path, ROOT)
    tree = ast.parse(io.open(path, encoding="utf-8").read())
    scopes = [(rel, "<module>", bound_in_scope(tree))]
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            scopes.append((rel, node.name, bound_in_scope(node)))
    for rel_, sname, names in scopes:
        low = {n.lower() for n in names}
        for cl, pats in CLUSTERS.items():
            found = sorted(n for n in names if n.lower() in pats)
            if len(found) >= 2:
                hits.append((rel_, sname, cl, found))

print(f"扫描 {len(files)} 文件\n")
print("=" * 76)
print("**同一作用域内同一含义 ≥2 个名字**（真实改漏危害面）")
print("=" * 76)
by_cl = defaultdict(list)
for rel_, sname, cl, found in hits:
    by_cl[cl].append((rel_, sname, found))
tot = 0
for cl in CLUSTERS:
    rows = by_cl.get(cl) or []
    if not rows:
        continue
    tot += len(rows)
    print(f"\n### {cl}：{len(rows)} 个作用域")
    for rel_, sname, found in rows[:10]:
        print(f"  {rel_}:{sname}  ->  {', '.join(found)}")
    if len(rows) > 10:
        print(f"  … 另 {len(rows) - 10} 处")
print(f"\n合计 {tot} 个作用域存在同义多名。")

print("\n=== 判读 ===")
print("  · 若计数小（几十）⟹ 可逐个收敛，风险可控。")
print("  · 若某簇占比高 ⟹ 优先收敛该簇；`pH`/`He` 是历史最久的两个。")
print("  · 与 namescan.py 的文件级数字对照：文件级是**风格总量**，")
print("    本脚本是**危害面**——两者相差越大，说明名字虽多但分布在各处、")
print("    并不集中，收敛的紧迫性越低。")
