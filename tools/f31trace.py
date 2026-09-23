# -*- coding: utf-8 -*-
"""第 209 轮 · 抓 F31 走步的**逐迭代轨迹**（为何 3~5 步就 no-cands）。

用引擎自带的 `CHEM_TRACE` 开关，把**全部** stdout 抓到 UTF-8 文件
（控制台是 GBK，中文会乱码；**不要**用 PowerShell 重定向，那会写 UTF-16）。

输出：`logs/f31trace.txt`
关注：每轮的 `[pick]`、`[ext]`、`[micro]`/`[micro-led]`（微步禁用点，
含禁用原因与 x_max）、以及最终的 `no-cands`。

用法： python tools/f31trace.py [用例前缀，默认 F31]
"""
import io
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

PRE = sys.argv[1] if len(sys.argv) > 1 else "F31"
runner = os.path.join(ROOT, "tools", "_f31trace_run.py")
code = '''# -*- coding: utf-8 -*-
import os, sys
sys.path.insert(0, r"%s")
import chemkit.engine as eng
from chemkit.data import load_tables
from chemkit.testsuit import load_cases, run_case
T = load_tables()
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith("%s ") or n == "%s"][0]
ok = run_case(c, T, verbose=False)
print("ASSERT:", "PASS" if ok else "FAIL")
''' % (ROOT, PRE, PRE)
with io.open(runner, "w", encoding="utf-8") as f:
    f.write(code)

out = os.path.join(ROOT, "logs", "f31trace.txt")
env = dict(os.environ)
env["CHEM_TRACE"] = "1"
env["PYTHONIOENCODING"] = "utf-8"
with io.open(out, "wb") as fh:
    r = subprocess.run([sys.executable, runner], stdout=fh,
                       stderr=subprocess.STDOUT, env=env, cwd=ROOT)
os.remove(runner)
print(f"exit={r.returncode}  ->  {os.path.relpath(out, ROOT)}")

txt = io.open(out, encoding="utf-8", errors="replace").read()
lines = txt.splitlines()
print(f"共 {len(lines)} 行")
print("\n--- 含 pick / ext / micro / no-cands / ASSERT 的行 ---")
keep = [ln for ln in lines
        if ("[pick]" in ln or "[ext]" in ln or "[micro" in ln
            or "no-cands" in ln or "ASSERT" in ln or "[joint" in ln
            or "[osc" in ln or "[freeze" in ln)]
for ln in keep[:60]:
    print("  " + ln[:200])
print(f"\n（过滤后 {len(keep)} 行；全文见 {os.path.relpath(out, ROOT)}）")
