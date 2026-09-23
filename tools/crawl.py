# -*- coding: utf-8 -*-
"""第 264 轮 · 抓「族内步 × `charge_pH` 再分配」爬行的**逐迭代轨迹**。

`tools/f31trace.py` 只看**开头** 60 行（它是为"几步就 no-cands"设计的）；
本轮要的是**尾部爬行区**：`N15` 跑满 3000 轮，得看第 2900 轮附近在做什么。

输出 `logs/crawl-<前缀>.txt`（UTF-8），并在终端打印尾部窗口。

用法：
    python tools/crawl.py N15            # 尾部 120 行
    python tools/crawl.py N15 --tail 400
    python tools/crawl.py N15 --grep micro
"""
from __future__ import annotations

import io
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

KEY = ("[pick]", "[ext]", "[micro", "[osc", "[freeze", "[joint",
       "no-cands", "max-iter", "ASSERT")


def main(argv: list[str]) -> int:
    pre = "N15"
    tail = 120
    only = None
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--tail":
            i += 1
            tail = int(argv[i])
        elif a == "--grep":
            i += 1
            only = argv[i]
        elif not a.startswith("-"):
            pre = a
        i += 1

    runner = os.path.join(ROOT, "tools", "_crawl_run.py")
    code = '''# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, r"%s")
from chemkit.data import load_tables
from chemkit.testsuit import load_cases, run_case
T = load_tables()
c = [v for n, v in {x["name"]: x for x in load_cases(None)}.items()
     if n.startswith("%s")][0]
print("CASE:", c["name"], flush=True)
ok = run_case(c, T, verbose=False)
print("ASSERT:", "PASS" if ok else "FAIL")
''' % (ROOT, pre)
    with io.open(runner, "w", encoding="utf-8") as fh:
        fh.write(code)

    out = os.path.join(ROOT, "logs", f"crawl-{pre}.txt")
    env = dict(os.environ)
    env["CHEM_TRACE"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    with io.open(out, "wb") as fh:
        r = subprocess.run([sys.executable, runner], stdout=fh,
                           stderr=subprocess.STDOUT, env=env, cwd=ROOT)
    os.remove(runner)

    txt = io.open(out, encoding="utf-8", errors="replace").read()
    lines = txt.splitlines()
    keep = [ln for ln in lines if any(k in ln for k in KEY)]
    print(f"exit={r.returncode}  共 {len(lines)} 行，关键行 {len(keep)} 行"
          f"  -> {os.path.relpath(out, ROOT)}")
    sel = keep
    if only:
        sel = [ln for ln in keep if only in ln]
        print(f"（过滤 '{only}'：{len(sel)} 行）")
    print(f"\n--- 尾部 {tail} 行 ---")
    for ln in sel[-tail:]:
        print("  " + ln[:190])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
