"""tools/dev.py —— **薄转发器**（第 303 轮）。

⚠️ 这里**不再有第二套测试实现**。唯一入口是：

    python -m chemkit.testsuit <verb> [args...]

本文件只是把它包成 `python tools/dev.py <verb>` 的老习惯写法，
把所有参数**逐字转发**过去，自身不含任何判定逻辑。

保留理由：历史命令与肌肉记忆（`dev.py suite` / `dev.py case ...`）；
删掉会让人不知所措，留着但明确它**不承载实现**。

    python tools/dev.py suite J05 NR23     ==  python -m chemkit.testsuit suite J05 NR23
    python tools/dev.py case 15            ==  python -m chemkit.testsuit case 15
    python tools/dev.py help               ==  python -m chemkit.testsuit help
"""
from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                    # pragma: no cover
        pass

from chemkit.testsuit import cli                          # noqa: E402

if __name__ == "__main__":
    sys.exit(cli(sys.argv[1:]))
