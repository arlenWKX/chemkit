# -*- coding: utf-8 -*-
"""`tools/jsondb.py` —— **标准化 JSON 读写**（数据表统一口径，第 147 轮）。

## 为什么（历史痛点）
`chemkit/data/*.json` 原先各自手工排版，`json.load → 改 → json.dump` 不能逐字节还原
⟹ 改数据只能靠**文本锚串**，于是反复踩：锚串跨用例撞车（第 137 轮 4 处）、
锚串尾部换行有无（第 147 轮 3 处）、JSON 里混入 ASCII 引号、Python 三引号结尾
连写四个引号造成歧义……这些坑与化学无关，纯粹是"排版不是数据"造成的。

**使用者已明确：可以接受数据库彻底格式化** ⟹ 从此按**统一口径**写：
`json.dumps(data, ensure_ascii=False, indent=1)` + 末尾换行；换行风格沿用原文件
（默认 CRLF，本仓库 Windows 检出即 CRLF）。所有数据表一次 `normalize` 之后，
后续按**键路径**改动都是可复现、可 diff 的。

## API
```python
doc = JsonDoc("chemkit/data/tests.json")
i = doc.case_index("F13 Ba+水")          # tests.json 专用：按用例名取下标
doc.get([i, "ph"])                        # 读（同 data[i]["ph"]）
doc.set([i, "ph"], [13.0, 14.5])          # 写
doc.pop([i, "eq"])                        # 删键
doc.del_value([i, "eq_has"], "…")         # 按值删数组元素
doc.add([i], "note", "…", after="changed")# 加键
doc.dump()                                # 落盘（写前 json.dumps 复验）
```

## CLI
```powershell
python tools/jsondb.py normalize chemkit/data/*.json      # 统一格式化（可多文件）
python tools/jsondb.py check chemkit/data/beta.json
python tools/jsondb.py get chemkit/data/tests.json 12 ph
python tools/jsondb.py set chemkit/data/tests.json 12 ph "[3.8, 5.5]"
```
"""
import argparse
import glob
import io
import json
import os
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

INDENT = 1


def dumps(data, nl: str = "\r\n") -> str:
    """统一口径：ensure_ascii=False + indent=1 + 末尾换行。"""
    out = json.dumps(data, ensure_ascii=False, indent=INDENT) + "\n"
    return out if nl == "\n" else out.replace("\n", "\r\n")


def load(path: str):
    with io.open(path, encoding="utf-8") as fh:
        return json.load(fh)


def save(path: str, data, nl: str | None = None) -> None:
    if nl is None:
        raw = io.open(path, encoding="utf-8", newline="").read()
        nl = "\r\n" if "\r\n" in raw else "\n"
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(dumps(data, nl))


class JsonDoc:
    """按路径改 JSON（内存树 + 统一口径落盘）。"""

    def __init__(self, path: str | None = None, data=None):
        self.path = path
        self.nl = "\n"
        if data is not None:
            self.data = data
            return
        raw = io.open(path, encoding="utf-8", newline="").read()
        self.nl = "\r\n" if "\r\n" in raw else "\n"
        self.data = json.loads(raw)

    # ---- 结构导航 ----
    def _at(self, path):
        v = self.data
        for k in path:
            v = v[k]
        return v

    def get(self, path):
        return self._at(path)

    def case_index(self, name: str) -> int:
        """tests.json 专用：按用例名取顶层下标（要求恰好命中一次）。"""
        if not isinstance(self.data, list):
            raise ValueError("根不是列表")
        hits = [i for i, c in enumerate(self.data)
                if isinstance(c, dict) and c.get("name") == name]
        if len(hits) != 1:
            raise KeyError(f"用例名 {name!r} 命中 {len(hits)} 次")
        return hits[0]

    # ---- 改 ----
    def set(self, path, value):
        *head, last = tuple(path)
        self._at(head)[last] = value

    def pop(self, path):
        *head, last = tuple(path)
        del self._at(head)[last]

    def del_value(self, path, value):
        arr = self._at(path)
        arr.remove(value)

    def add(self, path, key, value, after=None):
        obj = self._at(path)
        if key in obj:
            raise KeyError(f"键已存在：{key}")
        if after is None:
            obj[key] = value
            return
        items = list(obj.items())
        obj.clear()
        for k, v in items:
            obj[k] = v
            if k == after:
                obj[key] = value
        if key not in obj:
            raise KeyError(f"after 键不存在：{after}")

    def set_or_add(self, path, key, value):
        if key in self._at(path):
            self.set(list(path) + [key], value)
        else:
            self.add(path, key, value)

    # ---- 落盘 ----
    def dump(self, path: str | None = None):
        json.dumps(self.data, ensure_ascii=False)      # 复验
        save(path or self.path, self.data, self.nl)


def _path_from_args(args):
    out = []
    for a in args:
        out.append(int(a) if a.lstrip("-").isdigit() else a)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("cmd", choices=["check", "normalize", "get", "set", "del",
                                    "delval", "add"])
    ap.add_argument("file", nargs="+")
    ap.add_argument("--after", default=None)
    a = ap.parse_args()

    if a.cmd == "normalize":
        files = []
        for pat in a.file:
            files.extend(glob.glob(pat) or [pat])
        for f in files:
            data = load(f)
            before = os.path.getsize(f)
            save(f, data)
            print(f"[normalize] {f}: {before} -> {os.path.getsize(f)} B")
        return 0
    if a.cmd == "check":
        for f in a.file:
            data = load(f)
            print(f"[OK] {f}: JSON 合法（{len(json.dumps(data, ensure_ascii=False))} 字符）")
        return 0

    f = a.file[0]
    rest = a.file[1:]
    doc = JsonDoc(f)
    if a.cmd == "get":
        print(json.dumps(doc.get(_path_from_args(rest)), ensure_ascii=False,
                         indent=INDENT))
        return 0
    if a.cmd == "set":
        doc.set(_path_from_args(rest[:-1]), json.loads(rest[-1]))
        doc.dump()
        print(f"[写入] {f}")
        return 0
    if a.cmd == "del":
        doc.pop(_path_from_args(rest))
        doc.dump()
        print(f"[写入] {f}")
        return 0
    if a.cmd == "delval":
        doc.del_value(_path_from_args(rest[:-1]), json.loads(rest[-1]))
        doc.dump()
        print(f"[写入] {f}")
        return 0
    if a.cmd == "add":
        doc.add(_path_from_args(rest[:-2]), rest[-2], json.loads(rest[-1]),
                after=a.after)
        doc.dump()
        print(f"[写入] {f}")
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
