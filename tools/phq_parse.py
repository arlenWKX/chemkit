# -*- coding: utf-8 -*-
"""第 184 轮·修正 2：保留缩进的 PHREEQC 解析器（-log_k 无空格也能识别）。"""
import io
import os
import re
import sys

sys.path.insert(0, ".")
for s in (sys.stdout, sys.stderr):
    try:
        s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

ROOT = os.path.dirname(os.path.abspath(__file__))


def find(fname):
    for base in ("thirdparty", "docs"):
        b = os.path.join(ROOT, base)
        if not os.path.isdir(b):
            continue
        for d in os.listdir(b):
            for cand in (os.path.join(b, d, "database", fname),
                         os.path.join(b, d, fname)):
                if os.path.exists(cand):
                    return cand
    return None


def parse(path):
    """{section: {name: {'logk':v,'dh':v}}}；只用**行首不缩进**判块首。"""
    out, sec, name = {}, None, None
    for raw in io.open(path, encoding="utf-8", errors="replace"):
        line = raw.rstrip("\n").rstrip("\r")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indented = line[:1] in (" ", "\t")
        s = line.strip()
        if not indented:
            if re.fullmatch(r"[A-Z_]{4,}", s):
                sec = s
                out.setdefault(sec, {})
                name = None
                continue
        if sec is None:
            continue
        low = s.lstrip("-").lower()
        if low.startswith("log_k") or low.startswith("delta_h"):
            if name is None:
                continue
            t = s.split()
            val = None
            for tok in t[1:]:
                try:
                    val = float(tok)
                    break
                except ValueError:
                    continue
            if val is not None:
                key = "logk" if low.startswith("log_k") else "dh"
                out[sec].setdefault(name, {"logk": None, "dh": None})[key] = val
            continue
        if s.startswith("-"):
            continue
        if indented:
            # PHASES 段的相名（缩进一行、无 '='），或块内续行
            if "=" not in s and sec == "PHASES":
                name = s.split()[0]
                out[sec].setdefault(name, {"logk": None, "dh": None})
            continue
        # 非缩进、非段名 ⟹ 块首
        nm = s.split("=")[0].strip().split()[0] if "=" in s else s.split()[0]
        name = nm
        out[sec].setdefault(nm, {"logk": None, "dh": None})
        if "=" in s:
            out[sec][nm]["reaction"] = s      # 记下反应式（跨库匹配要用）
    return out


if __name__ == "__main__":
    for f in ("phreeqc.dat", "sit.dat", "minteq.v4.dat", "llnl.dat"):
        p = find(f)
        if not p:
            print(f"{f}: 未找到")
            continue
        d = parse(p)
        sp = d.get("SOLUTION_SPECIES", {})
        ph = d.get("PHASES", {})
        print(f"\n=== {f} ===")
        print(f"  SOLUTION_SPECIES {len(sp)}（带 logk "
              f"{sum(1 for v in sp.values() if v.get('logk') is not None)}）"
              f"   PHASES {len(ph)}（带 logk "
              f"{sum(1 for v in ph.values() if v.get('logk') is not None)}）")
        for nm, v in list(sp.items())[:4]:
            print(f"     {nm:16s} logk={v.get('logk')}")
        for nm, v in list(ph.items())[:3]:
            print(f"     [相] {nm:14s} logk={v.get('logk')}")
