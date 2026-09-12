"""tests.json 重复项审查/合并（修正标准前的第一步）。

三类重复，逐类不同处置：
  ① **同名**（name 完全相同）：几乎必是复制粘贴残留，合并为一个
     （断言取并集；冲突的数值断言并列报出交由人工裁决）
  ② **同实验**（subs + cond 相同、name 不同）：同一实验的多条用例——是
     真正的冗余还是有意分组（如 J06/J14 同物同温却给不同范围），必须
     人工判断，故只报不合
  ③ **同实验同断言**：完全冗余，可直接合并

用法：
    python tools/dedupe_tests.py            # 只报告
    python tools/dedupe_tests.py --write-同名   # 合并同名项并写回
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict

PATH = "chemkit/data/tests.json"
_ASSERT = ("degree", "eq", "eq_has", "has", "has_not", "has_range", "has_any",
           "has_initial", "ph", "changed", "reacted", "ann", "override")


def sig(c: dict) -> str:
    return json.dumps([c["subs"], c.get("cond")], sort_keys=True,
                      ensure_ascii=False)


def asig(c: dict) -> str:
    return json.dumps({k: c[k] for k in _ASSERT if k in c}, sort_keys=True,
                      ensure_ascii=False)


def _num_bounds(d) -> dict:
    return {k: v for k, v in (d or {}).items()}


def merge_pair(a: dict, b: dict) -> tuple[dict, list[str]]:
    """同实验合并：**取更严者**（合并后不弱于任一原条目）。

    实测 43 组同实验重复全是「严版 + 宽版」两套断言（不同批次生成）。
    故逐字段取更严：
      has      下界 → 逐物种取 max
      has_not  上界 → 逐物种取 min
      has_range/ph → 区间取交（空集则报冲突，保留前者）
      degree   允许集 → 取交（None = 任意，不参与）
      eq       → 任一为空取非空；都非空且不等则报冲突
      eq_has/has_initial/ann → 并集（都是"必须出现"清单）
      has_any  → 取交（空集退并集并报冲突）
      changed/reacted/override → 非空优先，取值不同报冲突
      note     → 拼接
    任何冲突都打印出来交人工裁决，不静默丢弃。
    """
    out = dict(a)
    bad: list[str] = []

    def _take(key):
        va, vb = a.get(key), b.get(key)
        if va is None:
            if vb is not None:
                out[key] = vb
        elif vb is not None and vb != va:
            bad.append(f"    {key}: {va!r} vs {vb!r}（保留前者）")
    for key in ("changed", "reacted", "override", "degree", "eq"):
        if key in ("changed", "reacted"):
            # 显式 null 会被 testsuit 当成"断言为 None"（`"reacted" in c` 判定），
            # 故两侧都为 None 时**不能写出该键**——首版合并把 4 例合并项写成
            # `reacted: null` 造成 "reacted=True 期望None" 假失败。
            v = a.get(key) if a.get(key) is not None else b.get(key)
            if v is None:
                out.pop(key, None)
            else:
                out[key] = v
            continue
        if key == "degree":
            da, db = a.get("degree"), b.get("degree")
            la = None if da is None else (da if isinstance(da, list) else [da])
            lb = None if db is None else (db if isinstance(db, list) else [db])
            if la is None:
                if lb is not None:
                    out["degree"] = db
            elif lb is not None:
                inter = [x for x in la if x in lb]
                if inter:
                    out["degree"] = inter
                else:
                    bad.append(f"    degree: {la} vs {lb} 无交集（保留前者）")
            continue
        _take(key)

    for key in ("has", "has_not"):
        va, vb = _num_bounds(a.get(key)), _num_bounds(b.get(key))
        if not va and not vb:
            continue
        m = dict(va)
        for sp, v in vb.items():
            if sp not in m:
                m[sp] = v
            else:
                m[sp] = max(m[sp], v) if key == "has" else min(m[sp], v)
        out[key] = m

    for key in ("has_range",):
        va, vb = a.get(key) or {}, b.get(key) or {}
        if not va or not vb:
            if vb and not va:
                out[key] = vb
            continue
        m = {}
        for sp in set(va) & set(vb):
            lo = max(va[sp][0], vb[sp][0])
            hi = min(va[sp][1], vb[sp][1])
            if lo > hi:
                bad.append(f"    {key}[{sp}]: {va[sp]} ∩ {vb[sp]} = 空"
                           f"（保留前者 {va[sp]}）")
                m[sp] = va[sp]
            else:
                m[sp] = [lo, hi]
        for sp in set(va) ^ set(vb):
            m[sp] = (va.get(sp) or vb.get(sp))
        out[key] = m

    pa, pb = a.get("ph"), b.get("ph")
    if pa and pb:
        lo, hi = max(pa[0], pb[0]), min(pa[1], pb[1])
        if lo > hi:
            bad.append(f"    ph: {pa} ∩ {pb} = 空（保留前者）")
        else:
            out["ph"] = [lo, hi]
    elif pb:
        out["ph"] = pb

    for key in ("eq_has", "has_initial", "ann"):
        va, vb = a.get(key) or [], b.get(key) or []
        if va or vb:
            seen, u = set(), []
            for x in list(va) + list(vb):
                k = json.dumps(x, sort_keys=True, ensure_ascii=False)
                if k not in seen:
                    seen.add(k)
                    u.append(x)
            out[key] = u
    ha, hb = a.get("has_any") or {}, b.get("has_any") or {}
    if ha and hb:
        inter = {k: max(ha[k], hb[k]) for k in set(ha) & set(hb)}
        if inter:
            out["has_any"] = inter
        else:
            bad.append(f"    has_any: {ha} vs {hb} 无交集（保留前者）")
    elif hb:
        out["has_any"] = hb

    na, nb = (a.get("note") or "").strip(), (b.get("note") or "").strip()
    if nb and nb not in na:
        out["note"] = (na + "；" + nb).strip("；") if na else nb
    return out, bad


def _detect_indent(text: str) -> int:
    for line in text.splitlines():
        if line.strip() == "{":
            return len(line) - len(line.lstrip(" "))
    return 1


def main(argv: list[str]) -> None:
    raw = open(PATH, encoding="utf-8", newline="").read()
    ind = _detect_indent(raw)
    rows = json.loads(raw)
    print(f"tests.json 共 {len(rows)} 条")

    by_name = defaultdict(list)
    for c in rows:
        by_name[c["name"]].append(c)
    dup_name = {k: v for k, v in by_name.items() if len(v) > 1}

    by_sig = defaultdict(list)
    for c in rows:
        by_sig[sig(c)].append(c)
    dup_exp = {k: v for k, v in by_sig.items() if len(v) > 1}

    by_asig = defaultdict(list)
    for c in rows:
        by_asig[(sig(c), asig(c))].append(c)
    dup_full = {k: v for k, v in by_asig.items() if len(v) > 1}

    print(f"\n① 同名重复 {len(dup_name)} 组：")
    for k, v in list(dup_name.items())[:30]:
        print(f"   ×{len(v)}  {k}")
    print(f"\n② 同实验（subs+cond 相同、名称不同）{len(dup_exp)} 组：")
    for k, v in list(dup_exp.items())[:30]:
        names = [x["name"] for x in v]
        same_as = len({asig(x) for x in v}) == 1
        print(f"   {'[断言也相同]' if same_as else '[断言不同]'} "
              + " | ".join(n[:38] for n in names))
    print(f"\n③ 同实验且同断言（完全冗余）{len(dup_full)} 组：")
    for _k, v in list(dup_full.items())[:20]:
        print("   " + " | ".join(x["name"][:40] for x in v))

    if "--write" not in argv:
        print("\n（--write 执行合并：同名与同实验两组都按「取更严者」合并）")
        return
    out, conflicts = [], []
    by_sig: dict = {}
    for c in rows:
        k = sig(c)
        if k in by_sig:
            i = by_sig[k]
            merged, bad = merge_pair(out[i], c)
            out[i] = merged
            conflicts += [f"  [{out[i]['name'][:44]}]" + (f" ← {c['name'][:30]}")
                          ] + bad
            continue
        by_sig[k] = len(out)
        out.append(dict(c))
    print(f"\n合并：{len(rows)} → {len(out)} 条（删去 {len(rows) - len(out)}）")
    print(f"需人工裁决的冲突 {sum(1 for m in conflicts if m.startswith('    '))} 条：")
    for m in conflicts:
        print(m)
    with open(PATH, "w", encoding="utf-8", newline="") as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=ind) + "\n")
    print(f"已写回 {PATH}（indent={ind}）")


if __name__ == "__main__":
    main(sys.argv[1:])
