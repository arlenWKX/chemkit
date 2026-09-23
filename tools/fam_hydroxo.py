# -*- coding: utf-8 -*-
"""fam_hydroxo.py -- FAMILY SCAN for the F31 (GaCl3+3NaOH) hydroxo-complex defect.

Question: is "strong base consumed as hydroxide LIGAND (hydroxo complex) instead
of precipitating the hydroxide SOLID" a general family defect across trivalent
metals, or is it specific to gallium?

Measurement only. Reads: chemkit.data.ksp / beta (solid + complex naming),
chemkit.testsuit (case list), chemkit.engine.judge(_probe=) (sanctioned
read-only diagnostic). Writes: logs/fam_hydroxo.json.

Buckets
  A solid_ok    : expected solid present in probe["ledger"] >= 0.5 * asserted lo
  B same_defect : solid ~0  AND  ledger dominated by a hydroxo complex of the
                  same metal  AND  some two-sided candidate with |S| > 3 whose
                  dis_why carries the engine's own "零推进" classification
  C other       : everything else (reason recorded)

Run:  python tools/fam_hydroxo.py
"""
import io
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import chemkit.engine as eng                                    # noqa: E402
from chemkit.core import charge_of, elements_of                 # noqa: E402
from chemkit.data import load_tables                            # noqa: E402
from chemkit.testsuit import load_cases                         # noqa: E402

X_MIN = 1e-6          # engine presence floor (mirrors chemkit.engine.X_MIN)
S_STRONG = 3.0        # "large |S|" threshold requested by the task
NEAR0_FRAC = 0.1      # solid below 10% of the asserted lower bound == "absent"
OK_FRAC = 0.5         # solid above 50% of the asserted lower bound == "formed"
TIME_BUDGET_S = 480.0

BASE_METALS = ("Na", "K", "Li", "Rb", "Cs", "Ba", "Sr", "Ca")
ZERO_PUSH = "\u96f6\u63a8\u8fdb"     # "零推进" -- engine's own classification tag


def _elems(name):
    try:
        return dict(elements_of(name))
    except Exception:
        return {}


def _metal_of_solid(solid):
    rest = sorted(e for e in _elems(solid) if e not in ("O", "H"))
    return rest[0] if len(rest) == 1 else None


def _oh_groups(name):
    """OH count for M(OH)_n / strong-base stoichiometry (O count == H count)."""
    els = _elems(name)
    if els.get("O") and els.get("O") == els.get("H"):
        return els["O"]
    return 0


def _is_strong_base(name):
    """NaOH / KOH / Ba(OH)2 ... but NOT Na[Al(OH)4] (contains a second metal)."""
    els = _elems(name)
    if not els.get("O") or els.get("O") != els.get("H"):
        return False
    if not any(b in els for b in BASE_METALS):
        return False
    return not (set(els) - set(BASE_METALS) - {"O", "H"})


def _why(a):
    w = a.get("dis_why")
    if not w:
        return ""
    if isinstance(w, (list, tuple)):
        return str(w[0]) if w else ""
    return str(w)


def _classify(probe, solid, metal, n_oh, lo):
    """Bucket + every task-4 field for one case."""
    led = probe.get("ledger") or {}
    m_solid = float(led.get(solid, 0.0))
    m_sp = {s: float(v) for s, v in led.items()
            if v > 0.0 and metal in _elems(s)}
    dissolved = {s: v for s, v in m_sp.items() if s != solid}
    free = {s: v for s, v in dissolved.items() if set(_elems(s)) == {metal}}
    hydroxo = {s: v for s, v in dissolved.items() if "(OH)" in s}
    tot_diss = sum(dissolved.values())
    dom_sp = max(hydroxo, key=lambda s: hydroxo[s]) if hydroxo else None
    dom_amt = hydroxo.get(dom_sp, 0.0) if dom_sp else 0.0
    free_sp = max(free, key=lambda s: free[s]) if free else None
    free_amt = free.get(free_sp, 0.0) if free_sp else 0.0

    active = list(probe.get("active") or [])
    tw = [a for a in active if a.get("two_sided")]
    tw_zero = [a for a in tw
               if abs(a.get("S", 0.0)) > S_STRONG and ZERO_PUSH in _why(a)]
    tw_sorted = sorted(tw, key=lambda a: -abs(a.get("S", 0.0)))
    mx = tw_sorted[0] if tw_sorted else None

    # candidates involving the expected solid; S normalised to the PRECIP direction
    precip = []
    for a in active:
        eq = a.get("eq") or ""
        if solid not in eq:
            continue
        left, sep, right = eq.partition(" -> ")
        if not sep:
            continue
        if solid in right:
            d, S = "precip", a.get("S")
        elif solid in left:
            d, S = "dissolve", -a.get("S")
        else:
            continue
        precip.append({"dir": d, "S": S, "eq": eq, "kind": a.get("kind"),
                       "two_sided": a.get("two_sided"),
                       "dis_why": _why(a), "ext_max": a.get("ext_max"),
                       "zero_push": ZERO_PUSH in _why(a)})
    undersat = [p for p in precip if p["dir"] == "precip" and p["S"] < 0.0]

    near0 = m_solid < max(NEAR0_FRAC * lo, X_MIN)
    dominated = (dom_sp is not None and dom_amt > free_amt and dom_amt > X_MIN
                 and dom_amt / max(tot_diss, X_MIN) >= 0.3)

    if not near0 and m_solid >= OK_FRAC * lo:
        bucket, why = "solid_ok", "solid formed"
    elif near0 and dominated and tw_zero:
        bucket = "same_defect"
        why = ("solid absent (%.6g); hydroxo complex %s = %.6g dominates "
               "dissolved metal (%.6g); %d two-sided |S|>%.1f candidate(s) "
               "tagged 零推进"
               % (m_solid, dom_sp, dom_amt, tot_diss, len(tw_zero), S_STRONG))
    else:
        bucket, bits = "other", []
        if not near0 and m_solid < OK_FRAC * lo:
            bits.append("solid only %.6g (%.1f%% of asserted lo %.6g)"
                        % (m_solid, 100.0 * m_solid / lo if lo else 0.0, lo))
        if near0:
            bits.append("solid ~0 (%.6g)" % m_solid)
            if not dominated:
                bits.append("no dominant hydroxo complex (dom=%s %.6g of "
                            "dissolved %.6g)" % (dom_sp, dom_amt, tot_diss))
            if not tw_zero:
                bits.append("no two-sided |S|>%.1f 零推进 candidate (max |S| = %s)"
                            % (S_STRONG,
                               ("%.3f" % abs(mx["S"])) if mx else "n/a"))
        why = "; ".join(bits) or "does not match either signature"

    return {
        "bucket": bucket, "bucket_why": why,
        "solid_expect": solid, "solid_n_oh": n_oh, "solid_assert_lo": lo,
        "solid_mol": round(m_solid, 9),
        "solid_frac_of_lo": (round(m_solid / lo, 6) if lo else None),
        "free_metal_sp": free_sp, "free_metal_mol": round(free_amt, 9),
        "dom_hydroxo_sp": dom_sp, "dom_hydroxo_mol": round(dom_amt, 9),
        "hydroxo_species": {s: round(v, 9) for s, v in
                            sorted(hydroxo.items(), key=lambda kv: -kv[1])},
        "dissolved_metal_total": round(tot_diss, 9),
        "dom_hydroxo_frac_of_dissolved": (round(dom_amt / tot_diss, 6)
                                          if tot_diss > 0 else None),
        "max_twosided_absS": (round(abs(mx["S"]), 3) if mx else None),
        "max_twosided": ({"kind": mx.get("kind"), "eq": mx.get("eq"),
                          "S": mx.get("S"), "ext_max": mx.get("ext_max"),
                          "dis_why": _why(mx), "zero_push": ZERO_PUSH in _why(mx)}
                         if mx else None),
        "n_twosided_absS_gt3_zero_push": len(tw_zero),
        "twosided_zero_push": [
            {"kind": a.get("kind"), "eq": a.get("eq"), "S": a.get("S"),
             "ext_max": a.get("ext_max")} for a in
            sorted(tw_zero, key=lambda a: -abs(a["S"]))],
        "precip_candidates": precip,
        "undersat_precip": undersat[0] if undersat else None,
        "n_undersat_precip": len(undersat),
    }


def main():
    t0 = time.time()
    T = load_tables()
    cases = load_cases(None)

    hydrox, oxide = {}, {}
    for e in T.ksp:
        s = e["solid"]
        m = _metal_of_solid(s)
        if m is None:
            continue
        els = set(_elems(s))
        if "(OH)" in s:
            hydrox[s] = {"metal": m, "n_oh": _oh_groups(s), "pKsp": e["pKsp"]}
        elif els == {m, "O"}:
            oxide[s] = {"metal": m, "pKsp": e["pKsp"]}
    print("solids in ksp.json: %d total | %d M(OH)n | %d MxOy"
          % (len(T.ksp), len(hydrox), len(oxide)))

    # --- data table: M(OH)n solid vs its hydroxo complexes -----------------
    print("\n--- DATA: M(OH)n solid vs its hydroxo complexes (beta.json) ---")
    print("%-10s %-5s %-7s %-22s %-4s %-8s %-11s %s" %
          ("solid", "metal", "pKsp", "complex", "nu", "logbeta", "logK_ampho",
           "OH/step"))
    ampho = []
    for s, info in sorted(hydrox.items(), key=lambda kv: kv[1]["metal"]):
        M, n = info["metal"], info["n_oh"]
        if not n:
            continue
        cands = []
        for b in T.beta:
            if b.get("ligand") != "OH^-":
                continue
            ce = _elems(b.get("center"))
            if set(ce) == {M} and charge_of(b["center"]) == n:
                cands.append(b)
        if not cands:
            continue
        hit = max(cands, key=lambda b: b.get("nu", 0))
        nu = int(hit.get("nu", 0))
        lg = float(hit["logb"])
        k = lg - float(info["pKsp"])          # per (nu-n) OH consumed
        step = nu - n
        ampho.append({"solid": s, "metal": M, "n": n, "pKsp": info["pKsp"],
                      "complex": hit["complex"], "nu": nu, "logb": lg,
                      "logK_ampho": round(k, 2), "oh_per_step": step,
                      "per_oh": round(k / step, 2) if step else None})
        print("%-10s %-5s %-7s %-22s %-4s %-8s %-11s %s" %
              (s, M, info["pKsp"], hit["complex"], nu, lg, round(k, 2),
               ("%d (logK/OH=%s)" % (step, round(k / step, 2))) if step else "-"))
    print("  logK_ampho = logbeta - pKsp = logK of M(OH)n(s) + (nu-n)OH^- ->"
          " [M(OH)nu]^(nu-n)-")

    # --- F31 data reconciliation ------------------------------------------
    try:
        from chemkit.core import pKw_of as _pKw
        from chemkit.speciation import ionic_strength, sit_fixpoint
        pk = _pKw(298.15)
        f31 = None
        for rec0 in cases:
            if rec0["name"].startswith("F31"):
                f31 = rec0
                break
        print("\n--- RECONCILIATION: F31 / [Ga(OH)_4]^- ---")
        print("  pKw(298.15) = %.6g" % pk)
        for b in T.beta:
            if b.get("center") == "Ga^{3+}" and b.get("ligand") == "OH^-":
                print("  beta.json: %-14s nu=%s logb=%s"
                      % (b["complex"], b["nu"], b["logb"]))
        print("  ksp.json : Ga(OH)_3 pKsp=%s"
              % T.ksp_by_solid["Ga(OH)_3"]["pKsp"])
        if f31 is not None:
            led = {"Ga^{3+}": 0.253447023, "[Ga(OH)_4]^-": 0.743146732,
                   "[Ga(OH)]^{2+}": 0.003406245, "Cl^-": 3.0, "Na^+": 3.0,
                   "H^+": 10.0 ** -1.62}
            I0 = ionic_strength(led, 1.0)
            Ie = sit_fixpoint(led, 1.0, ["Ga^{3+}", "[Ga(OH)_4]^-", "H^+",
                                         "H_2O"], I0)
            lk0 = 37.6 - 4.0 * pk
            print("  I(full ledger)=%s  I(reaction, fixpoint)=%s" % (I0, Ie))
            print("  logK(gamma=1, I->0) = logb4-4pKw = %.4f" % lk0)
            print("  SIT term (dz2=-4)   = %.4f"
                  % (-4.0 * 0.51 * Ie ** 0.5 / (1.0 + 1.5 * Ie ** 0.5)))
            print("  logK_eff            = %.4f"
                  % (lk0 - 4.0 * 0.51 * Ie ** 0.5 / (1.0 + 1.5 * Ie ** 0.5)))
            print("  logQ (ledger,pH=1.62) = %.6f"
                  % (__import__("math").log10(0.743146732 / 0.253447023)
                     - 4.0 * 1.62))
            print("  S = logK_eff - logQ = %.4f  (probe reports -13.305)"
                  % (lk0 - 4.0 * 0.51 * Ie ** 0.5 / (1.0 + 1.5 * Ie ** 0.5)
                     - (__import__("math").log10(0.743146732 / 0.253447023)
                        - 4.0 * 1.62)))
    except Exception as exc:                                    # noqa: BLE001
        print("  reconciliation skipped: %s" % exc)

    # ---------------------------------------------------------------- case scan
    selected, amphoteric, eq_mention = [], [], []
    for idx, c in enumerate(cases):
        hits = []
        for k in ("has", "has_range", "has_any"):
            for sp, v in (c.get(k) or {}).items():
                if sp in hydrox or sp in oxide:
                    hits.append({"solid": sp, "assert_key": k,
                                 "lo": v[0] if isinstance(v, (list, tuple)) else v,
                                 "spec": v,
                                 "class": "hydroxide" if sp in hydrox else "oxide",
                                 "metal": (hydrox.get(sp) or oxide.get(sp))["metal"],
                                 "n_oh": hydrox.get(sp, {}).get("n_oh")})
        if hits:
            selected.append({"idx": idx, "case": c, "hits": hits})
        for sp in (c.get("has_not") or {}):
            if sp in hydrox or sp in oxide:
                amphoteric.append({"idx": idx, "name": c["name"], "solid": sp,
                                   "cap": c["has_not"][sp]})
        txt = json.dumps({k: c.get(k) for k in ("eq", "eq_has") if k in c},
                         ensure_ascii=False)
        for sp in list(hydrox) + list(oxide):
            if sp in txt:
                eq_mention.append({"idx": idx, "name": c["name"], "solid": sp})

    print("\ncases in library: %d" % len(cases))
    print("cases asserting a metal hydroxide/oxide SOLID in has/has_range/"
          "has_any: %d" % len(selected))
    print("cases asserting has_not on such a solid (amphoteric ctx): %d"
          % len(amphoteric))
    print("cases mentioning such a solid only in eq/eq_has: %d"
          % len(eq_mention))

    for s in selected:
        c = s["case"]
        base_oh = 0.0
        for nm, mol in c["subs"]:
            if _is_strong_base(nm):
                base_oh += mol * _oh_groups(nm)
        s["base_oh_mol"] = base_oh
        s["metal_mol"] = {}
        for h in s["hits"]:
            M = h["metal"]
            tot = 0.0
            for nm, mol in c["subs"]:
                cnt = _elems(nm).get(M, 0)
                if cnt and not _is_strong_base(nm):
                    tot += cnt * mol
            s["metal_mol"][M] = tot

    print("\n--- task 1 inventory (asserted metal hydroxide/oxide solids) ---")
    print("%-34s %-10s %-4s %-8s %-8s %-7s" %
          ("case", "solid", "nOH", "base_OH", "metal", "OH/M"))
    for s in selected:
        c = s["case"]
        for h in s["hits"]:
            M = h["metal"]
            mm = s["metal_mol"].get(M, 0.0)
            print("%-34s %-10s %-4s %-8.4g %-8.4g %-7s"
                  % (c["name"][:34], h["solid"], h["n_oh"], s["base_oh_mol"],
                     mm, ("%.3f" % (s["base_oh_mol"] / mm)) if mm else "n/a"))

    # ---------------------------------------------------------------- run probes
    results = []
    t_used = 0.0
    for i, s in enumerate(selected):
        c = s["case"]
        if t_used > TIME_BUDGET_S:
            print("!! time budget hit, skipping %d remaining cases"
                  % (len(selected) - i))
            break
        probe = {}
        rec = {"name": c["name"], "idx": s["idx"], "subs": c["subs"],
               "cond": c.get("cond"), "note": c.get("note"),
               "asserts": {k: c[k] for k in
                           ("changed", "reacted", "degree", "has", "has_range",
                            "has_any", "has_not", "ph", "eq", "eq_has", "ann")
                           if k in c},
               "base_oh_mol": s["base_oh_mol"], "metal_mol": s["metal_mol"],
               "hits": s["hits"]}
        ts = time.time()
        try:
            r = eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                          c.get("cond") or {"V_L": 1.0}, T, _probe=probe)
            rec["judge_ok"] = True
            rec["override"] = r.get("override")
            rec["result"] = {k: r.get(k) for k in
                             ("changed", "reacted", "degree", "final_pH",
                              "annotations")}
            rec["result"]["production"] = [(p["name"], p["mol"])
                                           for p in r.get("production") or []]
            rec["result"]["consumption"] = [(p["name"], p["mol"])
                                            for p in r.get("consumption") or []]
            rec["result"]["steps"] = [(st.get("equation"), st.get("extent"))
                                      for st in (r.get("steps") or [])]
            rec["result"]["equations"] = [str(e) for e in
                                          (r.get("equations") or [])]
            # path signature: did a step that PRODUCES the asserted solid run?
            sig = {}
            for h in s["hits"]:
                sp = h["solid"]
                tot, n = 0.0, 0
                for eqn, ext in rec["result"]["steps"]:
                    if sp in (eqn or "").partition(" -> ")[2]:
                        n += 1
                        tot += float(ext or 0.0)
                sig[sp] = {"n_solid_forming_steps": n,
                           "solid_forming_extent_sum": round(tot, 9)}
            rec["solid_step_sig"] = sig
        except Exception as exc:                                # noqa: BLE001
            rec["judge_ok"] = False
            rec["error"] = "%s: %s" % (type(exc).__name__, exc)
        rec["ms"] = round((time.time() - ts) * 1000.0, 1)
        t_used += time.time() - ts
        rec["probe_empty"] = not probe

        if rec.get("judge_ok") and probe:
            cs = []
            for h in s["hits"]:
                a = _classify(probe, h["solid"], h["metal"], h["n_oh"], h["lo"])
                a["assert_key"] = h["assert_key"]
                a["assert_spec"] = h["spec"]
                cs.append(a)
            rec["classifications"] = cs
            rank = {"same_defect": 0, "other": 1, "solid_ok": 2}
            rec["bucket"] = sorted(cs, key=lambda a: rank[a["bucket"]])[0]["bucket"]
            rec["probe"] = {
                "exit": probe.get("exit"), "iters": probe.get("iters"),
                "hist": probe.get("hist"), "steps_n": probe.get("steps_n"),
                "pH": probe.get("pH"), "pH_solver": probe.get("pH_solver"),
                "H_excess": probe.get("H_excess"),
                "max_abs_S": probe.get("max_abs_S"),
                "frozen_n": probe.get("frozen_n"),
                "disabled_n": probe.get("disabled_n"),
                "ledger": probe.get("ledger"),
                "active": probe.get("active"),
            }
        else:
            rec["bucket"] = "other" if rec.get("judge_ok") else "run_error"
            rec["classifications"] = []
            rec["probe"] = {}
            if rec.get("judge_ok"):
                rec["classifications"] = [
                    dict(_classify({}, h["solid"], h["metal"], h["n_oh"], h["lo"]),
                         assert_key=h["assert_key"], assert_spec=h["spec"],
                         bucket="other",
                         bucket_why=("empty probe: judge returned before the "
                                     "walk (override=%s)" % rec.get("override")))
                    for h in s["hits"]]
        results.append(rec)
        print("[%3d/%3d] %-32s -> %-11s exit=%-9s it=%-5s pH=%-8s %6.0fms"
              % (i + 1, len(selected), c["name"][:32], rec["bucket"],
                 (rec.get("probe") or {}).get("exit"),
                 (rec.get("probe") or {}).get("iters"),
                 (rec.get("probe") or {}).get("pH"), rec["ms"]))

    # ---------------------------------------------------------------- summary
    counts = {}
    for rec in results:
        counts[rec["bucket"]] = counts.get(rec["bucket"], 0) + 1
    print("\n=== bucket counts ===")
    for k in ("solid_ok", "same_defect", "other", "run_error"):
        print("  %-12s %d" % (k, counts.get(k, 0)))

    print("\n=== B: same_defect (task 4 fields) ===")
    for rec in results:
        for a in rec["classifications"]:
            if a["bucket"] != "same_defect":
                continue
            h = next(x for x in rec["hits"] if x["solid"] == a["solid_expect"])
            mm = rec["metal_mol"].get(h["metal"], 0.0)
            p = rec.get("probe") or {}
            print("-" * 96)
            print("case            : %s" % rec["name"])
            print("subs            : %s   cond=%s" % (rec["subs"], rec.get("cond")))
            print("asserts         : %s" % json.dumps(rec["asserts"],
                                                      ensure_ascii=False))
            print("solid asserted  : %s %s=%s   (n_oh=%s)"
                  % (h["assert_key"], a["solid_expect"], h["spec"], h["n_oh"]))
            print("base_OH/metal   : %.6g / %.6g = %s"
                  % (rec["base_oh_mol"], mm,
                     ("%.3f" % (rec["base_oh_mol"] / mm)) if mm else "n/a"))
            print("probe           : exit=%s iters=%s hist=%s steps_n=%s "
                  "pH=%s pH_solver=%s H_excess=%s max_abs_S=%s frozen_n=%s "
                  "disabled_n=%s"
                  % (p.get("exit"), p.get("iters"), p.get("hist"),
                     p.get("steps_n"), p.get("pH"), p.get("pH_solver"),
                     p.get("H_excess"), p.get("max_abs_S"), p.get("frozen_n"),
                     p.get("disabled_n")))
            print("solid_mol       : %.9g  (%.6g%% of asserted lo %.6g)"
                  % (a["solid_mol"],
                     100.0 * a["solid_mol"] / h["lo"] if h["lo"] else 0.0, h["lo"]))
            print("dominant hydroxo: %s = %.9g  (%.3f%% of dissolved metal %.9g)"
                  % (a["dom_hydroxo_sp"], a["dom_hydroxo_mol"],
                     100.0 * (a["dom_hydroxo_frac_of_dissolved"] or 0.0),
                     a["dissolved_metal_total"]))
            print("hydroxo species : %s" % json.dumps(a["hydroxo_species"],
                                                      ensure_ascii=False))
            print("free metal ion  : %s = %.9g" % (a["free_metal_sp"],
                                                   a["free_metal_mol"]))
            print("max two-sided   : |S|=%.3f  kind=%s S=%s"
                  % (a["max_twosided_absS"], a["max_twosided"]["kind"],
                     a["max_twosided"]["S"]))
            print("                  eq=%s" % a["max_twosided"]["eq"])
            print("                  dis_why=%s  ext_max=%s  zero_push=%s"
                  % (json.dumps(a["max_twosided"]["dis_why"], ensure_ascii=False),
                     a["max_twosided"]["ext_max"], a["max_twosided"]["zero_push"]))
            print("n two-sided |S|>%.1f tagged 零推进 : %d"
                  % (S_STRONG, a["n_twosided_absS_gt3_zero_push"]))
            for z in a["twosided_zero_push"]:
                print("      S=%+9.3f ext_max=%-12.6g %-9s %s"
                      % (z["S"], z["ext_max"], z["kind"], z["eq"]))
            print("precip candidates involving %s (S normalised to precip):"
                  % a["solid_expect"])
            for pc in a["precip_candidates"]:
                print("      dir=%-8s S(precip)=%+9.3f kind=%-8s two_sided=%s "
                      "why=%s" % (pc["dir"], pc["S"], pc["kind"],
                                  pc["two_sided"],
                                  json.dumps(pc["dis_why"],
                                             ensure_ascii=False)[:44]))
                print("            %s" % pc["eq"])
            print("undersaturated precip candidate (S<0): %s"
                  % (json.dumps(a["undersat_precip"], ensure_ascii=False)
                     if a["undersat_precip"] else "NONE"))
            print("executed steps  :")
            for eqn, ext in (rec["result"].get("steps") or []):
                print("      ext=%-12s %s" % (ext, eqn))
            print("solid-forming steps: %s"
                  % json.dumps((rec.get("solid_step_sig") or {}).get(
                      a["solid_expect"]), ensure_ascii=False))
            print("net equations   : %s"
                  % json.dumps(rec["result"].get("equations"),
                               ensure_ascii=False))
            print("production      : %s"
                  % json.dumps(rec["result"].get("production"),
                               ensure_ascii=False))

    print("\n=== A: solid_ok, split by OH/M ratio and n_oh ===")
    triv = [r for r in results
            if r["bucket"] == "solid_ok"
            and any(a["solid_n_oh"] == 3 for a in r["classifications"])]
    print("-- n_oh == 3 (trivalent) : %d cases" % len(triv))
    print("%-34s %-6s %-8s %-10s %-8s %-9s %s"
          % ("case", "OH/M", "solid", "assert lo", "exit", "iters", "pH"))
    for r in triv:
        for a in r["classifications"]:
            if a["solid_n_oh"] != 3:
                continue
            h = next(x for x in r["hits"] if x["solid"] == a["solid_expect"])
            mm = r["metal_mol"].get(h["metal"], 0.0)
            print("%-34s %-6s %-8.5f %-10s %-8s %-9s %s"
                  % (r["name"][:34],
                     ("%.3f" % (r["base_oh_mol"] / mm)) if mm else "n/a",
                     a["solid_mol"], h["lo"], r["probe"]["exit"],
                     r["probe"]["iters"], r["probe"]["pH"]))
        continue

    print("\n=== C: other ===")
    for rec in results:
        if rec["bucket"] != "other":
            continue
        for a in rec["classifications"]:
            print("  %-32s %-10s %s" % (rec["name"][:32], a["solid_expect"],
                                        a["bucket_why"][:160]))
        if not rec["classifications"]:
            print("  %-32s (no solid assertion classified)" % rec["name"][:32])

    doc = {
        "meta": {
            "at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "cases_total": len(cases),
            "cases_selected": len(selected),
            "cases_probed": len(results),
            "elapsed_s": round(time.time() - t0, 1),
            "thresholds": {"X_MIN": X_MIN, "S_STRONG": S_STRONG,
                           "NEAR0_FRAC": NEAR0_FRAC, "OK_FRAC": OK_FRAC},
            "zero_push_tag": ZERO_PUSH,
            "solid_registry": {"hydroxides": hydrox, "oxides": oxide},
            "amphoteric_table": ampho,
        },
        "inventory": [{"name": s["case"]["name"], "subs": s["case"]["subs"],
                       "base_oh_mol": s["base_oh_mol"],
                       "metal_mol": s["metal_mol"], "hits": s["hits"]}
                      for s in selected],
        "amphoteric_has_not": amphoteric,
        "eq_mention_only": eq_mention,
        "bucket_counts": counts,
        "results": results,
    }
    out = os.path.join(ROOT, "logs", "fam_hydroxo.json")
    with io.open(out, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1, default=str)
    print("\nwrote %s  (%.1fs)" % (out, time.time() - t0))


if __name__ == "__main__":
    main()
