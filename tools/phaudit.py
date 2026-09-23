# -*- coding: utf-8 -*-
"""Round 204 - pH code-path divergence audit (READ-ONLY on chemkit/).

Audits the four coexisting pH producers for the same (ledger, H_excess) state:
  1. speciation.estimate_state / estimate_pH  (4-branch heuristic + 2 exact
     side paths: pinned block L973-993, non-pinned degenerate block L996-1010)
  2. speciation.closed_pH                     (fast path, None on weak species)
  3. speciation.exact_proton_pH               (B3 gate L1191, cold path)
  4. engine._presentation_He + presentation_pH (rewrites He to 0.0 on the
     phantom-base gate, then feeds that rewritten He onward)

Modes (argv[1]):
  recon    - environment reconnaissance, L08 identification
  sweep    - full-suite terminal-state divergence sweep  -> logs/phaudit_sweep.json
  gates    - L1191 gate-narrowing experiment             -> logs/phaudit_gates.json
  tags     - PH_TAGS branch-frequency instrumentation
  f31      - the proven F31 case, end to end (regression anchor)

Usage: python tools/phaudit.py <mode> [args]
"""
import io
import json
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                       # noqa: BLE001
        pass

import chemkit.acidbase as ab                                # noqa: E402
import chemkit.engine as eng                                 # noqa: E402
import chemkit.speciation as spec                            # noqa: E402
from chemkit.data import load_tables                         # noqa: E402
from chemkit.testsuit import load_cases                      # noqa: E402

LOGS = os.path.join(ROOT, "logs")
T = load_tables()
V, T_K = 1.0, 298.15
PKW = spec.pKw_of(T_K)


def wr(name: str, obj) -> str:
    os.makedirs(LOGS, exist_ok=True)
    p = os.path.join(LOGS, name)
    with io.open(p, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1, default=str)
    return p


def ledger_clean(led: dict) -> dict:
    """Ledger without engine metadata keys (__tot_*) - those are shadows."""
    return {s: float(m) for s, m in led.items()
            if not s.startswith("__") and s != spec.WATER and float(m) > 0.0}


def ledger_full(led: dict) -> dict:
    """Ledger as the engine rounded it, KEEPING __tot_* shadows (they are what
    _respeciate_strong_acids reads), dropping the solvent."""
    return {s: float(m) for s, m in led.items()
            if s != spec.WATER and float(m) >= 0.0}


def net_charge(led: dict) -> float:
    """Sum z_i * n_i over ledger species, excluding H+/OH- and __metadata."""
    tot = 0.0
    for s, m in led.items():
        if s.startswith("__") or m <= 0.0:
            continue
        if s in (spec.H_ION, "OH^-"):
            continue
        tot += spec._charge_cached(s) * m
    return tot


def truth_of(net: float, v: float, pkw: float):
    """Charge-consistent pH implied by a fixed-charge inventory `net` (mol).

    The solution's charge balance is net + V*(h - oh) = 0:
      net > 0  -> the ledger is cation-rich, free OH- must balance it,
      net < 0  -> free H+ must balance it,
      |net|/V <= X_MIN -> no constraint at all (the ledger is electroneutral
                          and pH is fixed by the weak-acid chemistry instead);
                          reported as pKw/2 with defined=False.
    Returns (pH, defined).
    """
    fixed = net / v
    if abs(fixed) <= spec.X_MIN:
        return pkw / 2.0, False
    if fixed > 0.0:
        return pkw + math.log10(fixed), True
    return -math.log10(-fixed), True


def safe(fn, *a, **kw):
    """(value, error) with exceptions turned into a marker string."""
    try:
        return fn(*a, **kw), None
    except Exception as exc:                                # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def ph_of(x):
    """closed_pH/exact_proton_pH return tuples or floats or None."""
    if x is None:
        return None
    if isinstance(x, tuple):
        return x[0]
    return x


# ---------------------------------------------------------------- recon

def mode_recon(argv=None):
    cases = load_cases(None)
    print(f"cases            = {len(cases)}")
    print(f"T.solids         = {len(T.solids)}")
    print(f"weak_species_set = {len(spec.weak_species_set(T))}")
    print(f"X_MIN            = {spec.X_MIN}")
    print(f"pKw(298.15)      = {PKW}")
    print(f"PINNED_TAKEOVER  = {spec.PINNED_TAKEOVER}")
    print(f"SIT_ALL          = {spec.SIT_ALL}")
    ksp_oh = [e["pair"][0] for e in T.ksp if e["pair"][1] == "OH^-"]
    print(f"Ksp-OH cations   = {len(set(ksp_oh))}: {sorted(set(ksp_oh))}")
    print(f"n Ksp entries    = {len(T.ksp)}")
    fam = spec.build_families(T)
    print(f"families         = {len(fam)}")
    # L08 identification
    for pre in ("L08", "L8"):
        hit = [c for c in cases if c["name"].startswith(pre + " ")]
        for c in hit:
            print(f"\nCASE {c['name']}: subs={c['subs']} cond={c.get('cond')}")
    # all cases whose name starts with L
    ln = [c["name"] for c in cases if c["name"].startswith("L")]
    print(f"\nL-family names ({len(ln)}): {ln}")
    return 0


# ------------------------------------------------------- B3 gate replica

def b3_info(led, he, res_set, fam_set):
    """Replica of exact_proton_pH's decision (speciation.py L1186-L1201),
    returning (result_or_None, reason, detail)."""
    net = 0.0
    nfam = 0
    fams = []
    for sp, m in led.items():
        if m <= 0.0 or sp == spec.WATER or sp.startswith("__"):
            continue
        if m > spec.X_MIN and sp in res_set:
            return None, "reservoir", sp
        net += spec._charge_cached(sp) * m
        if sp in fam_set:
            nfam += 1
            fams.append(sp)
    if not nfam:
        return None, "no_family", None
    if abs(net + he) > 1e-6:
        return None, "charge_mismatch", net + he
    return "ok", "ok", (net, nfam, fams)


def b3_state(T):
    b3 = getattr(T, "_b3_static", None)
    if b3 is None:
        b3 = T._b3_static = (
            frozenset(T.solids) | frozenset(e["pair"][0] for e in T.ksp
                                            if e["pair"][1] == "OH^-"),
            frozenset(spec.build_families(T)))
    return b3


def exact_proton_pH_solids_only(led, he, v, T, T_K, solids_set, fam_set):
    """The PROPOSED narrowing of the L1191 gate: reservoir = SOLID only."""
    r, why, det = b3_info(led, he, solids_set, fam_set)
    if r is None:
        return None, why, det
    return ab.charge_pH(led, v, T, T_K, fast=True), "ok", det


# ---------------------------------------------------------------- sweep

# Read-only instrumentation of the engine's phantom-base gate. `_presentation_He`
# is a module global of chemkit.engine, and engine.py resolves it at call time,
# so replacing the attribute only affects THIS process's measurement -- no file
# under chemkit/ is touched. Records (He_in, He_out) per call.
PRES_LOG = []
_PRES_ORIG = eng._presentation_He


def _pres_wrap(ledger, he, v, t, tk):
    out = _PRES_ORIG(ledger, he, v, t, tk)
    try:
        PRES_LOG.append((float(he), float(out)))
    except Exception:                                       # noqa: BLE001
        pass
    return out


eng._presentation_He = _pres_wrap

def measure_case(c, T, probe_out=None):
    """Run one test case and measure every pH producer on its terminal state."""
    subs = [{"name": n, "mol": m} for n, m in c["subs"]]
    cond = c.get("cond") or {"V_L": 1.0}
    v = float(cond.get("V_L", 1.0) or 1.0)
    tk = float(cond.get("T_K", 298.15) or 298.15)
    pkw = spec.pKw_of(tk)
    c_H = float(cond.get("c_H") or 0.0)
    c_OH = float(cond.get("c_OH") or 0.0)
    probe = {}
    PRES_LOG.clear()
    t0 = time.perf_counter()
    res, jerr = safe(eng.judge, subs, cond, T, _probe=probe)
    ms = (time.perf_counter() - t0) * 1000.0
    rec = {"name": c["name"], "ms_judge": round(ms, 1), "V_L": v, "T_K": tk,
           "c_H": c_H, "c_OH": c_OH, "c_pH": cond.get("pH"),
           "judge_err": jerr}
    if res is None or not probe:
        rec["no_probe"] = True
        rec["final_pH"] = (res or {}).get("final_pH")
        rec["H_excess_res"] = (res or {}).get("H_excess")
        return rec
    led = ledger_clean(probe.get("ledger") or {})
    led_raw = ledger_full(probe.get("ledger") or {})
    he_pres = float(probe.get("H_excess") or 0.0)
    rec["probe_pH"] = probe.get("pH")
    rec["probe_pH_solver"] = probe.get("pH_solver")
    rec["probe_H_excess"] = he_pres
    rec["final_pH"] = (res or {}).get("final_pH")
    rec["exit"] = probe.get("exit")
    rec["iters"] = probe.get("iters")
    rec["max_abs_S"] = probe.get("max_abs_S")
    rec["steps_n"] = probe.get("steps_n")
    # ---- the four producers, all on (probe ledger, probe H_excess)
    e, eerr = safe(eng.estimate_pH, led, he_pres, v, T, tk)
    rec["est"] = e
    rec["est_err"] = eerr
    cl, clerr = safe(spec.closed_pH, led, he_pres, v, T, tk)
    rec["cls"] = ph_of(cl)
    rec["cls_err"] = clerr
    ex, exerr = safe(spec.exact_proton_pH, led, he_pres, v, T, tk)
    rec["exact"] = ph_of(ex)
    rec["exact_err"] = exerr
    ch, cherr = safe(ab.charge_pH, led, v, T, tk)          # as prescribed
    rec["charge"] = ch
    rec["charge_err"] = cherr
    # the engine-internal convention: cond's strong-acid condition is a
    # proton-only charge that charge_pH must be told about (L316 docstring)
    if c_H:
        ch2, ch2err = safe(ab.charge_pH, led, v, T, tk, c_H=c_H)
        rec["charge_ch"] = ch2
        rec["charge_ch_err"] = ch2err
    else:
        rec["charge_ch"] = ch
        rec["charge_ch_err"] = cherr
    rec["net"] = net_charge(led)
    # truth A: literal prescription (net as the fixed charge)
    rec["truth"], rec["truth_defined"] = truth_of(rec["net"], v, pkw)
    # truth B: net minus the proton-only acid condition (engine invariant:
    #   sum(z*n) + He - c_H = 0  =>  the fixed charge the solution must
    #   balance with free H+/OH- is net - c_H)
    rec["truth_ch"], rec["truth_ch_defined"] = truth_of(
        rec["net"] - c_H, v, pkw)
    # ---- gate forensics + proposed narrowing
    res_set, fam_set = b3_state(T)
    _ok, why, det = b3_info(led, he_pres, res_set, fam_set)
    rec["gate_reason"] = why
    rec["gate_detail"] = det
    so, sowhy, sodet = exact_proton_pH_solids_only(
        led, he_pres, v, T, tk, frozenset(T.solids), fam_set)
    rec["exact_solid_only"] = so
    rec["exact_solid_reason"] = sowhy
    rec["exact_solid_detail"] = sodet
    # ---- solid / reservoir inventory
    solids = sorted(s for s in led if s in T.solids)
    cats = sorted(s for s in led
                  if s in {e2["pair"][0] for e2 in T.ksp
                           if e2["pair"][1] == "OH^-"})
    rec["solids"] = solids
    rec["ksp_cats"] = cats
    rec["nfam"] = det[1] if why == "ok" else None
    rec["fams"] = det[2] if why == "ok" else None
    # ---- presentation-He forensics (bucket iii)
    # The probe's H_excess is ALREADY the post-gate value, so recomputing
    # _presentation_He from it can never reveal a rewrite (it is idempotent).
    # The honest detector is a read-only wrapper installed around the engine's
    # module-global `_presentation_He` (see PRES_LOG / _pres_wrap): it records
    # (He_in, He_out) for every call the engine makes during this case.
    calls = list(PRES_LOG)
    rec["pres_calls"] = [list(x) for x in calls]
    rec["n_pres_calls"] = len(calls)
    if calls:
        # the LAST call is the one _finalize_result makes on the exit state,
        # i.e. the terminal presentation rewrite; earlier calls are the walk's
        # own phantom-base caps (engine L415).
        rec["He_walk"] = calls[-1][0]         # exit He, pre-gate
        rec["He_gate_out"] = calls[-1][1]
        rec["pres_rewrote"] = bool(calls[-1][0] != calls[-1][1])
        rec["pres_zeroed"] = bool(calls[-1][1] == 0.0 and calls[-1][0] != 0.0)
        rec["pres_rewrote_any"] = any(a != b for a, b in calls)
    else:
        rec["He_walk"] = None
        rec["He_gate_out"] = None
        rec["pres_rewrote"] = None
        rec["pres_zeroed"] = None
        rec["pres_rewrote_any"] = None
    if rec.get("He_walk") is not None:
        ew, ewerr = safe(eng.estimate_pH, led, rec["He_walk"], v, T, tk)
        rec["est_at_walk"] = ew
        rec["est_at_walk_err"] = ewerr
    else:
        rec["est_at_walk"] = None
        rec["est_at_walk_err"] = None
    # ---- branch tag of the heuristic AT THE TERMINAL STATE (PH_TAGS hook)
    _save = spec.PH_TAGS
    spec.PH_TAGS = []
    try:
        eng.estimate_pH(led, he_pres, v, T, tk)
        rec["tags"] = list(spec.PH_TAGS)
    finally:
        spec.PH_TAGS = _save
    rec["tag_last"] = rec["tags"][-1] if rec["tags"] else None
    # ---- divergence metrics over the pH-valued producers
    vals = {}
    for k in ("probe_pH", "probe_pH_solver", "est", "cls", "exact", "charge",
              "truth"):
        x = rec.get(k)
        if isinstance(x, bool):
            continue
        if isinstance(x, (int, float)) and math.isfinite(x):
            vals[k] = float(x)
    rec["n_vals"] = len(vals)
    for tag, keys, out in (
            ("", ("probe_pH", "probe_pH_solver", "est", "cls", "exact",
                  "charge", "truth"), "spread"),
            ("_ch", ("probe_pH", "probe_pH_solver", "est", "cls", "exact",
                     "charge_ch", "truth_ch"), "spread_ch")):
        vv = {k: vals[k] for k in keys if k in vals}
        if len(vv) >= 2:
            rec[out] = max(vv.values()) - min(vv.values())
            rec[out + "_lo"] = min(vv, key=vv.get)
            rec[out + "_hi"] = max(vv, key=vv.get)
        else:
            rec[out] = None
    # physical (guard-railed) spread: only meaningful when the ledger really
    # carries a free-acid/base inventory to balance (truth defined) and the
    # engine's charge_pH was told about c_H
    vv = {k: vals[k] for k in ("probe_pH", "probe_pH_solver", "est", "cls",
                               "exact", "charge_ch", "truth_ch") if k in vals}
    if not rec["truth_ch_defined"]:
        vv.pop("truth_ch", None)
    rec["spread_phys"] = (max(vv.values()) - min(vv.values())
                          if len(vv) >= 2 else None)
    if vv:
        rec["spread_phys_lo"] = min(vv, key=vv.get)
        rec["spread_phys_hi"] = max(vv, key=vv.get)
    rec["lo_key"] = rec.get("spread_lo")
    rec["hi_key"] = rec.get("spread_hi")
    if probe_out is not None:
        probe_out.append(rec)
    return rec


def _fmt_row(rec, keys):
    cells = []
    for k in keys:
        x = rec.get(k)
        if x is None:
            cells.append("     None")
        elif isinstance(x, float):
            cells.append(f"{x:9.3f}")
        else:
            cells.append(f"{x:>9}")
    return cells


SWEEP_KEYS = ("probe_pH", "probe_pH_solver", "est", "cls", "exact", "charge",
              "truth")


def mode_sweep(argv):
    """argv: [step] [--full]  step=1 -> every case, step=N -> every Nth."""
    step = 1
    if argv:
        try:
            step = max(1, int(argv[0]))
        except ValueError:
            step = 1
    cases = load_cases(None)
    sel = cases[::step]
    print(f"sweep: {len(sel)} of {len(cases)} cases (step={step})")
    recs = []
    t0 = time.perf_counter()
    for i, c in enumerate(sel):
        r = measure_case(c, T)
        recs.append(r)
        if (i + 1) % 25 == 0 or i + 1 == len(sel):
            el = time.perf_counter() - t0
            print(f"  {i + 1}/{len(sel)}  elapsed {el:.1f}s  "
                  f"eta {el / (i + 1) * (len(sel) - i - 1):.0f}s", flush=True)
    n1 = sum(1 for r in recs
             if r.get("spread") is not None and r["spread"] > 1.0)
    n3 = sum(1 for r in recs
             if r.get("spread") is not None and r["spread"] > 3.0)
    print(f"\n=== divergence (any two of {'/'.join(SWEEP_KEYS)}) ===")
    print(f"  cases with >=2 finite producers : "
          f"{sum(1 for r in recs if r.get('n_vals', 0) >= 2)}")
    print(f"  spread > 1 pH unit              : {n1}")
    print(f"  spread > 3 pH units             : {n3}")
    rank = sorted((r for r in recs if r.get("spread") is not None),
                  key=lambda r: -r["spread"])
    print(f"\n  {'case':28s} " + " ".join(f"{k:>9}" for k in SWEEP_KEYS)
          + "   spread  lo->hi")
    for r in rank[:30]:
        cells = _fmt_row(r, SWEEP_KEYS)
        print(f"  {r['name'][:28]:28s} " + " ".join(cells)
              + f"  {r['spread']:7.3f}  {r['lo_key']}->{r['hi_key']}")
    p = wr(f"phaudit_sweep_step{step}.json",
           {"step": step, "n_cases": len(recs),
            "n_spread_gt1": n1, "n_spread_gt3": n3, "records": recs})
    print(f"\nwrote {p}")
    return 0


# ---------------------------------------------------------------- one case

def _print_ledger(led, v, pkw, title):
    print(f"\n--- {title} ---")
    net = 0.0
    for sp, m in sorted(led.items(), key=lambda kv: -abs(kv[1])):
        if m <= 0.0 or sp == spec.WATER:
            continue
        if sp.startswith("__"):
            print(f"  {sp:24s} {m:14.9g}  (shadow)")
            continue
        z = spec._charge_cached(sp)
        net += z * m
        print(f"  {sp:24s} {m:14.9g}  z={z:+d}  z*n={z * m:+14.9g}")
    print(f"  {'SUM z*n':24s} {'':14s}          {net:+14.9g} mol")
    return net


def mode_one(argv):
    """One case, end to end: every producer + the terminal ledger."""
    pre = argv[0] if argv else "F31"
    cases = load_cases(None)
    hit = [c for c in cases if c["name"].startswith(pre + " ")
           or c["name"] == pre]
    if not hit:
        print(f"no case matching {pre!r}")
        return 1
    c = hit[0]
    print(f"CASE {c['name']}  subs={c['subs']}  cond={c.get('cond')}")
    v = float((c.get("cond") or {}).get("V_L", 1.0) or 1.0)
    tk = float((c.get("cond") or {}).get("T_K", 298.15) or 298.15)
    pkw = spec.pKw_of(tk)
    rec = measure_case(c, T)
    for k in sorted(rec):
        if k in ("tags", "fams", "solids", "ksp_cats", "gate_detail"):
            continue
        print(f"  {k:20s} = {rec[k]}")
    probe = {}
    eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
              c.get("cond") or {"V_L": 1.0}, T, _probe=probe)
    led = ledger_full(probe.get("ledger") or {})
    net = _print_ledger(led, v, pkw, "terminal ledger (as engine rounded it)")
    print(f"  probe H_excess (post-gate) = {probe.get('H_excess')}  "
          f"He_walk (pre-gate) = {rec.get('He_walk')}  "
          f"He_gate_out = {rec.get('He_gate_out')}  "
          f"pres_calls = {rec.get('pres_calls')}")
    print(f"  charge_pH(no c_H)     = {rec.get('charge')}")
    print(f"  charge_pH(c_H={rec.get('c_H')}) = {rec.get('charge_ch')}")
    print(f"  truth (literal formula) = {rec.get('truth')} "
          f"(defined={rec.get('truth_defined')})")
    print(f"  truth (c_H corrected)   = {rec.get('truth_ch')} "
          f"(defined={rec.get('truth_ch_defined')})")
    print(f"  gate_reason = {rec.get('gate_reason')} det={rec.get('gate_detail')}")
    print(f"  exact_solid_only = {rec.get('exact_solid_only')} "
          f"({rec.get('exact_solid_reason')})")
    print(f"  solids = {rec.get('solids')}")
    print(f"  Ksp-OH cations present = {rec.get('ksp_cats')}")
    print(f"  terminal branch tag = {rec.get('tag_last')} tags={rec.get('tags')}")
    return 0


# ---------------------------------------------------------------- report

SWEEP_KEYS = ("probe_pH", "probe_pH_solver", "est", "cls", "exact", "charge",
              "truth")
SWEEP_KEYS_CH = ("probe_pH", "probe_pH_solver", "est", "cls", "exact",
                 "charge_ch", "truth_ch")


def _fmt_row(rec, keys):
    cells = []
    for k in keys:
        x = rec.get(k)
        if x is None:
            cells.append("     None")
        elif isinstance(x, (int, float)):
            cells.append(f"{x:9.3f}")
        else:
            cells.append(f"{x:>9}")
    return cells


def _load_sweep():
    p = os.path.join(LOGS, "phaudit_sweep_step1.json")
    if not os.path.exists(p):
        raise SystemExit(f"missing {p}; run: python tools/phaudit.py sweep 1")
    with io.open(p, encoding="utf-8") as f:
        return json.load(f)


def _table(recs, keys, title, n=30):
    print(f"\n=== {title} ===")
    print(f"  {'case':30s} " + " ".join(f"{k:>9}" for k in keys)
          + "   spread  lo->hi")
    for r in recs[:n]:
        sp = r.get("spread" if keys is SWEEP_KEYS else "spread_ch")
        print(f"  {r['name'][:30]:30s} " + " ".join(_fmt_row(r, keys))
              + f"  {sp if sp is None else round(sp, 3):>7}  "
              + f"{r.get('lo_key')}->{r.get('hi_key')}")


def mode_report(argv):
    j = _load_sweep()
    recs = j["records"]
    ok = [r for r in recs if not r.get("no_probe")]
    print(f"records: {len(recs)}  with probe: {len(ok)}  "
          f"no-probe (override): {len(recs) - len(ok)}")
    n1 = sum(1 for r in ok if (r.get("spread") or 0) > 1.0)
    n3 = sum(1 for r in ok if (r.get("spread") or 0) > 3.0)
    p1 = sum(1 for r in ok if (r.get("spread_phys") or 0) > 1.0)
    p3 = sum(1 for r in ok if (r.get("spread_phys") or 0) > 3.0)
    print(f"\n### B. divergence counts (any two producers differ by > X)")
    print(f"  literal set  {SWEEP_KEYS}:")
    print(f"      > 1 pH unit : {n1}")
    print(f"      > 3 pH units: {n3}")
    print(f"  engine-consistent set {SWEEP_KEYS_CH} (c_H passed to charge_pH,"
          f" c_H-corrected truth, undefined-truth cases dropped):")
    print(f"      > 1 pH unit : {p1}")
    print(f"      > 3 pH units: {p3}")
    rank = sorted(ok, key=lambda r: -(r.get("spread") or 0))
    _table(rank, SWEEP_KEYS, "B1. top-30 by literal spread")
    rankp = sorted(ok, key=lambda r: -(r.get("spread_phys") or 0))
    _table(rankp, SWEEP_KEYS_CH, "B2. top-30 by engine-consistent spread")
    # ------------------------------------------------ C. buckets
    def cnt(pred):
        return [r for r in ok if pred(r)]

    undef = cnt(lambda r: not r.get("truth_defined"))
    chfix = cnt(lambda r: r.get("charge_ch") is not None
                and r.get("charge") is not None
                and abs(r["charge_ch"] - r["charge"]) > 1.0)
    b_i = cnt(lambda r: r.get("exact") is None and r.get("est") is not None
              and r.get("truth_defined")
              and abs(r["est"] - r["truth_ch"]) > 1.0)
    b_i_gate = cnt(lambda r: r.get("exact") is None
                   and r.get("exact_solid_only") is not None)
    b_ii = cnt(lambda r: r.get("truth_defined") and r.get("est") is not None
               and abs(r["est"] - r["truth_ch"]) > 2.0)
    b_iii = cnt(lambda r: r.get("pres_rewrote"))
    b_iii0 = cnt(lambda r: r.get("probe_H_excess") == 0.0
                 and abs(r.get("net") or 0.0) > spec.X_MIN * (r["V_L"]))
    b_iv = cnt(lambda r: r.get("cls") is None)
    b_v = cnt(lambda r: r.get("probe_pH") is not None
              and r.get("probe_pH_solver") is not None
              and abs(r["probe_pH"] - r["probe_pH_solver"]) > 0.5)
    b_vi = cnt(lambda r: r.get("exact") is not None)
    print("\n### C. buckets")
    print(f"  (0) truth undefined (electroneutral ledger, pH set by weak "
          f"chemistry): {len(undef)}")
    print(f"  (0b) charge_pH differs from charge_pH(c_H) by >1: {len(chfix)}")
    print(f"  (i) exact_proton_pH None while estimate_pH gives a value that "
          f"differs from truth by >1: {len(b_i)}")
    print(f"  (i') gate would PASS under solids-only narrowing: {len(b_i_gate)}")
    print(f"  (ii) |estimate_pH - truth| > 2: {len(b_ii)}")
    print(f"  (iii) _presentation_He rewrote He: {len(b_iii)}")
    print(f"  (iii') probe H_excess == 0.0 while ledger charge != 0: "
          f"{len(b_iii0)}")
    print(f"  (iv) closed_pH returned None: {len(b_iv)}")
    print(f"  (v) probe pH != probe pH_solver (>0.5): {len(b_v)}")
    print(f"  (vi) exact_proton_pH DID fire: {len(b_vi)}")
    for tag, lst in (("i", b_i), ("ii", b_ii), ("iii", b_iii),
                     ("iii'", b_iii0), ("iv", b_iv), ("v", b_v),
                     ("vi", b_vi), ("i'", b_i_gate)):
        print(f"\n  -- examples bucket {tag} (first 5) --")
        for r in lst[:5]:
            print(f"    {r['name'][:34]:34s} probe={r.get('probe_pH')} "
                  f"solver={r.get('probe_pH_solver')} est={r.get('est')} "
                  f"cls={r.get('cls')} exact={r.get('exact')} "
                  f"charge={r.get('charge')} truth={r.get('truth')} "
                  f"He={r.get('probe_H_excess')} He_walk={r.get('He_walk')} "
                  f"gate={r.get('gate_reason')}")
    # ------------------------------------------------ D. gate narrowing
    print("\n### D. L1191 gate narrowing experiment")
    solid_present = [r for r in ok if r.get("solids")]
    print(f"  cases with a solid in the terminal ledger: {len(solid_present)}")
    narrowed = [r for r in ok if r.get("exact") is None
                and r.get("exact_solid_only") is not None]
    print(f"  cases the narrowing would newly ADMIT: {len(narrowed)}")
    for r in narrowed[:40]:
        print(f"    {r['name'][:30]:30s} est={r.get('est')} "
              f"solids_only={r.get('exact_solid_only')} engine={r.get('probe_pH')}"
              f" charge_ch={r.get('charge_ch')} truth_ch={r.get('truth_ch')} "
              f"cats={r.get('ksp_cats')} gate={r.get('gate_reason')}")
    print(f"  cases where the gate rejects on a SOLID (narrowing keeps "
          f"rejecting): "
          f"{sum(1 for r in ok if r.get('gate_reason') == 'reservoir' and r.get('solids'))}")
    print(f"  cases where the gate rejects on a Ksp-OH CATION only: "
          f"{sum(1 for r in ok if r.get('gate_reason') == 'reservoir' and not r.get('solids'))}")
    print("  gate_reason histogram: " + json.dumps(
        {k: sum(1 for r in ok if r.get("gate_reason") == k)
         for k in {r.get("gate_reason") for r in ok}}, ensure_ascii=False))
    print("  exact_solid_reason histogram: " + json.dumps(
        {k: sum(1 for r in ok if r.get("exact_solid_reason") == k)
         for k in {r.get("exact_solid_reason") for r in ok}},
        ensure_ascii=False))
    print("  terminal branch tag histogram: " + json.dumps(
        {k: sum(1 for r in ok if r.get("tag_last") == k)
         for k in {r.get("tag_last") for r in ok}}, ensure_ascii=False))
    print("  closed_pH/measurement errors: " + json.dumps(
        {k: sum(1 for r in ok if r.get(k)) for k in
         ("est_err", "cls_err", "exact_err", "charge_err", "charge_ch_err",
          "judge_err")}, ensure_ascii=False))
    out = {"n1": n1, "n3": n3, "phys1": p1, "phys3": p3,
           "buckets": {"truth_undefined": len(undef), "charge_ch_fix":
                       len(chfix), "i": len(b_i), "i_gate": len(b_i_gate),
                       "ii": len(b_ii), "iii": len(b_iii), "iii0": len(b_iii0),
                       "iv": len(b_iv), "v": len(b_v), "vi": len(b_vi)},
           "top30_literal": rank[:30], "top30_phys": rankp[:30],
           "narrowed_admit": narrowed}
    print(f"\nwrote {wr('phaudit_report.json', out)}")
    return 0


# ---------------------------------------------------------------- tags

class TagCounter:
    """PH_TAGS stand-in that counts instead of accumulating."""

    __slots__ = ("n", "last")

    def __init__(self):
        self.n = 0
        self.last = None

    def append(self, why):
        self.n += 1
        self.last = why

    def __bool__(self):
        return self.n > 0

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        if i == -1:
            return self.last
        raise IndexError(i)

    def __iter__(self):
        return iter(())


def _wrap(fn, cnt, key):
    def w(*a, **kw):
        prev = spec.PH_TAGS
        spec.PH_TAGS = []
        try:
            r = fn(*a, **kw)
        finally:
            tags = spec.PH_TAGS
            if prev is not None and tags:
                prev.extend(tags)
            spec.PH_TAGS = prev
            t = tags[-1] if tags else "-none-"
            cnt[key][t] = cnt[key].get(t, 0) + 1
        return r
    return w


def mode_tags(argv):
    """Walk-level branch-frequency instrumentation via PH_TAGS (read-only)."""
    step = max(1, int(argv[0])) if argv else 1
    cases = load_cases(None)[::step]
    cnt = {"estimate_pH": {}, "estimate_state": {}, "closed_pH": {},
           "exact_proton_pH": {}}
    orig = (eng.estimate_pH, eng.estimate_state, eng.closed_pH,
            eng.exact_proton_pH)
    eng.estimate_pH = _wrap(orig[0], cnt, "estimate_pH")
    eng.estimate_state = _wrap(orig[1], cnt, "estimate_state")
    eng.closed_pH = _wrap(orig[2], cnt, "closed_pH")
    eng.exact_proton_pH = _wrap(orig[3], cnt, "exact_proton_pH")
    try:
        t0 = time.perf_counter()
        for i, c in enumerate(cases):
            try:
                eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                          c.get("cond") or {"V_L": 1.0}, T)
            except Exception as exc:                        # noqa: BLE001
                print(f"  ERR {c['name']}: {type(exc).__name__}: {exc}")
            if (i + 1) % 25 == 0:
                print(f"  {i + 1}/{len(cases)} {time.perf_counter() - t0:.0f}s",
                      flush=True)
    finally:
        (eng.estimate_pH, eng.estimate_state, eng.closed_pH,
         eng.exact_proton_pH) = orig
        spec.PH_TAGS = None
    print(f"\n=== branch tags over {len(cases)} cases "
          f"(step={step}), walk level ===")
    for key in ("estimate_pH", "estimate_state", "closed_pH",
                "exact_proton_pH"):
        d = dict(sorted(cnt[key].items(), key=lambda kv: -kv[1]))
        tot = sum(d.values()) or 1
        print(f"\n  {key}: total calls = {sum(d.values())}")
        for k, n in d.items():
            print(f"    {n:12d}  {100.0 * n / tot:6.2f}%  {k}")
    p = wr(f"phaudit_tags_step{step}.json", cnt)
    print(f"\nwrote {p}")
    return 0


# ---------------------------------------------------------------- F31 anchor

def mode_f31(argv):
    """The proven F31 case: every producer + the He axis step."""
    rec = mode_one(["F31"])
    print("\n=== He axis (same terminal ledger) ===")
    probe = {}
    c = [x for x in load_cases(None) if x["name"].startswith("F31 ")][0]
    eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
              c.get("cond") or {"V_L": 1.0}, T, _probe=probe)
    led = ledger_clean(probe.get("ledger") or {})
    print(f"  {'He':>12} {'estimate_pH':>12} {'closed_pH':>10} "
          f"{'presentation':>13}")
    for he in (-3.0, -1.0, -0.5, -0.1, -0.05, -0.024007, -0.02, -0.01,
               -0.001, 0.0, 1e-4, 0.001, 0.1, 1.0):
        e, _ = safe(eng.estimate_pH, led, he, 1.0, T, T_K)
        cl, _ = safe(spec.closed_pH, led, he, 1.0, T, T_K)
        pr, _ = safe(eng.presentation_pH, led, he, 1.0, T, T_K)
        print(f"  {he:+12.6g} {e if e is None else round(e, 4):>12} "
              f"{str(ph_of(cl)):>10} "
              f"{pr if pr is None else round(pr, 4):>13}")
    return rec


def mode_cliff(argv):
    """Locate the exact He cliff of the heuristic on ONE ledger.

    Uses both audit hooks: PH_TAGS (which branch of estimate_state fired) and
    PH_SRC (which species/value fed h_c / o_c).  Read-only: the hooks are
    module globals that production leaves at None.
    """
    pre = argv[0] if argv else "F31"
    which = argv[1] if len(argv) > 1 else "terminal"
    cases = load_cases(None)
    c = [x for x in cases if x["name"].startswith(pre + " ")][0]
    v = float((c.get("cond") or {}).get("V_L", 1.0) or 1.0)
    tk = float((c.get("cond") or {}).get("T_K", 298.15) or 298.15)
    if which == "terminal":
        probe = {}
        eng.judge([{"name": n, "mol": m} for n, m in c["subs"]],
                  c.get("cond") or {"V_L": 1.0}, T, _probe=probe)
        led = ledger_clean(probe.get("ledger") or {})
        he0 = float(probe.get("H_excess") or 0.0)
        print(f"terminal ledger of {c['name']}, probe He = {he0}, "
              f"He_walk = {PRES_LOG[-1][0] if PRES_LOG else None}")
    else:
        import chemkit.normalize as _norm
        led, he0, _st, _un = _norm.normalize(
            [{"name": n, "mol": m} for n, m in c["subs"]],
            {"V_L": v, "T_K": tk, "c_H": None, "c_OH": None, "pH": None,
             "p_kpa": eng.P_EXT_KPA}, T, {})
        led = ledger_clean(led)
        print(f"initial ledger of {c['name']}, He0 = {he0}")
    print(f"  net = {net_charge(led)!r}  charge_pH = {ab.charge_pH(led, v, T, tk)}"
          f"  truth = {truth_of(net_charge(led), v, spec.pKw_of(tk))}")
    base = -net_charge(led)
    grid = []
    for d in (2e-6, 1e-6, 5e-7, 2e-7, 1.2e-7, 1e-7, 8e-8, 6e-8, 4e-8, 2e-8,
              0.0, -2e-8, -1e-7, -1e-6, -1e-5, -1e-4):
        grid.append(base + d)
    print(f"  fine grid around He = -net = {base!r}")
    print(f"  {'He':>22} {'estimate_pH':>12} {'tag':>26}  h_c-src / o_c-src")
    for he in grid:
        save_t, save_s = spec.PH_TAGS, spec.PH_SRC
        spec.PH_TAGS, spec.PH_SRC = [], []
        try:
            e, _ = safe(eng.estimate_pH, led, he, v, T, tk)
            tags = list(spec.PH_TAGS)
            src = list(spec.PH_SRC)
        finally:
            spec.PH_TAGS, spec.PH_SRC = save_t, save_s
        hsrc = [x for x in src if not str(x[0]).startswith("__")][:2]
        osrc = [x for x in src if str(x[0]).startswith("__")]
        print(f"  {he:+22.15g} {e if e is None else round(e, 6):>12} "
              f"{(tags[-1] if tags else '-'):>26}  "
              f"{[(a, b, round(c2, 8)) for a, b, c2 in hsrc]} / "
              f"{[(a, b, round(c2, 10)) for a, b, c2 in osrc]}")
    return 0


MODES = {"recon": mode_recon, "sweep": mode_sweep, "one": mode_one,
         "report": mode_report, "tags": mode_tags, "f31": mode_f31,
         "cliff": mode_cliff}


def main(argv):
    if not argv or argv[0] not in MODES:
        print(__doc__)
        return 2
    t0 = time.perf_counter()
    rc = MODES[argv[0]](argv[1:])
    print(f"\n[phaudit {argv[0]}] elapsed {time.perf_counter() - t0:.2f} s")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
