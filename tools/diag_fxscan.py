import chemkit.engine as eng
from chemkit.data import load_tables
from chemkit.testsuit import load_cases, run_case

eng.ROOT_AUDIT = {"case": "B26", "n": 0, "every": 1, "rec": []}
T = load_tables()
cases = [c for c in load_cases(None) if c["name"].startswith("B26")]
run_case(cases[0], T, verbose=False)
rec = eng.ROOT_AUDIT["rec"]
eng.ROOT_AUDIT = None
for i, r in enumerate(rec[:8]):
    print(f"rec{i} kind={r['kind']} dir={r['dir']} eq={r['eq'][:40]!r} "
          f"f0={r['f0']:.4g} x_max={r['x_max']:.4g} x_bis={r['x_bis']:.4g} "
          f"n_bis={r['n_bis']} flips={r['flips']} n_ill={r['n_ill']}")
print("---- rec0 61-point sweep (x, pH, tag) ----")
for t in rec[0]["trace"]:
    print(f"  x={t[0]:.4f} pH={t[1]:.4f} solid={t[2]} tag={t[3]}")
