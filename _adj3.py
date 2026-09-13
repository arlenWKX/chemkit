"""锚定 pKw 后的失败例：逐条看净过程是"痕量平衡"还是"显著反应"。"""
import sys
sys.path.insert(0, '.')
from chemkit.data import load_tables
from chemkit.engine import judge
from chemkit.system import Reaction
from chemkit.testsuit import load_cases

T = load_tables()
names = ('16 AlCl3', '21 AgBr', 'E55', 'T64', 'NR18', 'NR19', 'NR24', 'NR48',
         'NR49', 'NR50', 'NR51', 'NR57', 'NR107', 'NR111', 'NR133', 'NR152',
         'NR153', 'NR168', 'V04', 'V05', 'V20', 'V21', 'Co33', 'EU01', 'RX12')
for c in load_cases(None):
    if not c['name'].startswith(names):
        continue
    r = judge([{'name': n, 'mol': m} for n, m in c['subs']],
              c.get('cond') or {'V_L': 1.0}, T)
    rxn = Reaction(r)
    big = [(s['kind'], s['extent'], s['equation']) for s in r.get('steps', [])
           if s['extent'] >= 1e-6]
    mx = max((e for _, e, _ in big), default=0.0)
    print(f"== {c['name'][:26]:28s} changed={int(r['changed'])} "
          f"reacted={int(bool(r.get('reacted')))} degree={r['degree']} "
          f"最大步={mx:.3g}")
    print(f"   净(精编) {rxn.net_equation.plain() if rxn.net_equation else None}")
    for k, e, q in sorted(big, key=lambda t: -t[1])[:4]:
        print(f"     [{k:9s}] {e:10.4g}  {q[:62]}")
