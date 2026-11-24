"""Build network v3: RegulonDB + A2 + Kuhlman 2007 + Ozbudak/A3, one edge table with evidence and numbers."""
import json, re, sys, openpyxl, pandas as pd

XLSX = sys.argv[1]
wb = openpyxl.load_workbook(XLSX)
rows = [r for r in wb['Combined network'].iter_rows(min_row=3, values_only=True) if r[0]]

PROT = {'crp':'CRP','fis':'Fis','hns':'H-NS','marA':'MarA','cra':'Cra','malT':'MalT','galS':'GalS','tdcA':'TdcA',
        'mhpR':'MhpR','idnR':'IdnR','fur':'Fur','araC':'AraC','uhpA':'UhpA','nac':'Nac','fliZ':'FliZ','qseB':'QseB',
        'fnr':'FNR','xylR':'XylR','lrp':'Lrp','rbsR_TF':'RbsR','flhDC (FlhDC)':'FlhDC','LacI':'LacI',
        'Inducer (TMG)':'inducer (TMG/IPTG)','LacI DNA loop':'LacI DNA loop','lac promoter':'lac promoter (lacZYA)',
        'LacY':'LacY','Glucose':'glucose','cyaA':'cyaA (makes cAMP)'}
def nm(x, as_target=False):
    if x == 'rbsR' and not as_target: return 'RbsR'
    return PROT.get(x, x)

A2 = {'hflC','tdcB','uhpT','cas3','rbsC','yahN','fecE','mhpF','flhC','flhD','araJ','rbsA','galK','rbsB','nadC','gyrA','rbsR','lamB','malK','preA','idnO'}
CONFLICT = {'fecE','mhpF','nadC','gyrA'}
edges = {}
def add(s, t, sign, ev, layer, strength='', params=None, in_a3=False, note=''):
    k = (s, t)
    if k in edges:
        e = edges[k]; e['evidence'] = sorted(set(e['evidence']) | set(ev))
        if strength and strength not in e['strength']: e['strength'] = (e['strength'] + '; ' + strength).strip('; ')
        if note and note not in e['note']: e['note'] = (e['note'] + ' ' + note).strip()
        if params: e['params'] = {**(e['params'] or {}), **params}
        e['in_a3_model'] |= in_a3
        return
    edges[k] = dict(source=s, target=t, sign=sign, layer=layer, evidence=sorted(set(ev)), strength=strength,
                    params=params, in_a3_model=in_a3, note=note)

SIGN = {'switches on (+)':'+','switches off (−)':'-','both (±)':'±'}
for f, t, sg, src, det in [r[:5] for r in rows]:
    s, tt, sign = nm(f), nm(t, True), SIGN[sg]
    det = det or ''
    if 'MINER TF–TF (LASSO)' in src and 'RegulonDB' not in src and 'A2' not in src:
        add(s, tt, sign, ['MINER_LASSO_association'], 'association', det, note='statistical, not causal'); continue
    if f in ('Glucose','cyaA','Inducer (TMG)','LacY','LacI DNA loop') or t in ('LacI DNA loop',) or (f=='LacI' and t=='lac promoter') \
       or (f=='crp' and t=='lac promoter') or (f=='lac promoter' and t=='LacY'):
        continue  # lac core rebuilt below with numbers
    ev = []
    if 'A2' in src: ev.append('A2_causal')
    m = re.search(r'RegulonDB.*?\|\s*(Confirmed|Strong|Weak)', f'{src} | {det}') or re.search(r'(Confirmed|Strong|Weak)', det)
    if 'RegulonDB' in src: ev.append(f'RegulonDB_{m.group(1) if m else "listed"}')
    if 'Metabolism' in src: ev.append('metabolism')
    if 'Gene -> protein' in src: ev.append('definition')
    if 'MINER' in src: ev.append('MINER_LASSO_association')
    layer = 'target' if (f == 'crp' and t in A2) else 'regulatory'
    note = ''
    if f == 'crp' and t == 'LacI':
        layer, note = 'regulatory', 'via mhpR-lacI read-through (mhpRp1); LacI kept constant in A3'
    if f == 'crp' and t in CONFLICT:
        note = 'SIGN CONFLICT: A2 says CRP lowers it, RegulonDB lists activation (Weak/Strong).'
    add(s, tt, sign, ev, layer, det.replace('A2 causal; ', '').replace('A2 causal', ''), note=note)

# ---------------- lac core, rebuilt with numbers -----------------------------------------------
K = 'Kuhlman2007_intervention'; O = 'Ozbudak2004'
add('glucose', 'cAMP', '-', [O, 'literature'], 'input', 'CRP activity at saturating glucose = 0.20 of max, half-effect 3.9 uM',
    {'c_min': 0.20, 'K_glc_uM': 3.9}, True, 'A3 fit to Ozbudak Fig 2c; PRECISE-1K: gat 2.2x lower on glucose')
add('cyaA (makes cAMP)', 'cAMP', '+', ['metabolism'], 'input', 'adenylate cyclase')
add('cAMP', 'CRP', '+', ['RegulonDB_effector', K], 'input', 'non-cooperative: Hill m = 1.0 +/- 0.1 (paper), 1.15 [1.03-1.28] (our refit)',
    {'C_cAMP_uM': 320, 'm': 1.0}, True, 'extracellular cAMP in delta-cyaA delta-cpdA cells')
add('CRP', 'lac promoter (lacZYA)', '+', ['RegulonDB_Confirmed_grade1', K, O], 'core',
    '201x rise from 0 to 1 mM cAMP with LacI off; crp-null: no response (fold 1.00); output ~ proportional to CRP occupancy',
    {'b0_MU_hr': 5.9, 'f_cAMP': 240, 'C_cAMP_uM': 320, 'm_cAMP': 1.0, 'basal_without_CRP_frac_of_max': 0.005}, True,
    'causal: graded intervention + crp-null control (Kuhlman Fig 1b). A3 currently sets CRP-free output to 0; Kuhlman says ~1/200 of max.')
add('LacI', 'lac promoter (lacZYA)', '-', ['RegulonDB_Confirmed_grade1', K, O], 'core',
    'no-loop repression R ~ 50; loop strength L0 = 3.7; ON/OFF ~250 without CRP',
    {'R': 50, 'L0': 3.7, 'K_IPTG_uM': 12.3, 'rho_A3': 170, 'n_A3': 2}, True,
    'Kuhlman Eq 3: activity = a0 / (1 + R/s^2 + R*L/s^4), s = 1 + I/K. A3 uses Ozbudak form (1+u^2)/(rho+u^2).')
add('inducer (TMG/IPTG)', 'LacI', '-', ['RegulonDB_effector', K, O], 'core', 'IPTG-LacI dissociation K = 12.3 uM; only weakly cooperative (Hill ~1.5 with DNA)',
    {'K_IPTG_uM': 12.3}, True, 'RegulonDB lists allolactose as the natural effector')
add('CRP', 'LacI DNA loop', '+', [K, 'literature_in_vitro'], 'core',
    'CRP bends DNA and helps LacI loop. Measured ON/OFF: 258 (0 cAMP), 228 (10 uM), 1,190 (100 uM), 1,650 (1 mM) = ~6x tighter OFF',
    {'Omega_paper': 10.3, 'Omega_refit_fig2a': 12.4, 'L0': 3.7}, False,
    'NEW EDGE. Model-free evidence: ON/OFF rises ~6x with cAMP. Omega rests on the paper route; our free joint fit of Fig 2b cannot pin it (2.8 [1.0-6.1], dAIC 0.9). Not yet in A3: needs rho(theta) and a glucose refit.')
add('LacI DNA loop', 'lac promoter (lacZYA)', '-', [K], 'core', 'looping raises IPTG Hill from ~1.5 to ~2.0 (no CRP) and 2.6-2.8 (with CRP)',
    {'L': 'L0 (1 + Omega x)/(1 + x), x = cAMP/C'}, False)
add('lac promoter (lacZYA)', 'LacY', '+', ['definition', O], 'core', 'max/basal LacY 1,506/9 molecules (A3 units)', None, True)
add('LacY', 'inducer inside cell', '+', [O, K], 'core',
    'uptake ~ TMG^0.6 (Ozbudak Fig 3d). LacY feedback steepens the IPTG response: Hill 4.5 (wild type) vs 2.6 (delta-lacY), Kuhlman Table 1',
    {'uptake_exponent': 0.6, 'k_up': 1.84}, True)
add('inducer (TMG/IPTG)', 'inducer inside cell', '+', [O], 'core', 'external -> internal', None, True)
add('inducer inside cell', 'LacI', '-', [O], 'core', 'A3: inducer-LacI Hill n = 2 (assumed; gives cusp at rho = 9)', {'n': 2}, True)
add('glucose', 'LacY', '-', [O, 'literature'], 'input', 'inducer exclusion: LacY activity 0.65 of max at saturating glucose, half-effect 7.4 uM',
    {'min': 0.65, 'K_glc_uM': 7.4}, True)
add('CRP', 'cyaA (makes cAMP)', '-', ['RegulonDB_Strong'], 'regulatory', 'negative feedback on cAMP production')

# possible indirect routes for the sign conflicts
for t, path in [('fecE', 'CRP -> Fur -| fecE gives a net negative sign, matching A2'),
                ('mhpF', 'CRP -> MhpR -> mhpF is positive, so does NOT explain A2; unresolved'),
                ('gyrA', 'CRP -| Fis -| gyrA (LASSO/RegulonDB) would be positive; unresolved'),
                ('nadC', 'Nac -| nadC; no CRP -> Nac edge in this set; unresolved')]:
    edges[('CRP', t)]['note'] += ' Indirect route: ' + path + '.'

E = pd.DataFrame(edges.values())
E['evidence'] = E.evidence.apply(lambda v: ';'.join(v))
E['params'] = E.params.apply(lambda p: json.dumps(p) if p else '')
def tier(ev):
    if 'Kuhlman2007_intervention' in ev or 'A2_causal' in ev or 'grade1' in ev: return 'causal'
    if 'RegulonDB' in ev or 'Ozbudak' in ev or 'metabolism' in ev or 'definition' in ev or 'literature' in ev: return 'curated'
    return 'association'
E['tier'] = E.evidence.apply(tier)
E.loc[E.note.str.contains('SIGN CONFLICT'), 'tier'] = 'causal_in_A2_sign_conflict'
order = {'input':0,'core':1,'regulatory':2,'target':3,'association':4}
E = E.sort_values(['layer','source','target'], key=lambda c: c.map(order) if c.name=='layer' else c).reset_index(drop=True)
E.to_csv('network_v3_edges.csv', index=False)
json.dump(json.loads(E.to_json(orient='records')), open('network_v3_edges.json','w'), indent=1, allow_nan=False)
print(E.groupby(['layer','tier']).size()); print(len(E), 'edges,', len(set(E.source)|set(E.target)), 'nodes')
