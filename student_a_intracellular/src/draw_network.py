import pandas as pd, subprocess, html
E = pd.read_csv('network_v3_edges.csv').fillna('')
ACT, REP, GREY, NEW, CONF = '#2e7d32', '#c62828', '#9aa0a6', '#1565c0', '#ef6c00'

def col(sign): return ACT if sign == '+' else REP if sign == '-' else '#6d4c41'
def head(sign): return 'normal' if sign == '+' else 'tee' if sign == '-' else 'odot'

# ---------- Figure 1: lac switch with numbers ----------
core = [
 ('glucose','cAMP','-','lowers cAMP'),
 ('cAMP','CRP','+','switches CRP on\n(1 cAMP per CRP, no cooperativity)'),
 ('CRP','lac','+','turns lac UP 200x (Kuhlman; no effect if crp deleted)\nNCBI RNA-seq: lacY tracks CRP (slope 1.07, r = 0.93)'),
 ('LacI','lac','-','blocks lac: model rho = 170\n(Ozbudak 170; NCBI RNA-seq 160; Kuhlman ~250 without CRP)'),
 ('CRP','loop','+','NEW: helps LacI loop\nOFF ~6x tighter'),
 ('loop','lac','-',''),
 ('lac','LacY','+','makes LacY'),
 ('LacY','ind','+','pumps inducer in\n(feedback: Hill 2.6 -> 4.5)'),
 ('ind','LacI','-','switches LacI off'),
 ('glucose','LacY','-','blocks the pump'),
]
lab = {'glucose':'Glucose','cAMP':'cAMP','CRP':'CRP','lac':'lac genes\n(lacZ lacY lacA)','LacI':'LacI\n(repressor)',
       'loop':'DNA loop','LacY':'LacY\n(pump)','ind':'Inducer\ninside cell'}
fill = {'glucose':'#fff3e0','cAMP':'#fff3e0','CRP':'#e3f2fd','lac':'#e8f5e9','LacI':'#fce4ec','loop':'#fce4ec','LacY':'#e8f5e9','ind':'#f3e5f5'}
d = ['digraph G { graph [rankdir=LR, bgcolor=white, pad=0.3, nodesep=0.6, ranksep=0.9, fontname="Helvetica"];',
     'node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=15, color="#555555", penwidth=1.3];',
     'edge [fontname="Helvetica", fontsize=11, penwidth=2.2];']
for k,v in lab.items(): d.append(f'"{k}" [label="{v}", fillcolor="{fill[k]}"];')
for s,t,sg,l in core:
    c = NEW if 'NEW' in l else col(sg)
    d.append(f'"{s}" -> "{t}" [color="{c}", fontcolor="{c}", arrowhead={head(sg)}, label="{l}"];')
d.append('{rank=same; "LacI"; "loop";} {rank=same; "lac";}')
d.append('labelloc=t; fontsize=18; label="The lac switch: every arrow has a number (Kuhlman 2007 + Ozbudak 2004 + RegulonDB + NCBI RNA-seq)\\nNCBI RNA-seq = 19 IPTG-induced E. coli MG1655 samples + 3 uninduced, NCBI SRA BioProject PRJEB36111";}')
open('fig_lac_switch_quantified.dot','w').write('\n'.join(d))
subprocess.run(['dot','-Tpng','-Gdpi=170','fig_lac_switch_quantified.dot','-o','fig_lac_switch_quantified.png'],check=True)

# ---------- Figure 2: full connected network (causal + curated; associations left out) ----------
F = E[E.tier != 'association']
nodes = sorted(set(F.source) | set(F.target))
A2 = {'hflC','tdcB','uhpT','cas3','rbsC','yahN','fecE','mhpF','flhC','flhD','araJ','rbsA','galK','rbsB','nadC','gyrA','rbsR','lamB','malK','preA','idnO'}
inputs = {'glucose','cAMP','cyaA (makes cAMP)','inducer (TMG/IPTG)','inducer inside cell'}
lacmod = {'LacI','LacI DNA loop','lac promoter (lacZYA)','LacY'}
def short(n): return {'lac promoter (lacZYA)':'lac genes','inducer (TMG/IPTG)':'inducer\n(TMG/IPTG)','inducer inside cell':'inducer\ninside','cyaA (makes cAMP)':'cyaA','LacI DNA loop':'DNA loop'}.get(n,n)
d = ['digraph G { graph [rankdir=LR, bgcolor=white, pad=0.3, nodesep=0.12, ranksep=1.1, fontname="Helvetica", newrank=true];',
     'node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=12, color="#777777", fillcolor=white];',
     'edge [penwidth=1.2, arrowsize=0.7];']
def cluster(name, label, members, fc):
    d.append(f'subgraph cluster_{name} {{ label="{label}"; fontsize=15; style="rounded,filled"; fillcolor="{fc}"; color="#cccccc";')
    for n in sorted(members): d.append(f'"{n}" [label="{short(n)}"' + (', fillcolor="#ffe0b2"' if n in {"fecE","mhpF","nadC","gyrA"} else '') + '];')
    d.append('}')
cluster('in','Signals',inputs & set(nodes),'#fff8e1')
cluster('lac','lac switch',lacmod & set(nodes),'#e8f5e9')
d.append('"CRP" [label="CRP\\n(hub)", fillcolor="#bbdefb", fontsize=16, penwidth=2];')
helpers = set(nodes) - inputs - lacmod - A2 - {'CRP','93 flagella genes','galK'}
cluster('tf','Other regulators (RegulonDB)',helpers,'#f5f5f5')
cluster('a2','21 genes CRP controls in our data (A2)',(A2 | {'93 flagella genes'}) & set(nodes),'#e3f2fd')
for _,e in F.iterrows():
    sg = e.sign
    if e.tier == 'causal_in_A2_sign_conflict': c, st, w = CONF, 'dashed', 1.6
    elif 'NEW EDGE' in e.note: c, st, w = NEW, 'bold', 2.6
    elif e.tier == 'causal': c, st, w = col(sg), 'solid', 2.0
    else: c, st, w = GREY if sg!='-' else '#e57373', 'solid', 1.0
    d.append(f'"{e.source}" -> "{e.target}" [color="{c}", style={st}, penwidth={w}, arrowhead={head(sg)}];')
d.append('legend [shape=plaintext, fillcolor=white, label=<<TABLE BORDER="0" CELLSPACING="2"><TR><TD ALIGN="LEFT"><FONT COLOR="%s"><B>green/red thick</B></FONT> = proven cause (our A2 test or Kuhlman experiment)</TD></TR><TR><TD ALIGN="LEFT"><FONT COLOR="#888888">thin</FONT> = curated database (RegulonDB) or literature</TD></TR><TR><TD ALIGN="LEFT"><FONT COLOR="%s"><B>blue</B></FONT> = new link: CRP helps LacI shut lac OFF</TD></TR><TR><TD ALIGN="LEFT"><FONT COLOR="%s"><B>orange dashed</B></FONT> = our data and the database disagree on sign</TD></TR><TR><TD ALIGN="LEFT">arrow = switches on, bar = switches off</TD></TR></TABLE>>];' % (ACT, NEW, CONF))
d.append('labelloc=t; fontsize=18; label="CRP/lac network v3: lac switch wired into the CRP network";}')
open('fig_network_v3_full.dot','w').write('\n'.join(d))
subprocess.run(['dot','-Tpng','-Gdpi=150','fig_network_v3_full.dot','-o','fig_network_v3_full.png'],check=True)
print('ok', len(nodes), 'nodes', len(F), 'edges drawn')
