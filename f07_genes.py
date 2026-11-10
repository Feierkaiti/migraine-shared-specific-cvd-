"""Figure 7: tissues, drug-target gene sets and genes of the two components. Usage: python -I f07_genes.py <root>"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import fig_style as FS
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
root = sys.argv[1]
FS.apply()
S, M, GREY, INK, MUTED, GRID, BAND = FS.S, FS.M, FS.GREY, FS.INK, FS.MUTED, FS.GRID, FS.BAND
G = os.path.join(root, 'results', 'magma')
rd = lambda n: pd.read_csv(os.path.join(G, n), sep=r'\s+', comment='#')
def tissue(n):
    t = rd(n); t['name'] = t.FULL_NAME if 'FULL_NAME' in t.columns else t.VARIABLE
    return t.set_index('name')
ts, tm = tissue('COMP_S_gtex.gsa.out'), tissue('COMP_M_gtex.gsa.out')
def grp(x):
    if x.startswith('Brain'): return 0, 'Brain'
    if x.startswith('Nerve') or x.startswith('Pituitary'): return 1, 'Nerve and pituitary'
    if x.startswith('Artery'): return 2, 'Arteries'
    if x.startswith('Heart'): return 3, 'Heart'
    return 4, 'Other'
t = pd.DataFrame({'pS': ts.P, 'pM': tm.reindex(ts.index).P})
t['g'], t['gn'] = zip(*[grp(i) for i in t.index]); t = t.sort_values(['g', 'pS'])
fig = plt.figure(figsize=(7.2, 6.5))
a = fig.add_axes([0.075, 0.675, 0.91, 0.305])
x = np.arange(len(t)); bon = -np.log10(0.05 / len(t))
cols = [FS.TISSUE[g] for g in t.gn]
ls, lm = -np.log10(t.pS), -np.log10(t.pM)
a.bar(x, ls, color=cols, width=0.72, lw=0); a.bar(x, -lm * 6, color=cols, width=0.72, lw=0, alpha=0.75)
a.axhline(0, color=INK, lw=0.5); a.axhline(bon, color=MUTED, lw=0.6, ls=(0, (3, 2))); a.axhline(-bon * 6, color=MUTED, lw=0.6, ls=(0, (3, 2)))
a.set_xlim(-0.8, len(t) - 0.2); yt = [-18, -12, -6, 0, 5, 10, 15, 20]; a.set_yticks(yt); a.set_yticklabels(['3', '2', '1', '0', '5', '10', '15', '20']); a.set_ylim(-22, 20)
a.set_ylabel('−log$_{10}$(P)')
a.set_xticks(x); SHORT = {'Brain_Anterior_cingulate_cortex_BA24': 'Anterior cingulate', 'Brain_Nucleus_accumbens_basal_ganglia': 'Nucleus accumbens', 'Brain_Caudate_basal_ganglia': 'Caudate',
         'Brain_Putamen_basal_ganglia': 'Putamen', 'Brain_Spinal_cord_cervical_c-1': 'Spinal cord', 'Brain_Frontal_Cortex_BA9': 'Frontal cortex', 'Brain_Cerebellar_Hemisphere': 'Cerebellum (hem.)',
         'Brain_Substantia_nigra': 'Substantia nigra', 'Cells_Cultured_fibroblasts': 'Fibroblasts', 'Cells_EBV-transformed_lymphocytes': 'Lymphocytes', 'Skin_Not_Sun_Exposed_Suprapubic': 'Skin (unexposed)',
         'Skin_Sun_Exposed_Lower_leg': 'Skin (exposed)', 'Esophagus_Gastroesophageal_Junction': 'GE junction', 'Esophagus_Muscularis': 'Muscularis', 'Esophagus_Mucosa': 'Mucosa',
         'Small_Intestine_Terminal_Ileum': 'Small intestine', 'Adipose_Visceral_Omentum': 'Adipose (visceral)', 'Adipose_Subcutaneous': 'Adipose (subcut.)', 'Breast_Mammary_Tissue': 'Breast',
         'Heart_Atrial_Appendage': 'Heart atrium', 'Heart_Left_Ventricle': 'Heart ventricle', 'Minor_Salivary_Gland': 'Salivary gland', 'Cervix_Ectocervix': 'Ectocervix', 'Cervix_Endocervix': 'Endocervix', 'Fallopian_Tube': 'Fallopian tube'}
a.set_xticklabels([SHORT.get(i, i.replace('Brain_', '').replace('_', ' ')) for i in t.index], rotation=90, fontsize=5.6); a.tick_params(axis='x', length=0, pad=1)
a.text(0.995, 0.97, 'Shared component', transform=a.transAxes, color=S, fontsize=7.5, fontweight='bold', ha='right', va='top')
a.text(0.995, 0.03, 'Migraine-specific component (axis magnified 6-fold)', transform=a.transAxes, color=M, fontsize=7.5, fontweight='bold', ha='right', va='bottom')
a.legend(handles=[Line2D([], [], marker='s', ls='', color=c, ms=5, label=l) for l, c in FS.TISSUE.items()], loc='upper center', ncol=5, fontsize=6.4, handletextpad=0.1, columnspacing=0.9, bbox_to_anchor=(0.56, 1.03), borderaxespad=0)
FS.panel(a, 'a', x=-0.045, y=0.98)
# b drug-target gene sets
b = fig.add_axes([0.255, 0.345, 0.385, 0.185])
def sets(tag):
    p = {}
    for k in ('drugsets', 'drugsplit'):
        d = rd(f'{tag}_{k}.gsa.out'); d['n'] = d.FULL_NAME if 'FULL_NAME' in d.columns else d.VARIABLE
        for _, r in d.iterrows():
            p[r.n] = r.P
    return p
P = {k: sets(k) for k in ('MIG', 'COMP_S', 'COMP_M')}
LAB = [('migraine_drug_targets', 'Migraine drug targets (all)'), ('migraine_targets_shared_with_depression', 'Shared with depression drugs'),
       ('migraine_targets_not_depression', 'Not shared (incl. CGRP pathway)'), ('cgrp_triptan', 'CGRP and triptan/ditan receptors'), ('depression_drug_targets', 'Depression drug targets')]
ys = np.arange(len(LAB))[::-1]
for (k, c_, dy, lab) in (('MIG', GREY, 0.26, 'Migraine'), ('COMP_S', S, 0, 'Shared'), ('COMP_M', M, -0.26, 'Migraine-specific')):
    v = [-np.log10(P[k][s]) for s, _ in LAB]
    b.barh(ys + dy, v, height=0.23, color=c_, lw=0, label=lab)
b.axvline(-np.log10(0.05), color=MUTED, lw=0.6, ls=(0, (3, 2))); b.axvline(-np.log10(0.05 / len(LAB)), color=INK, lw=0.6, ls=(0, (1, 1.5)))
b.set_yticks(ys); b.set_yticklabels([l for _, l in LAB], fontsize=6.5); b.tick_params(axis='y', length=0); b.spines['left'].set_visible(False); b.set_xlabel('−log$_{10}$(P)')
b.text(-np.log10(0.05), 4.62, 'P = 0.05', fontsize=5.8, color=MUTED, ha='center'); b.text(-np.log10(0.01), 4.62, '0.01', fontsize=5.8, color=INK, ha='center')
b.legend(loc='center right', bbox_to_anchor=(1.0, 0.4), fontsize=6.4, handlelength=1, handletextpad=0.3, borderaxespad=0)
FS.panel(b, 'b', x=-0.44, y=1.0)
# c overlap of Bonferroni-significant migraine genes
c = fig.add_axes([0.735, 0.345, 0.25, 0.185])
gl = pd.read_csv(os.path.join(root, 'data', 'ref', 'genes', 'NCBI37.3.gene.loc'), sep='\t', header=None, names=['id', 'chr', 's', 'e', 'str', 'sym']).set_index('id').sym
gg = {n: pd.read_csv(os.path.join(G, f'{n}.genes.out'), sep=r'\s+').set_index('GENE') for n in ('MIG', 'COMP_S', 'COMP_M')}
sig = {n: set(g[g.P < 0.05 / len(g)].index) for n, g in gg.items()}
mg = sig['MIG']; both = mg & sig['COMP_S'] & sig['COMP_M']; onlyM = (mg & sig['COMP_M']) - sig['COMP_S']; onlyS = (mg & sig['COMP_S']) - sig['COMP_M']; none = mg - sig['COMP_S'] - sig['COMP_M']
cats = [('M only', len(onlyM), M), ('S only', len(onlyS), S), ('Both', len(both), '#8C5A7F'), ('Neither', len(none), '#C4C4BE')]
xs = np.arange(len(cats)); c.bar(xs, [v for _, v, _ in cats], color=[k for _, _, k in cats], width=0.62, lw=0)
for x_, (_, v, _) in zip(xs, cats):
    c.text(x_, v + 0.6, str(v), ha='center', fontsize=7.2, color=INK, fontweight='bold')
c.set_xticks(xs); c.set_xticklabels([l for l, _, _ in cats], fontsize=6.6); c.set_xlabel('Also significant in component GWAS'); c.tick_params(axis='x', length=0); c.set_ylabel(f'Migraine genes (of {len(mg)})'); c.set_ylim(0, max(v for _, v, _ in cats) * 1.2)
FS.panel(c, 'c', x=-0.18, y=1.0)
# d gene heat strip
e = fig.add_axes([0.075, 0.095, 0.91, 0.125])
o = gg['MIG'].loc[list(mg)].sort_values('P').index; lab = [gl.get(i, str(i)) for i in o]
mat = np.array([[-np.log10(max(gg[n].P.get(i, 1), 1e-30)) for i in o] for n in ('MIG', 'COMP_M', 'COMP_S')])
cm = LinearSegmentedColormap.from_list('s', ['#FFFFFF', '#F4B2A8', '#E64B35', '#7A1F12'])
im = e.imshow(mat, aspect='auto', cmap=cm, vmin=0, vmax=20, interpolation='nearest')
for yy, n in enumerate(('MIG', 'COMP_M', 'COMP_S')):
    for xx, i in enumerate(o):
        if i in sig[n]:
            e.text(xx, yy, '•', ha='center', va='center', fontsize=6, color='white' if mat[yy, xx] > 9 else INK)
e.set_yticks([0, 1, 2]); e.set_yticklabels(['Migraine', 'M', 'S'], fontsize=6.7); e.set_xticks(range(len(lab))); e.set_xticklabels(lab, rotation=90, fontsize=5.6, fontstyle='italic'); e.tick_params(length=0, pad=1)
for s_ in e.spines.values():
    s_.set_visible(False)
e.set_xticks(np.arange(-0.5, len(lab)), minor=True); e.set_yticks([0.5, 1.5], minor=True); e.grid(which='minor', color='white', lw=0.6); e.tick_params(which='minor', length=0)
cb = fig.add_axes([0.83, 0.256, 0.13, 0.009]); plt.colorbar(im, cax=cb, orientation='horizontal'); cb.tick_params(labelsize=5.8, length=2); cb.set_xticks([0, 10, 20]); cb.set_title('−log$_{10}$(P)  • Bonferroni significant', fontsize=5.8, pad=2, color=INK)
FS.panel(e, 'd', x=-0.045, y=1.06)
FS.save(fig, root, 'Figure7_genes_tissues_drugs')
