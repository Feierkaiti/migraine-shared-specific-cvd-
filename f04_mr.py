"""Figure 6: component-stratified Mendelian randomization. Usage: python -I f04_mr.py <root>"""
import sys, os, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import fig_style as FS
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
root = sys.argv[1]
FS.apply()
S, M, GREY, INK, MUTED, GRID, BAND = FS.S, FS.M, FS.GREY, FS.INK, FS.MUTED, FS.GRID, FS.BAND
MR = os.path.join(root, 'results', 'mr'); EX = os.path.join(MR, 'robustness')
d = pd.read_csv(os.path.join(MR, 'MIG_instruments_decomp.csv'))
gl = pd.read_csv(os.path.join(root, 'data', 'ref', 'genes', 'NCBI37.3.gene.loc'), sep='\t', header=None, names=['id', 'chr', 'start', 'end', 'strand', 'symbol'], dtype={'chr': str})
gl = gl[~gl.symbol.str.match(r'^(LOC|LINC|MIR|SNORD|RNU|FAM\d)')]
def nearest(c, bp):
    g = gl[gl.chr == str(int(c))]; inside = g[(g.start <= bp) & (g.end >= bp)]
    if len(inside):
        return inside.symbol.iloc[0]
    return g.symbol.iloc[int(np.argmin(np.minimum(abs(g.start - bp), abs(g.end - bp)).values))]
d['gene'] = [nearest(c, b) for c, b in zip(d.CHR, d.BP)]
d = d.sort_values('frac_M').reset_index(drop=True)
THR = 0.75
inst = pd.read_csv(os.path.join(MR, 'MIG_instruments.tsv'), sep='\t').SNP.tolist()
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
def load(n):
    r = {}
    with gzip.open(os.path.join(root, 'data', 'std', n + '.tsv.gz'), 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); ix = {c: h.index(c) for c in h}
        for l in f:
            x = l.rstrip('\n').split('\t')
            if x[ix['SNP']] in S_:
                r[x[ix['SNP']]] = (x[ix['A1']], x[ix['A2']], float(x[ix['BETA']]), float(x[ix['SE']]))
    return r
S_ = set(inst); E, O = load('MIG'), load('CAD')
H = []
for s in inst:
    if s in E and s in O:
        a1, a2, b, se = E[s]; c1, c2, by, sy = O[s]
        if (c1, c2) == (a2, a1): by = -by
        elif (c1, c2) != (a1, a2) and {a1, a2} == {COMP[a1], COMP[a2]}: continue
        elif (COMP.get(c1), COMP.get(c2)) == (a2, a1): by = -by
        if b < 0: b, by = -b, -by
        H.append((s, b, se, by, sy))
H = pd.DataFrame(H, columns='SNP bx sx by sy'.split()).merge(d[['SNP', 'frac_M', 'group', 'gene']], on='SNP')
fig = plt.figure(figsize=(7.2, 7.6))
gs = GridSpec(3, 2, height_ratios=[1.32, 1, 0.92], width_ratios=[1, 1.32], hspace=0.4, wspace=0.42, left=0.155, right=0.985, top=0.955, bottom=0.075)
a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1])
y = np.arange(len(d)); col = np.where(d.frac_M >= THR, M, S)
a.axvspan(THR, 1.3, color='#FCE8E4', lw=0, zorder=0)
a.hlines(y, 0, d.frac_M.clip(-0.2, 1.3), color=GRID, lw=0.7, zorder=1)
a.scatter(d.frac_M, y, c=col, s=17, edgecolors='white', linewidths=0.5, zorder=3)
a.axvline(THR, color=MUTED, lw=0.6, ls=(0, (3, 2)))
a.set_yticks(y); a.set_yticklabels(d.gene, fontsize=6.2, fontstyle='italic'); a.set_xlim(0, 1.22); a.set_ylim(-0.7, len(d) - 0.3); a.tick_params(axis='y', length=0)
a.set_xlabel('Proportion of migraine association\ncarried by M'); a.spines['left'].set_visible(False)
a.text(1.0, 2.0, 'M-dominant', color=M, fontsize=6.5, ha='center', va='center', rotation=90); a.text(0.3, 2.0, 'Mixed', color=S, fontsize=6.5, ha='center', va='center')
FS.panel(a, 'a', x=-0.3, y=1.0)
# b forest: IVW and weighted median by group
sets = [('mr_MIG.csv', 'All', INK, 'D'), ('mr_MIG_inst_M.csv', 'M-dominant', M, 's'), ('mr_MIG_inst_mixed.csv', 'Mixed', S, 'o')]
outs = [('CAD', 'Coronary artery\ndisease'), ('HF', 'Heart failure'), ('IS', 'Ischaemic stroke'), ('AF', 'Atrial\nfibrillation')]
yy, ticks, labs = 0, [], []
for ii, (k, name) in enumerate(outs[::-1]):
    y0 = yy
    for f, lab, c, mk in sets[::-1]:
        r = pd.read_csv(os.path.join(MR, f)).set_index('outcome').loc[k]
        for est, se, filled in ((r.wm, r.wm_se, False), (r.ivw, r.ivw_se, True)):
            lo, hi = np.exp(est - 1.96 * se), np.exp(est + 1.96 * se)
            b.plot([lo, hi], [yy, yy], color=c, lw=1.0 if filled else 0.8, zorder=2, ls='-' if filled else (0, (2.5, 1.2)))
            b.plot(np.exp(est), yy, mk, color=c if filled else 'white', mec=c if not filled else 'white', mew=0.8 if not filled else 0.4, ms=3.8, zorder=3)
            yy += 1
    if ii % 2 == 0:
        b.axhspan(y0 - 0.65, yy - 0.35, color=BAND, zorder=0, lw=0)
    ticks.append((y0 + yy - 1) / 2); labs.append(name); yy += 0.6
b.axvline(1, color=MUTED, lw=0.6, zorder=1); b.set_xscale('log'); b.set_xticks([0.5, 0.7, 1, 1.4]); b.set_xticklabels(['0.5', '0.7', '1', '1.4']); b.xaxis.set_minor_formatter(plt.NullFormatter()); b.set_xlim(0.44, 1.55)
b.set_ylim(-0.8, yy - 0.3); b.set_yticks(ticks); b.set_yticklabels(labs, fontsize=6.8); b.tick_params(axis='y', length=0); b.spines['left'].set_visible(False)
b.set_xlabel('Odds ratio per unit higher log odds of migraine (95% CI)')
FS.panel(b, 'b', x=-0.2, y=1.0)
b.legend(handles=[Line2D([], [], marker=m, color=c, lw=1, ms=3.8, mec='white', label=l) for f, l, c, m in sets] +
         [Line2D([], [], color=MUTED, lw=1, label='IVW'), Line2D([], [], color=MUTED, lw=0.8, ls=(0, (2.5, 1.2)), marker='o', mfc='white', ms=3.5, label='Weighted median')],
         loc='lower center', bbox_to_anchor=(0.45, 1.0), ncol=5, fontsize=6.1, handlelength=1.5, columnspacing=0.8, handletextpad=0.3, borderaxespad=0.1)
# c scatter
c = fig.add_subplot(gs[1, 0])
for gname, col_, mk in (('mixed', S, 'o'), ('M-dominant', M, 's')):
    q = H[H.group == gname]
    c.errorbar(q.bx, q.by, xerr=1.96 * q.sx, yerr=1.96 * q.sy, fmt='none', ecolor=col_, elinewidth=0.4, alpha=0.5, zorder=2)
    c.scatter(q.bx, q.by, s=14, marker=mk, color=col_, edgecolors='white', linewidths=0.4, zorder=3)
    w = 1 / q.sy ** 2; beta = (w * q.bx * q.by).sum() / (w * q.bx ** 2).sum(); xs = np.array([0, 0.11]); c.plot(xs, beta * xs, color=col_, lw=0.9, zorder=2)
for g_, ox, oy, ha in (('PHACTR1', 0.0015, -0.012, 'left'), ('LRP1', 0.0025, 0.002, 'left'), ('FHL5', -0.002, -0.012, 'right')):
    r_ = H[H.gene == g_]
    if len(r_):
        c.text(r_.bx.iloc[0] + ox, r_.by.iloc[0] + oy, g_, fontsize=6, fontstyle='italic', color=INK, ha=ha, va='center', zorder=5)
c.axhline(0, color=MUTED, lw=0.5, zorder=1); c.set_xlim(0, 0.115); c.set_xlabel('Instrument effect on migraine (log odds)'); c.set_ylabel('Effect on coronary artery\ndisease (log odds)')
FS.panel(c, 'c', x=-0.3, y=1.0)
# d Wald ratio against the proportion carried by M
e = fig.add_subplot(gs[1, 1])
H['r'], H['rse'] = H.by / H.bx, H.sy / H.bx
sz = 8 + 42 * (1 / H.rse ** 2) / (1 / H.rse ** 2).max()
e.axvspan(THR, 1.05, color='#FCE8E4', lw=0, zorder=0)
e.vlines(H.frac_M.clip(0, 1), H.r - 1.96 * H.rse, H.r + 1.96 * H.rse, color=GRID, lw=0.6, zorder=1)
e.scatter(H.frac_M.clip(0, 1), H.r, s=sz, c=np.where(H.frac_M >= THR, M, S), edgecolors='white', linewidths=0.4, zorder=3)
ct = pd.read_csv(os.path.join(EX, 'continuous.csv')); x_ = np.array([0.38, 1.0])
for st, ls_, c_ in (('all', '-', INK), ('no_PHACTR1', (0, (3, 1.5)), MUTED)):
    q = ct[(ct.outcome == 'CAD') & (ct.set == st)].iloc[0]
    e.plot(x_, q.intercept + q.slope * x_, color=c_, lw=0.9, ls=ls_, zorder=2)
e.axhline(0, color=MUTED, lw=0.5, zorder=1); e.set_xlim(0.36, 1.05); e.set_ylim(-2.05, 0.75); e.set_xlabel('Proportion of migraine association carried by M'); e.set_ylabel('Wald ratio for coronary\nartery disease')
p1 = ct[(ct.outcome == 'CAD') & (ct.set == 'all')].iloc[0]; p2 = ct[(ct.outcome == 'CAD') & (ct.set == 'no_PHACTR1')].iloc[0]
e.legend(handles=[Line2D([], [], color=INK, lw=0.9, label=f'All instruments (P = {p1.p:.3f})'), Line2D([], [], color=MUTED, lw=0.9, ls=(0, (3, 1.5)), label=f'Without PHACTR1 (P = {p2.p:.3f})')], loc='lower left', fontsize=6.2, handlelength=2, borderaxespad=0.2)
r_ = H[H.gene == 'PHACTR1']
e.text(r_.frac_M.iloc[0] - 0.015, r_.r.iloc[0], 'PHACTR1', fontsize=6, fontstyle='italic', color=INK, ha='right', va='center')
FS.panel(e, 'd', x=-0.14, y=1.0)
# e thresholds
f = fig.add_subplot(gs[2, 0])
t = pd.read_csv(os.path.join(EX, 'thresholds.csv')); t = t[t.outcome == 'CAD']
xpos = np.arange(len(t))
for grp, col_, mk, dx in (('M', M, 's', -0.1), ('mixed', S, 'o', 0.1)):
    v = t[f'wm_{grp}'].values; se = t[f'wm_{grp}_se'].values
    f.vlines(xpos + dx, np.exp(v - 1.96 * se), np.exp(v + 1.96 * se), color=col_, lw=1, zorder=2); f.scatter(xpos + dx, np.exp(v), marker=mk, s=14, color=col_, edgecolors='white', linewidths=0.4, zorder=3)
f.axhline(1, color=MUTED, lw=0.5, zorder=1); f.set_yscale('log'); f.set_yticks([0.6, 0.8, 1, 1.25]); f.set_yticklabels(['0.6', '0.8', '1', '1.25']); f.yaxis.set_minor_formatter(plt.NullFormatter()); f.set_ylim(0.58, 1.4)
f.set_xticks(xpos); f.set_xticklabels([f'{int(th * 100)}%\n{int(a_)}/{int(b_)}' for th, a_, b_ in zip(t.threshold, t.k_M, t.k_mixed)], fontsize=6.2); f.tick_params(axis='x', length=0)
f.set_xlabel('Threshold for M-dominant (instruments M-dominant/mixed)'); f.set_ylabel('Weighted median odds ratio\nfor coronary artery disease')
FS.panel(f, 'e', x=-0.3, y=1.0)
# f leave-one-out for the M-dominant group
g = fig.add_subplot(gs[2, 1])
l = pd.read_csv(os.path.join(EX, 'leave_one_out.csv')); l = l[(l.outcome == 'CAD') & (l.group == 'M-dominant')].merge(d[['SNP', 'gene']], left_on='left_out', right_on='SNP').sort_values('ivw').reset_index(drop=True)
xp = np.arange(len(l))
g.vlines(xp - 0.14, np.exp(l.ivw - 1.96 * l.ivw_se), np.exp(l.ivw + 1.96 * l.ivw_se), color=M, lw=0.8, zorder=2); g.scatter(xp - 0.14, np.exp(l.ivw), marker='s', s=9, color=M, zorder=3)
g.vlines(xp + 0.14, np.exp(l.wm - 1.96 * l.wm_se), np.exp(l.wm + 1.96 * l.wm_se), color=M, lw=0.7, ls=(0, (2.5, 1.2)), zorder=2); g.scatter(xp + 0.14, np.exp(l.wm), marker='s', s=9, facecolors='white', edgecolors=M, linewidths=0.7, zorder=3)
g.axhline(1, color=MUTED, lw=0.5, zorder=1); g.set_yscale('log'); g.set_yticks([0.5, 0.7, 1]); g.set_yticklabels(['0.5', '0.7', '1']); g.yaxis.set_minor_formatter(plt.NullFormatter()); g.set_ylim(0.42, 1.12)
g.set_xticks(xp); g.set_xticklabels(l.gene, rotation=50, ha='right', fontsize=5.8, fontstyle='italic'); g.set_xlim(-0.6, len(l) - 0.4); g.tick_params(axis='x', length=0); g.set_xlabel('Instrument left out', labelpad=1)
g.set_ylabel('Odds ratio, M-dominant\ninstruments'); FS.panel(g, 'f', x=-0.14, y=1.0)
FS.save(fig, root, 'Figure6_mr_instruments')
