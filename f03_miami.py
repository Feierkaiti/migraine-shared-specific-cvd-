"""Figure 5: Miami plot of the shared (top) and migraine-specific (bottom) component GWAS, with lead-locus genes.
Usage: python -I f03_miami.py <root>"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root = sys.argv[1]
import fig_style as FS
FS.apply()
BLUE, BLUE2, ORANGE, ORANGE2, INK, MUTED, GRID = FS.S, FS.S_L, FS.M, FS.M_L, FS.INK, FS.MUTED, FS.GRID

g = os.path.join(root, 'results', 'gwas', 'primary')
D = {c: pd.read_csv(os.path.join(g, f'COMP_{c}_pos.tsv.gz'), sep='\t', usecols=['SNP', 'CHR', 'BP', 'P']) for c in 'SM'}
allpos = pd.concat([D['S'][['CHR', 'BP']], D['M'][['CHR', 'BP']]])
lens = allpos.groupby('CHR').BP.max().sort_index()
off = pd.Series(np.concatenate([[0], np.cumsum(lens.values[:-1] + 5e6)]), index=lens.index)
loci = pd.read_csv(os.path.join(root, 'results', 'loci', 'loci_table.csv'))
fig, ax = plt.subplots(figsize=(7.2, 4.4))
rng = np.random.default_rng(1)
ymax = 0
for c, sgn, cols in (('S', 1, (BLUE, BLUE2)), ('M', -1, (ORANGE, ORANGE2))):
    d = D[c].dropna()
    d = d[d.P > 0]
    lp = -np.log10(d.P.clip(lower=1e-300))
    keep = (lp > 2) | (rng.random(len(d)) < 0.08)        # thin the null bulk for plotting only
    d, lp = d[keep], lp[keep]
    x = d.BP.values + off.reindex(d.CHR.values).values
    col = np.where(d.CHR.values % 2 == 1, cols[0], cols[1])
    ax.scatter(x, sgn * lp.values, s=1.2, c=col, linewidths=0, rasterized=True)
    ymax = max(ymax, lp.max())
    t = loci[loci.set == f'COMP_{c}'].sort_values('P')
    t = t.drop_duplicates('gene')
    t = t.head(10) if c == 'S' else t                      # all migraine-specific loci, the ten strongest shared loci
    placed = []                                            # greedy de-overlap: push a label outwards until it is clear
    span = (off + lens).max()
    if c == 'M':                                           # migraine-specific loci: one row of labels with leader lines
        t = t.assign(x=[b + off[ch] for b, ch in zip(t.BP, t.CHR)]).sort_values('x')
        lane = -np.log10(t.P.min()) + 3.5
        xs, gap = list(t.x.astype(float)), 0.037 * span     # labels alternate between two rows, so same-row spacing is 2*gap
        for _ in range(500):                               # spread neighbours apart symmetrically around their loci
            for i in range(1, len(xs)):
                d = xs[i] - xs[i - 1]
                if d < gap:
                    xs[i - 1] -= (gap - d) / 2; xs[i] += (gap - d) / 2
            xs[0] = max(xs[0], 0.04 * span)
        for j, ((_, r), xl) in enumerate(zip(t.iterrows(), xs)):
            yy = -np.log10(max(r.P, 1e-300))
            new_locus = str(r.near_MIG) not in ('True', '1', '1.0')
            ly = lane + (2.4 if j % 2 else 0)
            ax.plot([r.x, xl], [-yy - 0.3, -ly], color=GRID, lw=0.4, zorder=2)      # leader lines under all labels
            ax.text(xl, -ly, r.gene.split(',')[0], ha='center', va='top', fontsize=6.5, fontstyle='italic',
                    color=MUTED if new_locus else INK, zorder=6, bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'))
        t = t.iloc[0:0]
    for _, r in t.iterrows():
        xx = r.BP + off[r.CHR]; yy = -np.log10(max(r.P, 1e-300))
        ly = yy + 1.0
        while any(abs(xx - px) < 0.055 * span and abs(ly - py) < 1.9 for px, py in placed):
            ly += 1.9
        placed.append((xx, ly))
        new_locus = c == 'M' and str(r.near_MIG) not in ('True', '1', '1.0')
        ax.annotate(r.gene.split(',')[0], (xx, sgn * yy), xytext=(xx, sgn * ly), textcoords='data',
                    ha='center', va='bottom' if sgn > 0 else 'top', fontsize=6.5, fontstyle='italic',
                    color=MUTED if new_locus else INK,
                    arrowprops=dict(arrowstyle='-', lw=0.4, color=MUTED, shrinkA=0, shrinkB=1) if ly - yy > 1.5 else None)
gw = -np.log10(5e-8)
for s in (1, -1):
    ax.axhline(s * gw, color=MUTED, lw=0.7, ls=(0, (3, 2)))
ax.axhline(0, color=INK, lw=0.6)
lim = min(ymax, 60) * 1.12
ax.set_ylim(-lim, lim)
yt = [t for t in np.arange(0, lim, 10 if lim > 30 else 5)]
ax.set_yticks([-t for t in yt[1:]][::-1] + yt)
ax.set_yticklabels([str(int(abs(t))) for t in [-t for t in yt[1:]][::-1] + yt])
mids = off + lens / 2
ax.set_xticks(mids.values); ax.set_xticklabels([str(c) if c <= 18 or c % 2 == 0 else '' for c in lens.index], fontsize=7)
ax.set_xlim(-1e7, (off + lens).max() + 1e7)
ax.set_xlabel('Chromosome'); ax.set_ylabel('−log$_{10}$(P)')
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.text(0.01, 0.98, 'Shared component (S)', transform=ax.transAxes, color=BLUE, fontweight='bold', va='top')
ax.text(0.01, 0.02, 'Migraine-specific component (M)', transform=ax.transAxes, color=ORANGE, fontweight='bold', va='bottom')
fig.tight_layout()
out = os.path.join(root, 'figures', 'Figure5_miami')
fig.savefig(out + '.png', dpi=300); fig.savefig(out + '.pdf', dpi=300)
print('saved', out)
