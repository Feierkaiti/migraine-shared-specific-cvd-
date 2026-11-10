"""Figure 1: decomposition model with estimates (a) and study workflow (b). Usage: python -I f02_figure1.py <root>"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import fig_style as FS
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Rectangle, Polygon
root = sys.argv[1]
FS.apply()
S, M, DEP, PAIN, GREY, INK, MUTED, GRID = FS.S, FS.M, FS.DEP, FS.PAIN, FS.GREY, FS.INK, FS.MUTED, FS.GRID
b0 = pd.read_csv(os.path.join(root, 'results', 'models', 'base_primary.csv')).iloc[0]
lo = pd.read_csv(os.path.join(root, 'results', 'loadings_primary.csv')).set_index('rhs').STD_Genotype
o = pd.read_csv(os.path.join(root, 'results', 'models', 'outcomes_primary.csv')).set_index('outcome')
fig = plt.figure(figsize=(7.2, 4.1))
A = fig.add_axes([0.0, 0.0, 0.53, 0.94]); B = fig.add_axes([0.545, 0.0, 0.455, 0.94])
for ax in (A, B):
    ax.axis('off'); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
ax_ratio = (0.53 * 7.2) / (0.94 * 4.1)
fig.text(0.008, 0.965, 'a', fontsize=10, fontweight='bold', color=INK, va='center'); fig.text(0.04, 0.965, 'Partitioning migraine liability', fontsize=7.5, va='center', color=INK)
fig.text(0.552, 0.965, 'b', fontsize=10, fontweight='bold', color=INK, va='center'); fig.text(0.584, 0.965, 'Study workflow', fontsize=7.5, va='center', color=INK)


def box(ax, x, y, w, h, t, fc='white', ec=INK, fs=7.2, bold=False, col=INK, lw=0.9):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0,rounding_size=1.4', fc=fc, ec=ec, lw=lw, zorder=3, mutation_aspect=ax_ratio))
    ax.text(x, y, t, ha='center', va='center', fontsize=fs, fontweight='bold' if bold else 'normal', color=col, zorder=4, linespacing=1.05)


def arr(ax, p0, p1, c, lw=1.3, style='-|>', cs=None):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=7, color=c, lw=lw, shrinkA=0, shrinkB=0, zorder=2, connectionstyle=cs or 'arc3'))


cS, cM, r = (28, 53), (73, 53), 8.2
for c_, col, t in ((cS, S, 'S'), (cM, M, 'M')):
    A.add_patch(plt.matplotlib.patches.Ellipse(c_, 2 * r, 2 * r * ax_ratio, fc=col, ec='none', zorder=3))
    A.text(*c_, t, color='white', fontsize=12, fontweight='bold', ha='center', va='center', zorder=4)
A.text(cS[0] - r - 1.5, cS[1] + 9, 'Shared\ncomponent', ha='right', va='center', fontsize=6.8, color=S, linespacing=1.05)
A.text(cM[0] + r + 1.5, cM[1] + 9, 'Migraine-\nspecific\ncomponent', ha='left', va='center', fontsize=6.8, color=M, linespacing=1.05)
A.annotate('', xy=(cM[0] - r, cM[1]), xytext=(cS[0] + r, cS[1]), arrowprops=dict(arrowstyle='<|-|>', color=GREY, lw=0.8, ls=(0, (3, 2)), mutation_scale=6))
A.text((cS[0] + cM[0]) / 2, cS[1] + 2, 'r = 0', ha='center', fontsize=6.6, color=MUTED)
ind = [('Depression', 12, lo['MDD25'], DEP), ('Back pain', 33, lo['BACK'], PAIN), ('Fibro-\nmyalgia', 54, lo['FIB'], PAIN), ('Migraine', 80, b0.loadS, S)]
yb, h = 16, 12
for t, x, l, c in ind:
    box(A, x, yb, 18, h, t, bold=t == 'Migraine', fc='#FAFAF8', lw=1.0)
    dx, dy = (x - cS[0]), (yb + h / 2 - cS[1]) / ax_ratio
    d = (dx ** 2 + dy ** 2) ** 0.5
    p0 = (cS[0] + r * dx / d, cS[1] + r * ax_ratio * dy / d)
    arr(A, p0, (x, yb + h / 2), S, lw=1.1 + 0.9 * abs(l))
    tt = 0.62 if t != 'Migraine' else 0.42
    mx, my = p0[0] + (x - p0[0]) * tt, p0[1] + (yb + h / 2 - p0[1]) * tt
    A.text(mx, my, f'{l:.2f}', ha='center', va='center', fontsize=6.6, color=S, zorder=5, fontweight='bold', bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none'))
    res = 1 - l ** 2 if t != 'Migraine' else 0
    A.text(x, yb - h / 2 - 3.2, f'residual {res:.2f}'.replace('0.', '0.') if t != 'Migraine' else 'residual 0', ha='center', va='center', fontsize=5.8, color=MUTED)
arr(A, (cM[0] + 1.5, cM[1] - r * ax_ratio), (cM[0] + 7, yb + h / 2), M, lw=1.1 + 0.9 * b0.loadM)
A.text(cM[0] + 7.5, 33, f'{b0.loadM:.2f}', ha='center', va='center', fontsize=6.6, color=M, fontweight='bold', zorder=5, bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none'))
box(A, 50.5, 90, 56, 10, 'Cardiovascular or cardiometabolic outcome', fc='#F1F1EE', ec=GREY, fs=6.8)
for c_, col, tgt, lab, ox in ((cS, S, (36, 85), f"bS = {o.loc['CAD','bS']:.2f}", -8.5), (cM, M, (65, 85), f"bM = {o.loc['CAD','bM']:.2f}".replace('-', '−'), 8.5)):
    dx, dy = (tgt[0] - c_[0]), (tgt[1] - c_[1]) / ax_ratio
    d = (dx ** 2 + dy ** 2) ** 0.5
    arr(A, (c_[0] + r * dx / d, c_[1] + r * ax_ratio * dy / d), tgt, col)
    A.text((c_[0] + tgt[0]) / 2 + ox, (c_[1] + tgt[1]) / 2 + 2, lab, color=col, fontsize=6.6, ha='center', va='center')
A.text(50.5, 79, 'values shown for\ncoronary artery disease', fontsize=5.8, color=MUTED, ha='center', va='center', linespacing=1.05)
A.text(50, 2.0, f"Migraine SNP heritability: shared {b0.share * 100:.0f}%, migraine-specific {(1 - b0.share) * 100:.0f}% (λS² / (λS² + λM²))", ha='center', va='center', fontsize=6.6, color=INK)

# Workflow with small pictograms
steps = [('Input GWAS', 'Migraine (58,730 cases), depression, back pain,\nfibromyalgia; 35 European-ancestry datasets', '#E8ECF4', S, 'data'),
         ('Genomic SEM', 'Share of migraine SNP heritability; S and M against\ncardiovascular outcomes; 13 sensitivity models', '#FCE8E4', M, 'sem'),
         ('Is S only depression?', 'Pain-only and depression-only models; three\ncomponents; adjustment for risk-factor liabilities', '#E0F1EF', PAIN, 'bars'),
         ('Component GWAS', '6.57 million variants: loci, genes (MAGMA),\n54 GTEx tissues, drug-target gene sets', '#FBF1DF', '#C98A1E', 'gwas'),
         ('Component-stratified MR', 'Instruments grouped by the proportion carried by M;\nrobust estimators, directionality, power', '#ECECF2', '#57628A', 'mr')]
ys = np.linspace(87, 11, len(steps))
for i, (t, d, fc, ec, ic) in enumerate(steps):
    y = ys[i]
    B.add_patch(FancyBboxPatch((2, y - 7.8), 95, 15.6, boxstyle='round,pad=0,rounding_size=1.6', fc=fc, ec='none', zorder=1, mutation_aspect=0.455 * 7.2 / (0.94 * 4.1)))
    B.add_patch(Rectangle((2, y - 7.8), 1.3, 15.6, fc=ec, ec='none', zorder=2))
    B.text(6, y + 3.8, t, fontsize=7.2, fontweight='bold', color=INK, va='center')
    B.text(6, y - 2.6, d, fontsize=5.9, color=MUTED, va='center', linespacing=1.15)
    x0, w0 = 82, 13
    if ic == 'data':
        for k in range(4):
            B.add_patch(Rectangle((x0 + k * 3.4, y - 4.5), 2.6, 9, fc='white', ec=ec, lw=0.6, zorder=3))
            for q in range(3):
                B.plot([x0 + k * 3.4 + 0.5, x0 + k * 3.4 + 2.1], [y + 2.6 - q * 2.5] * 2, color=ec, lw=0.5, zorder=4)
    elif ic == 'sem':
        B.add_patch(plt.matplotlib.patches.Ellipse((x0 + 3.5, y + 2.6), 5.2, 4.42, fc=S, ec='none', zorder=3)); B.add_patch(plt.matplotlib.patches.Ellipse((x0 + 10, y + 2.6), 5.2, 4.42, fc=M, ec='none', zorder=3))
        for xx in (x0 + 1, x0 + 6.5, x0 + 12):
            B.add_patch(Rectangle((xx - 1.3, y - 6), 2.6, 3, fc='white', ec=INK, lw=0.5, zorder=3))
        for xx in (x0 + 1, x0 + 6.5, x0 + 12):
            B.plot([x0 + 3.5, xx], [y + 1, y - 3], color=S, lw=0.5, zorder=2)
        B.plot([x0 + 10, x0 + 12], [y + 1, y - 3], color=M, lw=0.5, zorder=2)
    elif ic == 'bars':
        vals = [0.23, 0.16, 0.61]; xx = x0
        for v, c in zip(vals, (DEP, PAIN, M)):
            B.add_patch(Rectangle((xx, y - 2), w0 * v, 4.2, fc=c, ec='white', lw=0.4, zorder=3)); xx += w0 * v
    elif ic == 'gwas':
        rng = np.random.default_rng(3); xs = np.linspace(x0, x0 + w0, 44)
        hs = rng.exponential(1.1, 44); hs[[8, 21, 33]] = [5.5, 6.7, 5.0]
        B.vlines(xs, y - 5, y - 5 + hs + 0.3, color=ec, lw=0.6, zorder=3); B.plot([x0, x0 + w0], [y + 0.5] * 2, color=INK, lw=0.4, ls=(0, (2, 1.5)), zorder=4)
    else:
        for k, (m_, ci_) in enumerate(((0.3, 0.2), (0.75, 0.28), (0.55, 0.16))):
            yy = y + 3.5 - k * 3.5; c = M if k == 0 else (S if k == 2 else '#57628A')
            B.plot([x0 + w0 * (m_ - ci_), x0 + w0 * (m_ + ci_)], [yy, yy], color=c, lw=0.7, zorder=3); B.plot(x0 + w0 * m_, yy, 's', color=c, ms=2.6, zorder=4)
        B.plot([x0 + w0 * 0.52] * 2, [y - 5.5, y + 5.5], color=MUTED, lw=0.4, zorder=2)
FS.save(fig, root, 'Figure1_design')
