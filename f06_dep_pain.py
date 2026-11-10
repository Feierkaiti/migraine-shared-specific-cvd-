"""Figure 3: depression, chronic pain and migraine-specific liability, and adjustment for risk-factor liabilities.
Usage: python -I f06_dep_pain.py <root>"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import fig_style as FS
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
root = sys.argv[1]
FS.apply()
R = os.path.join(root, 'results', 'models')
S, M, DEP, PAIN, GREY, INK, MUTED, BAND = FS.S, FS.M, FS.DEP, FS.PAIN, FS.GREY, FS.INK, FS.MUTED, FS.BAND
base = {t: pd.read_csv(os.path.join(R, f'base_{t}.csv')).iloc[0] for t in ('primary', 'sh_dep', 'sh_pain')}
b3 = pd.read_csv(os.path.join(R, 'base_dep_pain.csv')).set_index('component')
o3 = pd.read_csv(os.path.join(R, 'outcomes_dep_pain.csv')).set_index('outcome')
fig = plt.figure(figsize=(7.2, 5.0))
gs = GridSpec(2, 6, height_ratios=[1, 1.16], hspace=0.52, wspace=1.15, left=0.135, right=0.985, top=0.945, bottom=0.095)
a = fig.add_subplot(gs[0, :3]); b = fig.add_subplot(gs[0, 3:])
# a: stacked partition of migraine SNP heritability
rows = [('Primary model', [(base['primary'].share, S, 'Shared'), (1 - base['primary'].share, M, 'Specific')]),
        ('Depression only', [(base['sh_dep'].share, DEP, 'Dep.'), (1 - base['sh_dep'].share, M, 'Specific')]),
        ('Pain only', [(base['sh_pain'].share, PAIN, 'Pain'), (1 - base['sh_pain'].share, M, 'Specific')]),
        ('Three components', [(b3.loc['D', 'share'], DEP, 'Dep.'), (b3.loc['P', 'share'], PAIN, 'Pain'), (b3.loc['M', 'share'], M, 'Specific')])]
for i, (lab, parts) in enumerate(rows):
    y = len(rows) - 1 - i; x = 0
    for w, c, t in parts:
        a.barh(y, w - 0.006, left=x + 0.003, color=c, height=0.62, lw=0)
        a.text(x + w / 2, y, f'{w * 100:.0f}%', ha='center', va='center', color='white', fontsize=7, fontweight='bold')
        x += w
a.set_yticks(range(len(rows))); a.set_yticklabels([r[0] for r in rows][::-1]); a.tick_params(axis='y', length=0)
a.set_xlim(0, 1); a.set_xticks([0, 0.25, 0.5, 0.75, 1]); a.set_xticklabels(['0', '25', '50', '75', '100'])
a.set_xlabel('Share of migraine SNP heritability (%)'); a.spines['left'].set_visible(False)
FS.panel(a, 'a', x=-0.28)
hl = [Line2D([], [], marker='s', ls='', color=c, ms=5, label=l) for c, l in ((S, 'Shared'), (DEP, 'Depression'), (PAIN, 'Chronic pain'), (M, 'Migraine-specific'))]
a.legend(handles=hl, loc='lower center', bbox_to_anchor=(0.42, 1.0), ncol=4, fontsize=6.4, handletextpad=0.1, columnspacing=0.7, borderaxespad=0.2)
# b: three-component coefficients
OUT = [('CAD', 'CAD'), ('HF', 'HF'), ('IS', 'IS'), ('AF', 'AF'), ('UKB_T2D', 'T2D'), ('BMI', 'BMI')]
for i, (k, lab) in enumerate(OUT):
    r = o3.loc[k]
    if i % 2 == 0:
        b.axvspan(i - 0.5, i + 0.5, color=BAND, lw=0, zorder=0)
    for v, se, c, mk, dx in ((r.bD, r.bD_se, DEP, 'o', -0.24), (r.bP, r.bP_se, PAIN, 'D', 0), (r.bM, r.bM_se, M, 's', 0.24)):
        b.plot([i + dx] * 2, [v - 1.96 * se, v + 1.96 * se], color=c, lw=1.1, zorder=2)
        b.plot(i + dx, v, mk, color=c, ms=3.9 if mk != 'D' else 3.5, mec='white', mew=0.4, zorder=3)
b.axhline(0, color=MUTED, lw=0.6, zorder=1)
b.set_xticks(range(len(OUT))); b.set_xticklabels([l for _, l in OUT]); b.set_xlim(-0.5, len(OUT) - 0.5); b.tick_params(axis='x', length=0)
b.set_ylabel('Genetic correlation (95% CI)'); b.set_ylim(-0.32, 0.42); b.set_yticks([-0.2, 0, 0.2, 0.4]); b.set_yticklabels(['−0.2', '0', '0.2', '0.4'])
FS.panel(b, 'b', x=-0.18)
b.legend(handles=[Line2D([], [], color=DEP, marker='o', lw=1, ms=3.9, mec='white', label='Depression'), Line2D([], [], color=PAIN, marker='D', lw=1, ms=3.5, mec='white', label='Pain beyond depression'),
                  Line2D([], [], color=M, marker='s', lw=1, ms=3.9, mec='white', label='Migraine-specific')], loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=6.4, handletextpad=0.3, columnspacing=0.8, handlelength=1.4, borderaxespad=0.2)
# c: effect of adjustment
TAGS = [('primary', 'None'), ('adj_edu_smk', 'Educ.,\nsmoking'), ('adj_bmi', 'BMI'), ('adj_bmi_t2d', 'BMI,\nT2D'), ('adj_all', 'All\nfour')]
O = {t: pd.read_csv(os.path.join(R, f'outcomes_{t}.csv')).set_index('outcome') for t, _ in TAGS}
for j, (k, name) in enumerate((('CAD', 'Coronary artery disease'), ('HF', 'Heart failure'), ('IS', 'Ischaemic stroke'))):
    c = fig.add_subplot(gs[1, 2 * j:2 * j + 2])
    xs = np.arange(len(TAGS))
    for key, col, mk, dx in (('S', S, 'o', -0.09), ('M', M, 's', 0.09)):
        v = np.array([O[t].loc[k, f'b{key}'] for t, _ in TAGS]); se = np.array([O[t].loc[k, f'b{key}_se'] for t, _ in TAGS])
        c.plot(xs[:4] + dx, v[:4], color=col, lw=0.9, zorder=2); c.plot(xs[3:] + dx, v[3:], color=col, lw=0.9, ls=(0, (2, 1.5)), zorder=2)
        c.vlines(xs + dx, v - 1.96 * se, v + 1.96 * se, color=col, lw=0.9, zorder=2)
        c.scatter(xs[:4] + dx, v[:4], marker=mk, s=15, color=col, edgecolors='white', linewidths=0.4, zorder=3)
        c.scatter(xs[4:] + dx, v[4:], marker=mk, s=15, facecolors='white', edgecolors=col, linewidths=0.9, zorder=3)
    c.axhline(0, color=MUTED, lw=0.6, zorder=1); c.set_ylim(-0.27, 0.46); c.set_yticks([-0.2, 0, 0.2, 0.4]); c.set_yticklabels(['−0.2', '0', '0.2', '0.4'])
    c.set_xticks(xs); c.set_xticklabels([l for _, l in TAGS], fontsize=6.2); c.set_xlim(-0.45, len(TAGS) - 0.55); c.tick_params(axis='x', length=0)
    c.set_title(name, fontsize=7.5, color=INK, pad=3)
    if j == 0:
        c.set_ylabel('Coefficient (95% CI)'); FS.panel(c, 'c', x=-0.2, y=1.05)
    if j == 1:
        c.set_xlabel('Risk-factor liabilities added to the model')
    if j == 2:
        c.legend(handles=[Line2D([], [], color=S, marker='o', lw=0.9, ms=3.5, mec='white', label='Shared'), Line2D([], [], color=M, marker='s', lw=0.9, ms=3.5, mec='white', label='Migraine-specific')],
                 loc='upper right', fontsize=6.4, handlelength=1.4, handletextpad=0.3, borderaxespad=0)
FS.save(fig, root, 'Figure3_depression_pain')
