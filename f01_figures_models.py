"""Model-level figures. Figure 2: genetic correlation map and component correlations with outcomes.
Figure 4: robustness across migraine definitions and shared indicators. Usage: python -I f01_figures_models.py <root>"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import fig_style as FS
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
root = sys.argv[1]
FS.apply()
R = os.path.join(root, 'results', 'models')
S, M, GREY, INK, MUTED, GRID, BAND = FS.S, FS.M, FS.GREY, FS.INK, FS.MUTED, FS.GRID, FS.BAND
NAMES = {'CAD': 'Coronary artery disease', 'IS': 'Ischaemic stroke', 'HF': 'Heart failure', 'AF': 'Atrial fibrillation',
         'SBP': 'Systolic blood pressure', 'AS': 'Any stroke', 'CES': 'Cardioembolic stroke', 'DBP': 'Diastolic blood pressure',
         'UKB_T2D': 'Type 2 diabetes (UKB)', 'HDL': 'HDL cholesterol', 'TG': 'Triglycerides', 'BMI': 'Body mass index',
         'UKB_IHD': 'Ischaemic heart disease (UKB)', 'FG_CHD': 'Coronary heart disease (FinnGen)', 'FG_STR': 'Stroke (FinnGen)',
         'FG_T2D': 'Type 2 diabetes (FinnGen)', 'FG_HTN': 'Hypertension (FinnGen)'}
GROUPS = [('Primary outcomes', ['CAD', 'IS', 'HF', 'AF', 'SBP']),
          ('Secondary outcomes', ['AS', 'CES', 'DBP', 'UKB_T2D', 'HDL', 'TG', 'BMI']),
          ('Other outcome sources', ['UKB_IHD', 'FG_CHD', 'FG_STR', 'FG_T2D', 'FG_HTN'])]
o = pd.read_csv(os.path.join(R, 'outcomes_primary.csv')).set_index('outcome')
se_mig = pd.read_csv(os.path.join(root, 'results', 'rg_mig_se.csv')).set_index('outcome')['se']

# Figure 2
fig = plt.figure(figsize=(7.2, 4.45))
a = fig.add_axes([0.0, 0.245, 0.475, 0.685]); b = fig.add_axes([0.665, 0.095, 0.20, 0.835]); tb = fig.add_axes([0.868, 0.095, 0.127, 0.835])
bg = fig.add_axes([0.495, 0.095, 0.50, 0.835], zorder=-1); bg.axis('off')
b.patch.set_alpha(0); tb.patch.set_alpha(0); tb.axis('off')
rg = pd.read_csv(os.path.join(root, 'results', 'rg_full.csv'))
TR = [('MIG', 'Migraine', S), ('MDD25', 'Depression', FS.DEP), ('BACK', 'Back pain', FS.PAIN), ('FIB', 'Fibromyalgia', FS.PAIN),
      ('CAD', 'CAD', M), ('IS', 'IS', M), ('HF', 'HF', M), ('AF', 'AF', M), ('SBP', 'SBP', M),
      ('UKB_T2D', 'T2D', '#E9A23B'), ('TG', 'TG', '#E9A23B'), ('HDL', 'HDL', '#E9A23B'), ('BMI', 'BMI', '#E9A23B')]
n = len(TR); mat = np.full((n, n), np.nan); pv = np.full((n, n), np.nan)
for i, (x, _, _) in enumerate(TR):
    for j, (y, _, _) in enumerate(TR):
        if i == j:
            continue
        q = rg[((rg.trait1 == x) & (rg.trait2 == y)) | ((rg.trait1 == y) & (rg.trait2 == x))]
        if len(q):
            mat[i, j], pv[i, j] = q.rg.iloc[0], q.p.iloc[0]
cm = LinearSegmentedColormap.from_list('d', ['#3C5488', '#9FB0D0', '#FFFFFF', '#F4B2A8', '#E64B35'])
nm = TwoSlopeNorm(vmin=-0.8, vcenter=0, vmax=0.8)
for i in range(n):
    for j in range(i + 1):
        y0 = n - 1 - i
        if i == j:
            a.add_patch(Rectangle((j, y0), 1, 1, fc=TR[i][2], ec='white', lw=0.8))
            a.text(j + 1.12, y0 + 0.5, TR[i][1], va='center', ha='left', fontsize=6.6, color=INK, fontweight='bold' if i == 0 else 'normal')
            continue
        v = mat[i, j]
        a.add_patch(Rectangle((j, y0), 1, 1, fc=cm(nm(v)), ec='white', lw=0.8))
        sig = pv[i, j] < 0.05 / (n * (n - 1) / 2)
        a.text(j + 0.5, y0 + 0.5, '0' if abs(v) < 0.005 else f'{v:.2f}'.replace('-', '−').replace('0.', '.', 1), ha='center', va='center',
               fontsize=5.3, color=('white' if abs(v) > 0.55 else INK) if sig else '#8c8c8c', fontweight='bold' if sig else 'normal')
a.set_xlim(-0.1, n + 1.8); a.set_ylim(-0.1, n + 0.1); a.set_aspect('equal'); a.axis('off')
fig.text(0.008, 0.965, 'a', fontsize=10, fontweight='bold', color=INK, va='center')
fig.text(0.04, 0.965, 'Genetic correlations between input traits and outcomes', fontsize=7.5, color=INK, va='center')
cax = fig.add_axes([0.30, 0.135, 0.15, 0.018])
g = np.linspace(-0.8, 0.8, 256)[None, :]
cax.imshow(g, aspect='auto', cmap=cm, norm=nm, extent=[-0.8, 0.8, 0, 1])
cax.set_yticks([]); cax.set_xticks([-0.8, 0, 0.8]); cax.set_xticklabels(['−0.8', '0', '0.8'], fontsize=6)
for s_ in cax.spines.values():
    s_.set_visible(False)
fig.text(0.375, 0.165, 'Genetic correlation (r$_g$)', fontsize=6.4, color=INK, ha='center')
lg = [('Migraine', S), ('Depression', FS.DEP), ('Chronic pain', FS.PAIN), ('Cardiovascular', M), ('Cardiometabolic', '#E9A23B')]
for i, (t, c) in enumerate(lg):
    x_ = 0.02 + (i % 3) * 0.085 if i < 3 else 0.02 + (i - 3) * 0.105
    y_ = 0.165 if i < 3 else 0.13
    fig.patches.append(Rectangle((x_, y_ - 0.008), 0.011, 0.016, transform=fig.transFigure, fc=c, ec='none'))
    fig.text(x_ + 0.015, y_, t, fontsize=6.2, va='center', color=INK)
fig.text(0.02, 0.085, 'Bold values: P < 0.05 after Bonferroni correction; grey values: not significant.', fontsize=6, color=MUTED, va='center')

rows = []
for g_, items in GROUPS:
    rows.append(('head', g_)); rows += [('item', k) for k in items]
yt, yl = [], []; y = len(rows)
band = 0
tr = b.get_yaxis_transform()
for kind, k in rows:
    y -= 1
    if kind == 'head':
        b.text(-0.83, y, k, transform=tr, ha='left', va='center', fontweight='bold', color=INK, fontsize=7.2, clip_on=False); band = 0
        continue
    r = o.loc[k]
    if band % 2 == 0:
        bg.axhspan(y - 0.5, y + 0.5, color=BAND, zorder=0, lw=0)
    band += 1
    for val, se, col, mk, dy in [(r.bS, r.bS_se, S, 'o', 0.2), (r.bM, r.bM_se, M, 's', -0.2)]:
        b.plot([val - 1.96 * se, val + 1.96 * se], [y + dy] * 2, color=col, lw=1.0, solid_capstyle='butt', zorder=2)
        b.plot(val, y + dy, mk, color=col, ms=3.5, mec='white', mew=0.4, zorder=3)
    b.plot([r.rg_mig - 1.96 * se_mig[k], r.rg_mig + 1.96 * se_mig[k]], [y, y], color=GREY, lw=0.6, zorder=2)
    b.plot(r.rg_mig, y, 'D', color=GREY, ms=2.5, zorder=3)
    pt = '<0.001' if r.diff_p < 0.001 else f'{r.diff_p:.2f}' if r.diff_p >= 0.01 else f'{r.diff_p:.3f}'
    sig_ = k in ('CAD', 'IS', 'HF', 'AF', 'SBP') and r.diff_p < 0.01
    tb.text(0.13, y, f'{r.bS:.2f}'.replace('-', '−'), ha='center', va='center', fontsize=6.2, color=S)
    tb.text(0.47, y, f'{r.bM:.2f}'.replace('-', '−'), ha='center', va='center', fontsize=6.2, color=M)
    tb.text(0.83, y, pt, ha='center', va='center', fontsize=6.2, color=INK, fontweight='bold' if sig_ else 'normal')
    yt.append(y); yl.append(NAMES[k])
top = len(rows) - 0.75
tb.text(0.13, top, 'S', ha='center', fontsize=6.6, color=S, fontweight='bold', va='center')
tb.text(0.47, top, 'M', ha='center', fontsize=6.6, color=M, fontweight='bold', va='center')
tb.text(0.83, top, 'P$_{diff}$', ha='center', fontsize=6.6, color=INK, va='center')
for ax_ in (b, tb, bg):
    ax_.set_ylim(-0.7, len(rows) - 0.3)
tb.set_xlim(0, 1); bg.set_xlim(0, 1)
b.axvline(0, color=MUTED, lw=0.6, zorder=1)
b.set_yticks(yt); b.set_yticklabels(yl, fontsize=6.7); b.tick_params(axis='y', length=0)
b.set_xlim(-0.32, 0.6); b.set_xticks([-0.2, 0, 0.2, 0.4]); b.set_xticklabels(['−0.2', '0', '0.2', '0.4'])
b.set_xlabel('Genetic correlation (95% CI)'); b.spines['left'].set_visible(False)
fig.text(0.492, 0.965, 'b', fontsize=10, fontweight='bold', color=INK, va='center')
h = [Line2D([], [], color=S, marker='o', lw=1, ms=3.5, mec='white', label='Shared'), Line2D([], [], color=M, marker='s', lw=1, ms=3.5, mec='white', label='Migraine-specific'),
     Line2D([], [], color=GREY, marker='D', lw=0.6, ms=2.5, label='Migraine (total)')]
fig.legend(handles=h, loc='center', bbox_to_anchor=(0.745, 0.965), ncol=3, fontsize=6.6, handlelength=1.4, columnspacing=1.0, handletextpad=0.4, borderaxespad=0)
FS.save(fig, root, 'Figure2_components_outcomes')
plt.close(fig)

# Figure 4
LAB = [('primary', 'Primary analysis'),
       ('def_FG_MIG', 'FinnGen diagnosis only'), ('def_MVP_MIG', 'MVP only'), ('def_FG_MIGTRIP', 'FinnGen, incl. triptan purchase'),
       ('def_CHOQ', 'UK Biobank + GERA'), ('def_FG_MA', 'With aura (FinnGen)'), ('def_FG_MO', 'Without aura (FinnGen)'), ('def_MVP_MA', 'With aura (MVP)'),
       ('sh_dep', 'Depression only'), ('sh_pain', 'Back pain + fibromyalgia only'), ('sh_dep_back', 'Depression + back pain'),
       ('sh_depclin', 'Clinical-cohort depression'), ('sh_nofib', 'No fibromyalgia; + anxiety, IBS'), ('sh_wide', 'Primary set + anxiety, IBS')]
base = {t: pd.read_csv(os.path.join(R, f'base_{t}.csv')).iloc[0] for t, _ in LAB}
out = {t: pd.read_csv(os.path.join(R, f'outcomes_{t}.csv')).set_index('outcome') for t, _ in LAB}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 4.25), sharey=True, gridspec_kw={'width_ratios': [1, 1.2, 1.2]})
ys = np.arange(len(LAB))[::-1].astype(float)
ys[1:] -= 0.5; ys[8:] -= 0.5
for ax in axs:
    ax.axhspan(ys[1] + 0.5, ys[7] - 0.5, color=BAND, zorder=0, lw=0)
    ax.axhspan(ys[0] - 0.42, ys[0] + 0.5, color='#E8ECF4', zorder=0, lw=0)
for i, (tag, lab) in enumerate(LAB):
    b_ = base[tag]; bold = tag == 'primary'
    axs[0].plot([b_.share - 1.96 * b_.share_se, b_.share + 1.96 * b_.share_se], [ys[i]] * 2, color=S, lw=1.1 if bold else 0.9, zorder=2)
    axs[0].plot(b_.share, ys[i], 'o', color=S, ms=4 if bold else 3.3, mec='white', mew=0.4, zorder=3)
    axs[0].text(0.66, ys[i], f'{b_.share * 100:.0f}%', ha='right', va='center', fontsize=6.2, color=INK)
    for j, yv in enumerate(['CAD', 'IS']):
        if yv not in out[tag].index:
            continue
        r = out[tag].loc[yv]
        for val, se, col, mk, dy in [(r.bS, r.bS_se, S, 'o', 0.19), (r.bM, r.bM_se, M, 's', -0.19)]:
            axs[j + 1].plot([val - 1.96 * se, val + 1.96 * se], [ys[i] + dy] * 2, color=col, lw=0.95, zorder=2)
            axs[j + 1].plot(val, ys[i] + dy, mk, color=col, ms=3.3, mec='white', mew=0.4, zorder=3)
axs[0].set_yticks(ys); axs[0].set_yticklabels([l for _, l in LAB], fontsize=6.7)
axs[0].get_yticklabels()[0].set_fontweight('bold')
axs[0].set_xlabel('Share of migraine SNP heritability\nin the shared component'); axs[0].set_xlim(0, 0.68); axs[0].set_xticks([0, 0.2, 0.4])
for j, t in enumerate(['coronary artery disease', 'ischaemic stroke']):
    axs[j + 1].axvline(0, color=MUTED, lw=0.6, zorder=1); axs[j + 1].set_xlim(-0.4, 0.55); axs[j + 1].set_xticks([-0.25, 0, 0.25, 0.5])
    axs[j + 1].set_xticklabels(['−0.25', '0', '0.25', '0.5']); axs[j + 1].set_xlabel(f'Genetic correlation with\n{t}')
for k, (t, ax) in enumerate(zip('abc', axs)):
    ax.tick_params(axis='y', length=0); ax.spines['left'].set_visible(False)
    ax.text(-0.02, 1.015, t, transform=ax.transAxes, fontweight='bold', fontsize=10, ha='right', va='bottom', color=INK)
fig.legend(handles=[Line2D([], [], color=S, marker='o', lw=1, ms=3.5, mec='white', label='Shared component'), Line2D([], [], color=M, marker='s', lw=1, ms=3.5, mec='white', label='Migraine-specific component')],
           loc='upper center', bbox_to_anchor=(0.62, 1.0), ncol=2, fontsize=7, borderaxespad=0.2)
fig.text(0.012, 0.615, 'Migraine definition', rotation=90, fontsize=6.8, color=MUTED, va='center')
fig.text(0.012, 0.275, 'Shared indicators', rotation=90, fontsize=6.8, color=MUTED, va='center')
fig.subplots_adjust(left=0.235, right=0.985, bottom=0.135, top=0.905, wspace=0.11)
FS.save(fig, root, 'Figure4_robustness')
