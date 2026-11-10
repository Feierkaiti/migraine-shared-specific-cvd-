"""Shared figure style: colours and matplotlib settings for all figures."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

S, M, DEP, PAIN, GREY = '#3C5488', '#E64B35', '#4DBBD5', '#00A087', '#7A7A7A'
S_L, M_L = '#A9B5CF', '#F4B2A8'
INK, MUTED, GRID, BAND = '#1A1A1A', '#595959', '#D9D9D6', '#F4F4F1'
TISSUE = {'Brain': '#3C5488', 'Nerve and pituitary': '#6F86B5', 'Arteries': '#E64B35', 'Heart': '#F39B7F', 'Other': '#B5B5B0'}


def apply():
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'axes.linewidth': 0.6,
                         'xtick.color': MUTED, 'ytick.color': INK, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                         'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'legend.frameon': False, 'axes.spines.top': False, 'axes.spines.right': False})


def panel(ax, letter, x=-0.02, y=1.03):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', color=INK, va='bottom', ha='right')


def save(fig, root, name, dpi=300):
    import os
    p = os.path.join(root, 'figures', name)
    fig.savefig(p + '.png', dpi=dpi); fig.savefig(p + '.pdf', dpi=dpi)
    print('saved', p)
