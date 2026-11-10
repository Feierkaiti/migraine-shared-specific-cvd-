"""Two-sample MR: IVW (multiplicative random effects), weighted median (bootstrap SE), MR-Egger.
Usage: python -I p08_mr.py <root> <exposure_std_name|path.tsv> <instruments.tsv> <out.csv> <outcome,...>
Exposure and outcome files are standardised layouts (SNP A1 A2 BETA SE ...)."""
import sys, os, gzip
import numpy as np
from scipy.stats import norm, t as tdist
root, expo, instf, outf, outs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5].split(',')
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
inst = [l.rstrip('\n').split('\t')[0] for i, l in enumerate(open(instf)) if i]
def load(name, snps):
    path = name if os.path.exists(name) else os.path.join(root, 'data', 'std', name + '.tsv.gz')
    op = gzip.open if path.endswith('.gz') else open
    d = {}
    with op(path, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t')
        ix = {c: h.index(c) for c in ('SNP', 'A1', 'A2', 'BETA', 'SE')}
        for l in f:
            x = l.rstrip('\n').split('\t')
            if x[ix['SNP']] in snps:
                d[x[ix['SNP']]] = (x[ix['A1']], x[ix['A2']], float(x[ix['BETA']]), float(x[ix['SE']]))
    return d
E = load(expo, set(inst))
def wmedian(b, se, nboot=2000, seed=1):
    def wm(b, w):
        o = np.argsort(b); b, w = b[o], w[o] / w.sum(); c = np.cumsum(w) - 0.5 * w
        return np.interp(0.5, c, b)
    w = 1 / se ** 2; est = wm(b, w); rng = np.random.default_rng(seed)
    bs = [wm(rng.normal(b, se), w) for _ in range(nboot)]
    return est, float(np.std(bs))
rows = []
for o in outs:
    O = load(o, set(inst))
    bx, sx, by, sy, used = [], [], [], [], []
    for s in inst:
        if s not in E or s not in O:
            continue
        a1, a2, b1, s1 = E[s]; c1, c2, b2, s2 = O[s]
        if (c1, c2) == (a1, a2): pass
        elif (c1, c2) == (a2, a1): b2 = -b2
        elif {a1, a2} == {COMP[a1], COMP[a2]}: continue      # ambiguous
        elif (COMP.get(c1), COMP.get(c2)) == (a1, a2): pass
        elif (COMP.get(c1), COMP.get(c2)) == (a2, a1): b2 = -b2
        else: continue
        if b1 < 0: b1, b2 = -b1, -b2                          # orient to the risk-increasing allele
        bx.append(b1); sx.append(s1); by.append(b2); sy.append(s2); used.append(s)
    bx, sx, by, sy = map(np.array, (bx, sx, by, sy)); k = len(bx)
    if k < 3:
        continue
    ratio, rse = by / bx, sy / np.abs(bx)
    w = 1 / rse ** 2
    ivw = np.sum(w * ratio) / np.sum(w); q = np.sum(w * (ratio - ivw) ** 2)
    phi = max(1.0, q / (k - 1)); ivw_se = np.sqrt(phi / np.sum(w))
    X = np.column_stack([np.ones(k), bx]); W = np.diag(1 / sy ** 2)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ by); res = by - X @ beta
    s2 = max(1.0, float(res @ W @ res) / (k - 2)); cov = np.linalg.inv(X.T @ W @ X) * s2
    wm_est, wm_se = wmedian(ratio, rse)
    rows.append(dict(outcome=o, n_snps=k, ivw=ivw, ivw_se=ivw_se, ivw_p=2 * norm.sf(abs(ivw / ivw_se)), Q=q, Q_p=float(1 - __import__('scipy').stats.chi2.cdf(q, k - 1)),
                     wm=wm_est, wm_se=wm_se, wm_p=2 * norm.sf(abs(wm_est / wm_se)),
                     egger=beta[1], egger_se=np.sqrt(cov[1, 1]), egger_p=2 * tdist.sf(abs(beta[1] / np.sqrt(cov[1, 1])), k - 2),
                     egger_int=beta[0], egger_int_p=2 * tdist.sf(abs(beta[0] / np.sqrt(cov[0, 0])), k - 2)))
import csv
with open(outf, 'w', newline='') as f:
    wtr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wtr.writeheader(); wtr.writerows(rows)
for r in rows:
    print(f"{r['outcome']:8s} k={r['n_snps']:3d} IVW {r['ivw']:.3f} ({r['ivw_se']:.3f}) p={r['ivw_p']:.2g} | WM {r['wm']:.3f} ({r['wm_se']:.3f}) | Egger {r['egger']:.3f} p={r['egger_p']:.2g} int p={r['egger_int_p']:.2g} | Q p={r['Q_p']:.2g}")
