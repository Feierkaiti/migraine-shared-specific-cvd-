"""Robustness of the Mendelian randomization by component: leave-one-out, alternative thresholds, continuous analysis,
uncertainty in the grouping, directionality (Steiger), power.
Usage: python -I p11_mr_robustness.py <root>"""
import sys, os, gzip
import numpy as np, pandas as pd
from scipy.stats import norm, chi2
root = sys.argv[1]
OUT = ['CAD', 'IS', 'HF', 'AF']
mr = os.path.join(root, 'results', 'mr'); ext = os.path.join(mr, 'robustness'); os.makedirs(ext, exist_ok=True)
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
inst = pd.read_csv(os.path.join(mr, 'MIG_instruments.tsv'), sep='\t').SNP.tolist()
dec = pd.read_csv(os.path.join(mr, 'MIG_instruments_decomp.csv'))
def load(n, snps):
    d = {}
    with gzip.open(os.path.join(root, 'data', 'std', n + '.tsv.gz'), 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); ix = {c: h.index(c) for c in h}
        for l in f:
            x = l.rstrip('\n').split('\t')
            if x[ix['SNP']] in snps:
                d[x[ix['SNP']]] = (x[ix['A1']], x[ix['A2']], float(x[ix['BETA']]), float(x[ix['SE']]), float(x[ix['N']]))
    return d
E = load('MIG', set(inst))
def harm(o):
    O = load(o, set(inst)); r = []
    for s in inst:
        if s not in E or s not in O:
            continue
        a1, a2, b1, s1, n1 = E[s]; c1, c2, b2, s2, n2 = O[s]
        if (c1, c2) == (a1, a2): pass
        elif (c1, c2) == (a2, a1): b2 = -b2
        elif {a1, a2} == {COMP[a1], COMP[a2]}: continue
        elif (COMP.get(c1), COMP.get(c2)) == (a1, a2): pass
        elif (COMP.get(c1), COMP.get(c2)) == (a2, a1): b2 = -b2
        else: continue
        if b1 < 0: b1, b2 = -b1, -b2
        r.append((s, b1, s1, b2, s2, n1, n2))
    return pd.DataFrame(r, columns='SNP bx sx by sy Nx Ny'.split())
H = {o: harm(o).merge(dec[['SNP', 'est_S', 'SE_S', 'est_M', 'SE_M', 'frac_M', 'group']], on='SNP', how='left') for o in OUT}
lS, lM = pd.read_csv(os.path.join(root, 'results', 'models', 'base_primary.csv')).iloc[0][['loadS', 'loadM']]
def ivw(d):
    r, se = (d.by / d.bx).values, (d.sy / d.bx.abs()).values; w = 1 / se ** 2; k = len(r)
    est = np.sum(w * r) / np.sum(w); q = np.sum(w * (r - est) ** 2); phi = max(1.0, q / (k - 1))
    return est, np.sqrt(phi / np.sum(w)), q, chi2.sf(q, k - 1) if k > 1 else np.nan, k, np.sqrt(1 / np.sum(w))
def wm(d, nboot=2000, seed=1):
    r, se = (d.by / d.bx).values, (d.sy / d.bx.abs()).values; w = 1 / se ** 2
    def f(b):
        o = np.argsort(b); bb, ww = b[o], w[o] / w.sum(); c = np.cumsum(ww) - 0.5 * ww
        return np.interp(0.5, c, bb)
    rng = np.random.default_rng(seed)
    return f(r), float(np.std([f(rng.normal(r, se)) for _ in range(nboot)]))
GROUP = {'all': lambda d: d, 'M-dominant': lambda d: d[d.group == 'M-dominant'], 'mixed': lambda d: d[d.group == 'mixed']}
# 1 leave-one-out
rows = []
for o in OUT:
    for g, fn in GROUP.items():
        d = fn(H[o])
        for s in d.SNP:
            x = d[d.SNP != s]; a = ivw(x); m = wm(x, nboot=500)
            rows.append(dict(outcome=o, group=g, left_out=s, k=a[4], ivw=a[0], ivw_se=a[1], Q=a[2], wm=m[0], wm_se=m[1]))
loo = pd.DataFrame(rows); loo.to_csv(os.path.join(ext, 'leave_one_out.csv'), index=False)
# contribution of each instrument to Cochran's Q, and Q with the two vascular loci removed
rows = []
for o in OUT:
    for g, fn in GROUP.items():
        d = fn(H[o]); r, se = (d.by / d.bx).values, (d.sy / d.bx.abs()).values; w = 1 / se ** 2; est = np.sum(w * r) / np.sum(w)
        for s, c in zip(d.SNP, w * (r - est) ** 2):
            rows.append(dict(outcome=o, group=g, SNP=s, ratio=float(r[list(d.SNP).index(s)]), Q_contrib=c))
pd.DataFrame(rows).to_csv(os.path.join(ext, 'q_contributions.csv'), index=False)
# 2 thresholds
rows = []
for thr in (0.60, 0.70, 0.75, 0.80, 0.90):
    for o in OUT:
        d = H[o].dropna(subset=['frac_M']); A, B = d[d.frac_M >= thr], d[d.frac_M < thr]
        rec = dict(threshold=thr, outcome=o, k_M=len(A), k_mixed=len(B))
        if len(A) >= 3 and len(B) >= 3:
            a, b = ivw(A), ivw(B); ma, mb = wm(A), wm(B)
            z = (a[0] - b[0]) / np.hypot(a[1], b[1]); zw = (ma[0] - mb[0]) / np.hypot(ma[1], mb[1])
            rec.update(ivw_M=a[0], ivw_M_se=a[1], ivw_mixed=b[0], ivw_mixed_se=b[1], p_diff_ivw=2 * norm.sf(abs(z)),
                       wm_M=ma[0], wm_M_se=ma[1], wm_mixed=mb[0], wm_mixed_se=mb[1], p_diff_wm=2 * norm.sf(abs(zw)))
        rows.append(rec)
pd.DataFrame(rows).to_csv(os.path.join(ext, 'thresholds.csv'), index=False)
# 3 continuous analysis: Wald ratio regressed on the proportion carried by M
def metareg(d):
    r, se = (d.by / d.bx).values, (d.sy / d.bx.abs()).values; x = d.frac_M.clip(0, 1).values
    X = np.column_stack([np.ones(len(r)), x]); W = np.diag(1 / se ** 2)
    C = np.linalg.inv(X.T @ W @ X); b = C @ X.T @ W @ r; res = r - X @ b
    phi = max(1.0, float(res @ W @ res) / (len(r) - 2)); s = np.sqrt(np.diag(C) * phi)
    return dict(k=len(r), intercept=b[0], slope=b[1], slope_se=s[1], p=2 * norm.sf(abs(b[1] / s[1])))
rows = []
for o in OUT:
    d = H[o].dropna(subset=['frac_M'])
    for lab, x in (('all', d), ('no_PHACTR1', d[d.SNP != 'rs9349379']), ('no_PHACTR1_LRP1', d[~d.SNP.isin(['rs9349379', 'rs11172113'])])):
        rows.append(dict(outcome=o, set=lab, **metareg(x)))
pd.DataFrame(rows).to_csv(os.path.join(ext, 'continuous.csv'), index=False)
# 4 uncertainty in the grouping: redraw component effects from their sampling distributions, regroup, re-estimate
rng = np.random.default_rng(13); rows = []
for o in OUT:
    d = H[o].dropna(subset=['frac_M']).copy(); diffs, nM, ps = [], [], []
    for _ in range(5000):
        s = rng.normal(d.est_S, d.SE_S); m = rng.normal(d.est_M, d.SE_M); f = (lM * m) / (lS * s + lM * m)
        A, B = d[f >= 0.75], d[f < 0.75]
        if len(A) < 3 or len(B) < 3:
            continue
        a, b = ivw(A), ivw(B); dd = a[0] - b[0]; diffs.append(dd); nM.append(len(A)); ps.append(2 * norm.sf(abs(dd / np.hypot(a[1], b[1]))))
    diffs = np.array(diffs); ps = np.array(ps)
    rows.append(dict(outcome=o, draws=len(diffs), median_nM=float(np.median(nM)), nM_lo=float(np.percentile(nM, 2.5)), nM_hi=float(np.percentile(nM, 97.5)),
                     diff_median=float(np.median(diffs)), diff_lo=float(np.percentile(diffs, 2.5)), diff_hi=float(np.percentile(diffs, 97.5)),
                     prop_neg=float((diffs < 0).mean()), prop_p05=float((ps < 0.05).mean())))
pd.DataFrame(rows).to_csv(os.path.join(ext, 'grouping_uncertainty.csv'), index=False)
# 5 Steiger: variance explained approximated by Z^2/(Z^2+Neff-2) on the effective-sample-size scale
rows = []; sm = []
for o in OUT:
    d = H[o].copy(); zx, zy = d.bx / d.sx, d.by / d.sy
    r2x, r2y = zx ** 2 / (zx ** 2 + d.Nx - 2), zy ** 2 / (zy ** 2 + d.Ny - 2)
    z = (np.arctanh(np.sqrt(r2x)) - np.arctanh(np.sqrt(r2y))) / np.sqrt(1 / (d.Nx - 3) + 1 / (d.Ny - 3))
    d['r2x'], d['r2y'], d['steiger_z'], d['correct'] = r2x, r2y, z, r2x > r2y
    d['outcome'] = o; rows.append(d[['outcome', 'SNP', 'r2x', 'r2y', 'steiger_z', 'correct', 'group']])
    for g, fn in GROUP.items():
        x = fn(d); y = x[x.correct]
        if len(y) >= 3:
            a, m = ivw(y), wm(y)
            sm.append(dict(outcome=o, group=g, k_all=len(x), k_pass=len(y), failed=';'.join(x[~x.correct].SNP), ivw=a[0], ivw_se=a[1], Q=a[2], wm=m[0], wm_se=m[1]))
pd.concat(rows).to_csv(os.path.join(ext, 'steiger_instruments.csv'), index=False)
pd.DataFrame(sm).to_csv(os.path.join(ext, 'steiger_mr.csv'), index=False)
# 6 power from the estimated precision of the IVW estimate, fixed-effect and allowing for the observed heterogeneity
rows = []
for o in OUT:
    for g, fn in GROUP.items():
        a = ivw(fn(H[o])); se_re, se_fe = a[1], a[5]; rec = dict(outcome=o, group=g, k=a[4], se_re=se_re, se_fe=se_fe)
        for lab, s in (('fe', se_fe), ('re', se_re)):
            for orr in (1.10, 1.15):
                rec[f'power_{lab}_{int(round(orr*100))}'] = float(norm.cdf(np.log(orr) / s - 1.959964))
            rec[f'mdor_{lab}'] = float(np.exp((1.959964 + 0.841621) * s))
        rows.append(rec)
pd.DataFrame(rows).to_csv(os.path.join(ext, 'power.csv'), index=False)
