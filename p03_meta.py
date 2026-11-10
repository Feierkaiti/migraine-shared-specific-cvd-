"""Fixed-effect inverse-variance meta-analysis of standardised files (same SNP layout as p02).
Usage: python -I p03_meta.py <root> <out_name> <in_name> <in_name> [...]
Effect alleles are aligned to the first file; strand-ambiguous SNPs are dropped when alleles must be flipped by strand."""
import sys, os, gzip, math
import numpy as np
from scipy.stats import norm, chi2

root, outname, ins = sys.argv[1], sys.argv[2], sys.argv[3:]
STD = os.path.join(root, 'data', 'std')
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
acc = {}
for k, nm in enumerate(ins):
    with gzip.open(os.path.join(STD, nm + '.tsv.gz'), 'rt') as f:
        next(f)
        for l in f:
            s, a1, a2, b, se, p, n, maf, inf = (l.rstrip('\n').split('\t') + [''])[:9]
            b, se, n = float(b), float(se), float(n)
            rec = acc.get(s)
            if rec is None:
                rec = acc[s] = [a1, a2, [], [], [], [], []]
            r1, r2 = rec[0], rec[1]
            if (a1, a2) == (r1, r2):
                pass
            elif (a1, a2) == (r2, r1):
                b = -b
            elif {a1, a2} != {COMP.get(r1), COMP.get(r2)} or {r1, r2} == {COMP.get(r1), COMP.get(r2)}:
                continue
            elif (COMP.get(a1), COMP.get(a2)) == (r2, r1):
                b = -b
            rec[2].append(b); rec[3].append(se); rec[4].append(n)
            rec[5].append(float(maf) if maf else np.nan); rec[6].append(float(inf) if inf else np.nan)
nk = len(ins)
n_out = n_all = 0
with gzip.open(os.path.join(STD, outname + '.tsv.gz'), 'wt') as o:
    o.write('SNP\tA1\tA2\tBETA\tSE\tP\tN\tMAF\tINFO\tNSTUDY\tQ_P\n')
    for s, rec in acc.items():
        if rec is None or len(rec[2]) == 0:
            continue
        b = np.array(rec[2]); se = np.array(rec[3]); w = 1 / se ** 2
        bm = np.sum(w * b) / np.sum(w); sem = math.sqrt(1 / np.sum(w))
        q = float(np.sum(w * (b - bm) ** 2)); qp = chi2.sf(q, len(b) - 1) if len(b) > 1 else 1.0
        p = 2 * norm.sf(abs(bm / sem))
        maf = np.nanmean(rec[5]) if np.isfinite(rec[5]).any() else np.nan
        inf = np.nanmin(rec[6]) if np.isfinite(rec[6]).any() else np.nan
        n = float(np.sum(rec[4]))
        o.write(f"{s}\t{rec[0]}\t{rec[1]}\t{bm:.6g}\t{sem:.6g}\t{p:.6g}\t{n:.0f}\t{'' if not np.isfinite(maf) else f'{maf:.5f}'}\t{'' if not np.isfinite(inf) else f'{inf:.3f}'}\t{len(b)}\t{qp:.4g}\n")
        n_out += 1
        n_all += len(b) == nk
print(outname, 'SNPs', n_out, 'in all inputs', n_all)
