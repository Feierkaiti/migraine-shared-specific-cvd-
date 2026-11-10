"""LD clumping against 1000 Genomes EUR (MAGMA g1000_eur PLINK files), reading genotypes in Python.
Usage: python -I p07_clump.py <root> <std_or_tsv_file> <out.tsv> [p1=5e-8] [r2=0.001] [kb=10000] [pcol=P]
The input needs SNP and P columns (tab-separated, optionally gzipped)."""
import sys, os, gzip
import numpy as np
root, src, out = sys.argv[1], sys.argv[2], sys.argv[3]
p1 = float(sys.argv[4]) if len(sys.argv) > 4 else 5e-8
r2max = float(sys.argv[5]) if len(sys.argv) > 5 else 0.001
kb = float(sys.argv[6]) if len(sys.argv) > 6 else 10000
pcol = sys.argv[7] if len(sys.argv) > 7 else 'P'
prefix = os.path.join(root, 'data', 'ref', 'g1000', 'g1000_eur')
op = gzip.open if src.endswith('.gz') else open
cand = {}
with op(src, 'rt') as f:
    h = f.readline().rstrip('\n').split('\t'); i_s, i_p = h.index('SNP'), h.index(pcol)
    for l in f:
        x = l.rstrip('\n').split('\t')
        try:
            p = float(x[i_p])
        except ValueError:
            continue
        if p < p1:
            cand[x[i_s]] = p
n = sum(1 for _ in open(prefix + '.fam')); nb = (n + 3) // 4
pos = {}
with open(prefix + '.bim') as f:
    for k, l in enumerate(f):
        x = l.split()
        if x[1] in cand:
            pos[x[1]] = (int(x[0]), int(x[3]), k)
bed = np.memmap(prefix + '.bed', dtype=np.uint8, mode='r', offset=3)
dec = np.array([[2, np.nan, 1, 0][(b >> (2 * j)) & 3] for b in range(256) for j in range(4)]).reshape(256, 4)
def geno(k):
    g = dec[np.asarray(bed[k * nb:(k + 1) * nb])].reshape(-1)[:n]
    m = np.nanmean(g); g = np.where(np.isnan(g), m, g); g = g - g.mean(); s = g.std()
    return g / s if s > 0 else None
snps = sorted([s for s in cand if s in pos], key=lambda s: cand[s])
G = {s: geno(pos[s][2]) for s in snps}
leads, clumped = [], set()
for s in snps:
    if s in clumped or G[s] is None:
        continue
    leads.append(s); c, bp, _ = pos[s]
    for t in snps:
        if t in clumped or t == s or G[t] is None:
            continue
        ct, bt, _ = pos[t]
        if ct == c and abs(bt - bp) <= kb * 1000:
            r = float(np.dot(G[s], G[t]) / n)
            if r * r > r2max:
                clumped.add(t)
    clumped.add(s)
with open(out, 'w') as o:
    o.write('SNP\tCHR\tBP\tP\n')
    for s in leads:
        o.write(f'{s}\t{pos[s][0]}\t{pos[s][1]}\t{cand[s]:.3g}\n')
print('candidates', len(cand), 'in reference', len(snps), 'leads', len(leads), 'not in 1000G', len(cand) - len(snps))
