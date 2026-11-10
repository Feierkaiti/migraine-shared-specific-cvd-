"""Independent loci of the component GWAS and of the migraine meta-analysis, with overlap and nearest genes.
Usage: python -I c02_loci.py <root>
Needs results/loci/{COMP_S,COMP_M,MIG}_clump.tsv from p07_clump.py (P < 5e-8, r2 < 0.1, 250 kb).
Leads within 250 kb are merged into one locus. Writes results/loci/loci_table.csv and results/loci_summary.csv."""
import sys, os
import numpy as np, pandas as pd
root = sys.argv[1]
L = os.path.join(root, 'results', 'loci')
gl = pd.read_csv(os.path.join(root, 'data', 'ref', 'genes', 'NCBI37.3.gene.loc'), sep='\t', header=None,
                 names=['id', 'chr', 'start', 'end', 'strand', 'symbol'], dtype={'chr': str})

def pos_table(name):
    if name == 'MIG':
        f = os.path.join(root, 'data', 'std', 'MIG.tsv.gz')
        ref = pd.read_csv(os.path.join(root, 'data', 'ref', 'ref_1000G_EUR_maf005.txt.gz'), sep='\t', usecols=['SNP', 'CHR', 'BP'])
        return ref
    return pd.read_csv(os.path.join(root, 'results', 'gwas', 'primary', f'{name}_pos.tsv.gz'), sep='\t', usecols=['SNP', 'CHR', 'BP'])

def merge(name, pos):
    c = pd.read_csv(os.path.join(L, f'{name}_clump.tsv'), sep='\t')
    c = c.merge(pos, on='SNP', how='left', suffixes=('', '_r'))
    if 'CHR_r' in c:
        c['CHR'] = c['CHR'].fillna(c['CHR_r']); c['BP'] = c['BP'].fillna(c['BP_r'])
    c = c.sort_values('P')
    keep = []
    for _, r in c.iterrows():
        if any(k['CHR'] == r.CHR and abs(k['BP'] - r.BP) < 250000 for k in keep):
            continue
        keep.append(r)
    return pd.DataFrame(keep)

def nearest(chr_, bp):
    g = gl[gl.chr == str(int(chr_))]
    inside = g[(g.start <= bp) & (g.end >= bp)]
    if len(inside):
        return ','.join(inside.symbol.head(2))
    d = np.minimum(abs(g.start - bp), abs(g.end - bp))
    return g.symbol.iloc[int(np.argmin(d.values))]

ref_pos = pos_table('MIG')
loci = {}
for nm in ('COMP_S', 'COMP_M', 'MIG'):
    if os.path.exists(os.path.join(L, f'{nm}_clump.tsv')):
        loci[nm] = merge(nm, ref_pos if nm == 'MIG' else pos_table(nm))
rows = []
for nm, d in loci.items():
    for _, r in d.iterrows():
        near = {o: any((x.CHR == r.CHR) and abs(x.BP - r.BP) < 500000 for _, x in od.iterrows()) for o, od in loci.items() if o != nm}
        rows.append(dict(set=nm, SNP=r.SNP, CHR=int(r.CHR), BP=int(r.BP), P=r.P, gene=nearest(r.CHR, r.BP),
                         **{f'near_{o}': v for o, v in near.items()}))
t = pd.DataFrame(rows)
for c in [c for c in t.columns if c.startswith('near_')]:
    t[c] = t[c].astype('boolean')   # own-set column is empty
# P values of each lead variant in the two component GWAS (lead variants missing from an input GWAS stay empty)
for nm, col in (('COMP_S', 'P_S'), ('COMP_M', 'P_M')):
    pp = pd.read_csv(os.path.join(root, 'results', 'gwas', 'primary', f'{nm}_pos.tsv.gz'), sep='\t', usecols=['SNP', 'P']).set_index('SNP').P
    t[col] = t.SNP.map(pp)
t.to_csv(os.path.join(L, 'loci_table.csv'), index=False)
s = {'S': len(loci.get('COMP_S', [])), 'M': len(loci.get('COMP_M', [])), 'MIG': len(loci.get('MIG', []))}
if 'COMP_S' in loci and 'COMP_M' in loci and 'MIG' in loci:
    mt = t[t.set == 'MIG']
    nS, nM = mt.near_COMP_S.fillna(False).astype(bool), mt.near_COMP_M.fillna(False).astype(bool)
    s['MIG_near_M'] = int(nM.sum()); s['MIG_near_S'] = int(nS.sum()); s['MIG_near_both'] = int((nS & nM).sum())
    s['M_near_MIG'] = int(t[t.set == 'COMP_M'].near_MIG.fillna(False).astype(bool).sum())
    s['M_new'] = s['M'] - s['M_near_MIG']
    s['S_near_MIG'] = int(t[t.set == 'COMP_S'].near_MIG.fillna(False).astype(bool).sum())
    s['MIG_lead_tested'] = int(mt.P_M.notna().sum())
    s['MIG_lead_M_gws'] = int((mt.P_M < 5e-8).sum()); s['MIG_lead_S_gws'] = int((mt.P_S < 5e-8).sum())
    s['MIG_lead_M_1e5'] = int((mt.P_M < 1e-5).sum()); s['MIG_lead_S_1e5'] = int((mt.P_S < 1e-5).sum())
pd.DataFrame([s]).to_csv(os.path.join(root, 'results', 'loci_summary.csv'), index=False)
print(s)
print(t.to_string(index=False))
