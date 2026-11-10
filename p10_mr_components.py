"""Decompose each migraine MR instrument into its shared (S) and migraine-specific (M) parts and regroup instruments.
Usage: python -I p10_mr_components.py <root> <stage,stage,...|full> [threshold=0.75]
For instrument j, the model-implied migraine effect is lS*bS_j + lM*bM_j; frac_M = lM*bM_j / (lS*bS_j + lM*bM_j)
(invariant to allele coding). Instruments with frac_M >= threshold are M-dominant; the others form the mixed group.
Writes results/mr/MIG_instruments_decomp.csv, inst_M.tsv, inst_mixed.tsv and results/mr/instrument_summary.csv."""
import sys, os
import numpy as np, pandas as pd
root, stages = sys.argv[1], sys.argv[2].split(',')
thr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.75
g = os.path.join(root, 'results', 'gwas', 'primary')
base = pd.read_csv(os.path.join(root, 'results', 'models', 'base_primary.csv')).iloc[0]
lS, lM = base.loadS, base.loadM
inst = pd.read_csv(os.path.join(root, 'results', 'mr', 'MIG_instruments.tsv'), sep='\t')
def comp(c):
    d = pd.concat([pd.read_csv(os.path.join(g, f"{c}_SNP{'' if s == 'full' else '_' + s}.tsv.gz"), sep='\t',
                               usecols=['SNP', 'A1', 'est', 'SE', 'Pval_Estimate']) for s in stages])
    return d.drop_duplicates('SNP').rename(columns={'A1': f'A1_{c}', 'est': f'est_{c}', 'SE': f'SE_{c}', 'Pval_Estimate': f'P_{c}'})
d = inst.merge(comp('S'), on='SNP').merge(comp('M'), on='SNP')
assert (d.A1_S == d.A1_M).all()
d['contrib_S'], d['contrib_M'] = lS * d.est_S, lM * d.est_M
tot = d.contrib_S + d.contrib_M
d['frac_M'] = d.contrib_M / tot
sgn = np.sign(tot)                      # orient each instrument to its migraine-risk allele
w_share = float((d.contrib_M * sgn).sum() / (tot * sgn).sum())
d['group'] = np.where(d.frac_M >= thr, 'M-dominant', 'mixed')
d.to_csv(os.path.join(root, 'results', 'mr', 'MIG_instruments_decomp.csv'), index=False)
d.loc[d.group == 'M-dominant', ['SNP']].to_csv(os.path.join(root, 'results', 'mr', 'inst_M.tsv'), sep='\t', index=False)
d.loc[d.group == 'mixed', ['SNP']].to_csv(os.path.join(root, 'results', 'mr', 'inst_mixed.tsv'), sep='\t', index=False)
s = dict(n_instruments=len(inst), n_decomposed=len(d), threshold=thr, weighted_share_M=round(w_share, 3),
         median_frac_M=round(float(d.frac_M.median()), 3), n_M_dominant=int((d.group == 'M-dominant').sum()),
         n_mixed=int((d.group == 'mixed').sum()), n_frac_M_gt_half=int((d.frac_M > 0.5).sum()),
         missing=';'.join(sorted(set(inst.SNP) - set(d.SNP))))
pd.DataFrame([s]).to_csv(os.path.join(root, 'results', 'mr', 'instrument_summary.csv'), index=False)
print(s)
print(d.sort_values('frac_M')[['SNP', 'CHR', 'BP', 'P', 'est_S', 'P_S', 'est_M', 'P_M', 'frac_M', 'group']].to_string(index=False))
