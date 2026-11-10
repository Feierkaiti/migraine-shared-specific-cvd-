"""Turn userGWAS output for S and M into standardised summary statistics.
Usage: python -I c01_components.py <root> <stage,stage,...|full>
Reads results/gwas/primary/{S,M}_SNP[_stage].tsv.gz, keeps SNPs that converged without warnings, computes the
effective sample size of each component as in GenomicSEM (mean of 1/(2*MAF*(1-MAF)*SE^2) over SNPs with MAF 0.10-0.40),
writes data/std/COMP_S.tsv.gz and COMP_M.tsv.gz (SNP A1 A2 BETA SE P N MAF INFO), the MAGMA input files
data/magma_in/COMP_{S,M}.pval and results/components_summary.csv."""
import sys, os
import numpy as np, pandas as pd
from scipy.stats import chi2
root, stages = sys.argv[1], sys.argv[2].split(',')
gdir = os.path.join(root, 'results', 'gwas', 'primary')
os.makedirs(os.path.join(root, 'data', 'magma_in'), exist_ok=True)
hm3 = set(l.split()[0] for i, l in enumerate(open(os.path.join(root, 'data', 'ref', 'w_hm3.snplist'))) if i)
summ = []
for comp in ('S', 'M'):
    parts = []
    for st in stages:
        f = os.path.join(gdir, f"{comp}_SNP{'' if st == 'full' else '_' + st}.tsv.gz")
        parts.append(pd.read_csv(f, sep='\t', usecols=['SNP', 'CHR', 'BP', 'MAF', 'A1', 'A2', 'est', 'SE', 'Pval_Estimate',
                                                       'error', 'warning'], dtype={'error': str, 'warning': str}))
    d = pd.concat(parts, ignore_index=True)
    n0 = len(d)
    ok = d.error.astype(str).isin(['0', '0.0']) & d.warning.astype(str).isin(['0', '0.0'])
    ok &= np.isfinite(d.est) & np.isfinite(d.SE) & (d.SE > 0)
    d = d[ok].drop_duplicates('SNP')
    m = d.MAF.between(0.10, 0.40)
    neff = float(np.mean(1 / ((2 * d.MAF[m] * (1 - d.MAF[m])) * d.SE[m] ** 2)))
    h = d.SNP.isin(hm3)   # HapMap3 variants were analysed regardless of P, so their median chi2 is unbiased
    z2 = (d.est[h] / d.SE[h]) ** 2
    lam = float(np.median(z2) / chi2.ppf(0.5, 1))
    out = pd.DataFrame({'SNP': d.SNP, 'A1': d.A1, 'A2': d.A2, 'BETA': d.est, 'SE': d.SE, 'P': d.Pval_Estimate,
                        'N': round(neff), 'MAF': d.MAF, 'INFO': ''})
    out.to_csv(os.path.join(root, 'data', 'std', f'COMP_{comp}.tsv.gz'), sep='\t', index=False, float_format='%.6g')
    out[['SNP', 'P', 'N']].to_csv(os.path.join(root, 'data', 'magma_in', f'COMP_{comp}.pval'), sep='\t', index=False)
    d[['SNP', 'CHR', 'BP', 'A1', 'A2', 'est', 'SE', 'Pval_Estimate']].rename(columns={'est': 'BETA', 'Pval_Estimate': 'P'}).to_csv(
        os.path.join(gdir, f'COMP_{comp}_pos.tsv.gz'), sep='\t', index=False, float_format='%.6g')
    summ.append(dict(component=comp, stages='+'.join(stages), snps_in=n0, snps_ok=len(d), neff=round(neff),
                     lambda_gc=round(lam, 3), n_hm3=int(h.sum()), n_gws=int((d.Pval_Estimate < 5e-8).sum())))
    print(summ[-1], flush=True)
pd.DataFrame(summ).to_csv(os.path.join(root, 'results', 'components_summary.csv'), index=False)
