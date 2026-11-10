"""Standardise raw summary statistics to one layout: SNP A1 A2 BETA SE P N MAF INFO
(A1 = effect allele; BETA = log-odds for binary traits, per-SD for continuous traits; N = effective N).
Variants without an rsID are mapped by GRCh37 position through the 1000 Genomes EUR reference.
Usage: python -I p02_standardize.py <root> <name> [<name> ...]"""
import sys, os, csv, gzip, math, re, subprocess
import numpy as np
from scipy.stats import norm

root = sys.argv[1]
names = sys.argv[2:]
RAW = os.path.join(root, 'data', 'raw')
STD = os.path.join(root, 'data', 'std'); os.makedirs(STD, exist_ok=True)
REF = os.path.join(root, 'data', 'ref', 'ref_1000G_EUR_maf005.txt.gz')
man = {r['name']: r for r in csv.DictReader(open(os.path.join(root, 'src', 'd00_manifest.csv')))}

_pos = None
def posmap():
    global _pos
    if _pos is None:
        _pos = {}
        with gzip.open(REF, 'rt') as f:
            next(f)
            for l in f:
                s, c, bp, maf, a1, a2 = l.split()
                _pos[(c, bp)] = s
    return _pos

ALIASES = dict(
    snp=['rsid', 'rs_id', 'hm_rsid', 'snp', 'markername', 'rsids', 'id', 'variant_id', 'rsnumber', 'marker'],
    a1=['hm_effect_allele', 'effect_allele', 'ea', 'alt', 'a1', 'allele1', 'tested_allele', 'effect_allele_'],
    a2=['hm_other_allele', 'other_allele', 'nea', 'ref', 'a2', 'allele2', 'non_effect_allele', 'other_allele_'],
    beta=['hm_beta', 'beta', 'beta_eur', 'effect', 'effect_size', 'b', 'logor'],
    orr=['hm_odds_ratio', 'odds_ratio', 'or'],
    se=['standard_error', 'se', 'sebeta', 'se_eur', 'stderr', 'se_beta'],
    lo=['hm_ci_lower', 'ci_lower'], hi=['hm_ci_upper', 'ci_upper'],
    p=['p_value', 'pval', 'p', 'pvalue', 'p-value', 'p_bolt_lmm'],
    nlp=['neglog10_pval_eur', 'mlogp', 'neg_log_10_p_value'],
    frq=['hm_effect_allele_frequency', 'effect_allele_frequency', 'af_alt', 'af_eur', 'eaf', 'freq1', 'frq', 'pooled_alt_af', 'freq_tested_allele_in_hrs', 'freq', 'fcon', 'af_controls_eur'],
    info=['imputation_quality', 'info', 'impinfo', 'r2', 'rsq', 'imputation'],
    n=['n_eff', 'neff', 'n', 'n_total', 'samplesize', 'totalsamplesize', 'sample_size'],
    ncas=['num_cases', 'n_cases', 'ncas', 'ncase'], ncon=['num_controls', 'n_controls', 'ncon', 'ncontrol'],
    chr=['chromosome', 'chr', '#chrom', 'chrom', 'hm_chrom'], bp=['base_pair_location', 'pos', 'bp', 'position', 'pos_b37'],
    lowconf=['low_confidence_eur'],
)
NA = {'', 'NA', 'nan', 'NaN', '.', '#NA', 'NULL', 'na', 'Inf', '-Inf'}

def pick(h, key, prefer_hm=True):
    for c in ALIASES[key]:
        if c in h:
            return h[c]
    return None

def opener(path):
    if path.endswith('.gz') or path.endswith('.bgz'):
        return subprocess.Popen(['gzip', '-dc', path], stdout=subprocess.PIPE, text=True, errors='ignore', bufsize=1 << 20).stdout
    return open(path, errors='ignore')

def neff_const(r):
    try:
        ca, co = float(r['cases']), float(r['controls'])
        return 4 / (1 / ca + 1 / co)
    except Exception:
        try:
            return float(r['cases']) if r['controls'] == '' else None   # continuous trait: N in cases column
        except Exception:
            return None

def run(name):
    r = man[name]
    raw = [os.path.join(RAW, f) for f in os.listdir(RAW) if f.split('.')[0] == name and not re.search(r'\.p\d+$|\.tmp$', f)]
    if not raw:
        return name, 'no raw file'
    raw = raw[0]
    out = os.path.join(STD, name + '.tsv.gz')
    f = opener(raw)
    header = None
    for line in f:
        if line.startswith('##'):
            continue
        header = line
        break
    sep = '\t' if '\t' in header else None
    cols = [c.strip().strip('"').lower() for c in header.rstrip('\n').split(sep)]
    h = {c: i for i, c in enumerate(cols)}
    ix = {k: pick(h, k) for k in ALIASES}
    # FinnGen and Pan-UKB report ref/alt; effect allele is alt
    if 'alt' in h and 'ref' in h and 'effect_allele' not in h and 'hm_effect_allele' not in h:
        ix['a1'], ix['a2'] = h['alt'], h['ref']
    # Nielsen 2018 AF: effects are for A2
    if 'effect_a2' in h:
        ix['beta'], ix['a1'], ix['a2'] = h['effect_a2'], h['a2'], h['a1']
        ix['frq'], ix['snp'] = h['freq_a2'], h['rs_dbsnp147']
    nconst = neff_const(r)
    ntot = (float(r['cases']) + float(r['controls'] or 0)) if nconst is not None else None
    pm = None
    if ix['snp'] is None or name.startswith(('BACK', 'SBP', 'DBP', 'EDU', 'SMK', 'UKB_')):
        pm = posmap()
    n_in = n_out = 0
    seen = set()
    with gzip.open(out + '.tmp', 'wt') as o:
        o.write('SNP\tA1\tA2\tBETA\tSE\tP\tN\tMAF\tINFO\n')
        for line in f:
            n_in += 1
            p_ = line.rstrip('\n').split(sep)
            if len(p_) < len(cols):
                continue
            g = lambda k: (p_[ix[k]].strip('"') if ix[k] is not None else '')
            if ix['lowconf'] is not None and g('lowconf') == 'true':
                continue
            snp = g('snp') if ix['snp'] is not None else ''
            if pm is not None and not snp.startswith('rs'):
                snp = pm.get((g('chr').replace('chr', ''), g('bp')), '')
            if ',' in snp:
                snp = snp.split(',')[0]
            m = re.match(r'(rs\d+)', snp)
            if not m:
                continue
            snp = m.group(1)
            if snp in seen:
                continue
            a1, a2 = g('a1').upper(), g('a2').upper()
            if len(a1) != 1 or len(a2) != 1:
                continue
            try:
                b = g('beta')
                if b in NA:
                    o_ = g('orr')
                    if o_ in NA:
                        continue
                    b = math.log(float(o_))
                else:
                    b = float(b)
                se = g('se')
                if se in NA or float(se) <= 0:
                    lo_, hi_ = g('lo'), g('hi')
                    if lo_ not in NA and hi_ not in NA and float(lo_) > 0:
                        se = (math.log(float(hi_)) - math.log(float(lo_))) / (2 * 1.959964)
                    else:
                        # no SE or CI (e.g. GCST90000016): back-calculate from the P value
                        pv0 = g('p')
                        if pv0 in NA or b == 0 or not 0 < float(pv0) < 1:
                            continue
                        se = abs(b) / norm.isf(float(pv0) / 2)
                else:
                    se = float(se)
                if not se > 0:
                    continue
                pv = g('p')
                if pv in NA:
                    nl = g('nlp')
                    pv = 10 ** (-float(nl)) if nl not in NA else 2 * norm.sf(abs(b / se))
                else:
                    pv = float(pv)
                if ix['n'] is not None and g('n') not in NA and name.startswith('MDD'):
                    n = float(g('n'))           # PGC: per-SNP NEFF = sum of cohort effective N (preferred over total-count formula)
                elif ix['ncas'] is not None and g('ncas') not in NA and g('ncon') not in NA:
                    n = 4 / (1 / float(g('ncas')) + 1 / float(g('ncon')))
                elif nconst is not None and ix['n'] is not None and g('n') not in NA and not name.startswith(('BACK', 'SBP', 'DBP', 'EDU', 'SMK', 'UKB_')):
                    n = nconst * float(g('n')) / ntot      # per-SNP N rescaled to effective N
                elif nconst is not None:
                    n = nconst
                elif ix['n'] is not None and g('n') not in NA:
                    n = float(g('n'))
                else:
                    n = float('nan')
                fq = g('frq')
                maf = '' if fq in NA else f"{min(float(fq), 1 - float(fq)):.5f}"
                inf = g('info') if ix['info'] is not None and g('info') not in NA else ''
            except ValueError:
                continue
            seen.add(snp)
            o.write(f"{snp}\t{a1}\t{a2}\t{b:.6g}\t{se:.6g}\t{pv:.6g}\t{n:.0f}\t{maf}\t{inf}\n")
            n_out += 1
    os.rename(out + '.tmp', out)
    return name, f'rows {n_in} kept {n_out} cols {cols[:14]}'

for nm in names:
    print(*run(nm), flush=True)
