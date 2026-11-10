"""Prepare inputs for GenomicSEM::sumstats: SNPs present in every indicator file and in the
1000 Genomes EUR reference with MAF >= 0.01, written uncompressed. SNPs with a reported INFO below 0.6 are dropped here, and INFO is not passed on.
Usage: python -I p06_prep_gwas.py <root> <out_tag> <name> [<name> ...]"""
import sys, os, gzip, subprocess
root, tag, names = sys.argv[1], sys.argv[2], sys.argv[3:]
out = os.path.join(root, 'data', 'gwas_in', tag); os.makedirs(out, exist_ok=True)
ref = set()
with gzip.open(os.path.join(root, 'data', 'ref', 'ref_1000G_EUR_maf005.txt.gz'), 'rt') as f:
    next(f)
    for l in f:
        x = l.split()
        if float(x[3]) >= 0.01:
            ref.add(x[0])
keep = ref
def read(nm):
    p = subprocess.Popen(['gzip', '-dc', os.path.join(root, 'data', 'std', nm + '.tsv.gz')], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    next(p.stdout)
    for l in p.stdout:
        yield l.rstrip('\n').split('\t')
for nm in names:
    s = set()
    for f in read(nm):
        if f[0] in keep and not (f[8] not in ('', 'nan') and float(f[8]) < 0.6):
            s.add(f[0])
    keep = keep & s
    print(nm, 'running intersection', len(keep), flush=True)
for nm in names:
    with open(os.path.join(out, nm + '.tsv'), 'w') as o:
        o.write('SNP\tA1\tA2\tBETA\tSE\tP\tN\tMAF\n')
        for f in read(nm):
            if f[0] in keep:
                o.write('\t'.join(f[:7] + [f[7] if f[7] else 'NA']) + '\n')
print('SNPs written', len(keep))
