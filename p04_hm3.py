"""Write HapMap3-only, uncompressed copies of standardised files for GenomicSEM munge.
An INFO column that is empty for every SNP is dropped. Usage: python -I p04_hm3.py <root> <name> [...]"""
import sys, os, gzip, subprocess
root = sys.argv[1]
hm = set(l.split()[0] for i, l in enumerate(open(os.path.join(root, 'data', 'ref', 'w_hm3.snplist'))) if i)
tmp = os.path.join(root, 'data', 'tmp'); os.makedirs(tmp, exist_ok=True)
refmaf = {}
with gzip.open(os.path.join(root, 'data', 'ref', 'ref_1000G_EUR_maf005.txt.gz'), 'rt') as fr:
    next(fr)
    for l in fr:
        x = l.split()
        if x[0] in hm:
            refmaf[x[0]] = x[3]
for nm in sys.argv[2:]:
    src = os.path.join(root, 'data', 'std', nm + '.tsv.gz')
    p = subprocess.Popen(['gzip', '-dc', src], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    head = p.stdout.readline().rstrip('\n').split('\t')[:9]
    rows = []; has_info = False; filled = 0
    for l in p.stdout:
        f = l.rstrip('\n').split('\t')
        if f[0] in hm:
            if f[7] == '' and f[0] in refmaf:      # no frequency in the file: use 1000 Genomes EUR
                f[7] = refmaf[f[0]]; filled += 1
            rows.append(f[:9]); has_info |= f[8] != ''
    keep = 9 if has_info else 8
    with open(os.path.join(tmp, nm + '.hm3.tsv'), 'w') as o:
        o.write('\t'.join(head[:keep]) + '\n')
        for f in rows:
            o.write('\t'.join((x if x != '' else 'NA') for x in f[:keep]) + '\n')
    print(nm, len(rows), 'INFO kept' if has_info else 'INFO dropped', 'MAF filled', filled, flush=True)
