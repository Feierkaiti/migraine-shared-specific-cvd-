"""Build a 1000 Genomes EUR allele-frequency reference (SNP, CHR, BP, MAF, A1, A2) from the
MAGMA g1000_eur PLINK files, in the format GenomicSEM::sumstats expects.
Usage: python -I p01_build_ref.py <prefix of .bed/.bim/.fam> <out.txt.gz>"""
import sys, gzip
import numpy as np

prefix, out = sys.argv[1], sys.argv[2]
n = sum(1 for _ in open(prefix + '.fam'))
nbytes = (n + 3) // 4
last = n - 4 * (nbytes - 1)          # samples coded in the last byte of each SNP

# per-byte lookup: A1 allele count and non-missing count for the 4 genotypes packed in one byte
code_a1 = {0: 2, 2: 1, 3: 0}         # 00 hom A1, 10 het, 11 hom A2, 01 missing
lut_a1 = np.zeros((256, 5), np.int32)
lut_nm = np.zeros((256, 5), np.int32)
for b in range(256):
    for k in range(1, 5):            # k = number of samples used in this byte
        a = m = 0
        for j in range(k):
            g = (b >> (2 * j)) & 3
            if g != 1:
                a += code_a1[g]; m += 1
        lut_a1[b, k], lut_nm[b, k] = a, m

bed = np.memmap(prefix + '.bed', dtype=np.uint8, mode='r', offset=3)
nsnp = bed.size // nbytes
bed = bed.reshape(nsnp, nbytes)

kept = 0
with open(prefix + '.bim') as bim, gzip.open(out, 'wt') as o:
    o.write('SNP\tCHR\tBP\tMAF\tA1\tA2\n')
    seen = set()
    for start in range(0, nsnp, 200000):
        blk = np.asarray(bed[start:start + 200000])
        a1 = lut_a1[blk[:, :-1], 4].sum(1) + lut_a1[blk[:, -1], last]
        nm = lut_nm[blk[:, :-1], 4].sum(1) + lut_nm[blk[:, -1], last]
        frq = np.where(nm > 0, a1 / (2.0 * np.maximum(nm, 1)), np.nan)
        for i in range(len(blk)):
            f = bim.readline().split()
            chrom, snp, bp, al1, al2 = f[0], f[1], f[3], f[4].upper(), f[5].upper()
            if not snp.startswith('rs') or snp in seen or not np.isfinite(frq[i]):
                continue
            if len(al1) != 1 or len(al2) != 1:
                continue
            maf = min(frq[i], 1 - frq[i])
            if maf < 0.005:
                continue
            seen.add(snp)
            o.write(f'{snp}\t{chrom}\t{bp}\t{maf:.5f}\t{al1}\t{al2}\n')
            kept += 1
print('samples', n, 'variants', nsnp, 'kept', kept)
