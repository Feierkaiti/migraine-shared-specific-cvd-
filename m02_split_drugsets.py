"""Split the migraine drug-target set into targets that are also targets of depression drugs and the remainder.
Writes data/ref/annot/genesets_split.txt (MAGMA set-annot format). Usage: python -I m02_split_drugsets.py <root>
Then: bash m01_magma.sh runs the split-set analysis for every GWAS when this file exists."""
import sys, os
root = sys.argv[1]
A = os.path.join(root, 'data', 'ref', 'annot')
sets = {l.split()[0]: set(l.split()[1:]) for l in open(os.path.join(A, 'genesets.txt'))}
m, d = sets['migraine_drug_targets'], sets['depression_drug_targets']
with open(os.path.join(A, 'genesets_split.txt'), 'w') as f:
    f.write('migraine_targets_shared_with_depression ' + ' '.join(sorted(m & d)) + '\n')
    f.write('migraine_targets_not_depression ' + ' '.join(sorted(m - d)) + '\n')
print('migraine targets', len(m), 'depression targets', len(d), 'shared', len(m & d), 'migraine only', len(m - d))
