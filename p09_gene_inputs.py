"""Gene-level inputs for MAGMA:
1. GTEx v8 tissue expression covariates (average log2(TPM+1) per tissue, 54 tissues plus their mean),
   keyed by Entrez ID (MAGMA NCBI37.3 gene IDs), as used in FUMA's MAGMA tissue analysis.
2. Gene sets: targets of drugs indicated for migraine and for depression (ChEMBL MeSH indications, phase >= 3),
   and a prespecified CGRP/triptan set.
Usage: python -I p09_gene_inputs.py <root>"""
import sys, os, gzip, json, time, urllib.request, urllib.parse
import numpy as np, pandas as pd
root = sys.argv[1]
A = os.path.join(root, 'data', 'ref', 'annot')
hg = pd.read_csv(os.path.join(A, 'hgnc_complete_set.txt'), sep='\t', dtype=str, usecols=['symbol', 'entrez_id', 'ensembl_gene_id', 'uniprot_ids', 'status'])
hg = hg[hg.status == 'Approved']
loc = pd.read_csv(os.path.join(root, 'data', 'ref', 'genes', 'NCBI37.3.gene.loc'), sep='\t', header=None, dtype=str)
magma_ids = set(loc[0])

# 1. GTEx tissue expression
g = pd.read_csv(os.path.join(A, 'gtex_v8_median_tpm.gct.gz'), sep='\t', skiprows=2)
g['ens'] = g['Name'].str.split('.').str[0]
m = hg.dropna(subset=['entrez_id', 'ensembl_gene_id']).drop_duplicates('ensembl_gene_id').set_index('ensembl_gene_id')['entrez_id']
g['GENE'] = g['ens'].map(m)
g = g.dropna(subset=['GENE'])
g = g[g.GENE.isin(magma_ids)].drop_duplicates('GENE')
tissues = [c for c in g.columns if c not in ('Name', 'Description', 'ens', 'GENE')]
x = np.log2(g[tissues].astype(float) + 1)
x.columns = [c.replace(' - ', '_').replace(' ', '_').replace('(', '').replace(')', '') for c in tissues]
x['Average'] = x.mean(axis=1)
x.insert(0, 'GENE', g['GENE'].values)
x.to_csv(os.path.join(A, 'gtex_v8_ts_avg_log2TPM_54.txt'), sep='\t', index=False)
print('GTEx genes', len(x), 'tissues', len(tissues))

# 2. ChEMBL drug targets
API = 'https://www.ebi.ac.uk/chembl/api/data/'
def get(path, params):
    out, off = [], 0
    while True:
        q = dict(params, limit=1000, offset=off, format='json')
        for k in range(5):
            try:
                d = json.load(urllib.request.urlopen(API + path + '?' + urllib.parse.urlencode(q), timeout=120)); break
            except Exception:
                time.sleep(5)
        key = [k for k in d if k != 'page_meta'][0]
        out += d[key]
        if not d['page_meta'].get('next'):
            return out
        off += 1000
def targets(efo_ids):
    mols = set()
    for e in efo_ids:
        for r in get('drug_indication', {'mesh_id': e, 'max_phase_for_ind__gte': 3}):
            mols.add(r['molecule_chembl_id'])
    genes, rows = set(), []
    mols = sorted(mols)
    for i in range(0, len(mols), 40):
        for r in get('mechanism', {'molecule_chembl_id__in': ','.join(mols[i:i + 40])}):
            t = r.get('target_chembl_id')
            if not t:
                continue
            td = get('target', {'target_chembl_id': t})
            for tt in td:
                if tt.get('organism') != 'Homo sapiens':
                    continue
                for comp in tt.get('target_components', []):
                    for syn in comp.get('target_component_synonyms', []):
                        if syn['syn_type'] == 'GENE_SYMBOL':
                            genes.add(syn['component_synonym']); rows.append((r['molecule_chembl_id'], t, syn['component_synonym']))
    return mols, genes, rows
sym2ent = hg.dropna(subset=['entrez_id']).set_index('symbol')['entrez_id'].to_dict()
sets = {}
# ChEMBL indications are keyed by MeSH (efo_id now holds MONDO terms)
for name, efos in [('migraine_drug_targets', ['D008881', 'D020325', 'D020326']), ('depression_drug_targets', ['D003865', 'D003866'])]:
    mols, genes, rows = targets(efos)
    ent = sorted({sym2ent[s] for s in genes if s in sym2ent and sym2ent[s] in magma_ids})
    sets[name] = ent
    pd.DataFrame(rows, columns=['molecule', 'target', 'gene']).drop_duplicates().to_csv(os.path.join(A, f'chembl_{name}.csv'), index=False)
    print(name, 'drugs', len(mols), 'gene symbols', len(genes), 'mapped', len(ent))
cgrp = ['CALCA', 'CALCB', 'CALCRL', 'RAMP1', 'HTR1B', 'HTR1D', 'HTR1F']
sets['cgrp_triptan'] = sorted({sym2ent[s] for s in cgrp if s in sym2ent})
with open(os.path.join(A, 'genesets.txt'), 'w') as o:
    for k, v in sets.items():
        o.write(k + ' ' + ' '.join(v) + '\n')
print({k: len(v) for k, v in sets.items()})
