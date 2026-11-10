"""File check for migraine inputs: direction and P at three established loci."""
import sys, os, gzip, subprocess
root = sys.argv[1]
LOCI = {'rs11172113': 'LRP1', 'rs9349379': 'PHACTR1', 'rs10166942': 'TRPM8'}
ref = {}
for nm in sys.argv[2:]:
    p = subprocess.Popen(['gzip', '-dc', os.path.join(root, 'data', 'std', nm + '.tsv.gz')], stdout=subprocess.PIPE, text=True)
    next(p.stdout); got = {}
    for l in p.stdout:
        f = l.split('\t')
        if f[0] in LOCI:
            got[f[0]] = (f[1], f[2], float(f[3]), float(f[5]))
            if len(got) == 3: p.kill(); break
    out = []
    for s, g in LOCI.items():
        if s not in got: out.append(f'{g}: missing'); continue
        a1, a2, b, pv = got[s]
        if nm == sys.argv[2]: ref[s] = (a1, b)
        sign = b if (s not in ref or a1 == ref[s][0]) else -b
        same = '' if s not in ref else ('same' if (sign > 0) == (ref[s][1] > 0) else 'OPPOSITE')
        out.append(f'{g}: P={pv:.1e} {same}')
    print(nm, '|', '; '.join(out), flush=True)
