"""Download raw GWAS summary statistics listed in d00_manifest.csv into data/raw/.
EBI files are fetched as byte ranges over a few connections (EBI caps each connection);
other hosts are fetched with a single stream. Usage: d01_fetch.py <root> <ebi|fast> [names]"""
import sys,os,csv,re,subprocess,time
from concurrent.futures import ThreadPoolExecutor
root=sys.argv[1];mode=sys.argv[2];only=set(sys.argv[3].split(',')) if len(sys.argv)>3 else None
RAW=os.path.join(root,'data','raw');os.makedirs(RAW,exist_ok=True)
rows=[r for r in csv.DictReader(open(os.path.join(root,'src','d00_manifest.csv')))]
def dest(r):
    u=r['url'];ext='.gz' if re.search(r'\.(gz|bgz)$',u) or 'ndownloader' in u else ('.tsv' if u.endswith('.tsv') else '.txt')
    if u.endswith('.bgz'): ext='.bgz'
    return os.path.join(RAW,r['name']+ext)
def size(u):
    h=subprocess.run(['curl','-s','-I','-L','-m','60',u],capture_output=True,text=True).stdout
    m=re.findall(r'content-length:\s*(\d+)',h,re.I);return int(m[-1]) if m else None
def ranged(u,d,ncon=10,chunk=8_000_000):
    n=size(u)
    parts=[(i,i*chunk,min(n-1,(i+1)*chunk-1)) for i in range(n//chunk+1) if i*chunk<n]
    def get(p):
        i,a,b=p;f=f"{d}.p{i}"
        for k in range(10):
            if os.path.exists(f) and os.path.getsize(f)==b-a+1: return
            subprocess.run(['curl','-s','-m','1200','-r',f'{a}-{b}','-o',f,u])
            if not (os.path.exists(f) and os.path.getsize(f)==b-a+1): time.sleep(15)
        if not (os.path.exists(f) and os.path.getsize(f)==b-a+1): raise RuntimeError('failed '+f)
    with ThreadPoolExecutor(ncon) as ex: list(ex.map(get,parts))
    with open(d+'.tmp','wb') as o:
        for i,_,_ in parts:
            with open(f"{d}.p{i}",'rb') as q: o.write(q.read())
    for i,_,_ in parts: os.remove(f"{d}.p{i}")
    os.rename(d+'.tmp',d)
    if os.path.getsize(d)!=n: raise RuntimeError('size mismatch '+d)
def plain(u,d):
    subprocess.run(['curl','-s','-L','--retry','5','-o',d+'.tmp',u],check=True)
    n=size(u)
    if n and os.path.getsize(d+'.tmp')!=n: raise RuntimeError('size mismatch '+d)
    os.rename(d+'.tmp',d)
def job(r):
    d=dest(r)
    if os.path.exists(d): return r['name'],'exists'
    t=time.time()
    try:
        if 'ftp.ebi.ac.uk' in r['url']: ranged(r['url'],d)
        else: plain(r['url'],d)
        return r['name'],f'ok {os.path.getsize(d)/1e6:.0f} MB in {time.time()-t:.0f}s'
    except Exception as e: return r['name'],'ERR '+repr(e)
todo=[r for r in rows if (only is None or r['name'] in only) and (('ftp.ebi.ac.uk' in r['url'])==(mode=='ebi'))]
with ThreadPoolExecutor(1 if mode=='ebi' else 4) as ex:
    for name,msg in ex.map(job,todo): print(time.strftime('%H:%M:%S'),name,msg,flush=True)
