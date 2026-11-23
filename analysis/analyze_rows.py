import math, warnings; warnings.filterwarnings("ignore")
import mpmath as mp, snappy
src=open(__import__('os').path.join(__import__('os').path.dirname(__file__),'trace_map.py')).read().split("for w, gref")[0]
exec(src.replace("random.seed(1)","pass"))
def vol(w,th):
    M=snappy.Manifold('b++'+w)
    if th>1e-9: M.dehn_fill((2*math.pi/th,0))
    return float(M.volume())
def logmu(w,th):
    best=None
    import random
    sols=[]
    for seed in range(1,40):
        random.seed(seed)
        sols=fixed(w, mp.mpf(th))
        if any(max(abs(mp.im(x)) for x in s)>1e-10 for s in sols): break
    for s in sols:
        if max(abs(mp.im(x)) for x in s) < 1e-10: continue
        v=float(mp.log(mu(w,s)[-1])); best=v if best is None else max(best,v)
    return best
out=open('analysis_out.txt','w')
import sys
for f in sys.argv[1:]:
    for l in open(f):
        p=l.split(); r,c,n,th,w=int(p[0]),int(p[1]),int(p[2]),float(p[3]),p[4]
        LL=[float(x) for x in p[6:]]
        S=4*math.pi/r*LL[0]; V=vol(w,th); lm=logmu(w,th)
        if lm is None: print("FAIL",w,th); lm=float("nan")
        g=[LL[i]-LL[i+1] for i in range(5)]
        s=f"{w} {r} {c} {th:.4f} | {S:.6f} {V:.6f} {S-V:.2e} | {' '.join(f'{x:.5f}' for x in g)} | {lm:.5f}"
        print(s); out.write(s+'\n')
