import sys,pickle,time
from reference_mprep import *
from multiprocessing import Pool
def job(r):
    C=MCat(r); out={}; t=time.time()
    for c in C.cols:
        ad,Ta,Tb=gens(C,c)
        Tbi=Tb**-1
        br=mp.norm(Ta*Tb*Ta-Tb*Ta*Tb)/mp.norm(Ta*Tb*Ta)
        # store high-prec eig for powers: store eigenvalues as mpc strings
        out[c]=dict(ad=ad,br=float(br),ev={w:[mp.nstr(e,60) for e in mp.eig(P2,left=False,right=False)] for w,P2 in (('LR',Tbi*Ta),('LLR',Tbi*Tbi*Ta))})
    return r,out,time.time()-t
if __name__=='__main__':
    rs=list(range(5,int(sys.argv[1])+1,2))
    rs.sort(key=lambda x:-x)
    res={}
    with Pool() as p:
        for r,o,t in p.imap_unordered(job,rs):
            res[r]=o; print(r,round(t,1),max(v['br'] for v in o.values()),flush=True)
            pickle.dump(res,open('eigs_mp.pkl','wb'))
