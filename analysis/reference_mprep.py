import mpmath as mp, sys
mp.mp.dps=120
class MCat:
    def __init__(self,r,s=2):
        self.r=r; self.s=s
        sn=lambda n: mp.sin(mp.pi*s*n/r)/mp.sin(mp.pi*s/r)
        self.qi=[sn(n) for n in range(0,3*r+5)]
        f=[mp.mpf(1)]
        for k in range(1,len(self.qi)): f.append(f[-1]*self.qi[k])
        self.fac=f
        self.cols=list(range(0,r-2,2))
        self.th={a:mp.expjpi(mp.mpf(2*s)/r*(a//2)*(a//2+1)) for a in self.cols}
        self.cache={}
        self.dl={}
    def d(self,a): return self.qi[a+1]
    def adm(self,a,b,c): return (a+b+c)%2==0 and abs(a-b)<=c<=a+b and a+b+c<=2*(self.r-2)
    def Delta(self,a,b,c):
        k=tuple(sorted((a,b,c)))
        if k not in self.dl:
            f=self.fac
            v=f[(a+b-c)//2]*f[(a-b+c)//2]*f[(-a+b+c)//2]/f[(a+b+c)//2+1]
            self.dl[k]=mp.sqrt(mp.mpc(v))
        return self.dl[k]
    def sixj(self,j1,j2,j3,j4,j5,j6):
        key=(j1,j2,j3,j4,j5,j6)
        if key in self.cache: return self.cache[key]
        T=[(j1,j2,j3),(j1,j5,j6),(j4,j2,j6),(j4,j5,j3)]
        if not all(self.adm(*t) for t in T): return mp.mpf(0)
        f=self.fac
        pref=1
        for t in T: pref*=self.Delta(*t)
        a1=[(j1+j2+j3)//2,(j1+j5+j6)//2,(j4+j2+j6)//2,(j4+j5+j3)//2]
        b1=[(j1+j2+j4+j5)//2,(j2+j3+j5+j6)//2,(j3+j1+j6+j4)//2]
        tot=0
        for z in range(max(a1),min(b1)+1):
            den=1
            for x in a1: den*=f[z-x]
            for y in b1: den*=f[y-z]
            tot+=(-1)**z*f[z+1]/den
        v=pref*tot; self.cache[key]=v; return v
def gens(C,c):
    ad=[a for a in C.cols if C.adm(a,a,c)]; n=len(ad)
    S=mp.zeros(n,n)
    for i,a in enumerate(ad):
        for j,b in enumerate(ad):
            tot=0
            for x in C.cols:
                if not C.adm(a,b,x): continue
                tot+=C.d(x)*(C.th[x]/(C.th[a]*C.th[b]))*C.sixj(a,a,c,b,b,x)*(-1)**((a+b+x)//2)
            S[i,j]=tot*mp.sqrt(mp.mpc(C.d(a)))*mp.sqrt(mp.mpc(C.d(b)))
    Ta=mp.diag([C.th[a] for a in ad])
    Tb=S*Ta*S**-1
    return ad,Ta,Tb
