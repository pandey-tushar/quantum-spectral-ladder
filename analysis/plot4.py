import re, numpy as np, matplotlib, matplotlib.ticker; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rows=[]
for l in list(open('../data/r201-1001_summary.txt'))+list(open('analysis_out.txt')):
    m=re.match(r'(LLR|LR) (\d+) (\d+) ([\d.]+) \| ([-\d.e+]+) ([-\d.e+]+) ([-\d.e+]+) \| (.*) \| ([\d.]+)',l.strip())
    if m:
        w,r,c,th,S,V,dS,gaps,mu=m.groups()
        rows.append(dict(w=w,r=int(r),th=float(th),dS=float(dS),g=[float(x) for x in gaps.split()],mu=float(mu)))
names={'LR':'LR','LLR':'L²R'}
qty=[('top eigenvalue: |S − Vol(θ)|',lambda x:abs(x['dS'])),('first gap: |log|λ₁/λ₂| − log|μ(θ)||',lambda x:abs(x['g'][0]-x['mu']))]
fig,ax=plt.subplots(2,2,figsize=(12,8.5))
for q,(title,f) in enumerate(qty):
    labs=[];vals=[];cols=[]
    for w,mk in [('LR','o'),('LLR','s')]:
        for t in range(4):
            R=sorted([x for x in rows if x['w']==w and round(x['th'])==t],key=lambda x:x['r'])
            r=np.array([x['r'] for x in R],float); e=np.array([f(x) for x in R])
            kw=dict(color=f'C{t}',ms=4,lw=0.9,mfc='white' if w=='LLR' else None)
            # local exponent between consecutive r
            p=np.diff(np.log(e))/np.diff(np.log(r)); rm=np.sqrt(r[1:]*r[:-1])
            ax[q,0].semilogx(rm,p,mk+'-',label=f'{names[w]}, θ≈{t}',**kw)
            ax[q,1].loglog(r,e,mk+'-',label=f'{names[w]}, θ≈{t}',**kw)
            labs.append(f'{names[w]} θ≈{t}'); vals.append(p[-1])
    ax[q,0].axhline(-1,color='k',ls='--',lw=1.2,label='exponent −1')
    ax[q,0].set_ylim(-1.1,-0.9); ax[q,0].set_xlabel('r'); ax[q,0].set_ylabel('local exponent  d log|err| / d log r')
    ax[q,0].set_title(f'({"ac"[q]}) {title}: approach to 1/r\n(at r=2001: exponent {min(vals):.3f} to {max(vals):.3f})',loc='left',fontsize=9)
    rr=np.array([200,2001.]); ax[q,1].loglog(rr,(0.045 if q==0 else 0.0013)*200/rr,'k--',lw=1.2,label='∝ 1/r')
    ax[q,1].set_xlabel('r'); ax[q,1].set_ylabel('|error|')
    ax[q,1].set_title(f'({"bd"[q]}) {title}: decay with r',loc='left',fontsize=9)
    for k in range(2):
        ax[q,k].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax[q,k].set_xticks([200,300,500,1000,2000]); ax[q,k].set_xticklabels(['200','300','500','1000','2000'])
    ax[q,0].legend(frameon=False,fontsize=7,ncol=2); ax[q,1].legend(frameon=False,fontsize=7,ncol=2)
plt.tight_layout(); plt.savefig('../figures/four_panel.png',dpi=130)
