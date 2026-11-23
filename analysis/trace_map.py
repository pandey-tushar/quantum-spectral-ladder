import mpmath as mp, random, itertools
mp.mp.dps = 40
R = lambda v: (v[0], v[2], v[0]*v[2]-v[1])
L = lambda v: (v[2], v[1], v[2]*v[1]-v[0])
def word(w):
    def f(v):
        for ch in reversed(w): v = (L if ch == 'L' else R)(v)
        return v
    return f
kap = lambda v: v[0]**2+v[1]**2+v[2]**2-v[0]*v[1]*v[2]-2
def fixed(w, theta):
    f = word(w); k0 = -2*mp.cos(theta/2)   # trace of commutator
    sols = []
    random.seed(1)
    for _ in range(60):
        g = [mp.mpc(random.uniform(-3,3), random.uniform(-3,3)) for _ in range(3)]
        try:
            s = mp.findroot(lambda a,b,c: [f((a,b,c))[0]-a, f((a,b,c))[1]-b, kap((a,b,c))-k0], g)
        except Exception: continue
        s = tuple(s)
        if max(abs(x) for x in f(s)) > 1e6: continue
        if abs(f(s)[2]-s[2]) > 1e-20: continue
        if all(max(abs(a-b) for a,b in zip(s,t)) > 1e-15 for t in sols): sols.append(s)
    return sols
def mu(w, s):
    f = word(w); J = mp.matrix(3,3)
    for j in range(3):
        J[:, j] = mp.matrix(mp.diff(lambda x: f(tuple(x if i==j else s[i] for i in range(3)))[k], s[j]) for k in range(3)) if False else 0
    h = mp.mpf('1e-15')
    for j in range(3):
        sp = list(s); sm = list(s); sp[j] += h; sm[j] -= h
        fp, fm = f(tuple(sp)), f(tuple(sm))
        for k in range(3): J[k, j] = (fp[k]-fm[k])/(2*h)
    ev = mp.eig(J)[0]
    return sorted([abs(e) for e in ev])
for w, gref in [('LR', {0:1.56678, 1:1.60091, 2:1.68172, 3:1.77227}), ('LLR', {0:2.22902, 1:2.24692, 2:2.30850, 3:2.40925})]:
    for th, gv in gref.items():
        best = None
        for s in fixed(w, mp.mpf(th)):
            if max(abs(mp.im(x)) for x in s) < 1e-10: continue  # want non-real (geometric) point
            m = mu(w, s); val = float(mp.log(m[-1]))
            best = val if best is None else max(best, val)
        print(w, th, 'log|mu_max| =', best, ' g =', gv)
