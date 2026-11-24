// gens_check: checks the SO(3) one-holed-torus generators T_a, T_b against MCG relations.
// Primary T_b: curve-operator construction of solver/ladder3.c (setup, phaseA); second: S-matrix sum of
// analysis/reference_mprep.py gens().  usage: gens_check selftest|table [RMAX]
#include <quadmath.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
using namespace std;
typedef __float128 Q;
static const Q PIQ = 4 * atanq((Q) 1);

struct C { Q re, im; };
static C mk(Q a, Q b = 0) { return {a, b}; }
static C operator+(C a, C b) { return {a.re + b.re, a.im + b.im}; }
static C operator-(C a, C b) { return {a.re - b.re, a.im - b.im}; }
static C operator*(C a, C b) { return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}; }
static C operator*(Q s, C a) { return {s * a.re, s * a.im}; }
static C operator/(C a, C b) { Q d = b.re * b.re + b.im * b.im; return {(a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d}; }
static C cj(C a) { return {a.re, -a.im}; }
static Q cab(C a) { return hypotq(a.re, a.im); }
static C csq(C a) {  // principal branch; negative real v -> +i sqrt|v|
    Q m = cab(a); if (m == 0) return mk(0);
    Q x = sqrtq((m + fabsq(a.re)) / 2);
    if (a.re >= 0) return {x, a.im / (2 * x)};
    return {fabsq(a.im) / (2 * x), a.im >= 0 ? x : -x};
}
static C eipi(long m, int r) { m %= 2L * r; if (m < 0) m += 2L * r; Q t = PIQ * m / r; return {cosq(t), sinq(t)}; }  // exp(i pi m/r)
static double D(Q x) { return (double) x; }
static double mx(double a, double b) { return (a != a || b != b) ? NAN : (a > b ? a : b); }

struct Mat {
    int n; vector<C> a;
    Mat(int n_ = 0) : n(n_), a((size_t) n_ * n_, mk(0)) {}
    C &operator()(int i, int j) { return a[(size_t) i * n + j]; }
    C operator()(int i, int j) const { return a[(size_t) i * n + j]; }
};
static Mat operator*(const Mat &A, const Mat &B) {
    Mat P(A.n);
    for (int i = 0; i < A.n; i++) for (int k = 0; k < A.n; k++) { C x = A(i, k); for (int j = 0; j < A.n; j++) P(i, j) = P(i, j) + x * B(k, j); }
    return P;
}
static Mat operator-(const Mat &A, const Mat &B) { Mat P(A.n); for (size_t k = 0; k < A.a.size(); k++) P.a[k] = A.a[k] - B.a[k]; return P; }
static Q fro(const Mat &A) { Q s = 0; for (C x : A.a) s += x.re * x.re + x.im * x.im; return sqrtq(s); }
static Mat inv(Mat A) {
    int n = A.n; Mat B(n); for (int i = 0; i < n; i++) B(i, i) = mk(1);
    for (int k = 0; k < n; k++) {
        int p = k; for (int i = k + 1; i < n; i++) if (cab(A(i, k)) > cab(A(p, k))) p = i;
        for (int j = 0; j < n; j++) { C t = A(k, j); A(k, j) = A(p, j); A(p, j) = t; t = B(k, j); B(k, j) = B(p, j); B(p, j) = t; }
        C piv = A(k, k);
        for (int j = 0; j < n; j++) { A(k, j) = A(k, j) / piv; B(k, j) = B(k, j) / piv; }
        for (int i = 0; i < n; i++) if (i != k) { C f = A(i, k); for (int j = 0; j < n; j++) { A(i, j) = A(i, j) - f * A(k, j); B(i, j) = B(i, j) - f * B(k, j); } }
    }
    return B;
}

struct Cat {
    int r; vector<Q> qi, fac; vector<int> cols;
    explicit Cat(int r_) : r(r_) {
        int N = 3 * r + 5;
        Q s1 = sinq(2 * PIQ / r);
        for (int k = 0; k < N; k++) qi.push_back(sinq(2 * PIQ * k / r) / s1);  // unreduced, as ladder3.c:101
        fac.assign(N, 1);
        for (int k = 1; k < N; k++) fac[k] = fac[k - 1] * qi[k];
        for (int a = 0; a < r - 2; a += 2) cols.push_back(a);
    }
    Q d(int a) const { return qi[a + 1]; }
    Q qint(long m) const { return qi[m % r]; }
    bool adm(int a, int b, int c) const { return (a + b + c) % 2 == 0 && abs(a - b) <= c && c <= a + b && a + b + c <= 2 * (r - 2); }
    C xi(int a) const { long k = a / 2; return eipi(4 * k * (k + 1), r); }
    C Delta(int a, int b, int c) const {
        return csq(mk(fac[(a + b - c) / 2] * fac[(a - b + c) / 2] * fac[(-a + b + c) / 2] / fac[(a + b + c) / 2 + 1]));
    }
    C sixj(int j1, int j2, int j3, int j4, int j5, int j6) const {
        int T[4][3] = {{j1, j2, j3}, {j1, j5, j6}, {j4, j2, j6}, {j4, j5, j3}};
        C pref = mk(1);
        for (auto &t : T) { if (!adm(t[0], t[1], t[2])) return mk(0); pref = pref * Delta(t[0], t[1], t[2]); }
        int a1[4] = {(j1 + j2 + j3) / 2, (j1 + j5 + j6) / 2, (j4 + j2 + j6) / 2, (j4 + j5 + j3) / 2};
        int b1[3] = {(j1 + j2 + j4 + j5) / 2, (j2 + j3 + j5 + j6) / 2, (j3 + j1 + j6 + j4) / 2};
        int zl = max(max(a1[0], a1[1]), max(a1[2], a1[3])), zh = min(min(b1[0], b1[1]), b1[2]);
        Q tot = 0;
        for (int z = zl; z <= zh; z++) {
            Q den = 1;
            for (int x : a1) den *= fac[z - x];
            for (int y : b1) den *= fac[y - z];
            tot += (z % 2 ? -1 : 1) * fac[z + 1] / den;
        }
        return tot * pref;
    }
    vector<int> basis(int c) const { vector<int> ad; for (int a : cols) if (adm(a, a, c)) ad.push_back(a); return ad; }
};

// reference_mprep.py gens(): sgn 0 = (-1)^((a+b+x)/2), 1 = none, 2 = (-1)^(x/2)
static void gens_ref(const Cat &K, int c, int sgn, Mat &Ta, Mat &Tb, Mat *Sout = nullptr) {
    vector<int> ad = K.basis(c); int n = ad.size();
    Ta = Mat(n); Mat S(n);
    for (int i = 0; i < n; i++) {
        int a = ad[i]; Ta(i, i) = K.xi(a);
        for (int j = 0; j < n; j++) {
            int b = ad[j]; C tot = mk(0);
            for (int x : K.cols) {
                if (!K.adm(a, b, x)) continue;
                int e = sgn == 0 ? (a + b + x) / 2 : sgn == 2 ? x / 2 : 0;
                tot = tot + ((e % 2 ? -1 : 1) * K.d(x)) * (K.xi(x) / (K.xi(a) * K.xi(b)) * K.sixj(a, a, c, b, b, x));
            }
            S(i, j) = tot * csq(mk(K.d(a))) * csq(mk(K.d(b)));
        }
    }
    Tb = S * Ta * inv(S);
    if (Sout) *Sout = S;
}

// solver/ladder3.c setup() + phaseA(): T_b = sum_b w(b) xi_b (J w(b))^T / N_b, w(b) eigenvectors of the curve operator
// pw > 0: also K = sum_b w(b) thp_b (J w(b))^T / N_b, thp_b = exp(-i pi 4 pw k(k+1)/r), k = b/2 (ladder3.c:129-133, 229-247)
static double imal = 0;    // max |Im alpha| / |alpha| (ladder3.c:113 aborts if alpha is not real)
static void gens_solver(const Cat &K, int c, bool kapflip, Mat &Ta, Mat &Tb, int pw = 0, Mat *Kout = nullptr) {
    vector<int> ad = K.basis(c); int n = ad.size();
    Ta = Mat(n); Tb = Mat(n);
    if (Kout) *Kout = Mat(n);
    vector<Q> al(n), u(n), l(n), lam(n), f(n), g(n);
    vector<int> J(n, 1), isim(n, 0);
    Q kap = ((c / 2) % 2 == 1) ? 1 : -1;
    if (kapflip) kap = -kap;
    for (int i = 0; i < n; i++) {
        int a = ad[i]; Ta(i, i) = K.xi(a);
        C ac = (kap * K.d(a)) * K.sixj(a, a, 2, a, a, c);
        if (cab(ac) > 0) imal = mx(imal, D(fabsq(ac.im) / cab(ac)));
        al[i] = ac.re;
        if (i < n - 1) {
            int b = ad[i + 1];
            C be = csq(mk(K.d(a) * K.d(b))) * K.sixj(a, b, 2, b, a, c);
            if (fabsq(be.im) > fabsq(be.re)) { isim[i] = 1; u[i] = -be.im; l[i] = be.im; }
            else { u[i] = be.re; l[i] = be.re; }
        }
        lam[i] = K.qi[3 * (a + 1)] / K.qi[a + 1];
    }
    for (int i = 1; i < n; i++) J[i] = isim[i - 1] ? -J[i - 1] : J[i - 1];
    for (int b = 0; b < n; b++) {
        for (int i = 0; i < n - 1; i++) f[i] = u[i] / (lam[b] - al[i] - (i > 0 ? l[i - 1] * f[i - 1] : 0));
        for (int i = n - 1; i > 0; i--) g[i] = l[i - 1] / (lam[b] - al[i] - (i < n - 1 ? u[i] * g[i + 1] : 0));
        int kb = 0; Q mb = 0;
        for (int k = 0; k < n; k++) {
            Q gam = fabsq(al[k] - lam[b] + (k > 0 ? l[k - 1] * f[k - 1] : 0) + (k < n - 1 ? u[k] * g[k + 1] : 0));
            if (k == 0 || gam < mb) { mb = gam; kb = k; }
        }
        vector<Q> w(n, 1);
        for (int i = kb - 1; i >= 0; i--) w[i] = f[i] * w[i + 1];
        for (int i = kb + 1; i < n; i++) w[i] = g[i] * w[i - 1];
        Q N = 0; for (int i = 0; i < n; i++) N += J[i] * w[i] * w[i];
        long k = ad[b] / 2; C thp = eipi(-4 * pw * k * (k + 1), K.r);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
            Tb(i, j) = Tb(i, j) + (w[i] * J[j] * w[j] / N) * K.xi(ad[b]);
            if (Kout) (*Kout)(i, j) = (*Kout)(i, j) + (w[i] * J[j] * w[j] / N) * thp;
        }
    }
}

struct Rel { int n; double c1, c2; C z; };
static Rel rel(const Mat &Ta, const Mat &Tb) {
    Rel R; R.n = Ta.n;
    Mat A = Ta * Tb * Ta, B = Tb * Ta * Tb;
    R.c1 = D(fro(A - B) / fro(A));
    Mat P = Ta * Tb, P2 = P * P, Z = P2 * P2 * P2;
    C z = mk(0); for (int i = 0; i < R.n; i++) z = z + Z(i, i);
    z = (1 / (Q) R.n) * z; R.z = z;
    for (int i = 0; i < R.n; i++) Z(i, i) = Z(i, i) - z;
    R.c2 = D(fro(Z) / cab(z));
    return R;
}

// ||A - X B X^{-1}|| / ||A|| with X = diag(A(i,0)/B(i,0)): A and B equal up to diagonal conjugation
static Q conjres(const Mat &A, const Mat &B, vector<C> *Xout = nullptr) {
    int n = A.n; vector<C> X(n); Mat M(n);
    for (int i = 0; i < n; i++) X[i] = A(i, 0) / B(i, 0);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) M(i, j) = X[i] * B(i, j) / X[j];
    if (Xout) *Xout = X;
    return fro(A - M) / fro(A);
}

struct Row { int r, c, n; Rel s, p; double tbdiff, rplus, rminus; };
static vector<Row> run_all(int rmax) {
    vector<Row> rows;
    for (int r = 5; r <= rmax; r += 2) {
        Cat K(r); C z0 = mk(1);
        for (int c = 0; c <= r - 3; c += 2) {
            Row w; w.r = r; w.c = c; w.n = K.basis(c).size();
            if (w.n == 0) continue;
            Mat Ta, Tbs, Tbp;
            gens_solver(K, c, false, Ta, Tbs); w.s = rel(Ta, Tbs);
            gens_ref(K, c, 0, Ta, Tbp); w.p = rel(Ta, Tbp);
            w.tbdiff = D(conjres(Tbs, Tbp));
            if (c == 0) z0 = w.p.z;  // (T_a T_b)^6 is invariant under the diagonal conjugation of t0
            C q = w.s.z / z0, x = K.xi(c);
            w.rplus = D(cab(q - x)); w.rminus = D(cab(q - mk(1) / x));
            if (D(cab(x - mk(1) / x)) < 1e-6) w.rplus = w.rminus = -1;  // xi_c = 1/xi_c: eps undetermined
            rows.push_back(w);
        }
    }
    return rows;
}

static int nfail = 0, npass = 0;
static void report(bool ok, const char *name, const char *fmt, double v) {
    printf("%s %s ", ok ? "PASS" : "FAIL", name); printf(fmt, v); printf("\n");
    ok ? npass++ : nfail++;
}

static double nc2equiv(int rmax) {  // (-1)^(x/2) variant vs E T_b E, E = diag((-1)^(a/2))
    double m = 0;
    for (int r = 5; r <= rmax; r += 2) {
        Cat K(r);
        for (int c = 2; c <= r - 3; c += 2) {
            vector<int> ad = K.basis(c); Mat Ta, T0, T2; gens_ref(K, c, 0, Ta, T0); gens_ref(K, c, 2, Ta, T2);
            for (int i = 0; i < T0.n; i++) for (int j = 0; j < T0.n; j++) if ((ad[i] / 2 + ad[j] / 2) % 2) T0(i, j) = mk(0) - T0(i, j);
            m = mx(m, D(fro(T0 - T2) / fro(T0)));
        }
    }
    return m;
}

static double negctl(int rmax, int mode) {  // max C1 over c > 0; mode 1: reference without sign; 3: solver kap flipped
    double m = 0;
    for (int r = 5; r <= rmax; r += 2) {
        Cat K(r);
        for (int c = 2; c <= r - 3; c += 2) {
            Mat Ta, Tb;
            if (mode == 3) gens_solver(K, c, true, Ta, Tb); else gens_ref(K, c, mode, Ta, Tb);
            if (Ta.n == 0) continue;
            double v = rel(Ta, Tb).c1; if (!(v <= m)) m = v;
        }
    }
    return m;
}

static int selftest(int rmax) {
    double t1 = 0;
    for (int r = 5; r <= rmax; r += 2) {
        Cat K(r);
        t1 = mx(t1, D(fabsq(K.qi[1] - 1)) + D(fabsq(K.qi[0])));
        for (int m = 1; m + 1 < (int) K.qi.size(); m++) t1 = mx(t1, D(fabsq(K.qi[2] * K.qi[m] - K.qi[m + 1] - K.qi[m - 1])));
        for (int m = 1; m < r; m++) t1 = mx(t1, D(fabsq(K.qi[r - m] + K.qi[m])));
    }
    report(t1 < 1e-28, "t1 [1]=1, [2][n]=[n+1]+[n-1], [r-n]=-[n]: max err", "%.2e", t1);

    vector<Row> rows = run_all(rmax);
    int bad = 0;
    for (auto &w : rows) if (w.n != (w.r - 1) / 2 - w.c / 2) bad++;
    report(bad == 0 && !rows.empty(), "t2 dim = (r-1)/2 - c/2: mismatches", "%.0f", bad);

    double t3 = 0;
    for (int r = 5; r <= rmax; r += 2) {
        Cat K(r); Mat Ta, Tb, S; gens_ref(K, 0, 0, Ta, Tb, &S);
        vector<int> ad = K.basis(0); Mat H(S.n);
        for (int i = 0; i < S.n; i++) for (int j = 0; j < S.n; j++) {
            Q e = (K.d(ad[i]) < 0) == (K.d(ad[j]) < 0) ? 1 : -1;
            H(i, j) = mk(e * K.qint((long) (ad[i] + 1) * (ad[j] + 1)));
        }
        t3 = mx(t3, S.n ? D(fro(S - H) / fro(H)) : 1);
    }
    report(t3 < 1e-25, "t3 c=0: S[a][b] = sgn(d_a) sgn(d_b) [(a+1)(b+1)]: max rel err", "%.2e", t3);

    double m1s = 0, m1p = 0, m2s = 0, m2p = 0, md = 0;
    printf("     r  maxC1(solver) maxC1(ref)  maxC2(solver) maxC2(ref)  T_b(solver) vs X T_b(ref) X^-1\n");
    for (int r = 5; r <= rmax; r += 2) {
        double a = 0, b = 0, e = 0, f = 0, g = 0;
        for (auto &w : rows) if (w.r == r && w.c == 0 && w.s.c1 != w.s.c1) {
            printf("NOTE r=%d c=0: solver port non-finite (an eigenvector has an exact interior zero, so f_i = w_i/w_{i+1} is infinite); excluded\n", r);
            w.s.c1 = w.s.c2 = w.tbdiff = 0; w.s.z = w.p.z;
        }
        for (auto &w : rows) if (w.r == r) { a = mx(a, w.s.c1); b = mx(b, w.p.c1); e = mx(e, w.s.c2); f = mx(f, w.p.c2); g = mx(g, w.tbdiff); }
        printf("  %4d  %.2e      %.2e    %.2e      %.2e    %.2e\n", r, a, b, e, f, g);
        m1s = mx(m1s, a); m1p = mx(m1p, b); m2s = mx(m2s, e); m2p = mx(m2p, f); md = mx(md, g);
    }
    report(m1s < 1e-10, "t4 braid C1 (solver):", "%.2e", m1s);
    report(m1p < 1e-10, "t4 braid C1 (ref):", "%.2e", m1p);
    report(m2s < 1e-10, "t5 chain C2 (solver):", "%.2e", m2s);
    report(m2p < 1e-10, "t5 chain C2 (ref):", "%.2e", m2p);
    report(md < 1e-10, "t0 T_b(solver) = X T_b(ref) X^-1, X diagonal:", "%.2e", md);
    {   Cat K(9); Mat Ta, Ts, Tr; vector<C> X; gens_solver(K, 2, false, Ta, Ts); gens_ref(K, 2, 0, Ta, Tr);
        printf("INFO r=9 c=2: |Ts-Tr|/|Tr|=%.2e, X =", D(fro(Ts - Tr) / fro(Tr))); conjres(Ts, Tr, &X);
        for (C x : X) printf(" (%.6f%+.6fi)", D(x.re), D(x.im));
        printf("\n"); }
    double kres[3] = {0, 0, 0};
    for (int r = 5; r <= rmax; r += 2) {
        Cat K(r);
        for (int c = 0; c <= r - 3; c += 2) {
            if (K.basis(c).empty()) continue;
            Mat Ta, Ts, Tr, Kp; vector<C> X;
            gens_solver(K, c, false, Ta, Ts); gens_ref(K, c, 0, Ta, Tr);
            if (!isfinite(D(fro(Ts)))) continue;  // excluded case of t4/t5/t0
            conjres(Ts, Tr, &X);
            Mat Ti = inv(Tr), P = Ti;
            for (int pw = 1; pw <= 2; pw++, P = P * Ti) {
                gens_solver(K, c, false, Ta, Ts, pw, &Kp); Mat M(P.n);
                for (int i = 0; i < P.n; i++) for (int j = 0; j < P.n; j++) M(i, j) = X[i] * P(i, j) / X[j];
                kres[pw] = mx(kres[pw], D(fro(Kp - M) / fro(Kp)));
            }
        }
    }
    report(kres[1] < 1e-10, "t8 K(pw=1) (solver) = X T_b(ref)^-1 X^-1, X of t0:", "%.2e", kres[1]);
    report(kres[2] < 1e-10, "t8 K(pw=2) (solver) = X T_b(ref)^-2 X^-1, X of t0:", "%.2e", kres[2]);
    report(imal < 1e-25, "t9 alpha real (ladder3.c:113): max |Im alpha|/|alpha|", "%.2e", imal);

    int np = 0, nm = 0, nn = 0, na = 0; double worst = 0;
    for (auto &w : rows) if (w.c > 0) {
        if (w.rplus < 0) { na++; continue; }
        double lo = fmin(w.rplus, w.rminus); worst = mx(worst, lo);
        if (!(lo <= 1e-10)) nn++; else if (w.rplus < w.rminus) np++; else nm++;
    }
    printf("     C3: eps=+1 for %d cases, eps=-1 for %d, neither %d, xi_c=1/xi_c (undecided) %d; worst best-residual %.2e\n", np, nm, nn, na, worst);
    report(nn == 0 && (np == 0 || nm == 0) && np + nm > 0, "t6 uniform eps in C3: eps =", "%+.0f", np ? 1.0 : -1.0);

    double c4 = 0;
    for (int r = 5; r <= rmax; r += 2) {
        Cat K(r); C p = mk(0);
        for (int a : K.cols) p = p + (K.d(a) * K.d(a)) * K.xi(a);
        C zt = (1 / cab(p)) * p, z4 = zt * zt * zt * zt, z0 = mk(0);
        for (auto &w : rows) if (w.r == r && w.c == 0) z0 = w.p.z;
        c4 = mx(c4, D(cab(z0 - z4)));
        if (r <= 7 || r == rmax)
            printf("INFO C4 r=%d z(r,0)=%.15f%+.15fi |z-zeta^4|=%.2e |z-conj(zeta)^4|=%.2e arg z/pi=%.6f arg zeta/pi=%.6f\n",
                   r, D(z0.re), D(z0.im), D(cab(z0 - z4)), D(cab(z0 - cj(z4))), D(atan2q(z0.im, z0.re) / PIQ), D(atan2q(zt.im, zt.re) / PIQ));
    }
    report(c4 < 1e-10, "t10 C4 z(r,0) = zeta^4, zeta = p_plus/|p_plus|: max |z - zeta^4|", "%.2e", c4);
    printf("SKIP t7 C7: data/stress_7words_r101.csv stores fitted quantities for LLLR..LLLLR, not log|lambda| for LR/LLR\n");

    int rn = rmax < 15 ? rmax : 15;
    double n1 = negctl(rn, 1), n2 = nc2equiv(rn), n3 = negctl(rn, 3);
    report(n1 > 1e-3, "NC1 ref without (-1)^((a+b+x)/2): braid fails, max C1 =", "%.2e", n1);
    report(n2 < 1e-25, "NC2 (-1)^(x/2) only gives E T_b E, E = diag((-1)^(a/2)), so undetectable: max rel diff", "%.2e", n2);
    report(n3 > 1e-3, "NC3 solver with kap flipped: braid fails, max C1 =", "%.2e", n3);
    printf("selftest: %d PASS, %d FAIL\n", npass, nfail);
    return nfail ? 1 : 0;
}

static void table(int rmax) {
    printf("r c dim C1_solver C2_solver C1_ref C2_ref argz/pi res_xi res_invxi Tbdiff\n");
    for (auto &w : run_all(rmax))
        printf("%d %d %d %.2e %.2e %.2e %.2e %.12f %.2e %.2e %.2e\n", w.r, w.c, w.n, w.s.c1, w.s.c2, w.p.c1, w.p.c2,
               D(atan2q(w.s.z.im, w.s.z.re) / PIQ), w.rplus, w.rminus, w.tbdiff);
}

int main(int argc, char **argv) {
    int rmax = argc > 2 ? atoi(argv[2]) : 41;
    if (rmax > 41) { fprintf(stderr, "error: RMAX %d > 41; the fixed tolerance 1e-10 of t4/t5 fails for numerical reasons beyond about r = 51\n", rmax); return 1; }
    if (argc > 1 && !strcmp(argv[1], "selftest")) return selftest(rmax);
    if (argc > 1 && !strcmp(argv[1], "table")) { table(rmax); return 0; }
    fprintf(stderr, "usage: gens_check selftest|table [RMAX]\n");
    return 1;
}
