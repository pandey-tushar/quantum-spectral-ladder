/* ladder2: top eigenvalues (log|lambda|) of SO(3) quantum representations of the once-punctured torus
   mapping class group at q = e^{2 pi i / r}, phi = LR (A = T_b^{-1} T_a) or LLR (A = T_b^{-2} T_a),
   puncture color c.  Same numbers and output line as ladder.c; see README_ladder2.md.

   Method
   1. The loop operator Y (tridiagonal: alpha real, beta real except a few purely imaginary entries) is made
      real by the diagonal similarity S = diag(i^{e_i}):  Y' = S^-1 Y S, J = S^2 = diag(+-1).
   2. Eigenvectors w_b of Y' (eigenvalue [3(b+1)]/[b+1]) by twisted factorization in real arithmetic at
      precision d1.  They are never stored: only rows 0 and n-1 of K = S^-1 T_b^-p S
      (K = sum_b theta_b^-p w_b w_b^T J / (w_b^T J w_b)) are accumulated.  O(n^2), parallel over b.
   3. K commutes with Y', so every row of K follows from rows 0 and n-1 by a three-term recurrence, run from
      both ends and meeting at the largest |alpha| (precision d2).  The mismatch of the two sweeps at the
      meeting row measures the accuracy of K.  O(n^2).
   4. Block subspace iteration with A' = K Theta (similar to A) at precision d3.  K is stored with one
      fixed-point exponent per row, as residues modulo ~59-bit primes ("mm" kernel, default): each
      product K*Y is a set of independent nmod_mat_mul calls (one per prime, run in parallel) followed by CRT.
      The "acb" kernel stores K as acb_mat (about half the memory, several times slower).

   usage: ladder2 r c word [options]
     -t T        threads (default: OMP_NUM_THREADS or all cores)
     -kernel K   mm (default) | acb
     -d1/-d2/-d3 D   override the precisions (decimal digits) of the three stages (default: automatic)
     -m M  -ntop K  -tol E  -maxit N   block size (12), eigenvalues printed (8), stop when the max change of
                 log|lambda_1..K| between iterations is < 1e-E (22), iteration cap (100)
     -ck FILE    checkpoint: X is written after every iteration and read at start if FILE matches
     -q          quiet
   stdout: r c n theta word it log|lambda_1..ntop| (30 digits), exactly as ladder.c.
   stderr: precisions, accuracy self-checks, timings, peak memory. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/nmod_mat.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>

#define LOG10_2 0.30102999566398120
#define DIG2BITS(d) ((slong)((d) * 3.3219280948873623) + 32)

static slong r, cc, n, pw;
static slong *ad;
static int quiet = 0;

static double wall(void) { struct timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + 1e-6 * tv.tv_usec; }
static double maxrss_gb(void) { struct rusage u; getrusage(RUSAGE_SELF, &u); return u.ru_maxrss / 1048576.0; }
static inline void amid(arb_t x) { mag_zero(arb_radref(x)); }
static inline void cmid(acb_t z) { mag_zero(arb_radref(acb_realref(z))); mag_zero(arb_radref(acb_imagref(z))); }
/* upper bound for log2|x| (WORD_MIN for 0) */
static inline slong lexp(const arb_t x) { return arf_is_zero(arb_midref(x)) ? WORD_MIN : arf_abs_bound_lt_2exp_si(arb_midref(x)); }
static inline slong lexpc(const acb_t z) { slong a = lexp(acb_realref(z)), b = lexp(acb_imagref(z)); return a > b ? a : b; }

/* ------------------------------------------------------------------ 6j symbols */
static int adm(slong a, slong b, slong c) {
    return ((a + b + c) % 2 == 0) && labs(a - b) <= c && c <= a + b && a + b + c <= 2 * (r - 2);
}
static void Delta(acb_t res, arb_srcptr fac, slong a, slong b, slong c, slong prec) {
    arb_t v; arb_init(v);
    arb_mul(v, fac + (a + b - c) / 2, fac + (a - b + c) / 2, prec);
    arb_mul(v, v, fac + (-a + b + c) / 2, prec);
    arb_div(v, v, fac + (a + b + c) / 2 + 1, prec);
    acb_set_arb(res, v); acb_sqrt(res, res, prec);
    arb_clear(v);
}
static void sixj(acb_t res, arb_srcptr fac, slong j1, slong j2, slong j3, slong j4, slong j5, slong j6, slong prec) {
    slong T[4][3] = {{j1,j2,j3},{j1,j5,j6},{j4,j2,j6},{j4,j5,j3}};
    for (int k = 0; k < 4; k++) if (!adm(T[k][0],T[k][1],T[k][2])) { acb_zero(res); return; }
    acb_t pref, t; acb_init(pref); acb_init(t); acb_one(pref);
    for (int k = 0; k < 4; k++) { Delta(t, fac, T[k][0],T[k][1],T[k][2], prec); acb_mul(pref, pref, t, prec); }
    slong a1[4] = {(j1+j2+j3)/2,(j1+j5+j6)/2,(j4+j2+j6)/2,(j4+j5+j3)/2};
    slong b1[3] = {(j1+j2+j4+j5)/2,(j2+j3+j5+j6)/2,(j3+j1+j6+j4)/2};
    slong zl = a1[0], zh = b1[0];
    for (int k = 1; k < 4; k++) if (a1[k] > zl) zl = a1[k];
    for (int k = 1; k < 3; k++) if (b1[k] < zh) zh = b1[k];
    arb_t tot, den, term; arb_init(tot); arb_init(den); arb_init(term);
    for (slong z = zl; z <= zh; z++) {
        arb_one(den);
        for (int k = 0; k < 4; k++) arb_mul(den, den, fac + (z - a1[k]), prec);
        for (int k = 0; k < 3; k++) arb_mul(den, den, fac + (b1[k] - z), prec);
        arb_div(term, fac + z + 1, den, prec);
        if (z % 2) arb_sub(tot, tot, term, prec); else arb_add(tot, tot, term, prec);
    }
    acb_mul_arb(res, pref, tot, prec);
    arb_clear(tot); arb_clear(den); arb_clear(term); acb_clear(pref); acb_clear(t);
}

/* ------------------------------------------------------------------ problem data at precision p1 */
static arb_ptr al, u, l, lam;   /* real tridiagonal Y' (diag al, upper u, lower l), eigenvalues lam */
static int *J;
static acb_ptr thp;             /* theta_b^{-pw} */

static void setup(slong prec) {
    slong nqi = 3 * r + 5;
    arb_ptr qi = _arb_vec_init(nqi), fac = _arb_vec_init(nqi);
    arb_t pi, x, s1; arb_init(pi); arb_init(x); arb_init(s1);
    arb_const_pi(pi, prec);
    arb_mul_si(x, pi, 2, prec); arb_div_si(x, x, r, prec); arb_sin(s1, x, prec);
    for (slong k = 0; k < nqi; k++) { arb_mul_si(x, pi, 2 * k, prec); arb_div_si(x, x, r, prec); arb_sin(qi + k, x, prec); arb_div(qi + k, qi + k, s1, prec); }
    arb_one(fac);
    for (slong k = 1; k < nqi; k++) arb_mul(fac + k, fac + k - 1, qi + k, prec);
    al = _arb_vec_init(n); u = _arb_vec_init(n); l = _arb_vec_init(n); lam = _arb_vec_init(n);
    J = malloc(sizeof(int) * n);
    int *isim = calloc(n, sizeof(int));
    slong kap = ((cc / 2) % 2 == 1) ? 1 : -1;
    acb_t t, uu, ac, be; acb_init(t); acb_init(uu); acb_init(ac); acb_init(be);
    for (slong i = 0; i < n; i++) {
        slong a = ad[i];
        sixj(t, fac, a, a, 2, a, a, cc, prec);
        acb_mul_arb(ac, t, qi + a + 1, prec); acb_mul_si(ac, ac, kap, prec);
        if (arf_cmpabs(arb_midref(acb_imagref(ac)), arb_midref(acb_realref(ac))) > 0) { fprintf(stderr, "error: alpha not real at %ld\n", (long)i); exit(2); }
        arb_set(al + i, acb_realref(ac)); amid(al + i);
        if (i < n - 1) {
            slong b = ad[i + 1];
            acb_set_arb(uu, qi + a + 1); acb_mul_arb(uu, uu, qi + b + 1, prec); acb_sqrt(uu, uu, prec);
            sixj(t, fac, a, b, 2, b, a, cc, prec); acb_mul(be, uu, t, prec);
            if (arf_cmpabs(arb_midref(acb_imagref(be)), arb_midref(acb_realref(be))) > 0) {
                isim[i] = 1;   /* beta = i b: upper -> -b, lower -> b */
                arb_neg(u + i, acb_imagref(be)); arb_set(l + i, acb_imagref(be));
            } else { arb_set(u + i, acb_realref(be)); arb_set(l + i, acb_realref(be)); }
            amid(u + i); amid(l + i);
        }
        arb_div(lam + i, qi + 3 * (a + 1), qi + a + 1, prec); amid(lam + i);
    }
    J[0] = 1;
    for (slong i = 1; i < n; i++) J[i] = isim[i - 1] ? -J[i - 1] : J[i - 1];
    thp = _acb_vec_init(n);
    for (slong i = 0; i < n; i++) {
        slong k = ad[i] / 2;
        arb_mul_si(x, pi, -4 * pw * k * (k + 1), prec); arb_div_si(x, x, r, prec);
        acb_zero(thp + i); arb_set(acb_imagref(thp + i), x); acb_exp(thp + i, thp + i, prec); cmid(thp + i);
    }
    acb_clear(t); acb_clear(uu); acb_clear(ac); acb_clear(be);
    arb_clear(pi); arb_clear(x); arb_clear(s1);
    _arb_vec_clear(qi, nqi); _arb_vec_clear(fac, nqi); free(isim);
}
static void setup_clear(void) {
    _arb_vec_clear(al, n); _arb_vec_clear(u, n); _arb_vec_clear(l, n); _arb_vec_clear(lam, n);
    _acb_vec_clear(thp, n); free(J);
}

/* ------------------------------------------------------------------ phase A: eigenvectors -> rows 0, n-1 of K */
static double lcomp;   /* log10 of the completeness residual  max_j |sum_b w_b0 w_bj J_j / N_b - delta_0j| (rows 0 and n-1) */
static void phaseA(acb_ptr R0, acb_ptr R1, double *lvmax, slong prec) {
    double lvm = -1e300;
    arb_ptr E0 = _arb_vec_init(n), E1 = _arb_vec_init(n);
    _acb_vec_zero(R0, n); _acb_vec_zero(R1, n);
    #pragma omp parallel
    {
        arb_ptr f = _arb_vec_init(n), g = _arb_vec_init(n), w = _arb_vec_init(n);
        acb_ptr r0 = _acb_vec_init(n), r1 = _acb_vec_init(n);
        arb_ptr e0 = _arb_vec_init(n), e1 = _arb_vec_init(n);
        arb_t t0, t1; arb_init(t0); arb_init(t1);
        arb_t t, gam, N; arb_init(t); arb_init(gam); arb_init(N);
        acb_t c0, c1; acb_init(c0); acb_init(c1);
        mag_t mb, mc; mag_init(mb); mag_init(mc);
        double lv = -1e300;
        #pragma omp for schedule(dynamic, 1)
        for (slong b = 0; b < n; b++) {
            /* forward:  f_i = u_i / (lam - al_i - l_{i-1} f_{i-1}),  w_i = f_i w_{i+1} */
            for (slong i = 0; i < n - 1; i++) {
                arb_sub(t, lam + b, al + i, prec);
                if (i > 0) arb_submul(t, l + i - 1, f + i - 1, prec);
                arb_div(f + i, u + i, t, prec); amid(f + i);
            }
            /* backward: g_i = l_{i-1} / (lam - al_i - u_i g_{i+1}),  w_i = g_i w_{i-1} */
            for (slong i = n - 1; i > 0; i--) {
                arb_sub(t, lam + b, al + i, prec);
                if (i < n - 1) arb_submul(t, u + i, g + i + 1, prec);
                arb_div(g + i, l + i - 1, t, prec); amid(g + i);
            }
            /* twist index: argmin |gamma_k| */
            slong kb = 0;
            for (slong k = 0; k < n; k++) {
                arb_sub(gam, al + k, lam + b, prec);
                if (k > 0) arb_addmul(gam, l + k - 1, f + k - 1, prec);
                if (k < n - 1) arb_addmul(gam, u + k, g + k + 1, prec);
                arb_get_mag(mc, gam);
                if (k == 0 || mag_cmp(mc, mb) < 0) { mag_set(mb, mc); kb = k; }
            }
            arb_one(w + kb);
            for (slong i = kb - 1; i >= 0; i--) { arb_mul(w + i, f + i, w + i + 1, prec); amid(w + i); }
            for (slong i = kb + 1; i < n; i++) { arb_mul(w + i, g + i, w + i - 1, prec); amid(w + i); }
            arb_zero(N);
            for (slong i = 0; i < n; i++) { if (J[i] > 0) arb_addmul(N, w + i, w + i, prec); else arb_submul(N, w + i, w + i, prec); }
            amid(N);
            {   /* log10 of max_i |w_i| / sqrt|N|  = max entry of the normalized eigenvector matrix V */
                slong em = WORD_MIN;
                for (slong i = 0; i < n; i++) { slong e = lexp(w + i); if (e > em) em = e; }
                double e2 = (em - 0.5 * lexp(N)) * LOG10_2;
                if (e2 > lv) lv = e2;
            }
            arb_div(t0, w + 0, N, prec); acb_mul_arb(c0, thp + b, t0, prec);
            arb_div(t1, w + n - 1, N, prec); acb_mul_arb(c1, thp + b, t1, prec);
            for (slong j = 0; j < n; j++) {
                if (J[j] > 0) { arb_addmul(e0 + j, t0, w + j, prec); arb_addmul(e1 + j, t1, w + j, prec); }
                else { arb_submul(e0 + j, t0, w + j, prec); arb_submul(e1 + j, t1, w + j, prec); }
                amid(e0 + j); amid(e1 + j);
                if (J[j] > 0) { acb_addmul_arb(r0 + j, c0, w + j, prec); acb_addmul_arb(r1 + j, c1, w + j, prec); }
                else { acb_submul_arb(r0 + j, c0, w + j, prec); acb_submul_arb(r1 + j, c1, w + j, prec); }
                cmid(r0 + j); cmid(r1 + j);
            }
        }
        #pragma omp critical
        {
            if (lv > lvm) lvm = lv;
            for (slong j = 0; j < n; j++) { acb_add(R0 + j, R0 + j, r0 + j, prec); acb_add(R1 + j, R1 + j, r1 + j, prec); cmid(R0 + j); cmid(R1 + j);
                                            arb_add(E0 + j, E0 + j, e0 + j, prec); arb_add(E1 + j, E1 + j, e1 + j, prec); }
        }
        _arb_vec_clear(e0, n); _arb_vec_clear(e1, n); arb_clear(t0); arb_clear(t1);
        _arb_vec_clear(f, n); _arb_vec_clear(g, n); _arb_vec_clear(w, n);
        _acb_vec_clear(r0, n); _acb_vec_clear(r1, n);
        arb_clear(t); arb_clear(gam); arb_clear(N); acb_clear(c0); acb_clear(c1); mag_clear(mb); mag_clear(mc);
    }
    *lvmax = lvm;
    arb_sub_ui(E0 + 0, E0 + 0, 1, prec); arb_sub_ui(E1 + n - 1, E1 + n - 1, 1, prec);
    slong em = WORD_MIN;
    for (slong j = 0; j < n; j++) { slong e = lexp(E0 + j); if (e > em) em = e; e = lexp(E1 + j); if (e > em) em = e; }
    lcomp = em == WORD_MIN ? -1e300 : em * LOG10_2;
    _arb_vec_clear(E0, n); _arb_vec_clear(E1, n);
}

/* ------------------------------------------------------------------ K storage */
static int kern_mm = 2;                /* 0: acb, 1: mmsym (symmetric block storage), 2: mm (full) */
static acb_mat_t Kacb;                 /* acb kernel */
/* mm (symmetric) kernel: M' = K J is symmetric; its lower block triangle (block size BS) is stored, one global
   fixed-point exponent Eglob: M'_ij = k_ij 2^(Eglob - BK) */
static slong BS = 256;
static slong nbk, nblk, Eglob = 0, Emaxseen = WORD_MIN;
static nmod_mat_struct *Bm;            /* [(t*2 + part)*nblk + I(I+1)/2 + J] */
#define KSLACK 16                      /* bits of headroom above Eglob */
static slong np = 0, BK, BY;           /* mm kernel: number of primes, fixed-point bits of K and Y */
static mp_limb_t *primes;
static nmod_mat_struct *Kr, *Ki, *Ym;  /* per prime */
static slong *Erow;                    /* row exponents: K_ij = k_ij 2^(Erow_i - BK) */
static fmpz_comb_t comb;
static int kalloc = 0;

static slong nprimes_for(slong p3) {
    slong lg = 1; while ((WORD(1) << lg) < n) lg++;
    return (2 * p3 + lg + 6 + (kern_mm == 1 ? KSLACK : 0) + 57) / 58;
}
static void kstore_alloc(slong p3) {
    if (!kern_mm) { acb_mat_init(Kacb, n, n); kalloc = 1; return; }
    BK = p3; BY = p3;
    np = nprimes_for(p3);
    primes = malloc(sizeof(mp_limb_t) * np);
    mp_limb_t p = UWORD(1) << 58;
    for (slong t = 0; t < np; t++) { p = n_nextprime(p, 1); primes[t] = p; }
    fmpz_comb_init(comb, primes, np);
    Ym = NULL; Emaxseen = WORD_MIN;
    if (kern_mm == 2) {
        Kr = malloc(sizeof(nmod_mat_struct) * np); Ki = malloc(sizeof(nmod_mat_struct) * np);
        for (slong t = 0; t < np; t++) { nmod_mat_init(Kr + t, n, n, primes[t]); nmod_mat_init(Ki + t, n, n, primes[t]); }
        Erow = malloc(sizeof(slong) * n);
    } else {
        nbk = (n + BS - 1) / BS; nblk = nbk * (nbk + 1) / 2;
        Bm = malloc(sizeof(nmod_mat_struct) * np * 2 * nblk);
        for (slong t = 0; t < np; t++) for (int part = 0; part < 2; part++)
            for (slong I = 0; I < nbk; I++) for (slong Jb = 0; Jb <= I; Jb++) {
                slong ri = (I + 1) * BS < n ? BS : n - I * BS, cj = (Jb + 1) * BS < n ? BS : n - Jb * BS;
                nmod_mat_init(Bm + (t * 2 + part) * nblk + I * (I + 1) / 2 + Jb, ri, cj, primes[t]);
            }
    }
    kalloc = 1;
}
static void kstore_free(void) {
    if (!kalloc) return;
    if (!kern_mm) acb_mat_clear(Kacb);
    else {
        if (kern_mm == 2) {
            for (slong t = 0; t < np; t++) { nmod_mat_clear(Kr + t); nmod_mat_clear(Ki + t); }
            free(Kr); free(Ki); free(Erow);
        } else {
            for (slong q = 0; q < np * 2 * nblk; q++) nmod_mat_clear(Bm + q);
            free(Bm);
        }
        if (Ym != NULL) { for (slong t = 0; t < np; t++) nmod_mat_clear(Ym + t); free(Ym); Ym = NULL; }
        fmpz_comb_clear(comb); free(primes);
    }
    kalloc = 0;
}
/* store row i of K (given at higher precision); parallel over j inside */
static void kstore_row(slong i, acb_srcptr row, slong p3, slong *rowexp) {
    slong em = WORD_MIN;
    for (slong j = 0; j < n; j++) { slong e = lexpc(row + j); if (e > em) em = e; }
    *rowexp = em;
    if (!kern_mm) {
        #pragma omp parallel for schedule(static)
        for (slong j = 0; j < n; j++) acb_set_round(acb_mat_entry(Kacb, i, j), row + j, p3);
        return;
    }
    if (kern_mm == 1) {
        if (em > Emaxseen) Emaxseen = em;
        slong I = i / BS, jend = (I + 1) * BS < n ? (I + 1) * BS : n;
        #pragma omp parallel
        {
            fmpz_t z; fmpz_init(z);
            fmpz_comb_temp_t ct; fmpz_comb_temp_init(ct, comb);
            mp_limb_t *res = malloc(sizeof(mp_limb_t) * np);
            #pragma omp for schedule(static)
            for (slong j = 0; j < jend; j++) {
                slong Jb = j / BS, bi = I * (I + 1) / 2 + Jb;
                for (int part = 0; part < 2; part++) {
                    arf_get_fmpz_fixed_si(z, arb_midref(part ? acb_imagref(row + j) : acb_realref(row + j)), Eglob - BK);
                    if (J[j] < 0) fmpz_neg(z, z);          /* M' = K J */
                    fmpz_multi_mod_ui(res, z, comb, ct);
                    for (slong t = 0; t < np; t++) nmod_mat_entry(Bm + (t * 2 + part) * nblk + bi, i - I * BS, j - Jb * BS) = res[t];
                }
            }
            free(res); fmpz_comb_temp_clear(ct); fmpz_clear(z);
        }
        return;
    }
    if (em == WORD_MIN) em = 0;
    Erow[i] = em;
    #pragma omp parallel
    {
        fmpz_t z; fmpz_init(z);
        fmpz_comb_temp_t ct; fmpz_comb_temp_init(ct, comb);
        mp_limb_t *res = malloc(sizeof(mp_limb_t) * np);
        #pragma omp for schedule(static)
        for (slong j = 0; j < n; j++) {
            arf_get_fmpz_fixed_si(z, arb_midref(acb_realref(row + j)), em - BK);
            fmpz_multi_mod_ui(res, z, comb, ct);
            for (slong t = 0; t < np; t++) nmod_mat_entry(Kr + t, i, j) = res[t];
            arf_get_fmpz_fixed_si(z, arb_midref(acb_imagref(row + j)), em - BK);
            fmpz_multi_mod_ui(res, z, comb, ct);
            for (slong t = 0; t < np; t++) nmod_mat_entry(Ki + t, i, j) = res[t];
        }
        free(res); fmpz_comb_temp_clear(ct); fmpz_clear(z);
    }
}

/* ------------------------------------------------------------------ phase B: rows of K by the commuting recurrence */
static void phaseB(acb_srcptr R0, acb_srcptr R1, slong p2, slong p3, double *lkmax, double *qK) {
    kstore_alloc(p3);
    slong em_all = WORD_MIN, e;
    double lmis = -1e300;
    if (n < 3) {
        kstore_row(0, R0, p3, &e); if (e > em_all) em_all = e;
        if (n == 2) { kstore_row(1, R1, p3, &e); if (e > em_all) em_all = e; }
        *lkmax = em_all * LOG10_2; *qK = 1e9; return;
    }
    slong h = 0;
    for (slong i = 1; i < n; i++) if (arf_cmpabs(arb_midref(al + i), arb_midref(al + h)) > 0) h = i;
    if (h > n - 2) h = n - 2;
    arb_ptr al2 = _arb_vec_init(n), u2 = _arb_vec_init(n), l2 = _arb_vec_init(n), iu = _arb_vec_init(n), il = _arb_vec_init(n);
    for (slong i = 0; i < n; i++) {
        arb_set_round(al2 + i, al + i, p2);
        if (i < n - 1) {
            arb_set_round(u2 + i, u + i, p2); arb_set_round(l2 + i, l + i, p2);
            arb_inv(iu + i, u + i, p2); arb_inv(il + i, l + i, p2); amid(iu + i); amid(il + i);
        }
    }
    acb_ptr chk = _acb_vec_init(n);
    for (int dir = 0; dir < 2; dir++) {
        acb_ptr pr = _acb_vec_init(n), cu = _acb_vec_init(n), nx = _acb_vec_init(n);
        acb_srcptr R = dir == 0 ? R0 : R1;
        /* dir 0: rows 0 -> h (stored) and h+1 (check);  dir 1: rows n-1 -> h+1 (stored) */
        slong i0 = dir == 0 ? 0 : n - 1, st = dir == 0 ? 1 : -1;
        for (slong j = 0; j < n; j++) acb_set_round(cu + j, R + j, p2);
        kstore_row(i0, cu, p3, &e); if (e > em_all) em_all = e;
        for (slong i = i0; i != h + 1; i += st) {
            #pragma omp parallel
            {
                acb_t s; arb_t da; acb_init(s); arb_init(da);
                #pragma omp for schedule(static)
                for (slong j = 0; j < n; j++) {
                    /* (Y'K)_ij = (KY')_ij:  l_{i-1}K_{i-1,j} + al_i K_ij + u_i K_{i+1,j} = u_{j-1}K_{i,j-1} + al_j K_ij + l_j K_{i,j+1} */
                    arb_sub(da, al2 + j, al2 + i, p2);
                    acb_mul_arb(s, cu + j, da, p2);
                    if (j > 0) acb_addmul_arb(s, cu + j - 1, u2 + j - 1, p2);
                    if (j < n - 1) acb_addmul_arb(s, cu + j + 1, l2 + j, p2);
                    if (dir == 0) { if (i > 0) acb_submul_arb(s, pr + j, l2 + i - 1, p2); acb_mul_arb(nx + j, s, iu + i, p2); }
                    else { if (i < n - 1) acb_submul_arb(s, pr + j, u2 + i, p2); acb_mul_arb(nx + j, s, il + i - 1, p2); }
                    cmid(nx + j);
                }
                acb_clear(s); arb_clear(da);
            }
            if (dir == 1 || i + st <= h) { kstore_row(i + st, nx, p3, &e); if (e > em_all) em_all = e; }
            acb_ptr tmp = pr; pr = cu; cu = nx; nx = tmp;
        }
        if (dir == 0) _acb_vec_set(chk, cu, n);
        else {
            acb_t d; acb_init(d); slong em = WORD_MIN;
            for (slong j = 0; j < n; j++) { acb_sub(d, chk + j, cu + j, p2); slong ee = lexpc(d); if (ee > em) em = ee; }
            lmis = em == WORD_MIN ? -1e300 : em * LOG10_2;
            acb_clear(d);
        }
        _acb_vec_clear(pr, n); _acb_vec_clear(cu, n); _acb_vec_clear(nx, n);
    }
    _acb_vec_clear(chk, n);
    _arb_vec_clear(al2, n); _arb_vec_clear(u2, n); _arb_vec_clear(l2, n); _arb_vec_clear(iu, n); _arb_vec_clear(il, n);
    *lkmax = em_all * LOG10_2;
    *qK = *lkmax - lmis;
}

/* ------------------------------------------------------------------ products Z = K Y */
static void kmul(acb_mat_t Z, const acb_mat_t Y, slong prec) {
    slong m = acb_mat_ncols(Y);
    if (!kern_mm) {
        slong bs = 16, nb = (n + bs - 1) / bs;
        #pragma omp parallel for schedule(dynamic, 1)
        for (slong b = 0; b < nb; b++) {
            slong i0 = b * bs, i1 = i0 + bs < n ? i0 + bs : n;
            acb_mat_t Kw, Zw;
            acb_mat_window_init(Kw, Kacb, i0, 0, i1, n);
            acb_mat_window_init(Zw, Z, i0, 0, i1, m);
            acb_mat_approx_mul(Zw, Kw, Y, prec);
            acb_mat_window_clear(Kw); acb_mat_window_clear(Zw);
        }
        return;
    }
    if (Ym == NULL) {
        Ym = malloc(sizeof(nmod_mat_struct) * np);
        for (slong t = 0; t < np; t++) nmod_mat_init(Ym + t, n, 2 * m, primes[t]);
    }
    /* column exponents of Y (for mm: W = J Y is what multiplies M') */
    slong *F = malloc(sizeof(slong) * m);
    for (slong k = 0; k < m; k++) {
        slong em = WORD_MIN;
        for (slong j = 0; j < n; j++) { slong e = lexpc(acb_mat_entry(Y, j, k)); if (e > em) em = e; }
        F[k] = em == WORD_MIN ? 0 : em;
    }
    #pragma omp parallel
    {
        fmpz_t z; fmpz_init(z);
        fmpz_comb_temp_t ct; fmpz_comb_temp_init(ct, comb);
        mp_limb_t *res = malloc(sizeof(mp_limb_t) * np);
        #pragma omp for schedule(static)
        for (slong j = 0; j < n; j++) for (slong k = 0; k < m; k++) for (int part = 0; part < 2; part++) {
            arf_get_fmpz_fixed_si(z, arb_midref(part ? acb_imagref(acb_mat_entry(Y, j, k)) : acb_realref(acb_mat_entry(Y, j, k))), F[k] - BY);
            if (kern_mm == 1 && J[j] < 0) fmpz_neg(z, z);
            fmpz_multi_mod_ui(res, z, comb, ct);
            for (slong t = 0; t < np; t++) nmod_mat_entry(Ym + t, j, k + part * m) = res[t];
        }
        free(res); fmpz_comb_temp_clear(ct); fmpz_clear(z);
    }
    /* residues of Z: layout [(i*m + k)*2 + part][t] */
    mp_limb_t *RZ = malloc(sizeof(mp_limb_t) * n * m * 2 * np);
    #pragma omp parallel
    {
        nmod_mat_t ZR, ZI;
        nmod_mat_init(ZR, n, 2 * m, 3); nmod_mat_init(ZI, n, 2 * m, 3);
        nmod_mat_t WT; nmod_mat_init(WT, 2 * m, n, 3);
        #pragma omp for schedule(dynamic, 1)
        for (slong t = 0; t < np; t++) {
            mp_limb_t p = primes[t];
            nmod_mat_set_mod(ZR, p); nmod_mat_set_mod(ZI, p);
            if (kern_mm == 2) {
                nmod_mat_mul(ZR, Kr + t, Ym + t);
                nmod_mat_mul(ZI, Ki + t, Ym + t);
            } else {
                nmod_mat_set_mod(WT, p);
                nmod_mat_transpose(WT, Ym + t);
                for (int part = 0; part < 2; part++) {
                    nmod_mat_struct *Zp = part ? ZI : ZR;
                    nmod_mat_zero(Zp);
                    for (slong I = 0; I < nbk; I++) for (slong Jb = 0; Jb <= I; Jb++) {
                        nmod_mat_struct *B = Bm + (t * 2 + part) * nblk + I * (I + 1) / 2 + Jb;
                        slong ri = B->r, cj = B->c, i0 = I * BS, j0 = Jb * BS;
                        nmod_mat_t Wj, T;
                        /* Z[I] += B W[J] */
                        nmod_mat_window_init(Wj, Ym + t, j0, 0, j0 + cj, 2 * m);
                        nmod_mat_init(T, ri, 2 * m, p);
                        nmod_mat_mul(T, B, Wj);
                        for (slong a = 0; a < ri; a++) for (slong k = 0; k < 2 * m; k++)
                            nmod_mat_entry(Zp, i0 + a, k) = n_addmod(nmod_mat_entry(Zp, i0 + a, k), nmod_mat_entry(T, a, k), p);
                        nmod_mat_clear(T); nmod_mat_window_clear(Wj);
                        if (Jb < I) {   /* Z[J] += B^T W[I] = (W[I]^T B)^T */
                            nmod_mat_t WiT, T2;
                            nmod_mat_window_init(WiT, WT, 0, i0, 2 * m, i0 + ri);
                            nmod_mat_init(T2, 2 * m, cj, p);
                            nmod_mat_mul(T2, WiT, B);
                            for (slong b = 0; b < cj; b++) for (slong k = 0; k < 2 * m; k++)
                                nmod_mat_entry(Zp, j0 + b, k) = n_addmod(nmod_mat_entry(Zp, j0 + b, k), nmod_mat_entry(T2, k, b), p);
                            nmod_mat_clear(T2); nmod_mat_window_clear(WiT);
                        }
                    }
                }
            }
            for (slong i = 0; i < n; i++) for (slong k = 0; k < m; k++) {
                RZ[((i * m + k) * 2 + 0) * np + t] = n_submod(nmod_mat_entry(ZR, i, k), nmod_mat_entry(ZI, i, k + m), p);
                RZ[((i * m + k) * 2 + 1) * np + t] = n_addmod(nmod_mat_entry(ZR, i, k + m), nmod_mat_entry(ZI, i, k), p);
            }
        }
        nmod_mat_clear(ZR); nmod_mat_clear(ZI); nmod_mat_clear(WT);
    }
    #pragma omp parallel
    {
        fmpz_t z; fmpz_init(z);
        fmpz_comb_temp_t ct; fmpz_comb_temp_init(ct, comb);
        #pragma omp for schedule(static)
        for (slong idx = 0; idx < n * m * 2; idx++) {
            slong i = idx / (2 * m), k = (idx / 2) % m; int part = idx % 2;
            fmpz_multi_CRT_ui(z, RZ + idx * np, comb, ct, 1);
            arb_ptr x = part ? acb_imagref(acb_mat_entry(Z, i, k)) : acb_realref(acb_mat_entry(Z, i, k));
            arf_set_round_fmpz(arb_midref(x), z, prec, ARF_RND_NEAR);
            arf_mul_2exp_si(arb_midref(x), arb_midref(x), (kern_mm == 1 ? Eglob : Erow[i]) - BK + F[k] - BY);
            mag_zero(arb_radref(x));
        }
        fmpz_comb_temp_clear(ct); fmpz_clear(z);
    }
    free(RZ); free(F);
}

/* ------------------------------------------------------------------ orthonormalization (CGS2, parallel over rows) */
static void orth(acb_mat_t X, slong prec) {
    slong m = acb_mat_ncols(X);
    int nt = omp_get_max_threads();
    acb_ptr part = _acb_vec_init(nt * (m + 1)), s = _acb_vec_init(m + 1);
    for (slong j = 0; j < m; j++) {
        for (int pass = 0; pass < 3; pass++) {   /* passes 0,1: project out columns < j; pass 2: normalize */
            slong nk = pass < 2 ? j : 1;
            if (nk == 0) continue;
            _acb_vec_zero(part, nt * (m + 1));
            #pragma omp parallel
            {
                int tid = omp_get_thread_num(); acb_t c; acb_init(c);
                #pragma omp for schedule(static)
                for (slong i = 0; i < n; i++) {
                    if (pass < 2) for (slong k = 0; k < j; k++) { acb_conj(c, acb_mat_entry(X, i, k)); acb_addmul(part + tid * (m + 1) + k, c, acb_mat_entry(X, i, j), prec); }
                    else { acb_conj(c, acb_mat_entry(X, i, j)); acb_addmul(part + tid * (m + 1), c, acb_mat_entry(X, i, j), prec); }
                }
                acb_clear(c);
            }
            _acb_vec_zero(s, m + 1);
            for (int t = 0; t < nt; t++) for (slong k = 0; k < nk; k++) acb_add(s + k, s + k, part + t * (m + 1) + k, prec);
            if (pass < 2) {
                #pragma omp parallel for schedule(static)
                for (slong i = 0; i < n; i++) {
                    for (slong k = 0; k < j; k++) acb_submul(acb_mat_entry(X, i, j), s + k, acb_mat_entry(X, i, k), prec);
                    cmid(acb_mat_entry(X, i, j));
                }
            } else {
                arb_t nr; arb_init(nr);
                arb_rsqrt(nr, acb_realref(s), prec); amid(nr);
                #pragma omp parallel for schedule(static)
                for (slong i = 0; i < n; i++) { acb_mul_arb(acb_mat_entry(X, i, j), acb_mat_entry(X, i, j), nr, prec); cmid(acb_mat_entry(X, i, j)); }
                arb_clear(nr);
            }
        }
    }
    _acb_vec_clear(part, nt * (m + 1)); _acb_vec_clear(s, m + 1);
}

/* twist eigenvalues theta_a = exp(4 pi i k(k+1)/r), k = a/2 */
static void compute_tha(acb_ptr tha, slong prec) {
    arb_t pi, x; arb_init(pi); arb_init(x); arb_const_pi(pi, prec);
    for (slong i = 0; i < n; i++) {
        slong k = ad[i] / 2;
        arb_mul_si(x, pi, 4 * k * (k + 1), prec); arb_div_si(x, x, r, prec);
        acb_zero(tha + i); arb_set(acb_imagref(tha + i), x); acb_exp(tha + i, tha + i, prec); cmid(tha + i);
    }
    arb_clear(pi); arb_clear(x);
}

static void usage(void) {
    fprintf(stderr, "usage: ladder2 r c LR|LLR [-t threads] [-kernel mm|mmsym|acb] [-d1 D] [-d2 D] [-d3 D] [-m 12] [-ntop 8] [-tol 22] [-maxit 100] [-ck file] [-q]\n");
    exit(1);
}

int main(int argc, char **argv) {
    if (argc < 4) usage();
    r = atol(argv[1]); cc = atol(argv[2]); const char *word = argv[3];
    if (!strcmp(word, "LR")) pw = 1; else if (!strcmp(word, "LLR")) pw = 2; else usage();
    slong m = 12, ntop = 8, tolexp = 22, maxit = 100, d1 = -1, d2 = -1, d3 = -1;
    int nth = 0;
    const char *ck = NULL;
    for (int a = 4; a < argc; a++) {
        if (!strcmp(argv[a], "-d1") && a + 1 < argc) d1 = atol(argv[++a]);
        else if (!strcmp(argv[a], "-d2") && a + 1 < argc) d2 = atol(argv[++a]);
        else if (!strcmp(argv[a], "-d3") && a + 1 < argc) d3 = atol(argv[++a]);
        else if (!strcmp(argv[a], "-m") && a + 1 < argc) m = atol(argv[++a]);
        else if (!strcmp(argv[a], "-ntop") && a + 1 < argc) ntop = atol(argv[++a]);
        else if (!strcmp(argv[a], "-tol") && a + 1 < argc) tolexp = atol(argv[++a]);
        else if (!strcmp(argv[a], "-maxit") && a + 1 < argc) maxit = atol(argv[++a]);
        else if (!strcmp(argv[a], "-t") && a + 1 < argc) nth = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-ck") && a + 1 < argc) ck = argv[++a];
        else if (!strcmp(argv[a], "-kernel") && a + 1 < argc) { a++; if (!strcmp(argv[a], "acb")) kern_mm = 0; else if (!strcmp(argv[a], "mmsym")) kern_mm = 1; else if (!strcmp(argv[a], "mm")) kern_mm = 2; else usage(); }
        else if (!strcmp(argv[a], "-q")) quiet = 1;
        else usage();
    }
    if (nth > 0) omp_set_num_threads(nth);
    nth = omp_get_max_threads();
    flint_set_num_threads(1);
    if (getenv("LADDER2_BS")) BS = atol(getenv("LADDER2_BS"));
    int user1 = d1 > 0, user2 = d2 > 0, user3 = d3 > 0;
    double T0 = wall();

    ad = malloc(sizeof(slong) * (r + 1)); n = 0;
    for (slong a = 0; a < r - 2; a += 2) if (adm(a, a, cc)) ad[n++] = a;
    if (n == 0) { fprintf(stderr, "error: no admissible colors\n"); return 1; }
    if (m > n) m = n;
    if (ntop > m) ntop = m;

    /* ---- precision control (decimal digits), see README ----
       maxV  = log10 max|V_ij|, largest entry of the normalized T_b eigenvector matrix; a priori estimate
               r (0.128 - 0.3 (c/r - 1/2)^2), measured exactly in stage 1.
       stage 1 (d1): completeness residual  max_j |sum_b w_b0 w_bj J_j / N_b - delta_0j|  must be
               <= 10^(0.02 r - 30); it behaves like 10^(L1 - d1) with L1 ~ (3.52 - 1.13 c/r) maxV.
       stage 2 (d2): two-sided mismatch of the K recurrence gives the absolute error of K; required
               maxV - log10(error) >= qreq = 0.06 r + 30; error ~ 10^(L2 - d2), L2 ~ (1.67 - 1.43 c/r) maxV
               (+ (0.95 - 0.8 c/r) maxV for LLR).
       stage 3 (d3): needs log10 max|K| - log10|lambda_1| + a r + 30, a = 0.068 (LR), 0.098 (LLR); K is first
               stored with d3 = log10 max|K| + 35 (log10 max|K| = maxV for LR, ~1.4 maxV for LLR), lambda_1 is
               read off after 3 iterations, and K is rebuilt at higher precision if needed.
       Each stage is re-run with more digits if its check fails. */
    double xr = (double) cc / r, mvest = r * (0.128 - 0.3 * (xr - 0.5) * (xr - 0.5));
    if (mvest < 0.02 * r) mvest = 0.02 * r;
    double qreq = 0.06 * r + 30, ctarget = 0.02 * r - 30;
    if (!user1) d1 = (slong)((3.62 - 1.13 * xr) * mvest - ctarget + 10);

    acb_mat_t X; acb_mat_init(X, n, m);
    slong ckd3 = -1;
    int loaded = 0; slong it0 = 0;
    if (ck) {
        FILE *fp = fopen(ck, "r");
        if (fp) {
            long rr, c2, nn, mm, itf, d3f = -1; char wd[8], line[256];
            if (fgets(line, sizeof line, fp) && sscanf(line, "%ld %ld %7s %ld %ld %ld %ld", &rr, &c2, wd, &nn, &mm, &itf, &d3f) >= 6
                && rr == r && c2 == cc && !strcmp(wd, word) && nn == n && mm == m) {
                size_t cap = 1 << 22; char *buf = malloc(cap);
                int ok = 1;
                for (slong i = 0; i < n && ok; i++) for (slong j = 0; j < m && ok; j++) for (int part = 0; part < 2; part++) {
                    if (!fgets(buf, cap, fp)) { ok = 0; break; }
                    buf[strcspn(buf, "\n")] = 0;
                    if (arb_load_str(part ? acb_imagref(acb_mat_entry(X,i,j)) : acb_realref(acb_mat_entry(X,i,j)), buf)) { ok = 0; break; }
                }
                free(buf);
                if (ok) { loaded = 1; it0 = itf; ckd3 = d3f; fprintf(stderr, "# resumed from checkpoint %s (iteration %ld)\n", ck, itf); }
            }
            fclose(fp);
        }
    }
    /* a checkpoint written after the stage-3 precision check carries the final d3: reuse it */
    int d3fromck = loaded && ckd3 > 0 && !user3;
    if (d3fromck) d3 = ckd3;
    double lvmax = 0, lkmax = 0, qK = 0, qabs = 0, tA = 0, tB = 0;
    acb_ptr R0 = _acb_vec_init(n), R1 = _acb_vec_init(n);
    int d2set = 0;
    for (int outer = 0; ; outer++) {
        /* ---- stage 1 ---- */
        for (int attempt = 0; ; attempt++) {
            double ta = wall();
            setup(DIG2BITS(d1));
            phaseA(R0, R1, &lvmax, DIG2BITS(d1));
            tA += wall() - ta;
            if (!quiet) fprintf(stderr, "# stage 1: d1 %ld | log10 max|V| %.1f | completeness residual 1e%.1f (need <= 1e%.0f) | %.1fs\n",
                                (long)d1, lvmax, lcomp, ctarget, wall() - ta);
            if (lcomp <= ctarget || user1 || attempt >= 3) {
                if (lcomp > ctarget) fprintf(stderr, "# WARNING: stage 1 completeness residual too large; results may be inaccurate\n");
                break;
            }
            d1 += (slong)(lcomp - ctarget) + 15;
            setup_clear();
        }
        /* ---- stage 2 ---- */
        double lkest = (pw == 1 ? 1.02 : 1.4) * lvmax;
        if (!user3 && !d3fromck && outer == 0) d3 = (slong)(lkest + 35);
        /* recurrence loss ~ (1.67 - 1.43 c/r) maxV (LR); LLR adds ~ (0.95 - 0.8 c/r) maxV */
        if (!user2 && !d2set) {
            d2 = (slong)(qreq + (1.77 - 1.43 * xr) * lvmax + (pw == 2 ? (1.0 - 0.8 * xr) * lvmax : 0) + 20);
            if (kern_mm == 1) d2 = (slong)(d3 + (1.77 - 1.43 * xr) * lvmax + (pw == 2 ? (1.0 - 0.8 * xr) * lvmax : 0) + 10);
            d2set = 1;
        }
        if (outer == 0) Eglob = (slong) ceil(lkest / LOG10_2) + 4;
        int d1limited = 0; double prevq = -1e300; slong prevraise = 0;
        for (int attempt = 0; ; attempt++) {
            if (d2 < d3) d2 = d3;
            if (!quiet && kern_mm) {
                slong npr = nprimes_for(DIG2BITS(d3));
                double ent = kern_mm == 1 ? 0.5 * n * n + 0.5 * n * BS : (double) n * n;
                fprintf(stderr, "# K storage (%s kernel): %ld primes, %.2f GB\n", kern_mm == 1 ? "mmsym" : "mm", (long)npr, 2.0 * npr * ent * 8 / 1e9);
            }
            double tb = wall();
            phaseB(R0, R1, DIG2BITS(d2), DIG2BITS(d3), &lkmax, &qK);
            tB += wall() - tb;
            qabs = lvmax - (lkmax - qK);
            /* required accuracy of K (digits relative to max|K|): absolute error <= 10^(maxV - qreq) for the full
               kernels; symmetric storage mixes the two recurrence sweeps and needs ~ d3 */
            double qneed = kern_mm == 1 ? d3 - 3 : qreq + lkmax - lvmax;
            if (!quiet) fprintf(stderr, "# stage 2: d2 %ld d3 %ld | log10 max|K| %.1f | K accuracy %.1f digits (need %.1f) | %.1fs\n",
                                (long)d2, (long)d3, lkmax, qK, qneed, wall() - tb);
            int redo = 0;
            if (attempt < 4) {
                if (!user3 && lkmax + 33 > d3) { d3 = (slong)(lkmax + 35); redo = 1; if (kern_mm == 1) qneed = d3 - 3; }
                if (!user2 && qK < qneed) {
                    if (prevraise > 0 && qK - prevq < 0.5 * prevraise && !user1) { d1limited = 1; d1 += (slong)(qneed - qK) + 20; break; }
                    prevraise = (slong)(qneed - qK) + 20; d2 += prevraise; redo = 1;
                } else prevraise = 0;
                if (kern_mm == 1 && Emaxseen > Eglob + KSLACK - 1) { Eglob = Emaxseen + 2; redo = 1; }   /* fixed-point overflow */
            }
            prevq = qK;
            if (!redo) {
                if (qK < qneed) fprintf(stderr, "# WARNING: K accuracy too low; results may be inaccurate\n");
                break;
            }
            kstore_free();
        }
        if (d1limited && outer < 2) {
            if (!quiet) fprintf(stderr, "# stage 2 limited by stage 1 precision: d1 -> %ld\n", (long)d1);
            kstore_free(); setup_clear(); continue;
        }
        break;
    }
    double T2 = wall();

    /* ---------------- subspace iteration with A' = K Theta ---------------- */
    slong prec = DIG2BITS(d3);
    acb_ptr tha = _acb_vec_init(n);
    int d3checked = user3 || d3fromck;
    slong minit = 0;   /* after a rebuild of K, iterate at least this far before accepting convergence */
    double arate = pw == 1 ? 0.068 : 0.098;
    compute_tha(tha, prec);
    acb_mat_t Y, Z, H; acb_mat_init(Y, n, m); acb_mat_init(Z, n, m); acb_mat_init(H, m, m);
    if (!loaded) {
        flint_rand_t st; flint_randinit(st);
        for (slong i = 0; i < n; i++) for (slong j = 0; j < m; j++) {
            arb_set_d(acb_realref(acb_mat_entry(X,i,j)), (double)n_randint(st, 1000000) / 1e6 - 0.5);
            arb_set_d(acb_imagref(acb_mat_entry(X,i,j)), (double)n_randint(st, 1000000) / 1e6 - 0.5);
        }
        flint_randclear(st);
    }
    orth(X, prec);
    acb_ptr ev = _acb_vec_init(m);
    arb_ptr logs = _arb_vec_init(m), plogs = _arb_vec_init(m);
    arb_t tol, dif, tmp; arb_init(tol); arb_init(dif); arb_init(tmp);
    arb_set_si(tol, 10); arb_pow_ui(tol, tol, tolexp, prec); arb_inv(tol, tol, prec);
    slong it;
    int conv = 0;
    for (it = it0; it < it0 + maxit; it++) {
        double ta = wall();
        #pragma omp parallel for schedule(static)
        for (slong i = 0; i < n; i++) for (slong j = 0; j < m; j++) { acb_mul(acb_mat_entry(Y,i,j), tha + i, acb_mat_entry(X,i,j), prec); cmid(acb_mat_entry(Y,i,j)); }
        double tm = wall();
        kmul(Z, Y, prec);
        tm = wall() - tm;
        /* Rayleigh-Ritz: H = X^H Z */
        #pragma omp parallel
        {
            acb_ptr a = _acb_vec_init(n), b = _acb_vec_init(n);
            #pragma omp for schedule(dynamic, 1)
            for (slong q = 0; q < m * m; q++) {
                slong ia = q / m, ib = q % m;
                for (slong i = 0; i < n; i++) { acb_conj(a + i, acb_mat_entry(X, i, ia)); acb_set(b + i, acb_mat_entry(Z, i, ib)); }
                acb_dot(acb_mat_entry(H, ia, ib), NULL, 0, a, 1, b, 1, n, prec);
                cmid(acb_mat_entry(H, ia, ib));
            }
            _acb_vec_clear(a, n); _acb_vec_clear(b, n);
        }
        acb_mat_approx_eig_qr(ev, NULL, NULL, H, NULL, 0, prec);
        for (slong k = 0; k < m; k++) { acb_abs(logs + k, ev + k, prec); arb_log(logs + k, logs + k, prec); amid(logs + k); }
        for (slong a = 0; a < m; a++) for (slong b = a + 1; b < m; b++) if (arf_cmp(arb_midref(logs + b), arb_midref(logs + a)) > 0) arb_swap(logs + a, logs + b);
        int done = 0;
        if (it > it0) {
            arb_zero(dif);
            for (slong k = 0; k < ntop; k++) { arb_sub(tmp, logs + k, plogs + k, prec); arb_abs(tmp, tmp); if (arf_cmp(arb_midref(tmp), arb_midref(dif)) > 0) arb_set(dif, tmp); }
            if (arf_cmp(arb_midref(dif), arb_midref(tol)) < 0) done = 1;
        }
        for (slong k = 0; k < m; k++) arb_set(plogs + k, logs + k);
        if (!quiet) fprintf(stderr, "it %ld dif %.3e  (%.2fs, product %.2fs)\n", (long)it, it > it0 ? arf_get_d(arb_midref(dif), ARF_RND_NEAR) : 0.0, wall() - ta, tm);
        if (!d3checked && it >= it0 + 2) {
            /* iteration stage needs d3 >= log10 max|K| - log10|lambda_1| + a r + 30 (a = 0.068 LR, 0.098 LLR) */
            d3checked = 1;
            double l1 = arf_get_d(arb_midref(logs + 0), ARF_RND_NEAR) / log(10.0);
            slong need = (slong) ceil(lkmax - (l1 - 2) + arate * r + 30);
            /* effective iteration digits: the global fixed-point exponent of the mm kernel may exceed max|K| */
            slong d3eff = d3 - (kern_mm == 1 && Eglob > Emaxseen ? (slong) ceil((Eglob - Emaxseen) * LOG10_2) : 0);
            if (!quiet) fprintf(stderr, "# log10|lambda_1| ~ %.1f: stage 3 needs d3 >= %ld (have %ld)\n", l1, (long)need, (long)d3eff);
            if (need > d3eff) {
                double tb = wall();
                slong d3old = d3;
                if (need > d3) d3 = need;
                if (kern_mm == 1 && !user2) d2 += d3 - d3old;
                if (d2 < d3) d2 = d3;
                Eglob = Emaxseen + 2;
                kstore_free();
                phaseB(R0, R1, DIG2BITS(d2), DIG2BITS(d3), &lkmax, &qK);
                if (kern_mm == 1 && qK < d3 - 3 && !user2) {
                    d2 += (slong)(d3 - 3 - qK) + 20; kstore_free();
                    phaseB(R0, R1, DIG2BITS(d2), DIG2BITS(d3), &lkmax, &qK);
                }
                if (kern_mm == 1 && qK < d3 - 3) fprintf(stderr, "# WARNING: K accuracy %.1f < %ld; results may be inaccurate\n", qK, (long)(d3 - 3));
                tB += wall() - tb;
                if (!quiet) fprintf(stderr, "# rebuilt K at d3 %ld (%.1fs)\n", (long)d3, wall() - tb);
                prec = DIG2BITS(d3);
                compute_tha(tha, prec);
                minit = it + 7;
                acb_mat_swap(X, Z); orth(X, prec);
                continue;
            }
        }
        if (done && it >= minit) { conv = 1; break; }
        acb_mat_swap(X, Z); orth(X, prec);
        if (ck) {
            char tmpn[4096]; snprintf(tmpn, sizeof tmpn, "%s.tmp", ck);
            FILE *fp = fopen(tmpn, "w");
            if (fp) {
                fprintf(fp, "%ld %ld %s %ld %ld %ld %ld\n", (long)r, (long)cc, word, (long)n, (long)m, (long)(it + 1), (long)(d3checked ? d3 : -1));
                for (slong i = 0; i < n; i++) for (slong j = 0; j < m; j++) for (int part = 0; part < 2; part++) {
                    char *s = arb_dump_str(part ? acb_imagref(acb_mat_entry(X,i,j)) : acb_realref(acb_mat_entry(X,i,j)));
                    fprintf(fp, "%s\n", s); flint_free(s);
                }
                fclose(fp); rename(tmpn, ck);
            }
        }
    }
    if (!conv) { fprintf(stderr, "# WARNING: not converged after %ld iterations\n", (long)maxit); if (it > it0) it--; }
    double T4 = wall();
    double thv = fabs(2 * M_PI - 4 * M_PI * (cc + 1) / (double)r);
    printf("%ld %ld %ld %.10f %s %ld", (long)r, (long)cc, (long)n, thv, word, (long)it);
    for (slong k = 0; k < ntop; k++) { char *s = arb_get_str(logs + k, 30, ARB_STR_NO_RADIUS); printf(" %s", s); flint_free(s); }
    printf("\n"); fflush(stdout);
    fprintf(stderr, "# r %ld c %ld %s n %ld threads %d kernel %s | digits d1 %ld d2 %ld d3 %ld | maxV %.1f maxK %.1f Kacc %.1f comp %.1f | eigvec %.1fs recur %.1fs iter %.1fs (%ld its) total %.1fs | peak RSS %.2f GB\n",
            (long)r, (long)cc, word, (long)n, nth, kern_mm == 1 ? "mmsym" : kern_mm == 2 ? "mm" : "acb", (long)d1, (long)d2, (long)d3, lvmax, lkmax, qabs, lcomp, tA, tB, T4 - T2, (long)(it - it0 + 1), T4 - T0, maxrss_gb());
    return 0;
}
