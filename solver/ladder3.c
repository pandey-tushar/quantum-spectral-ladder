/* ladder3: same computation and output as ladder2 (see ladder2.c, README_ladder2.md), with bounded memory.

   Differences from ladder2 (mm kernel):
   - K is stored with a uniform absolute fixed-point tolerance 2^Etol, Etol = log2 max|K| - BK, in tiles of
     RB rows x CB columns.  A tile whose entries need b bits is stored as residues modulo only as many 59-bit
     primes as a b-bit by b-bit product needs (NL precision levels), and Y is truncated to the same level.
     Tiles below the tolerance are dropped.  This cuts storage and product work (K decays away from its
     central peak), and needs one CRT per row, column and level.
   - With -scratch DIR, K is written to DIR/K.bin and streamed through a bounded buffer each iteration
     (prefetched by an I/O thread); -mem GB bounds the RAM used for K (resident part + stream buffers).
     The build is checkpointed (stage 1 rows, stage 2 recurrence state, finished K) and resumes; the
     iteration checkpoint defaults to DIR/iter.ck.
   - Stage 1 accumulates into shared, chunk-locked sums: 2 vectors per thread instead of 9.

   usage: ladder3 r c LR|LLR [-t T] [-scratch DIR] [-mem GB] [-nocache] [-ckint SEC] [-levels NL] [-rb RB]
                             [-cb CB] [-d1/-d2/-d3 D] [-m 12] [-ntop 8] [-tol 22] [-maxit 100] [-ck FILE] [-q] */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/nmod_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>

#define LOG10_2 0.30102999566398120
#define DIG2BITS(d) ((slong)((d) * 3.3219280948873623) + 32)
#define NLMAX 32
#define KSLACK 48       /* bits of headroom of the top level above BK (max|K| estimate is a sampled lower bound) */
#define YM 8            /* extra bits of Y beyond the tile's bits */

static slong r, cc, n, pw;
static slong *ad;
static int quiet = 0;

static double wall(void) { struct timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + 1e-6 * tv.tv_usec; }
static double maxrss_gb(void) { struct rusage u; getrusage(RUSAGE_SELF, &u); return u.ru_maxrss / 1048576.0; }
static inline void amid(arb_t x) { mag_zero(arb_radref(x)); }
static inline void cmid(acb_t z) { mag_zero(arb_radref(acb_realref(z))); mag_zero(arb_radref(acb_imagref(z))); }
static inline slong lexp(const arb_t x) { return arf_is_zero(arb_midref(x)) ? WORD_MIN : arf_abs_bound_lt_2exp_si(arb_midref(x)); }
static inline slong lexpc(const acb_t z) { slong a = lexp(acb_realref(z)), b = lexp(acb_imagref(z)); return a > b ? a : b; }

/* ------------------------------------------------------------------ 6j symbols (as ladder2) */
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

/* ------------------------------------------------------------------ problem data at precision p1 (as ladder2) */
static arb_ptr al, u, l, lam;
static int *J;
static acb_ptr thp;
static slong hrow;   /* row of the largest |alpha| (meeting row of stage 2) */

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
                isim[i] = 1;
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
    hrow = 0;
    for (slong i = 1; i < n; i++) if (arf_cmpabs(arb_midref(al + i), arb_midref(al + hrow)) > 0) hrow = i;
    if (hrow > n - 2) hrow = n - 2;
    acb_clear(t); acb_clear(uu); acb_clear(ac); acb_clear(be);
    arb_clear(pi); arb_clear(x); arb_clear(s1);
    _arb_vec_clear(qi, nqi); _arb_vec_clear(fac, nqi); free(isim);
}
static void setup_clear(void) {
    _arb_vec_clear(al, n); _arb_vec_clear(u, n); _arb_vec_clear(l, n); _arb_vec_clear(lam, n);
    _acb_vec_clear(thp, n); free(J);
}

/* ------------------------------------------------------------------ text dump helpers */
static int dump_acbvec(FILE *fp, acb_srcptr v, slong len) {
    for (slong j = 0; j < len; j++) for (int part = 0; part < 2; part++) {
        char *s = arb_dump_str(part ? acb_imagref(v + j) : acb_realref(v + j));
        if (fprintf(fp, "%s\n", s) < 0) { flint_free(s); return 1; }
        flint_free(s);
    }
    return 0;
}
static int load_acbvec(FILE *fp, acb_ptr v, slong len) {
    size_t cap = 1 << 22; char *buf = malloc(cap); int ok = 1;
    for (slong j = 0; j < len && ok; j++) for (int part = 0; part < 2 && ok; part++) {
        if (!fgets(buf, cap, fp)) { ok = 0; break; }
        buf[strcspn(buf, "\n")] = 0;
        if (arb_load_str(part ? acb_imagref(v + j) : acb_realref(v + j), buf)) ok = 0;
    }
    free(buf);
    return !ok;
}
static void atomic_close(FILE *fp, const char *tmp, const char *fin) {
    fflush(fp); fsync(fileno(fp)); fclose(fp); rename(tmp, fin);
}

/* ------------------------------------------------------------------ stage 1: eigenvectors -> rows 0, n-1 of K */
static double lcomp;
static slong Esample;          /* upper bound of log2 of the largest sampled |K_{h,j}| */
#define NHS 65
static void phaseA(acb_ptr R0, acb_ptr R1, double *lvmax, slong prec) {
    double lvm = -1e300;
    arb_ptr E0 = _arb_vec_init(n), E1 = _arb_vec_init(n);
    slong hs_j[NHS], nhs = 0, sp = n / 64 > 1 ? n / 64 : 1;
    for (slong q = -(NHS / 2); q <= NHS / 2; q++) { slong j = hrow + q * sp; if (j >= 0 && j < n) hs_j[nhs++] = j; }
    acb_ptr HS = _acb_vec_init(nhs);
    _acb_vec_zero(R0, n); _acb_vec_zero(R1, n);
    int nt = omp_get_max_threads();
    slong nch = 4 * nt; if (nch > n) nch = n;
    omp_lock_t *lk = malloc(sizeof(omp_lock_t) * nch);
    for (slong q = 0; q < nch; q++) omp_init_lock(lk + q);
    #pragma omp parallel
    {
        arb_ptr f = _arb_vec_init(n), g = _arb_vec_init(n);
        acb_ptr hs = _acb_vec_init(nhs);
        char *done = malloc(nch);
        arb_t t0, t1, t, gam, N, one; arb_init(t0); arb_init(t1); arb_init(t); arb_init(gam); arb_init(N); arb_init(one); arb_one(one);
        acb_t c0, c1, ch; acb_init(c0); acb_init(c1); acb_init(ch);
        mag_t mb, mc; mag_init(mb); mag_init(mc);
        double lv = -1e300;
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic, 1)
        for (slong b = 0; b < n; b++) {
            for (slong i = 0; i < n - 1; i++) {
                arb_sub(t, lam + b, al + i, prec);
                if (i > 0) arb_submul(t, l + i - 1, f + i - 1, prec);
                arb_div(f + i, u + i, t, prec); amid(f + i);
            }
            for (slong i = n - 1; i > 0; i--) {
                arb_sub(t, lam + b, al + i, prec);
                if (i < n - 1) arb_submul(t, u + i, g + i + 1, prec);
                arb_div(g + i, l + i - 1, t, prec); amid(g + i);
            }
            slong kb = 0;
            for (slong k = 0; k < n; k++) {
                arb_sub(gam, al + k, lam + b, prec);
                if (k > 0) arb_addmul(gam, l + k - 1, f + k - 1, prec);
                if (k < n - 1) arb_addmul(gam, u + k, g + k + 1, prec);
                arb_get_mag(mc, gam);
                if (k == 0 || mag_cmp(mc, mb) < 0) { mag_set(mb, mc); kb = k; }
            }
            /* eigenvector in place: w_i = f_i (i < kb), 1 (i = kb), g_i (i > kb) */
            for (slong i = kb - 2; i >= 0; i--) { arb_mul(f + i, f + i, f + i + 1, prec); amid(f + i); }
            for (slong i = kb + 2; i < n; i++) { arb_mul(g + i, g + i, g + i - 1, prec); amid(g + i); }
            arb_one(f + kb);
            #define W(i) ((i) <= kb ? f + (i) : g + (i))
            arb_zero(N);
            for (slong i = 0; i < n; i++) { if (J[i] > 0) arb_addmul(N, W(i), W(i), prec); else arb_submul(N, W(i), W(i), prec); }
            amid(N);
            {
                slong em = WORD_MIN;
                for (slong i = 0; i < n; i++) { slong e = lexp(W(i)); if (e > em) em = e; }
                double e2 = (em - 0.5 * lexp(N)) * LOG10_2;
                if (e2 > lv) lv = e2;
            }
            arb_div(t0, W(0), N, prec); acb_mul_arb(c0, thp + b, t0, prec);
            arb_div(t1, W(n - 1), N, prec); acb_mul_arb(c1, thp + b, t1, prec);
            arb_div(t, W(hrow), N, prec); acb_mul_arb(ch, thp + b, t, prec);
            for (slong q = 0; q < nhs; q++) {
                slong j = hs_j[q];
                if (J[j] > 0) acb_addmul_arb(hs + q, ch, W(j), prec); else acb_submul_arb(hs + q, ch, W(j), prec);
                cmid(hs + q);
            }
            /* accumulate into the shared sums, one locked chunk of j at a time */
            memset(done, 0, nch);
            slong left = nch, q = (tid * 7 + b) % nch;
            while (left > 0) {
                if (!done[q] && omp_test_lock(lk + q)) {
                    slong j0 = q * n / nch, j1 = (q + 1) * n / nch;
                    for (slong j = j0; j < j1; j++) {
                        if (J[j] > 0) { arb_addmul(E0 + j, t0, W(j), prec); arb_addmul(E1 + j, t1, W(j), prec);
                                        acb_addmul_arb(R0 + j, c0, W(j), prec); acb_addmul_arb(R1 + j, c1, W(j), prec); }
                        else { arb_submul(E0 + j, t0, W(j), prec); arb_submul(E1 + j, t1, W(j), prec);
                               acb_submul_arb(R0 + j, c0, W(j), prec); acb_submul_arb(R1 + j, c1, W(j), prec); }
                        amid(E0 + j); amid(E1 + j); cmid(R0 + j); cmid(R1 + j);
                    }
                    omp_unset_lock(lk + q);
                    done[q] = 1; left--;
                }
                q = q + 1 == nch ? 0 : q + 1;
            }
            #undef W
        }
        #pragma omp critical
        {
            if (lv > lvm) lvm = lv;
            for (slong q = 0; q < nhs; q++) { acb_add(HS + q, HS + q, hs + q, prec); cmid(HS + q); }
        }
        free(done);
        _arb_vec_clear(f, n); _arb_vec_clear(g, n); _acb_vec_clear(hs, nhs);
        arb_clear(t0); arb_clear(t1); arb_clear(t); arb_clear(gam); arb_clear(N); arb_clear(one);
        acb_clear(c0); acb_clear(c1); acb_clear(ch); mag_clear(mb); mag_clear(mc);
    }
    for (slong q = 0; q < nch; q++) omp_destroy_lock(lk + q);
    free(lk);
    *lvmax = lvm;
    arb_sub_ui(E0 + 0, E0 + 0, 1, prec); arb_sub_ui(E1 + n - 1, E1 + n - 1, 1, prec);
    slong em = WORD_MIN;
    for (slong j = 0; j < n; j++) { slong e = lexp(E0 + j); if (e > em) em = e; e = lexp(E1 + j); if (e > em) em = e; }
    lcomp = em == WORD_MIN ? -1e300 : em * LOG10_2;
    Esample = WORD_MIN;
    for (slong q = 0; q < nhs; q++) { slong e = lexpc(HS + q); if (e > Esample) Esample = e; }
    _arb_vec_clear(E0, n); _arb_vec_clear(E1, n); _acb_vec_clear(HS, nhs);
}

/* ------------------------------------------------------------------ K store (tiled, levelled residues; RAM or disk) */
static slong RB = 32, CB = 32;
static int NL = 6;
static double segmb = 64;                  /* target segment size (MB) */
static const char *scratch = NULL;
static double memgb = -1;
static int nocache = 0;
static double ckint = 600;

static slong BK, Eref, Etol, Emaxseen, npmax, lgn;
static int overflow;
static mp_limb_t *primes;
static fmpz_comb_struct combs[NLMAX];
static slong npl[NLMAX], bl[NLMAX];
static int kalloc = 0;

typedef struct {
    int lev; slong B, r0, nr;  /* level, row block, its first row and row count */
    slong w, ntile; slong *tiles;
    off_t off; size_t bytes;   /* in K.bin */
    mp_limb_t *mem;            /* resident data or NULL */
    mp_limb_t *RZ;             /* product residues (during kmul) */
} seg_t;
static seg_t *segs; static slong nseg, segcap;
static slong nblk, *brow0, *bnr, *blk_of, *bfill;
static fmpz *rowbuf;                       /* RB x n x 2 fixed-point integers of the current block(s) */
static slong rb_dir0blk;                   /* blocks [0, rb_dir0blk) belong to sweep 0 */
static int kfd = -1; static off_t kend = 0;
static size_t res_bytes = 0, res_cap = (size_t)-1, tot_bytes = 0, bufsz = 0;
static double ref_bytes = 0;               /* ladder2 mm storage for the same d3 */
static int plan_ok = 0;
static char sigbuf[512];

static char *spath(const char *name) { static char p[4][4096]; static int k = 0; k = (k + 1) & 3; snprintf(p[k], 4096, "%s/%s", scratch, name); return p[k]; }

static slong P2bits = 0;
static void store_init(slong p3) {
    BK = p3;
    lgn = 1; while ((WORD(1) << lgn) < n) lgn++;
    npmax = (2 * (BK + KSLACK) + YM + lgn + 4 + 57) / 58;
    if (NL > npmax) NL = (int) npmax;
    primes = malloc(sizeof(mp_limb_t) * npmax);
    mp_limb_t p = UWORD(1) << 58;
    for (slong t = 0; t < npmax; t++) { p = n_nextprime(p, 1); primes[t] = p; }
    for (int v = 0; v < NL; v++) {
        npl[v] = ((v + 1) * npmax + NL - 1) / NL;
        bl[v] = (58 * npl[v] - lgn - 4 - YM) / 2;
        fmpz_comb_init(combs + v, primes, npl[v]);
    }
    Etol = Eref - BK;
    overflow = 0; Emaxseen = WORD_MIN;
    nseg = 0; segcap = 1024; segs = malloc(sizeof(seg_t) * segcap);
    res_bytes = 0; tot_bytes = 0; plan_ok = 0;
    {   slong np0 = (2 * p3 + lgn + 6 + 57) / 58; ref_bytes = 2.0 * np0 * (double) n * n * 8; }
    if (scratch) {
        double m = memgb > 0 ? memgb : 1.5;
        bufsz = (size_t)(m * 1e9 / 8);
        if (bufsz > ((size_t)256 << 20)) bufsz = (size_t)256 << 20;
        if (bufsz < ((size_t)16 << 20)) bufsz = (size_t)16 << 20;
        double rc = m * 1e9 - 2.0 * bufsz;
        res_cap = rc > 0 ? (size_t) rc : 0;
        if (kfd < 0) { kfd = open(spath("K.bin"), O_RDWR | O_CREAT, 0644); if (kfd < 0) { perror("K.bin"); exit(3); } }
    } else { res_cap = memgb > 0 ? (size_t)(memgb * 1e9) : (size_t)-1; bufsz = 0; }
    snprintf(sigbuf, sizeof sigbuf, "%ld %ld %s %ld P2 %ld BK %ld Eref %ld NL %d RB %ld CB %ld seg %.0f np %ld",
             (long)r, (long)cc, pw == 1 ? "LR" : "LLR", (long)n, (long)P2bits, (long)BK, (long)Eref, NL, (long)RB, (long)CB, segmb, (long)npmax);
    kalloc = 1;
}
static void plan_free(void);
static void store_free(void) {
    if (!kalloc) return;
    plan_free();
    for (slong s = 0; s < nseg; s++) { free(segs[s].tiles); free(segs[s].mem); free(segs[s].RZ); }
    free(segs); segs = NULL; nseg = 0;
    for (int v = 0; v < NL; v++) fmpz_comb_clear(combs + v);
    free(primes);
    kalloc = 0; plan_ok = 0;
}

static void seg_commit(seg_t *S) {
    if (nseg == segcap) { segcap *= 2; segs = realloc(segs, sizeof(seg_t) * segcap); }
    S->RZ = NULL;
    tot_bytes += S->bytes;
    if (scratch) {
        S->off = kend;
        size_t done = 0;
        while (done < S->bytes) {
            ssize_t w = pwrite(kfd, (char *) S->mem + done, S->bytes - done, kend + done);
            if (w <= 0) { if (errno == EINTR) continue; perror("write K.bin"); exit(3); }
            done += w;
        }
        kend += S->bytes;
    } else S->off = -1;
    if (res_bytes + S->bytes <= res_cap) res_bytes += S->bytes;
    else { if (!scratch) { fprintf(stderr, "error: K exceeds -mem and no -scratch given\n"); exit(3); } free(S->mem); S->mem = NULL; }
    segs[nseg++] = *S;
}

/* block B complete in rowbuf: choose tile levels, write segments */
static void finalize_block(slong B) {
    slong nr = bnr[B], nct = (n + CB - 1) / CB, base = (B < rb_dir0blk) ? 0 : RB;   /* two rowbuf halves: sweeps 0 and 1 */
    fmpz *rb = rowbuf + base * n * 2;
    slong *bt = malloc(sizeof(slong) * nct);
    #pragma omp parallel for schedule(dynamic, 1)
    for (slong T = 0; T < nct; T++) {
        slong mb = 0, j1 = (T + 1) * CB < n ? (T + 1) * CB : n;
        for (slong a = 0; a < nr; a++) for (slong j = T * CB; j < j1; j++) for (int part = 0; part < 2; part++) {
            slong b = fmpz_bits(rb + (a * n + j) * 2 + part); if (b > mb) mb = b;
        }
        bt[T] = mb;
    }
    int *lv = malloc(sizeof(int) * nct);
    for (slong T = 0; T < nct; T++) {
        if (bt[T] == 0) { lv[T] = -1; continue; }
        int v = 0; while (v < NL - 1 && bl[v] < bt[T]) v++;
        if (bl[v] < bt[T]) overflow = 1;
        lv[T] = v;
    }
    /* merge a thin level into the next level present in this block when the extra product work is
       cheaper than the CRT it saves (measured: one CRT ~ 0.0007 np^2 + 0.5 us, one product word ~ 0.9 ns) */
    {
        slong wl[NLMAX] = {0};
        for (slong T = 0; T < nct; T++) if (lv[T] >= 0) wl[lv[T]] += (T + 1) * CB < n ? CB : n - T * CB;
        for (int v = 0; v < NL - 1; v++) {
            if (wl[v] == 0) continue;
            int v2 = v + 1; while (v2 < NL && wl[v2] == 0) v2++;
            if (v2 == NL) break;
            double extra = 0.9e-3 * wl[v] * (npl[v2] - npl[v]) * 2, saved = 0.0007 * npl[v] * npl[v] + 0.5;   /* us per row and column of Y */
            if (extra < saved) { for (slong T = 0; T < nct; T++) if (lv[T] == v) lv[T] = v2; wl[v2] += wl[v]; wl[v] = 0; }
        }
    }
    slong *tl = malloc(sizeof(slong) * nct);
    for (int v = 0; v < NL; v++) {
        slong np = npl[v];
        slong gw = (slong)(segmb * 1048576.0 / (16.0 * nr * np)); gw = gw / CB * CB; if (gw < CB) gw = CB;
        slong T = 0;
        while (T < nct) {
            slong nt = 0, w = 0;
            for (; T < nct; T++) {
                if (lv[T] != v) continue;
                slong tw = (T + 1) * CB < n ? CB : n - T * CB;
                if (w + tw > gw && nt > 0) break;
                tl[nt++] = T; w += tw;
            }
            if (nt == 0) break;
            seg_t S; S.lev = v; S.B = B; S.r0 = brow0[B]; S.nr = nr; S.w = w; S.ntile = nt;
            S.tiles = malloc(sizeof(slong) * nt); memcpy(S.tiles, tl, sizeof(slong) * nt);
            S.bytes = (size_t) np * 2 * nr * w * sizeof(mp_limb_t);
            S.mem = malloc(S.bytes);
            slong *col = malloc(sizeof(slong) * w);
            for (slong q = 0, c = 0; q < nt; q++) { slong j1 = (tl[q] + 1) * CB < n ? (tl[q] + 1) * CB : n; for (slong j = tl[q] * CB; j < j1; j++) col[c++] = j; }
            mp_limb_t *D = S.mem;
            #pragma omp parallel
            {
                fmpz_comb_temp_t ct; fmpz_comb_temp_init(ct, combs + v);
                mp_limb_t *res = malloc(sizeof(mp_limb_t) * np);
                #pragma omp for schedule(static)
                for (slong e = 0; e < nr * w; e++) {
                    slong a = e / w, c = e % w;
                    for (int part = 0; part < 2; part++) {
                        fmpz_multi_mod_ui(res, rb + (a * n + col[c]) * 2 + part, combs + v, ct);
                        for (slong t = 0; t < np; t++) D[((t * 2 + part) * nr + a) * w + c] = res[t];
                    }
                }
                free(res); fmpz_comb_temp_clear(ct);
            }
            free(col);
            seg_commit(&S);
        }
    }
    free(tl); free(lv); free(bt);
    #pragma omp parallel for schedule(static)
    for (slong e = 0; e < nr * n * 2; e++) fmpz_zero(rb + e);
}

/* convert row i (precision >= p3) to fixed point at 2^Etol into rowbuf; finalize its block when complete */
static int blockdone_flag = 0;
static void rowstore(slong i, acb_srcptr row, slong *rowexp) {
    slong em = WORD_MIN;
    for (slong j = 0; j < n; j++) { slong e = lexpc(row + j); if (e > em) em = e; }
    *rowexp = em;
    slong B = blk_of[i], a = i - brow0[B], base = (B < rb_dir0blk) ? 0 : RB;
    fmpz *rb = rowbuf + (base + a) * n * 2;
    #pragma omp parallel for schedule(static)
    for (slong j = 0; j < n; j++) {
        arf_get_fmpz_fixed_si(rb + j * 2 + 0, arb_midref(acb_realref(row + j)), Etol);
        arf_get_fmpz_fixed_si(rb + j * 2 + 1, arb_midref(acb_imagref(row + j)), Etol);
    }
    if (++bfill[B] == bnr[B]) { finalize_block(B); blockdone_flag = 1; }
}

static void blocks_init(void) {
    slong h = hrow, cap = n / RB + 4;
    brow0 = malloc(sizeof(slong) * cap); bnr = malloc(sizeof(slong) * cap); bfill = calloc(cap, sizeof(slong));
    blk_of = malloc(sizeof(slong) * n);
    nblk = 0;
    for (slong i = 0; i <= h; i += RB) { brow0[nblk] = i; bnr[nblk] = (i + RB <= h + 1) ? RB : h + 1 - i; nblk++; }
    rb_dir0blk = nblk;
    for (slong e = n; e > h + 1; e -= RB) { slong s = e - RB > h + 1 ? e - RB : h + 1; brow0[nblk] = s; bnr[nblk] = e - s; nblk++; }
    for (slong B = 0; B < nblk; B++) for (slong i = brow0[B]; i < brow0[B] + bnr[B]; i++) blk_of[i] = B;
    rowbuf = _fmpz_vec_init(2 * RB * n * 2);
}
static void blocks_clear(void) {
    free(brow0); free(bnr); free(bfill); free(blk_of); _fmpz_vec_clear(rowbuf, 2 * RB * n * 2); rowbuf = NULL;
}

/* index save/load (text) */
static void save_index(FILE *fp) {
    fprintf(fp, "%ld\n", (long)nseg);
    for (slong s = 0; s < nseg; s++) {
        fprintf(fp, "%d %ld %ld %ld %ld %ld %lld %zu", segs[s].lev, (long)segs[s].B, (long)segs[s].r0, (long)segs[s].nr, (long)segs[s].w, (long)segs[s].ntile, (long long)segs[s].off, segs[s].bytes);
        for (slong q = 0; q < segs[s].ntile; q++) fprintf(fp, " %ld", (long)segs[s].tiles[q]);
        fprintf(fp, "\n");
    }
}
static int load_index(FILE *fp) {
    long ns; if (fscanf(fp, "%ld", &ns) != 1) return 1;
    for (long s = 0; s < ns; s++) {
        seg_t S; long B, r0, nr, w, nt; long long off; int lev;
        if (fscanf(fp, "%d %ld %ld %ld %ld %ld %lld %zu", &lev, &B, &r0, &nr, &w, &nt, &off, &S.bytes) != 8) return 1;
        S.lev = lev; S.B = B; S.r0 = r0; S.nr = nr; S.w = w; S.ntile = nt; S.off = off; S.mem = NULL; S.RZ = NULL;
        S.tiles = malloc(sizeof(slong) * nt);
        for (long q = 0; q < nt; q++) { long x; if (fscanf(fp, "%ld", &x) != 1) return 1; S.tiles[q] = x; }
        if (nseg == segcap) { segcap *= 2; segs = realloc(segs, sizeof(seg_t) * segcap); }
        segs[nseg++] = S; tot_bytes += S.bytes;
    }
    int ch; while ((ch = fgetc(fp)) != '\n' && ch != EOF) ;
    return 0;
}
static void pread_all(void *buf, size_t len, off_t off) {
    size_t done = 0;
    while (done < len) {
        ssize_t k = pread(kfd, (char *) buf + done, len - done, off + done);
        if (k <= 0) { if (k < 0 && errno == EINTR) continue; perror("read K.bin"); exit(3); }
        done += k;
    }
    if (nocache) posix_fadvise(kfd, off, len, POSIX_FADV_DONTNEED);
}
/* make segments resident while the budget allows (after a resume) */
static void load_resident(void) {
    for (slong s = 0; s < nseg; s++) {
        if (segs[s].mem) continue;
        if (res_bytes + segs[s].bytes > res_cap) continue;
        segs[s].mem = malloc(segs[s].bytes);
        pread_all(segs[s].mem, segs[s].bytes, segs[s].off);
        res_bytes += segs[s].bytes;
    }
}

/* ------------------------------------------------------------------ stage 2: rows of K by the commuting recurrence */
static void save_build_state(int dir, slong inext, acb_srcptr pr, acb_srcptr cu, acb_srcptr chk, slong em_all) {
    fdatasync(kfd);
    char *tmp = spath("build.state.tmp"), *fin = spath("build.state");
    FILE *fp = fopen(tmp, "w"); if (!fp) return;
    fprintf(fp, "L3B %s\n%d %ld %ld %ld %d %lld\n", sigbuf, dir, (long)inext, (long)em_all, (long)Emaxseen, overflow, (long long)kend);
    save_index(fp);
    int bad = 0;
    if (inext >= 0) bad |= dump_acbvec(fp, pr, n) | dump_acbvec(fp, cu, n);
    if (dir == 1) bad |= dump_acbvec(fp, chk, n);
    if (bad) { fclose(fp); return; }
    atomic_close(fp, tmp, fin);
}
static double tbuild_io = 0;
static void phaseB(acb_srcptr R0, acb_srcptr R1, slong p2, slong p3, double *lkmax, double *qK) {
    P2bits = p2; store_init(p3);
    slong h = hrow;
    char line[1024];
    if (scratch) {   /* finished K for these parameters? */
        FILE *fp = fopen(spath("K.meta"), "r");
        if (fp) {
            double lk, q; long em, ke; int ov;
            if (fgets(line, sizeof line, fp) && !strncmp(line, "L3K ", 4) && !strncmp(line + 4, sigbuf, strlen(sigbuf)) && line[4 + strlen(sigbuf)] == '\n'
                && fscanf(fp, "%lf %lf %ld %d %ld", &lk, &q, &em, &ov, &ke) == 5 && !load_index(fp)) {
                *lkmax = lk; *qK = q; Emaxseen = em; overflow = ov; kend = ke;
                fclose(fp);
                load_resident();
                if (!quiet) fprintf(stderr, "# K loaded from %s (%zu segments, %.2f GB)\n", spath("K.bin"), (size_t)nseg, tot_bytes / 1e9);
                return;
            }
            fclose(fp);
            for (slong s = 0; s < nseg; s++) free(segs[s].tiles);
            nseg = 0; tot_bytes = 0;
        }
    }
    blocks_init();
    slong em_all = WORD_MIN, e;
    double lmis = -1e300;
    arb_ptr al2 = _arb_vec_init(n), u2 = _arb_vec_init(n), l2 = _arb_vec_init(n), iu = _arb_vec_init(n), il = _arb_vec_init(n);
    for (slong i = 0; i < n; i++) {
        arb_set_round(al2 + i, al + i, p2);
        if (i < n - 1) {
            arb_set_round(u2 + i, u + i, p2); arb_set_round(l2 + i, l + i, p2);
            arb_inv(iu + i, u + i, p2); arb_inv(il + i, l + i, p2); amid(iu + i); amid(il + i);
        }
    }
    acb_ptr chk = _acb_vec_init(n), pr = _acb_vec_init(n), cu = _acb_vec_init(n), nx = _acb_vec_init(n);
    int dir0 = 0; slong ires = -1;
    if (scratch) {   /* resume a partial build */
        FILE *fp = fopen(spath("build.state"), "r");
        int ok = 0;
        if (fp) {
            int d, ov; long in, ea, em; long long ke;
            if (fgets(line, sizeof line, fp) && !strncmp(line, "L3B ", 4) && !strncmp(line + 4, sigbuf, strlen(sigbuf)) && line[4 + strlen(sigbuf)] == '\n'
                && fscanf(fp, "%d %ld %ld %ld %d %lld", &d, &in, &ea, &em, &ov, &ke) == 6 && !load_index(fp)) {
                ok = 1;
                if (in >= 0 && (load_acbvec(fp, pr, n) || load_acbvec(fp, cu, n))) ok = 0;
                if (ok && d == 1 && load_acbvec(fp, chk, n)) ok = 0;
                if (ok) {
                    dir0 = d; ires = in; em_all = ea; Emaxseen = em; overflow = ov; kend = ke;
                    if (ftruncate(kfd, kend)) perror("ftruncate");
                    for (slong s = 0; s < nseg; s++) bfill[segs[s].B] = bnr[segs[s].B];
                    if (!quiet) fprintf(stderr, "# resumed stage 2 from %s (sweep %d, row %ld, %ld segments)\n", spath("build.state"), d, in, (long)nseg);
                }
            }
            fclose(fp);
        }
        if (!ok) {
            for (slong s = 0; s < nseg; s++) free(segs[s].tiles);
            nseg = 0; tot_bytes = 0; kend = 0;
            if (ftruncate(kfd, 0)) perror("ftruncate");
        }
    }
    double lastck = wall();
    for (int dir = dir0; dir < 2; dir++) {
        acb_srcptr R = dir == 0 ? R0 : R1;
        slong i0 = dir == 0 ? 0 : n - 1, st = dir == 0 ? 1 : -1, i;
        if (dir == dir0 && ires >= 0) i = ires;
        else {
            for (slong j = 0; j < n; j++) acb_set_round(cu + j, R + j, p2);
            rowstore(i0, cu, &e); if (e > em_all) em_all = e;
            i = i0;
        }
        for (; i != h + 1; i += st) {
            #pragma omp parallel
            {
                acb_t s; arb_t da; acb_init(s); arb_init(da);
                #pragma omp for schedule(static)
                for (slong j = 0; j < n; j++) {
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
            blockdone_flag = 0;
            if (dir == 1 || i + st <= h) { rowstore(i + st, nx, &e); if (e > em_all) em_all = e; }
            acb_ptr tmp = pr; pr = cu; cu = nx; nx = tmp;
            if (scratch && blockdone_flag && wall() - lastck > ckint && i + st != h + 1) {
                double tc = wall();
                if (em_all > Emaxseen) Emaxseen = em_all;
                save_build_state(dir, i + st, pr, cu, chk, em_all);
                tbuild_io += wall() - tc; lastck = wall();
            }
        }
        if (dir == 0) {
            _acb_vec_set(chk, cu, n);
            if (scratch && wall() - lastck > ckint) { if (em_all > Emaxseen) Emaxseen = em_all; save_build_state(1, -1, pr, cu, chk, em_all); lastck = wall(); }
        } else {
            acb_t d; acb_init(d); slong em = WORD_MIN;
            for (slong j = 0; j < n; j++) { acb_sub(d, chk + j, cu + j, p2); slong ee = lexpc(d); if (ee > em) em = ee; }
            lmis = em == WORD_MIN ? -1e300 : em * LOG10_2;
            acb_clear(d);
        }
    }
    _acb_vec_clear(chk, n); _acb_vec_clear(pr, n); _acb_vec_clear(cu, n); _acb_vec_clear(nx, n);
    _arb_vec_clear(al2, n); _arb_vec_clear(u2, n); _arb_vec_clear(l2, n); _arb_vec_clear(iu, n); _arb_vec_clear(il, n);
    blocks_clear();
    if (em_all > Emaxseen) Emaxseen = em_all;
    *lkmax = em_all * LOG10_2;
    *qK = *lkmax - lmis;
    if (scratch) {
        fdatasync(kfd);
        char *tmp = spath("K.meta.tmp");
        FILE *fp = fopen(tmp, "w");
        if (fp) {
            fprintf(fp, "L3K %s\n%.17g %.17g %ld %d %lld\n", sigbuf, *lkmax, *qK, (long)Emaxseen, overflow, (long long)kend);
            save_index(fp);
            atomic_close(fp, tmp, spath("K.meta"));
        }
        unlink(spath("build.state"));
        load_resident();
    }
}

/* ------------------------------------------------------------------ products Z = K Y */
typedef struct { slong s, t0, t1; int last; size_t boff; } piece_t;
typedef struct { slong p0, p1; int lev; size_t sbytes; } chunk_t;
static piece_t *pcs; static slong npcs;
static chunk_t *chs; static slong nchs;
static mp_limb_t *sbuf[2];
static slong *keyfirst, *keynext, *keyleft;   /* segments of key (lev * nblk + B): linked list */
static slong nkeys;
static mp_limb_t *Yr = NULL; static slong Yr_m = 0, Yr_np = 0;
static double t_iowait = 0, t_crt = 0, t_mul = 0; static double io_bytes = 0;

static void plan_free(void) {
    if (!plan_ok) return;
    free(pcs); free(chs); free(keyfirst); free(keynext); free(keyleft);
    free(sbuf[0]); free(sbuf[1]); sbuf[0] = sbuf[1] = NULL;
    plan_ok = 0;
}
static void plan_build(void) {
    plan_free();
    slong maxB = 0; for (slong s = 0; s < nseg; s++) if (segs[s].B > maxB) maxB = segs[s].B;
    nblk = maxB + 1; nkeys = NL * nblk;
    keyfirst = malloc(sizeof(slong) * nkeys); keyleft = calloc(nkeys, sizeof(slong)); keynext = malloc(sizeof(slong) * (nseg + 1));
    for (slong k = 0; k < nkeys; k++) keyfirst[k] = -1;
    for (slong s = nseg - 1; s >= 0; s--) { slong k = segs[s].lev * nblk + segs[s].B; keynext[s] = keyfirst[k]; keyfirst[k] = s; }
    slong cap = 64; pcs = malloc(sizeof(piece_t) * cap); npcs = 0;
    slong ccap = 64; chs = malloc(sizeof(chunk_t) * ccap); nchs = 0;
    for (int v = 0; v < NL; v++) {
        chunk_t C = { npcs, npcs, v, 0 };
        for (slong s = 0; s < nseg; s++) {
            if (segs[s].lev != v) continue;
            slong np = npl[v];
            size_t pp = (size_t) 2 * segs[s].w * segs[s].nr * sizeof(mp_limb_t);
            slong tper = segs[s].mem ? np : (slong)(bufsz / pp);
            if (tper < 1) tper = 1;
            for (slong t0 = 0; t0 < np; t0 += tper) {
                slong t1 = t0 + tper < np ? t0 + tper : np;
                size_t by = segs[s].mem ? 0 : pp * (t1 - t0);
                if (C.sbytes + by > bufsz && C.p1 > C.p0) {
                    if (nchs == ccap) { ccap *= 2; chs = realloc(chs, sizeof(chunk_t) * ccap); }
                    chs[nchs++] = C; C.p0 = C.p1 = npcs; C.sbytes = 0;
                }
                if (npcs == cap) { cap *= 2; pcs = realloc(pcs, sizeof(piece_t) * cap); }
                piece_t P = { s, t0, t1, t1 == np, C.sbytes };
                pcs[npcs++] = P; C.p1 = npcs; C.sbytes += by;
            }
        }
        if (C.p1 > C.p0) { if (nchs == ccap) { ccap *= 2; chs = realloc(chs, sizeof(chunk_t) * ccap); } chs[nchs++] = C; }
    }
    size_t mx = 0; for (slong c = 0; c < nchs; c++) if (chs[c].sbytes > mx) mx = chs[c].sbytes;
    if (mx) { sbuf[0] = malloc(mx); sbuf[1] = malloc(mx); }
    plan_ok = 1;
}

typedef struct { slong c; } rdarg_t;
static void *reader(void *arg) {
    slong c = ((rdarg_t *) arg)->c;
    mp_limb_t *buf = sbuf[c & 1];
    for (slong p = chs[c].p0; p < chs[c].p1; p++) {
        seg_t *S = segs + pcs[p].s;
        if (S->mem) continue;
        size_t pp = (size_t) 2 * S->w * segs[pcs[p].s].nr * sizeof(mp_limb_t);
        pread_all((char *) buf + pcs[p].boff, pp * (pcs[p].t1 - pcs[p].t0), S->off + pp * pcs[p].t0);
    }
    return NULL;
}

static void kmul(acb_mat_t Z, const acb_mat_t Y, slong prec) {
    slong m = acb_mat_ncols(Y), m2 = 2 * m;
    if (!plan_ok) plan_build();
    if (Yr == NULL || Yr_m != m || Yr_np < npmax) { free(Yr); Yr = malloc(sizeof(mp_limb_t) * npmax * n * m2); Yr_m = m; Yr_np = npmax; }
    slong *F = malloc(sizeof(slong) * m);
    for (slong k = 0; k < m; k++) {
        slong em = WORD_MIN;
        for (slong j = 0; j < n; j++) { slong e = lexpc(acb_mat_entry(Y, j, k)); if (e > em) em = e; }
        F[k] = em == WORD_MIN ? 0 : em;
    }
    acb_mat_zero(Z);
    for (slong k = 0; k < nkeys; k++) { keyleft[k] = 0; for (slong s = keyfirst[k]; s >= 0; s = keynext[s]) keyleft[k]++; }
    int cur = -1;
    pthread_t th; rdarg_t ra[2]; int running = 0;
    if (nchs > 0 && chs[0].sbytes) { double tw = wall(); ra[0].c = 0; reader(ra); t_iowait += wall() - tw; io_bytes += chs[0].sbytes; }
    slong *done = malloc(sizeof(slong) * (nkeys + 1));
    for (slong c = 0; c < nchs; c++) {
        if (c + 1 < nchs && chs[c + 1].sbytes) { ra[(c + 1) & 1].c = c + 1; pthread_create(&th, NULL, reader, ra + ((c + 1) & 1)); running = 1; io_bytes += chs[c + 1].sbytes; }
        int v = chs[c].lev; slong np = npl[v], bet = bl[v] + YM;
        if (v != cur) {   /* residues of Y truncated to beta_v bits */
            cur = v;
            #pragma omp parallel
            {
                fmpz_t z; fmpz_init(z);
                fmpz_comb_temp_t ct; fmpz_comb_temp_init(ct, combs + v);
                mp_limb_t *res = malloc(sizeof(mp_limb_t) * np);
                #pragma omp for schedule(static)
                for (slong j = 0; j < n; j++) for (slong k = 0; k < m; k++) for (int part = 0; part < 2; part++) {
                    arf_get_fmpz_fixed_si(z, arb_midref(part ? acb_imagref(acb_mat_entry(Y, j, k)) : acb_realref(acb_mat_entry(Y, j, k))), F[k] - bet);
                    fmpz_multi_mod_ui(res, z, combs + v, ct);
                    for (slong t = 0; t < np; t++) Yr[(t * n + j) * m2 + k + part * m] = res[t];
                }
                free(res); fmpz_comb_temp_clear(ct); fmpz_clear(z);
            }
        }
        /* tasks (piece, prime) */
        slong ntask = 0;
        for (slong p = chs[c].p0; p < chs[c].p1; p++) {
            seg_t *S = segs + pcs[p].s;
            if (S->RZ == NULL) S->RZ = malloc(sizeof(mp_limb_t) * segs[pcs[p].s].nr * m2 * np);
            ntask += pcs[p].t1 - pcs[p].t0;
        }
        slong *tp = malloc(sizeof(slong) * 2 * ntask);
        for (slong p = chs[c].p0, q = 0; p < chs[c].p1; p++) for (slong t = pcs[p].t0; t < pcs[p].t1; t++) { tp[2 * q] = p; tp[2 * q + 1] = t; q++; }
        double tm = wall();
        mp_limb_t *cbuf = sbuf[c & 1];
        #pragma omp parallel
        {
            mp_limb_t *W = malloc(sizeof(mp_limb_t) * n * m2), *Zb = malloc(sizeof(mp_limb_t) * 2 * RB * m2);
            mp_limb_t **rw = malloc(sizeof(mp_limb_t *) * (n + 4 * RB));
            #pragma omp for schedule(dynamic, 1)
            for (slong q = 0; q < ntask; q++) {
                slong p = tp[2 * q], t = tp[2 * q + 1];
                seg_t *S = segs + pcs[p].s;
                slong nr = segs[pcs[p].s].nr, w = S->w;
                mp_limb_t pt = primes[t];
                mp_limb_t *base = S->mem ? S->mem + (size_t) t * 2 * nr * w : cbuf + pcs[p].boff / sizeof(mp_limb_t) + (size_t)(t - pcs[p].t0) * 2 * nr * w;
                slong c0 = 0;
                for (slong g = 0; g < S->ntile; g++) {
                    slong j0 = S->tiles[g] * CB, j1 = j0 + CB < n ? j0 + CB : n;
                    memcpy(W + c0 * m2, Yr + ((size_t) t * n + j0) * m2, sizeof(mp_limb_t) * (j1 - j0) * m2);
                    c0 += j1 - j0;
                }
                nmod_mat_struct A, Bv, C1, C2;
                nmod_t mod; nmod_init(&mod, pt);
                mp_limb_t **rA = rw, **rB = rw + 2 * RB, **rC = rw + 2 * RB + w;
                for (slong a = 0; a < w; a++) rB[a] = W + a * m2;
                for (slong a = 0; a < nr; a++) { rC[a] = Zb + a * m2; rC[RB + a] = Zb + (RB + a) * m2; }
                Bv.entries = W; Bv.r = w; Bv.c = m2; Bv.rows = rB; Bv.mod = mod;
                C1.entries = Zb; C1.r = nr; C1.c = m2; C1.rows = rC; C1.mod = mod;
                C2.entries = Zb + RB * m2; C2.r = nr; C2.c = m2; C2.rows = rC + RB; C2.mod = mod;
                for (int part = 0; part < 2; part++) {
                    mp_limb_t *Kp = base + (size_t) part * nr * w;
                    for (slong a = 0; a < nr; a++) rA[part * RB + a] = Kp + a * w;
                }
                A.entries = base; A.r = nr; A.c = w; A.rows = rA; A.mod = mod;
                nmod_mat_mul(&C1, &A, &Bv);
                A.entries = base + (size_t) nr * w; A.rows = rA + RB;
                nmod_mat_mul(&C2, &A, &Bv);
                mp_limb_t *RZ = S->RZ;
                for (slong a = 0; a < nr; a++) for (slong k = 0; k < m; k++) {
                    RZ[((a * m + k) * 2 + 0) * np + t] = n_submod(rC[a][k], rC[RB + a][k + m], pt);
                    RZ[((a * m + k) * 2 + 1) * np + t] = n_addmod(rC[a][k + m], rC[RB + a][k], pt);
                }
            }
            free(W); free(Zb); free(rw);
        }
        free(tp);
        t_mul += wall() - tm;
        /* keys completed in this chunk: CRT, scale, add into Z */
        slong nd = 0;
        for (slong p = chs[c].p0; p < chs[c].p1; p++) if (pcs[p].last) {
            slong k = segs[pcs[p].s].lev * nblk + segs[pcs[p].s].B;
            if (--keyleft[k] == 0) done[nd++] = k;
        }
        if (nd) {
            double tc = wall();
            slong per = RB * m * 2;
            #pragma omp parallel
            {
                fmpz_t z; fmpz_init(z); arf_t x; arf_init(x);
                fmpz_comb_temp_t ct; fmpz_comb_temp_init(ct, combs + v);
                mp_limb_t *res = malloc(sizeof(mp_limb_t) * np);
                #pragma omp for schedule(dynamic, 16)
                for (slong q = 0; q < nd * per; q++) {
                    slong k = done[q / per], e = q % per, s = keyfirst[k];
                    slong a = e / (2 * m), kk = (e / 2) % m; int part = e % 2;
                    if (a >= segs[s].nr) continue;
                    for (slong t = 0; t < np; t++) res[t] = segs[s].RZ[e * np + t];
                    for (slong s2 = keynext[s]; s2 >= 0; s2 = keynext[s2])
                        for (slong t = 0; t < np; t++) res[t] = n_addmod(res[t], segs[s2].RZ[e * np + t], primes[t]);
                    fmpz_multi_CRT_ui(z, res, combs + v, ct, 1);
                    arf_set_round_fmpz(x, z, prec, ARF_RND_NEAR);
                    arf_mul_2exp_si(x, x, Etol + F[kk] - bet);
                    arb_ptr y = part ? acb_imagref(acb_mat_entry(Z, segs[s].r0 + a, kk)) : acb_realref(acb_mat_entry(Z, segs[s].r0 + a, kk));
                    arf_add(arb_midref(y), arb_midref(y), x, prec, ARF_RND_NEAR);
                }
                free(res); fmpz_comb_temp_clear(ct); fmpz_clear(z); arf_clear(x);
            }
            for (slong d = 0; d < nd; d++) for (slong s = keyfirst[done[d]]; s >= 0; s = keynext[s]) { free(segs[s].RZ); segs[s].RZ = NULL; }
            t_crt += wall() - tc;
        }
        if (running) { double tw = wall(); pthread_join(th, NULL); running = 0; t_iowait += wall() - tw; }
    }
    free(done); free(F);
}

/* ------------------------------------------------------------------ orthonormalization (as ladder2) */
static void orth(acb_mat_t X, slong prec) {
    slong m = acb_mat_ncols(X);
    int nt = omp_get_max_threads();
    acb_ptr part = _acb_vec_init(nt * (m + 1)), s = _acb_vec_init(m + 1);
    for (slong j = 0; j < m; j++) {
        for (int pass = 0; pass < 3; pass++) {
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
    fprintf(stderr, "usage: ladder3 r c LR|LLR [-t threads] [-scratch DIR] [-mem GB] [-nocache] [-ckint SEC] [-levels NL] [-rb RB] [-cb CB]\n"
                    "               [-d1 D] [-d2 D] [-d3 D] [-m 12] [-ntop 8] [-tol 22] [-maxit 100] [-ck file] [-q]\n");
    exit(1);
}

/* stage 1 state (rows 0, n-1 of K) */
static void save_stage1(slong d1, double lvmax, acb_srcptr R0, acb_srcptr R1) {
    char *tmp = spath("stage1.state.tmp");
    FILE *fp = fopen(tmp, "w"); if (!fp) return;
    fprintf(fp, "L3S %ld %ld %s %ld\n%ld %.17g %.17g %ld\n", (long)r, (long)cc, pw == 1 ? "LR" : "LLR", (long)n, (long)d1, lvmax, lcomp, (long)Esample);
    if (dump_acbvec(fp, R0, n) || dump_acbvec(fp, R1, n)) { fclose(fp); return; }
    atomic_close(fp, tmp, spath("stage1.state"));
}
static int load_stage1(slong *d1, double *lvmax, acb_ptr R0, acb_ptr R1) {
    FILE *fp = fopen(spath("stage1.state"), "r"); if (!fp) return 0;
    long rr, c2, nn, dd, es; char wd[8]; double lv, lc; int ok = 0;
    if (fscanf(fp, "L3S %ld %ld %7s %ld %ld %lf %lf %ld", &rr, &c2, wd, &nn, &dd, &lv, &lc, &es) == 8 && rr == r && c2 == cc && nn == n && !strcmp(wd, pw == 1 ? "LR" : "LLR")) {
        int ch; while ((ch = fgetc(fp)) != '\n' && ch != EOF) ;
        if (!load_acbvec(fp, R0, n) && !load_acbvec(fp, R1, n)) { ok = 1; *d1 = dd; *lvmax = lv; lcomp = lc; Esample = es; }
    }
    fclose(fp);
    return ok;
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
        else if (!strcmp(argv[a], "-scratch") && a + 1 < argc) scratch = argv[++a];
        else if (!strcmp(argv[a], "-mem") && a + 1 < argc) memgb = atof(argv[++a]);
        else if (!strcmp(argv[a], "-ckint") && a + 1 < argc) ckint = atof(argv[++a]);
        else if (!strcmp(argv[a], "-levels") && a + 1 < argc) { NL = atoi(argv[++a]); if (NL < 1) NL = 1; if (NL > NLMAX) NL = NLMAX; }
        else if (!strcmp(argv[a], "-rb") && a + 1 < argc) RB = atol(argv[++a]);
        else if (!strcmp(argv[a], "-cb") && a + 1 < argc) CB = atol(argv[++a]);
        else if (!strcmp(argv[a], "-segmb") && a + 1 < argc) segmb = atof(argv[++a]);
        else if (!strcmp(argv[a], "-nocache")) nocache = 1;
        else if (!strcmp(argv[a], "-q")) quiet = 1;
        else usage();
    }
    if (nth > 0) omp_set_num_threads(nth);
    nth = omp_get_max_threads();
    flint_set_num_threads(1);
    int user1 = d1 > 0, user2 = d2 > 0, user3 = d3 > 0;
    double T0 = wall();
    char ckdef[4096];
    if (scratch && !ck) { snprintf(ckdef, sizeof ckdef, "%s/iter.ck", scratch); ck = ckdef; }

    ad = malloc(sizeof(slong) * (r + 1)); n = 0;
    for (slong a = 0; a < r - 2; a += 2) if (adm(a, a, cc)) ad[n++] = a;
    if (n < 3) { fprintf(stderr, "error: n = %ld < 3; use ladder2\n", (long)n); return 1; }
    if (m > n) m = n;
    if (ntop > m) ntop = m;

    /* precision control: as ladder2 (see README_ladder2.md) */
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
    int d3fromck = loaded && ckd3 > 0 && !user3;
    if (d3fromck) d3 = ckd3;
    double lvmax = 0, lkmax = 0, qK = 0, qabs = 0, tA = 0, tB = 0, rssA = 0, rssB = 0;
    acb_ptr R0 = _acb_vec_init(n), R1 = _acb_vec_init(n);
    int d2set = 0;
    for (int outer = 0; ; outer++) {
        /* ---- stage 1 ---- */
        for (int attempt = 0; ; attempt++) {
            double ta = wall();
            if (scratch && outer == 0 && attempt == 0 && !user1) {
                slong dd; double lv;
                if (load_stage1(&dd, &lv, R0, R1)) {
                    d1 = dd; lvmax = lv;
                    setup(DIG2BITS(d1));
                    if (!quiet) fprintf(stderr, "# stage 1 loaded from %s: d1 %ld | log10 max|V| %.1f | completeness residual 1e%.1f\n", spath("stage1.state"), (long)d1, lvmax, lcomp);
                    break;
                }
            }
            setup(DIG2BITS(d1));
            phaseA(R0, R1, &lvmax, DIG2BITS(d1));
            tA += wall() - ta;
            if (!quiet) fprintf(stderr, "# stage 1: d1 %ld | log10 max|V| %.1f | completeness residual 1e%.1f (need <= 1e%.0f) | %.1fs | RSS %.2f GB\n",
                                (long)d1, lvmax, lcomp, ctarget, wall() - ta, maxrss_gb());
            if (lcomp <= ctarget || user1 || attempt >= 3) {
                if (lcomp > ctarget) fprintf(stderr, "# WARNING: stage 1 completeness residual too large; results may be inaccurate\n");
                if (scratch) save_stage1(d1, lvmax, R0, R1);
                break;
            }
            d1 += (slong)(lcomp - ctarget) + 15;
            setup_clear();
        }
        rssA = maxrss_gb();
        /* ---- stage 2 ---- */
        double lkest = (pw == 1 ? 1.02 : 1.4) * lvmax;
        if (!user3 && !d3fromck && outer == 0) d3 = (slong)(lkest + 35);
        if (!user2 && !d2set) {
            d2 = (slong)(qreq + (1.77 - 1.43 * xr) * lvmax + (pw == 2 ? (1.0 - 0.8 * xr) * lvmax : 0) + 20);
            d2set = 1;
        }
        if (outer == 0) Eref = Esample;   /* fixed-point reference: sampled max|K| (row h) */
        int d1limited = 0; double prevq = -1e300; slong prevraise = 0;
        for (int attempt = 0; ; attempt++) {
            if (d2 < d3) d2 = d3;
            double tb = wall();
            phaseB(R0, R1, DIG2BITS(d2), DIG2BITS(d3), &lkmax, &qK);
            tB += wall() - tb;
            qabs = lvmax - (lkmax - qK);
            double qneed = qreq + lkmax - lvmax;
            if (!quiet) fprintf(stderr, "# stage 2: d2 %ld d3 %ld | log10 max|K| %.1f | K accuracy %.1f digits (need %.1f) | %.1fs | log2 max|K| %ld (sampled %ld) | K %.2f GB (ladder2 %.2f GB), %ld segments, resident %.2f GB | RSS %.2f GB\n",
                                (long)d2, (long)d3, lkmax, qK, qneed, wall() - tb, (long)Emaxseen, (long)Esample, tot_bytes / 1e9, ref_bytes / 1e9, (long)nseg, res_bytes / 1e9, maxrss_gb());
            int redo = 0;
            if (attempt < 4) {
                if (!user3 && lkmax + 33 > d3) { d3 = (slong)(lkmax + 35); redo = 1; }
                if (!user2 && qK < qneed) {
                    if (prevraise > 0 && qK - prevq < 0.5 * prevraise && !user1) { d1limited = 1; d1 += (slong)(qneed - qK) + 20; break; }
                    prevraise = (slong)(qneed - qK) + 20; d2 += prevraise; redo = 1;
                } else prevraise = 0;
                if (overflow) { Eref = Emaxseen + 2; redo = 1; if (!quiet) fprintf(stderr, "# fixed-point range exceeded: Eref -> %ld\n", (long)Eref); }
            }
            prevq = qK;
            if (!redo) {
                if (qK < qneed) fprintf(stderr, "# WARNING: K accuracy too low; results may be inaccurate\n");
                if (overflow) fprintf(stderr, "# WARNING: K fixed-point overflow; results may be inaccurate\n");
                break;
            }
            store_free();
        }
        if (d1limited && outer < 2) {
            if (!quiet) fprintf(stderr, "# stage 2 limited by stage 1 precision: d1 -> %ld\n", (long)d1);
            store_free(); setup_clear(); continue;
        }
        break;
    }
    rssB = maxrss_gb();
    double T2 = wall();

    /* ---------------- subspace iteration with A' = K Theta ---------------- */
    slong prec = DIG2BITS(d3);
    acb_ptr tha = _acb_vec_init(n);
    int d3checked = user3 || d3fromck;
    slong minit = 0;
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
        /* Rayleigh-Ritz: H = X^H Z (Y is free after the product: it holds conj(X); no per-thread copies) */
        #pragma omp parallel for schedule(static)
        for (slong i = 0; i < n; i++) for (slong j = 0; j < m; j++) acb_conj(acb_mat_entry(Y, i, j), acb_mat_entry(X, i, j));
        #pragma omp parallel for schedule(dynamic, 1)
        for (slong q = 0; q < m * m; q++) {
            slong ia = q / m, ib = q % m;
            acb_dot(acb_mat_entry(H, ia, ib), NULL, 0, acb_mat_entry(Y, 0, ia), m, acb_mat_entry(Z, 0, ib), m, n, prec);
            cmid(acb_mat_entry(H, ia, ib));
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
            d3checked = 1;
            double l1 = arf_get_d(arb_midref(logs + 0), ARF_RND_NEAR) / log(10.0);
            slong need = (slong) ceil(lkmax - (l1 - 2) + arate * r + 30);
            /* effective digits: the fixed-point tolerance is 2^(Eref - BK); the true max|K| is 2^Emaxseen */
            slong d3eff = d3 + (slong) floor((Emaxseen - 1 - Eref) * LOG10_2);
            if (!quiet) fprintf(stderr, "# log10|lambda_1| ~ %.1f: stage 3 needs d3 >= %ld (have %ld)\n", l1, (long)need, (long)d3eff);
            if (need > d3eff) {
                double tb = wall();
                if (need > d3) d3 = need;
                if (d2 < d3) d2 = d3;
                store_free();
                phaseB(R0, R1, DIG2BITS(d2), DIG2BITS(d3), &lkmax, &qK);
                if (overflow) { Eref = Emaxseen + 2; store_free(); phaseB(R0, R1, DIG2BITS(d2), DIG2BITS(d3), &lkmax, &qK); }
                tB += wall() - tb;
                if (!quiet) fprintf(stderr, "# rebuilt K at d3 %ld (%.1fs) | K %.2f GB (ladder2 %.2f GB) | RSS %.2f GB\n", (long)d3, wall() - tb, tot_bytes / 1e9, ref_bytes / 1e9, maxrss_gb());
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
    fprintf(stderr, "# r %ld c %ld %s n %ld threads %d | digits d1 %ld d2 %ld d3 %ld | maxV %.1f maxK %.1f Kacc %.1f comp %.1f | eigvec %.1fs recur %.1fs iter %.1fs (%ld its) total %.1fs | "
                    "K %.2f GB (%.2f of ladder2), resident %.2f GB, streamed %.1f GB, product %.1fs crt %.1fs io-wait %.1fs | peak RSS: stage1 %.2f stage2 %.2f final %.2f GB\n",
            (long)r, (long)cc, word, (long)n, nth, (long)d1, (long)d2, (long)d3, lvmax, lkmax, qabs, lcomp, tA, tB, T4 - T2, (long)(it - it0 + 1), T4 - T0,
            tot_bytes / 1e9, tot_bytes / ref_bytes, res_bytes / 1e9, io_bytes / 1e9, t_mul, t_crt, t_iowait, rssA, rssB, maxrss_gb());
    return 0;
}
