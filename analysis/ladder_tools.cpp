// Trace-map fixed points, multipliers, volumes, twist coefficients and ladder fits.
// robust: fit variants (drop one r, level ranges, linear term, 1/r^2) against pi*Im(b), and the direct test
// gamma_row - pi*Im b(theta_row) -> 0 per raw row; gaps: the same from the r = 201..1001 level gaps;
// stress: reference check and predictions for longer words, conjugate chosen by the signed volume sum (rule S).
// order2: third/fourth differences and the D(n) fit against d3_pred = (4 pi^2/3) Re c, c the next Birkhoff coefficient;
// alpha: the theta = 0 constant alpha, alpha0 = alpha + gamma/4, S2 and continued-fraction convergents.
// Least-squares designs with rcond < 1e-8 print n/a, tagged singular.
// Build: g++ -O2 -std=c++17 -Wall -Wextra analysis/ladder_tools.cpp -o /tmp/ladder_tools
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

typedef long double R;
typedef complex<R> C;
typedef array<C, 3> V;
typedef array<array<C, 3>, 3> M3;
const R PI = acosl(-1.0L);
const C I1(0, 1);

static void require(bool ok, const string& msg) { if (!ok) throw runtime_error(msg); }

// ---------- trace map ----------
V letter(char ch, const V& v) {
    return ch == 'L' ? V{v[2], v[1], v[2] * v[1] - v[0]} : V{v[0], v[2], v[0] * v[2] - v[1]};
}
V word_apply(const string& w, V v) {
    for (int i = (int)w.size() - 1; i >= 0; i--) v = letter(w[i], v);
    return v;
}
C kappa(const V& v) { return v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - v[0] * v[1] * v[2] - (R)2; }
M3 mmul(const M3& a, const M3& b) {
    M3 c{};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) c[i][j] += a[i][k] * b[k][j];
    return c;
}
M3 dletter(char ch, const V& v) {
    C x = v[0], y = v[1], z = v[2];
    if (ch == 'R') return M3{{{1, 0, 0}, {0, 0, 1}, {z, -1, x}}};
    return M3{{{0, 0, 1}, {0, 1, 0}, {-1, z, y}}};
}
M3 word_jac(const string& w, V v) {
    M3 J{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    for (int i = (int)w.size() - 1; i >= 0; i--) { J = mmul(dletter(w[i], v), J); v = letter(w[i], v); }
    return J;
}

bool solve3(M3 a, V b, V& x) {
    for (int c = 0; c < 3; c++) {
        int p = c;
        for (int r = c + 1; r < 3; r++) if (abs(a[r][c]) > abs(a[p][c])) p = r;
        if (abs(a[p][c]) == 0) return false;
        swap(a[p], a[c]); swap(b[p], b[c]);
        for (int r = c + 1; r < 3; r++) {
            C f = a[r][c] / a[c][c];
            for (int k = c; k < 3; k++) a[r][k] -= f * a[c][k];
            b[r] -= f * b[c];
        }
    }
    for (int c = 2; c >= 0; c--) {
        C s = b[c];
        for (int k = c + 1; k < 3; k++) s -= a[c][k] * x[k];
        x[c] = s / a[c][c];
    }
    return true;
}

R vdist(const V& a, const V& b) { return max({abs(a[0] - b[0]), abs(a[1] - b[1]), abs(a[2] - b[2])}); }

vector<V> fixed_points(const string& w, R theta) {
    C k0 = -2 * cosl(theta / 2);
    uint64_t seed = 12345;
    auto rnd = [&] {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        return (R)(seed >> 11) / 9007199254740992.0L * 12 - 6;
    };
    auto Fn = [&](const V& v) { V fv = word_apply(w, v); return V{fv[0] - v[0], fv[1] - v[1], kappa(v) - k0}; };
    auto nrm = [](const V& f) { return max({abs(f[0]), abs(f[1]), abs(f[2])}); };
    // starts 0..399: plain Newton (the original search; its candidates and their order are kept);
    // starts 400..799: the same 400 starting points, damped Newton (halve the step while the residual does not decrease)
    vector<V> sols;
    for (int start = 0; start < 800; start++) {
        bool damped = start >= 400;
        if (start == 400) seed = 12345;
        V v;
        for (auto& c : v) { R re = rnd(); c = C(re, rnd()); }
        V F = Fn(v);
        for (int it = 0; it < (damped ? 400 : 200); it++) {
            M3 J = word_jac(w, v), A;
            A[0] = J[0]; A[0][0] -= (R)1;
            A[1] = J[1]; A[1][1] -= (R)1;
            A[2] = {(R)2 * v[0] - v[1] * v[2], (R)2 * v[1] - v[0] * v[2], (R)2 * v[2] - v[0] * v[1]};
            V dv, nv, nF;
            if (!solve3(A, F, dv)) break;
            R lam = 1;
            for (int h = 0; h < 40; h++) {
                for (int i = 0; i < 3; i++) nv[i] = v[i] - lam * dv[i];
                nF = Fn(nv);
                if (!damped || nrm(nF) < nrm(F)) break;
                lam /= 2;
            }
            v = nv; F = nF;
            R n = max({abs(v[0]), abs(v[1]), abs(v[2])});
            if (!(n < 1e6L)) break;
            if (lam * max({abs(dv[0]), abs(dv[1]), abs(dv[2])}) < 1e-18L * (1 + n)) break;
        }
        if (!(vdist(word_apply(w, v), v) < 1e-15L)) continue;
        bool dup = false;
        for (auto& s : sols) if (vdist(s, v) < 1e-9L) dup = true;
        if (!dup) sols.push_back(v);
    }
    return sols;
}

C multiplier(const string& w, const V& v) {
    M3 J = word_jac(w, v);
    C s = J[0][0] + J[1][1] + J[2][2] - (R)1;
    C mu = (s + sqrt(s * s - (R)4)) / (R)2;
    return abs(mu) < 1 ? (R)1 / mu : mu;
}

// ---------- volume ----------
R cl2(R x) {
    static vector<R> zeta;
    if (zeta.empty()) {
        zeta.assign(61, 0);
        zeta[1] = PI * PI / 6;
        for (int n = 2; n <= 60; n++) {
            R s = 0;
            for (int k = 1; k < n; k++) s += zeta[k] * zeta[n - k];
            zeta[n] = s / (n + 0.5L);
        }
    }
    x = fmodl(x, 2 * PI);
    if (x > PI) x -= 2 * PI;
    if (x < -PI) x += 2 * PI;
    if (x == 0) return 0;
    R u = (x / (2 * PI)) * (x / (2 * PI)), p = 1, s = 1 - logl(fabsl(x));
    for (int n = 1; n <= 60; n++) { p *= u; s += zeta[n] / (n * (2 * n + 1)) * p; }
    return x * s;
}
R lob(R t) { return cl2(2 * t) / 2; }
R bloch_wigner(C w) {
    if (w.imag() == 0) return 0;
    if (w.imag() < 0) return -bloch_wigner(conj(w));
    return lob(arg(w)) + lob(arg((R)1 / ((R)1 - w))) + lob(arg((R)1 - (R)1 / w));
}

struct Shape { C w; R D; };
// signed Bloch-Wigner sum; tflip takes the other root t -> 1/t, pmode 1/2 forces eigenvector p1/p2 (0: larger)
R volume_signed(const string& w, V v, R theta, int sgn, vector<Shape>* shapes = nullptr, bool tflip = false, int pmode = 0) {
    typedef array<array<C, 2>, 2> M2;
    auto mm = [](const M2& a, const M2& b) {
        M2 c{};
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j];
        return c;
    };
    auto mv = [](const M2& a, const array<C, 2>& p) { return array<C, 2>{a[0][0] * p[0] + a[0][1] * p[1], a[1][0] * p[0] + a[1][1] * p[1]}; };
    auto det = [](const array<C, 2>& a, const array<C, 2>& b) { return a[0] * b[1] - a[1] * b[0]; };
    C lam0 = -exp(I1 * (R)(sgn * theta / 2));
    R sum = 0;
    for (int i = (int)w.size() - 1; i >= 0; i--) {
        C x = v[0], y = v[1], z = v[2];
        C t = (z + (tflip ? (R)-1 : (R)1) * sqrt(z * z - (R)4)) / (R)2;
        M2 A{{{x, -1}, {1, 0}}}, B{{{0, t}, {-(R)1 / t, y}}};
        M2 Ai{{{0, 1}, {-1, x}}}, Bi{{{y, -t}, {(R)1 / t, 0}}};
        M2 K = mm(mm(Ai, Bi), mm(A, B));
        array<C, 2> p1{K[0][1], lam0 - K[0][0]}, p2{lam0 - K[1][1], K[1][0]};
        array<C, 2> p = pmode == 1 ? p1 : pmode == 2 ? p2 : norm(p1[0]) + norm(p1[1]) >= norm(p2[0]) + norm(p2[1]) ? p1 : p2;
        array<C, 2> q0 = p, q1 = mv(A, p), q2 = mv(mm(A, B), p), q3 = mv(B, p);
        C sh = det(q0, q2) * det(q1, q3) / (det(q0, q3) * det(q1, q2));
        R D = bloch_wigner(sh);
        if (shapes) shapes->push_back({sh, D});
        sum += D;
        v = letter(w[i], v);
    }
    return sum;
}
R volume(const string& w, V v, R theta, int sgn, vector<Shape>* shapes = nullptr) { return fabsl(volume_signed(w, v, theta, sgn, shapes)); }

// ---------- geometric fixed point ----------
struct Cand { V v; C mu; R vol; };
struct Geo { Cand best; bool maxmu; vector<Cand> all; };
bool nonreal(const V& v) { return max({fabsl(v[0].imag()), fabsl(v[1].imag()), fabsl(v[2].imag())}) > 1e-10L; }
Geo geometric(const string& w, R theta) {
    Geo g;
    for (auto& v : fixed_points(w, theta)) g.all.push_back({v, multiplier(w, v), volume(w, v, theta, 1)});
    int ib = -1, im = -1;
    for (int i = 0; i < (int)g.all.size(); i++) {
        if (!nonreal(g.all[i].v)) continue;
        // conjugate pairs tie in Vol: prefer Im x > 0
        if (ib < 0 || g.all[i].vol > g.all[ib].vol + 1e-12L ||
            (fabsl(g.all[i].vol - g.all[ib].vol) <= 1e-12L && g.all[i].v[0].imag() > g.all[ib].v[0].imag())) ib = i;
        if (im < 0 || abs(g.all[i].mu) > abs(g.all[im].mu)) im = i;
    }
    require(ib >= 0, "no non-real fixed point for " + w);
    g.best = g.all[ib];
    g.maxmu = abs(g.all[ib].mu) >= abs(g.all[im].mu) * (1 - 1e-12L);
    return g;
}

// ---------- twist coefficient ----------
struct Twist { C mu, c1, c2, b, c2p, j1, b2, c; R detlin, quadmax, cubmax, jmax, consist; };
// bivariate series truncated at total degree 5: c[i][j] multiplies s^i t^j
struct S {
    C c[6][6] = {};
    S() {}
    S(C a) { c[0][0] = a; }
};
const int DEG = 5;
S operator+(const S& a, const S& b) { S r; for (int i = 0; i <= DEG; i++) for (int j = 0; i + j <= DEG; j++) r.c[i][j] = a.c[i][j] + b.c[i][j]; return r; }
S operator-(const S& a, const S& b) { S r; for (int i = 0; i <= DEG; i++) for (int j = 0; i + j <= DEG; j++) r.c[i][j] = a.c[i][j] - b.c[i][j]; return r; }
S operator*(C k, const S& a) { S r; for (int i = 0; i <= DEG; i++) for (int j = 0; i + j <= DEG; j++) r.c[i][j] = k * a.c[i][j]; return r; }
S operator*(const S& a, const S& b) {
    S r;
    for (int i = 0; i <= DEG; i++) for (int j = 0; i + j <= DEG; j++)
        for (int k = 0; i + j + k <= DEG; k++) for (int l = 0; i + j + k + l <= DEG; l++)
            r.c[i + k][j + l] += a.c[i][j] * b.c[k][l];
    return r;
}
S var(int k) { S r; r.c[k == 0][k == 1] = 1; return r; }
S compose(const S& f, const S& a, const S& b) {  // f(a, b), a and b without constant term
    S pa[DEG + 1], pb[DEG + 1], r;
    pa[0] = pb[0] = S((R)1);
    for (int k = 1; k <= DEG; k++) { pa[k] = pa[k - 1] * a; pb[k] = pb[k - 1] * b; }
    for (int i = 0; i <= DEG; i++) for (int j = 0; i + j <= DEG; j++) r = r + f.c[i][j] * (pa[i] * pb[j]);
    return r;
}
typedef array<S, 2> Map;
Map compose(const Map& f, const Map& g) { return {compose(f[0], g[0], g[1]), compose(f[1], g[0], g[1])}; }
C coef(const S& a, int i, int j) { return a.c[i][j]; }
S deriv(const S& a, int k) {
    S r;
    for (int i = 0; i < DEG; i++) for (int j = 0; i + j < DEG; j++) r.c[i][j] = k == 0 ? (R)(i + 1) * a.c[i + 1][j] : (R)(j + 1) * a.c[i][j + 1];
    return r;
}
Twist normal_form(const Map& G);

Twist twist_at(const string& w, const V& p0) {
    C x0 = p0[0], y0 = p0[1], z0 = p0[2], k0 = kappa(p0);
    S s = var(0), t = var(1);
    C dk = (R)2 * z0 - x0 * y0;
    require(abs(dk) > 1e-8L, "dkappa/dz = 0");
    // a. z(s,t) on the level set
    S X = S(x0) + s, Y = S(y0) + t, Z = S(z0);
    for (int it = 0; it < 8; it++) {
        S kap = X * X + Y * Y + Z * Z - X * Y * Z - S((R)2 + k0);
        Z = Z - ((R)1 / dk) * kap;
    }
    // b. F(s,t)
    array<S, 3> v{X, Y, Z};
    for (int i = (int)w.size() - 1; i >= 0; i--)
        v = w[i] == 'L' ? array<S, 3>{v[2], v[1], v[2] * v[1] - v[0]} : array<S, 3>{v[0], v[2], v[0] * v[2] - v[1]};
    Map F{v[0] - S(x0), v[1] - S(y0)};
    F[0].c[0][0] = F[1].c[0][0] = 0;
    // c. g = 1/(2z - xy), U = int_0^s g, inverse s = Sinv(U,t)
    S P = (R)2 * Z - X * Y, Q = ((R)-1 / P.c[0][0]) * (P - S(P.c[0][0])), gs((R)1), pw((R)1);
    for (int k = 1; k <= DEG; k++) { pw = pw * Q; gs = gs + pw; }
    gs = ((R)1 / P.c[0][0]) * gs;
    S U;
    for (int i = 0; i < DEG; i++) for (int j = 0; i + 1 + j <= DEG; j++) U.c[i + 1][j] = gs.c[i][j] / (R)(i + 1);
    C g00 = U.c[1][0];
    S rest = U - g00 * s, Sinv = ((R)1 / g00) * s;
    for (int it = 0; it < 6; it++) Sinv = ((R)1 / g00) * (s - compose(rest, Sinv, t));
    // d. G = Phi o F o Phi^-1
    Map Fi = compose(F, Map{Sinv, t});
    Map G{compose(U, Fi[0], Fi[1]), Fi[1]};
    return normal_form(G);
}
// G: planar map with det-1 linear part, no constant term. mu, b as before; c from H(I) = I log mu + b I^2/2 + c I^3/3
Twist normal_form(const Map& G) {
    C a = G[0].c[1][0], b = G[0].c[0][1], c = G[1].c[1][0], d = G[1].c[0][1];
    Twist out;
    out.detlin = abs(a * d - b * c - (R)1);
    require(out.detlin < 1e-12L, "det of linear part != 1");
    // e. diagonalise
    C tr = a + d, l1 = (tr + sqrt(tr * tr - (R)4)) / (R)2;
    if (abs(l1) < 1) l1 = (R)1 / l1;
    C l2 = (R)1 / l1;
    auto ev = [&](C l) {
        array<C, 2> p{b, l - a}, q{l - d, c};
        return norm(p[0]) + norm(p[1]) >= norm(q[0]) + norm(q[1]) ? p : q;
    };
    auto e1 = ev(l1), e2 = ev(l2);
    C sq = sqrt(e1[0] * e2[1] - e2[0] * e1[1]);
    C M00 = e1[0] / sq, M10 = e1[1] / sq, M01 = e2[0] / sq, M11 = e2[1] / sq;
    S u = var(0), vv = var(1);
    Map Mm{M00 * u + M01 * vv, M10 * u + M11 * vv}, Mi{M11 * u - M01 * vv, -M10 * u + M00 * vv};
    Map Gp = compose(Mi, compose(G, Mm));
    // f. remove quadratic terms
    C lam[2] = {l1, l2};
    Map h2;
    for (int k = 0; k < 2; k++)
        for (int i = 0; i <= 2; i++) {
            int j = 2 - i;
            h2[k].c[i][j] = Gp[k].c[i][j] / (pow(l1, i) * pow(l2, j) - lam[k]);
        }
    Map id{u, vv}, h{u + h2[0], vv + h2[1]}, hi = id;
    for (int it = 0; it < 6; it++) { Map hk = compose(h2, hi); hi = {u - hk[0], vv - hk[1]}; }
    Map N = compose(hi, compose(Gp, h));
    out.quadmax = 0;
    for (int k = 0; k < 2; k++) for (int i = 0; i <= 2; i++) out.quadmax = max(out.quadmax, abs(N[k].c[i][2 - i]));
    require(out.quadmax < 1e-10L, "quadratic terms of N do not vanish");
    out.mu = l1;
    out.c1 = coef(N[0], 2, 1);
    out.c2 = coef(N[1], 1, 2);
    out.b = out.c1 / l1;
    // g. remove non-resonant cubic terms (u^2 v in N[0], u v^2 in N[1] stay); quartic terms do not reach u^3 v^2
    Map h3;
    for (int k = 0; k < 2; k++)
        for (int i = 0; i <= 3; i++) {
            int j = 3 - i;
            if (i - j == (k == 0 ? 1 : -1)) continue;
            h3[k].c[i][j] = N[k].c[i][j] / (pow(l1, i) * pow(l2, j) - lam[k]);
        }
    Map g3{u + h3[0], vv + h3[1]}, g3i = id;
    for (int it = 0; it < 3; it++) { Map hk = compose(h3, g3i); g3i = {u - hk[0], vv - hk[1]}; }
    Map N3 = compose(g3i, compose(N, g3));
    out.cubmax = 0;
    for (int k = 0; k < 2; k++) for (int i = 0; i <= 3; i++) if (i - (3 - i) != (k == 0 ? 1 : -1)) out.cubmax = max(out.cubmax, abs(N3[k].c[i][3 - i]));
    require(out.cubmax < 1e-10L, "non-resonant cubic terms of N do not vanish");
    // h. h o g3 is not symplectic: density j = det D(h o g3), action I = s + j1 s^2/2, s = uv
    Map ht = compose(h, g3);
    S jd = deriv(ht[0], 0) * deriv(ht[1], 1) - deriv(ht[0], 1) * deriv(ht[1], 0);
    out.jmax = max({abs(jd.c[1][0]), abs(jd.c[0][1]), abs(jd.c[2][0]), abs(jd.c[0][2])});
    require(out.jmax < 1e-9L, "density has degree-1 or u^2, v^2 terms");
    out.j1 = jd.c[1][1];
    out.c2p = coef(N3[0], 3, 2);
    out.b2 = (out.c2p - out.c1 * out.j1 / (R)2) / l1;
    out.c = out.b2 - out.b * out.b / (R)2;
    // N[1] = v / a(s), a(s) = mu + c1 s + c2p s^2: u^2 v^3 coefficient (b^2 - c2p/mu)/mu.
    // Rounding: observed rel <= 1.3e-10 for LR, LLR (selftest NF5 holds those to 1e-9), up to 2.2e-7 for LRLRR
    C want = (out.b * out.b - out.c2p / l1) / l1;
    out.consist = abs(coef(N3[1], 2, 3) - want) / abs(want);
    require(out.consist < 1e-6L, "u^2 v^3 coefficient of N[1] != (b^2 - c2p/mu)/mu");
    return out;
}
Twist twist(const string& w, R theta) { return twist_at(w, geometric(w, theta).best.v); }

// ---------- data ----------
struct Row { int r, c; R theta_file; string word; vector<R> L; };
R row_theta(int r, int c) { return fabsl(2 * PI - 4 * PI * (c + 1) / (R)r); }
vector<Row> read_raw(const vector<string>& files) {
    vector<Row> rows;
    for (auto& f : files) {
        ifstream in(f);
        require((bool)in, "cannot open " + f);
        string line;
        while (getline(in, line)) {
            istringstream ss(line);
            Row row; int n, it; string th;
            if (!(ss >> row.r >> row.c >> n >> th >> row.word >> it)) continue;
            row.theta_file = strtold(th.c_str(), nullptr);
            string s;
            while (ss >> s) row.L.push_back(strtold(s.c_str(), nullptr));
            require(row.L.size() == 8, "bad row in " + f);
            rows.push_back(row);
        }
    }
    return rows;
}
vector<vector<string>> read_csv(const string& f) {
    ifstream in(f);
    require((bool)in, "cannot open " + f);
    vector<vector<string>> out;
    string line;
    getline(in, line);
    while (getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.find_first_not_of(" \t") == string::npos) continue;
        vector<string> cells; string cell; istringstream ss(line);
        while (getline(ss, cell, ',')) cells.push_back(cell);
        out.push_back(cells);
    }
    return out;
}
R num_cell(const string& s, const string& f) {
    char* end;
    R x = strtold(s.c_str(), &end);
    require(end != s.c_str(), "bad number '" + s + "' in " + f);
    return x;
}
struct Ref { string word; int th; R alpha, gamma, logmu, vol; };
vector<Ref> read_fits() {
    const string f = "data/ladder_subleading_fits.csv";
    vector<Ref> out;
    for (auto& c : read_csv(f)) {
        require(c.size() == 18, "expected 18 cells per row in " + f);
        out.push_back({c[0], (int)num_cell(c[1], f), num_cell(c[2], f), num_cell(c[3], f), num_cell(c[12], f), num_cell(c[13], f)});
    }
    return out;
}

// least squares via normal equations (small, real)
vector<R> lsq(const vector<vector<R>>& X, const vector<R>& y) {
    int p = X[0].size();
    vector<vector<R>> A(p, vector<R>(p + 1, 0));
    for (size_t k = 0; k < X.size(); k++)
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) A[i][j] += X[k][i] * X[k][j];
            A[i][p] += X[k][i] * y[k];
        }
    for (int c = 0; c < p; c++) {
        int q = c;
        for (int r = c + 1; r < p; r++) if (fabsl(A[r][c]) > fabsl(A[q][c])) q = r;
        swap(A[q], A[c]);
        for (int r = 0; r < p; r++) if (r != c) {
            R f = A[r][c] / A[c][c];
            for (int k = c; k <= p; k++) A[r][k] -= f * A[c][k];
        }
    }
    vector<R> b(p);
    for (int i = 0; i < p; i++) b[i] = A[i][p] / A[i][i];
    return b;
}

// reciprocal condition estimate of a design: columns scaled to unit norm, column-pivoted Gram-Schmidt, min/max |R_kk|
R rcond(const vector<vector<R>>& X) {
    int m = X.size(), p = X[0].size();
    if (m < p) return 0;
    vector<vector<R>> Q(p, vector<R>(m));
    for (int j = 0; j < p; j++) {
        R s = 0;
        for (int i = 0; i < m; i++) s += X[i][j] * X[i][j];
        if (s == 0) return 0;
        for (int i = 0; i < m; i++) Q[j][i] = X[i][j] / sqrtl(s);
    }
    vector<bool> used(p, false);
    R dmax = 0, dmin = INFINITY;
    for (int k = 0; k < p; k++) {
        int q = -1; R best = -1;
        for (int j = 0; j < p; j++) if (!used[j]) {
            R n = 0;
            for (int i = 0; i < m; i++) n += Q[j][i] * Q[j][i];
            if (n > best) { best = n; q = j; }
        }
        used[q] = true;
        R d = sqrtl(best);
        dmax = max(dmax, d); dmin = min(dmin, d);
        if (d == 0) return 0;
        for (int i = 0; i < m; i++) Q[q][i] /= d;
        for (int j = 0; j < p; j++) if (!used[j]) {
            R dot = 0;
            for (int i = 0; i < m; i++) dot += Q[q][i] * Q[j][i];
            for (int i = 0; i < m; i++) Q[j][i] -= dot * Q[q][i];
        }
    }
    return dmin / dmax;
}

struct Resid { Row row; R theta, vol, logmu; vector<R> res; bool maxmu; };
vector<R> resid_vec(int r, const vector<R>& L, R vol, R logmu) {
    vector<R> res;
    for (int n = 0; n < 8; n++) res.push_back(r * (L[n] - r * vol / (4 * PI) + (n + 0.5L) * logmu));
    return res;
}
Resid resid_row(const Row& row) {
    Resid o; o.row = row; o.theta = row_theta(row.r, row.c);
    Geo g = geometric(row.word, o.theta);
    o.vol = g.best.vol; o.logmu = logl(abs(g.best.mu)); o.maxmu = g.maxmu;
    o.res = resid_vec(row.r, row.L, o.vol, o.logmu);
    return o;
}

struct Group { string word; int th; vector<Resid> rs; };
vector<Group> make_groups(const vector<string>& files) {
    map<pair<string, int>, vector<Resid>> m;
    for (auto& row : read_raw(files)) {
        Resid o = resid_row(row);
        m[{row.word, (int)lroundl(o.theta)}].push_back(o);
    }
    vector<Group> out;
    for (auto& [key, rs] : m) out.push_back({key.first, key.second, rs});
    return out;
}
// first stage regressors {1, 1/r} (+ 1/r^2) (+ theta - theta_nom unless theta_nom = 0, where theta = 2 pi/r)
vector<R> regs(int r, R theta, int th, bool rsq) {
    vector<R> x{1, 1 / (R)r};
    if (rsq) x.push_back(1 / ((R)r * r));
    if (th != 0) x.push_back(theta - th);
    return x;
}
// second stage b_n = alpha (+ lin (n+1/2)) + gamma (n+1/2)^2 over n = 0..nlev-1
struct Stage2 { R alpha, lin, gamma; };
vector<vector<R>> stage2_X(int nlev, bool lin) {
    vector<vector<R>> X;
    for (int n = 0; n < nlev; n++) {
        R h = n + 0.5L;
        X.push_back(lin ? vector<R>{1, h, h * h} : vector<R>{1, h * h});
    }
    return X;
}
Stage2 stage2(const vector<R>& b, int nlev, bool lin) {
    auto p = lsq(stage2_X(nlev, lin), vector<R>(b.begin(), b.begin() + nlev));
    return lin ? Stage2{p[0], p[1], p[2]} : Stage2{p[0], NAN, p[1]};
}
// a design with rcond() below this is reported as singular (observed: singular designs <= 3.3e-17, all others >= 1.6e-2)
const R RCOND_MIN = 1e-8L;
struct Variant { string name; int nrows, nreg; bool ok, exact, singular = false; Stage2 s; vector<R> b; };
Variant fit_variant(const Group& g, const string& name, int drop_r, bool rsq, int nlev, bool lin) {
    Variant v; v.name = name;
    vector<const Resid*> rs;
    for (auto& o : g.rs) if (o.row.r != drop_r) rs.push_back(&o);
    v.nrows = rs.size(); v.nreg = regs(1, 0, g.th, rsq).size();
    v.ok = v.nrows >= v.nreg; v.exact = v.nrows == v.nreg;
    v.s = {NAN, NAN, NAN};
    if (!v.ok) { v.b.assign(8, NAN); return v; }
    vector<vector<R>> X;
    for (auto o : rs) X.push_back(regs(o->row.r, o->theta, g.th, rsq));
    v.singular = rcond(X) < RCOND_MIN || rcond(stage2_X(nlev, lin)) < RCOND_MIN;
    for (int n = 0; n < 8; n++) {
        vector<R> y;
        for (auto o : rs) y.push_back(o->res[n]);
        v.b.push_back(lsq(X, y)[0]);
    }
    v.s = stage2(v.b, nlev, lin);
    return v;
}

struct FitOut { string word; int th; R alpha, gamma, rms; vector<R> b; int nrows; };
vector<FitOut> fit(const vector<string>& files) {
    vector<FitOut> out;
    for (auto& g : make_groups(files)) {
        Variant v = fit_variant(g, "base", 0, false, 8, false);
        FitOut f; f.word = g.word; f.th = g.th; f.nrows = v.nrows; f.b = v.b; f.alpha = v.s.alpha; f.gamma = v.s.gamma;
        R ss = 0;
        for (int n = 0; n < 8; n++) { R h = n + 0.5L, e = f.b[n] - f.alpha - f.gamma * (h * h); ss += e * e; }
        f.rms = sqrtl(ss / 8);
        out.push_back(f);
    }
    return out;
}

vector<Variant> robust_variants(const Group& g) {
    vector<int> rs;
    for (auto& o : g.rs) if (find(rs.begin(), rs.end(), o.row.r) == rs.end()) rs.push_back(o.row.r);
    sort(rs.begin(), rs.end());
    vector<Variant> vs{fit_variant(g, "base", 0, false, 8, false)};
    for (int r : rs) vs.push_back(fit_variant(g, "drop_r" + to_string(r), r, false, 8, false));
    vs.push_back(fit_variant(g, "n0-3", 0, false, 4, false));
    vs.push_back(fit_variant(g, "n0-5", 0, false, 6, false));
    vs.push_back(fit_variant(g, "lin", 0, false, 8, true));
    vs.push_back(fit_variant(g, "rsq", 0, true, 8, false));
    return vs;
}
string fmt(const char* f, ...) __attribute__((format(printf, 1, 2)));
string robust_line(const Group& g, const Variant& v, R imb) {
    string s = fmt("%s,%d,%s%s,%d,", g.word.c_str(), g.th, v.name.c_str(), v.exact ? "*" : "", v.nrows);
    if (!v.ok || v.singular) return s + fmt("n/a,n/a,n/a,%.10Lf,n/a,n/a,%s", imb, v.ok ? "singular" : "rows<regressors");
    s += fmt("%.8Lf,%.10Lf,", v.s.alpha, v.s.gamma);
    if (!isnan(v.s.lin)) s += fmt("%.8Lf", v.s.lin);
    return s + fmt(",%.10Lf,%.8Lf,%+.3Le,", imb, v.s.gamma / imb, v.s.gamma / (PI * imb) - 1);
}

// ---------- gaps r = 201..1001 ----------
struct GapRow { string word; int r, c; R theta, logmu, pimb; vector<R> d; };
// row: word r c theta | S Vol S-Vol | d1..d5 | logmu ; theta and logmu are recomputed, not read; pimb = pi Im b(theta)
vector<GapRow> read_gaps(const string& f) {
    ifstream in(f);
    require((bool)in, "cannot open " + f);
    vector<GapRow> out;
    string line;
    getline(in, line);
    while (getline(in, line)) {
        istringstream ss(line);
        vector<string> t; string s;
        while (ss >> s) if (s != "|") t.push_back(s);
        if (t.empty()) continue;
        require(t.size() == 13, "bad row in " + f);
        GapRow g; g.word = t[0]; g.r = atoi(t[1].c_str()); g.c = atoi(t[2].c_str());
        for (int k = 7; k < 12; k++) g.d.push_back(strtold(t[k].c_str(), nullptr));
        g.theta = row_theta(g.r, g.c);
        Geo geo = geometric(g.word, g.theta);
        g.logmu = logl(abs(geo.best.mu));
        g.pimb = PI * twist_at(g.word, geo.best.v).b.imag();
        out.push_back(g);
    }
    return out;
}
// y_n = r (log|mu| - delta_n) / (2 (n+1)); intercept of the fit over r is gamma_n
// direct: w_n = y_n - pi Im b(theta_row), regressors {1, 1/r} (+ 1/r^2) only; the intercept should be 0
// bound_n = sum_k |w_k| r_k 5e-6 / (2 (n+1)), w = intercept weights: worst case of the 5-decimal rounding
struct GapFit { string word; int th, nrows, rmax; vector<R> g, bound; R mean, spread, bmean; bool singular = false; };
vector<GapFit> gap_fits(const vector<GapRow>& rows, bool rsq, bool direct) {
    map<pair<string, int>, vector<const GapRow*>> m;
    for (auto& g : rows) m[{g.word, (int)lroundl(g.theta)}].push_back(&g);
    vector<GapFit> out;
    for (auto& [key, rs] : m) {
        GapFit f; f.word = key.first; f.th = key.second; f.nrows = rs.size(); f.rmax = 0;
        vector<vector<R>> X;
        for (auto p : rs) { X.push_back(regs(p->r, p->theta, direct ? 0 : f.th, rsq)); f.rmax = max(f.rmax, p->r); }
        f.singular = rcond(X) < RCOND_MIN;
        vector<R> w;
        for (size_t k = 0; k < rs.size(); k++) { vector<R> e(rs.size(), 0); e[k] = 1; w.push_back(lsq(X, e)[0]); }
        for (int n = 0; n < 5; n++) {
            vector<R> y; R bd = 0;
            for (size_t k = 0; k < rs.size(); k++) {
                y.push_back(rs[k]->r * (rs[k]->logmu - rs[k]->d[n]) / (2 * (n + 1)) - (direct ? rs[k]->pimb : 0));
                bd += fabsl(w[k]) * rs[k]->r * 5e-6L / (2 * (n + 1));
            }
            f.g.push_back(lsq(X, y)[0]); f.bound.push_back(bd);
        }
        R s = 0, sb = 0;
        for (int n = 0; n < 5; n++) { s += f.g[n]; sb += f.bound[n]; }
        f.mean = s / 5; f.bmean = sb / 5;
        f.spread = *max_element(f.g.begin(), f.g.end()) - *min_element(f.g.begin(), f.g.end());
        out.push_back(f);
    }
    return out;
}

// ---------- direct fits ----------
// per raw row: gamma_row from the second stage over n of res[n] alone (no fit over r, no theta regressor)
const vector<tuple<string, int, bool>> STAGES{{"n0-7", 8, false}, {"n0-3", 4, false}, {"n0-5", 6, false}, {"lin", 8, true}};
struct DirectRow { string word; int th, r, c; R theta, pimb; vector<R> gam; };
vector<DirectRow> direct_rows(const Group& g) {
    vector<DirectRow> out;
    for (auto& o : g.rs) {
        DirectRow d{g.word, g.th, o.row.r, o.row.c, o.theta, PI * twist(g.word, o.theta).b.imag(), {}};
        for (auto& [name, nlev, lin] : STAGES) d.gam.push_back(stage2(o.res, nlev, lin).gamma);
        out.push_back(d);
    }
    return out;
}
// z_row = gamma_row - pi Im b(theta_row) fitted over the rows with {1, 1/r} (+ 1/r^2); claim: intercept z0 = 0
struct DirectFit { string word; int th; string stage; bool rsq; int nrows; bool singular; R z0, rel, rms; };
vector<DirectFit> direct_fits(const Group& g, const vector<DirectRow>& rows) {
    R pnom = PI * twist(g.word, g.th).b.imag();
    vector<DirectFit> out;
    for (size_t s = 0; s < STAGES.size(); s++) for (bool rsq : {false, true}) {
        if (rsq && rows.size() <= 3) continue;
        vector<vector<R>> X; vector<R> y;
        for (auto& d : rows) { X.push_back(regs(d.r, 0, 0, rsq)); y.push_back(d.gam[s] - d.pimb); }
        auto p = lsq(X, y);
        R ss = 0;
        for (size_t k = 0; k < X.size(); k++) {
            R e = y[k];
            for (size_t j = 0; j < p.size(); j++) e -= p[j] * X[k][j];
            ss += e * e;
        }
        out.push_back({g.word, g.th, get<0>(STAGES[s]), rsq, (int)rows.size(), rcond(X) < RCOND_MIN, p[0], p[0] / pnom, sqrtl(ss / X.size())});
    }
    return out;
}

// ---------- conjugate choice by the signed sum ----------
// rule S: among the max-Vol non-real points take one whose signed Bloch-Wigner sum has the sign it has at the
// present LR theta = 0 point (the one calibration case); ties among those by the larger Im x
int sign_S() {
    static int s = volume_signed("LR", geometric("LR", 0).best.v, 0, 1) > 0 ? 1 : -1;
    return s;
}
Cand select_S(const string& w, R theta) {
    Geo g = geometric(w, theta);
    int ib = -1;
    for (int i = 0; i < (int)g.all.size(); i++) {
        auto& c = g.all[i];
        if (!nonreal(c.v) || c.vol < g.best.vol - 1e-12L || (volume_signed(w, c.v, theta, 1) > 0 ? 1 : -1) != sign_S()) continue;
        if (ib < 0 || c.v[0].imag() > g.all[ib].v[0].imag()) ib = i;
    }
    require(ib >= 0, "rule S: no point with the calibrated sign for " + w);
    return g.all[ib];
}

// ---------- stress words ----------
struct StressRow { string word; R theta, vol_f, logmu_f, S, S_err, d0, d0_err, d1, d1_err; Geo g; };
vector<StressRow> read_stress(const string& f) {
    vector<StressRow> out;
    for (auto& c : read_csv(f)) {
        require(c.size() == 22, "expected 22 cells per row in " + f);
        auto num = [&](int i) { return num_cell(c[i], f); };
        StressRow s{c[0], num(1), num(2), num(3), num(4), num(5), num(6), num(7), num(8), num(9), {}};
        s.g = geometric(s.word, s.theta);
        out.push_back(s);
    }
    return out;
}

// ---------- order 2 ----------
// Bohr-Sommerfeld with hbar = 2 pi i/r: the c I^3/3 term of H gives d3_pred m^3/r^2 in L_n, d3_pred = (4 pi^2/3) Re c
R d3_pred(const Twist& t) { return 4 * PI * PI / 3 * t.c.real(); }
// per row: pimb = pi Im b, d3p = d3_pred, both at the row's exact theta and the point used by direct_rows
struct O2Row { string word; int th, r; R theta, pimb, d3p; vector<R> L, res; };
O2Row o2_row(const Resid& o, int th) {
    Twist t = twist(o.row.word, o.theta);
    return {o.row.word, th, o.row.r, o.theta, PI * t.b.imag(), d3_pred(t), o.row.L, o.res};
}
vector<O2Row> o2_rows(const Group& g) {
    vector<O2Row> out;
    for (auto& o : g.rs) out.push_back(o2_row(o, g.th));
    return out;
}
// k-th forward differences of L, n = 0..7-k
vector<R> fdiff(vector<R> L, int k) {
    for (int i = 0; i < k; i++) for (size_t n = 0; n + 1 < L.size() - i; n++) L[n] = L[n + 1] - L[n];
    L.resize(L.size() - k);
    return L;
}
// T_n = r^2 Delta^3 L_n -> 6 d3, n = 0..4
vector<R> o2_T(const O2Row& o) { auto d = fdiff(o.L, 3); for (auto& x : d) x *= (R)o.r * o.r; return d; }
// W_n = r (r Delta^2 L_n - 2 pi Im b) -> 2 d2 + 6 d3 (n + 3/2), n = 0..5
vector<R> o2_W(const O2Row& o) { auto d = fdiff(o.L, 2); for (auto& x : d) x = o.r * (o.r * x - 2 * o.pimb); return d; }
// fit y over the rows by {1, 1/r} (+ 1/r^2)
struct RFit { bool ok, exact, singular; vector<R> p; R rms; };
RFit fit_r(const vector<O2Row>& rows, const vector<R>& y, bool rsq) {
    RFit f;
    size_t nreg = rsq ? 3 : 2;
    f.ok = rows.size() >= nreg; f.exact = rows.size() == nreg; f.singular = false; f.rms = NAN;
    if (!f.ok) { f.p.assign(nreg, NAN); return f; }
    vector<vector<R>> X;
    for (auto& o : rows) X.push_back(regs(o.r, 0, 0, rsq));
    f.singular = rcond(X) < RCOND_MIN;
    f.p = lsq(X, y);
    R ss = 0;
    for (size_t k = 0; k < X.size(); k++) {
        R e = y[k];
        for (size_t j = 0; j < nreg; j++) e -= f.p[j] * X[k][j];
        ss += e * e;
    }
    f.rms = sqrtl(ss / X.size());
    return f;
}
// B4: y_n = res[n] - pi Im b m^2 over the rows; [rsq][n], intercept alpha, 1/r slope D(n)
array<vector<RFit>, 2> b4_fits(const vector<O2Row>& rows) {
    array<vector<RFit>, 2> out;
    for (int rsq = 0; rsq < 2; rsq++)
        for (int n = 0; n < 8; n++) {
            vector<R> y;
            for (auto& o : rows) y.push_back(o.res[n] - o.pimb * (n + 0.5L) * (n + 0.5L));
            out[rsq].push_back(fit_r(rows, y, rsq));
        }
    return out;
}
// least squares of D(n) in powers of m = n + 1/2
struct PolyFit { vector<R> p; R rms; };
PolyFit poly_m(const vector<R>& D, const vector<int>& pows) {
    vector<vector<R>> X;
    for (size_t n = 0; n < D.size(); n++) {
        vector<R> x;
        for (int k : pows) x.push_back(powl(n + 0.5L, k));
        X.push_back(x);
    }
    PolyFit f{lsq(X, D), 0};
    for (size_t n = 0; n < D.size(); n++) {
        R e = D[n];
        for (size_t j = 0; j < pows.size(); j++) e -= f.p[j] * X[n][j];
        f.rms += e * e;
    }
    f.rms = sqrtl(f.rms / D.size());
    return f;
}
struct Spread { R mean, spread; };
Spread spread_of(const vector<R>& x) {
    R s = 0;
    for (R v : x) s += v;
    return {s / x.size(), *max_element(x.begin(), x.end()) - *min_element(x.begin(), x.end())};
}

// ---------- commands ----------
string cs(C z) { char b[96]; snprintf(b, sizeof b, "%.15Lg%+.15Lgi", z.real(), z.imag()); return b; }

int fails = 0;
void check(const string& name, bool ok, const string& detail) {
    printf("%s %s %s\n", ok ? "PASS" : "FAIL", name.c_str(), detail.c_str());
    if (!ok) fails++;
}
string fmt(const char* f, ...) __attribute__((format(printf, 1, 2)));
string fmt(const char* f, ...) {
    char b[512]; va_list a; va_start(a, f); vsnprintf(b, sizeof b, f, a); va_end(a); return b;
}

int selftest() {
    const vector<string> words{"LR", "LLR"};
    auto guard = [](const string& name, auto fn) {
        try { fn(); } catch (exception& e) { check(name, false, string("exception: ") + e.what()); }
    };
    check("Lob", fabsl(6 * lob(PI / 3) - 2.029883212819307L) < 1e-15L, fmt("6Lob(pi/3)=%.18Lf", 6 * lob(PI / 3)));
    guard("T1", [&] {
        Geo g = geometric("LR", 0);
        C s = g.best.mu + (R)1 / g.best.mu;
        C k = kappa(g.best.v);
        check("T1", abs(s - (R)5) < 1e-15L && abs(k + (R)2) < 1e-15L,
              "mu+1/mu=" + cs(s) + " kappa=" + cs(k));
    });
    auto refs = read_fits();
    for (auto& rf : refs) guard("T2/T4 " + rf.word, [&] {
        Geo g = geometric(rf.word, rf.th);
        R lm = logl(abs(g.best.mu));
        check("T2", fabsl(lm - rf.logmu) < 2e-5L, fmt("%s th=%d logmu=%.10Lf ref=%.5Lf diff=%.2Le", rf.word.c_str(), rf.th, lm, rf.logmu, lm - rf.logmu));
        check("T4", fabsl(g.best.vol - rf.vol) < 6e-6L, fmt("%s th=%d Vol=%.10Lf ref=%.5Lf diff=%.2Le", rf.word.c_str(), rf.th, g.best.vol, rf.vol, g.best.vol - rf.vol));
        printf("INFO %s th=%d max-Vol point is max-|mu| point: %s\n", rf.word.c_str(), rf.th, g.maxmu ? "yes" : "no");
    });
    for (auto& c : read_csv("data/extrap_r121.csv")) guard("T3 " + c[0], [&] {
        R th = strtold(c[1].c_str(), nullptr), ref = strtold(c[2].c_str(), nullptr);
        R v = geometric(c[0], th).best.vol;
        check("T3", fabsl(v - ref) < 1e-9L, fmt("%s th=%.2Lf Vol=%.13Lf ref=%.13Lf diff=%.2Le", c[0].c_str(), th, v, ref, v - ref));
    });
    for (auto& row : read_raw({"data/r2001_raw.txt"})) {
        R th = row_theta(row.r, row.c);
        check("T5", fabsl(th - row.theta_file) < 1e-9L, fmt("r=%d c=%d theta=%.12Lf file=%.10Lf diff=%.2Le", row.r, row.c, th, row.theta_file, th - row.theta_file));
    }
    for (auto [w, th] : vector<pair<string, int>>{{"LR", 1}, {"LLR", 2}}) guard("T6", [&] {
        Geo g = geometric(w, th);
        R a = volume(w, g.best.v, th, 1), b = volume(w, g.best.v, th, -1);
        check("T6", fabsl(a - b) < 1e-10L, fmt("%s th=%d Vol(+)=%.15Lf Vol(-)=%.15Lf", w.c_str(), th, a, b));
    });
    for (auto& w : words) for (int th = 0; th <= 3; th++) guard("T7", [&] {
        Twist t = twist(w, th);
        C s = t.c1 / t.mu + t.c2 * t.mu;
        R rel = abs(s) / max(abs(t.c1 / t.mu), abs(t.c2 * t.mu));
        check("T7", rel < 1e-9L, fmt("%s th=%d c1/mu+c2*mu rel=%.2Le b=", w.c_str(), th, rel) + cs(t.b));
    });
    guard("T8", [&] {
        auto fs = fit({"data/r2001_raw.txt", "data/r3001-5001_raw.txt", "data/r10001_raw.txt"});
        for (auto& rf : refs) {
            bool found = false;
            for (auto& f : fs) if (f.word == rf.word && f.th == rf.th) {
                found = true;
                check("T8", fabsl(f.alpha - rf.alpha) < 5e-3L && fabsl(f.gamma - rf.gamma) < 5e-3L,
                      fmt("%s th=%d alpha=%.5Lf ref=%.5Lf gamma=%.5Lf ref=%.5Lf", rf.word.c_str(), rf.th, f.alpha, rf.alpha, f.gamma, rf.gamma));
            }
            if (!found) check("T8", false, rf.word + " group missing");
        }
    });
    const vector<string> raw3{"data/r2001_raw.txt", "data/r3001-5001_raw.txt", "data/r10001_raw.txt"};
    vector<Group> gs;
    guard("R1", [&] {
        gs = make_groups(raw3);
        auto fs = fit(raw3);
        check("R1", gs.size() == 8 && fs.size() == 8, fmt("groups=%zu fits=%zu", gs.size(), fs.size()));
        for (size_t i = 0; i < gs.size() && i < fs.size(); i++) {
            auto vs = robust_variants(gs[i]);
            bool ok = !vs.empty() && vs[0].name == "base" && vs[0].s.alpha == fs[i].alpha && vs[0].s.gamma == fs[i].gamma;
            check("R1", ok, fmt("%s th=%d base alpha=%.12Lf gamma=%.12Lf fit alpha=%.12Lf gamma=%.12Lf", gs[i].word.c_str(), gs[i].th,
                                vs.empty() ? NAN : vs[0].s.alpha, vs.empty() ? NAN : vs[0].s.gamma, fs[i].alpha, fs[i].gamma));
        }
    });
    guard("R2", [&] {
        check("R2", gs.size() == 8, fmt("groups=%zu", gs.size()));
        for (auto& g : gs) {
            auto vs = robust_variants(g);
            vector<int> rs;
            for (auto& o : g.rs) if (find(rs.begin(), rs.end(), o.row.r) == rs.end()) rs.push_back(o.row.r);
            int ndrop = 0; bool ok = !vs.empty();
            for (auto& v : vs) {
                if (v.name.rfind("drop_r", 0) == 0) { ndrop++; ok = ok && v.nrows == vs[0].nrows - 1; }
                ok = ok && v.ok == (v.nrows >= v.nreg) && v.exact == (v.nrows == v.nreg);
                if (!v.ok) ok = ok && robust_line(g, v, 1).find(",n/a,n/a,n/a,") != string::npos;
            }
            check("R2", ok && ndrop == (int)rs.size() && vs.size() == rs.size() + 5,
                  fmt("%s th=%d variants=%zu drop=%d distinct r=%zu", g.word.c_str(), g.th, vs.size(), ndrop, rs.size()));
            if (g.word == "LLR" && g.th == 3) {
                // 3 rows, 3 regressors: base, n0-3, n0-5, lin interpolate exactly; drop-one (2 rows) and rsq (4 regressors) are n/a
                string got, want = "base*:3 drop_r2001:n/a drop_r3001:n/a drop_r10001:n/a n0-3*:3 n0-5*:3 lin*:3 rsq:n/a ";
                for (auto& v : vs) got += v.name + (v.exact ? "*" : "") + ":" + (v.ok ? to_string(v.nrows) : "n/a") + " ";
                check("R2", got == want, "LLR th=3 " + got);
            }
        }
    });
    guard("R3", [&] {
        const R A = 0.44621L, Lc = -0.0371L, G = 1.38195L;
        vector<R> b, b4, b6;
        for (int n = 0; n < 8; n++) {
            R h = n + 0.5L;
            b.push_back(A + Lc * h + G * h * h);
            b4.push_back(n < 4 ? A + G * h * h : 1000 * (n + 1));
            b6.push_back(n < 6 ? A + G * h * h : -1000 * (n + 1));
        }
        Stage2 s = stage2(b, 8, true), s4 = stage2(b4, 4, false), s6 = stage2(b6, 6, false);
        check("R3", fabsl(s.alpha - A) < 1e-12L && fabsl(s.lin - Lc) < 1e-12L && fabsl(s.gamma - G) < 1e-12L,
              fmt("lin: alpha=%.15Lf lin=%.15Lf gamma=%.15Lf", s.alpha, s.lin, s.gamma));
        check("R3", fabsl(s4.alpha - A) < 1e-12L && fabsl(s4.gamma - G) < 1e-12L, fmt("n0-3: alpha=%.15Lf gamma=%.15Lf", s4.alpha, s4.gamma));
        check("R3", fabsl(s6.alpha - A) < 1e-12L && fabsl(s6.gamma - G) < 1e-12L, fmt("n0-5: alpha=%.15Lf gamma=%.15Lf", s6.alpha, s6.gamma));
    });
    guard("R4", [&] {
        // file prints Vol at full double precision (17 digits): tolerance 1e-12; logmu at 8 decimals: half unit 5e-9 (+1e-13)
        auto ss = read_stress("data/stress_7words_r101.csv");
        check("R4", ss.size() == 20, fmt("rows=%zu", ss.size()));
        for (auto& s : ss) {
            R dv = s.g.best.vol - s.vol_f, dl = logl(abs(s.g.best.mu)) - s.logmu_f;
            check("R4", fabsl(dv) < 1e-12L && fabsl(dl) < 5e-9L + 1e-13L,
                  fmt("%s th=%.0Lf Vol=%.13Lf file=%.13Lf diff=%.2Le logmu=%.10Lf file=%.8Lf diff=%.2Le", s.word.c_str(), s.theta,
                      s.g.best.vol, s.vol_f, dv, logl(abs(s.g.best.mu)), s.logmu_f, dl));
        }
    });
    vector<GapRow> gaprows;
    guard("R5", [&] {
        // direct block {1, 1/r, 1/r^2}: |mean intercept of w_n| <= worst-case propagation of the 5e-6 gap rounding.
        // The bound covers rounding of the gaps only, not truncation; the two-term {1, 1/r} block exceeds its own
        // bound in 4 of 8 cases (LLR th=1, LLR th=3, LR th=0, LR th=2).
        gaprows = read_gaps("data/r201-1001_summary.txt");
        auto gf = gap_fits(gaprows, true, true);
        check("R5", gf.size() == 8, fmt("cases=%zu", gf.size()));
        for (auto& f : gf)
            check("R5", !f.singular && fabsl(f.mean) <= f.bmean, fmt("%s th=%d mean w0=%+.3Le bound=%.3Le", f.word.c_str(), f.th, f.mean, f.bmean));
    });
    guard("A1", [&] {
        // rcond < 1e-8 is singular: exactly those two drop_r10001 rows in robust; theta-regressor gaps blocks at theta_nom 1, 2
        int nsing = 0;
        for (auto& g : gs) for (auto& v : robust_variants(g)) if (v.singular) {
            nsing++;
            check("A1", g.th == 2 && v.name == "drop_r10001", fmt("%s th=%d %s singular", g.word.c_str(), g.th, v.name.c_str()));
        }
        check("A1", nsing == 2, fmt("singular robust variants=%d", nsing));
        for (bool rsq : {false, true}) {
            for (auto& f : gap_fits(gaprows, rsq, false))
                check("A1", f.singular == (f.th == 1 || f.th == 2), fmt("gaps theta block rsq=%d %s th=%d singular=%d", rsq, f.word.c_str(), f.th, f.singular));
            for (auto& f : gap_fits(gaprows, rsq, true))
                check("A1", !f.singular, fmt("gaps direct block rsq=%d %s th=%d singular=%d", rsq, f.word.c_str(), f.th, f.singular));
        }
    });
    guard("R6", [&] {
        // twist at each row's exact (non-integer) theta: same identity as T7; direct z fits present for every case
        for (auto& g : gs) {
            R worst = 0;
            for (auto& o : g.rs) {
                Twist t = twist(g.word, o.theta);
                worst = max(worst, abs(t.c1 / t.mu + t.c2 * t.mu) / max(abs(t.c1 / t.mu), abs(t.c2 * t.mu)));
            }
            auto dr = direct_rows(g);
            auto df = direct_fits(g, dr);
            size_t want = 4 * (g.rs.size() > 3 ? 2 : 1);
            bool ok = dr.size() == g.rs.size() && df.size() == want;
            for (auto& f : df) ok = ok && !f.singular && !isnan(f.z0) && f.nrows == (int)g.rs.size();
            check("R6", worst < 1e-9L && ok, fmt("%s th=%d rows=%zu fits=%zu (want %zu) twist identity worst rel=%.2Le", g.word.c_str(), g.th, dr.size(), df.size(), want, worst));
        }
    });
    guard("C1", [&] {
        // LR/LLR selected points from before the damped search (pass 1 output, 19 digits)
        struct P { string w; int th; V v; };
        vector<P> old{
            {"LR", 0, {C(1.500000000000000000e+00L, 8.660254037844386467e-01L), C(1.500000000000000000e+00L, -8.660254037844386467e-01L), C(1.500000000000000000e+00L, -8.660254037844386467e-01L)}},
            {"LR", 1, {C(1.539752315004402910e+00L, 8.418238761447597556e-01L), C(1.539752315004402910e+00L, -8.418238761447597556e-01L), C(1.539752315004402910e+00L, -8.418238761447597556e-01L)}},
            {"LR", 2, {C(1.640139790744088769e+00L, 7.682584515032127150e-01L), C(1.640139790744088769e+00L, -7.682584515032127150e-01L), C(1.640139790744088769e+00L, -7.682584515032127150e-01L)}},
            {"LR", 3, {C(1.763474912943654736e+00L, 6.458374852125562577e-01L), C(1.763474912943654736e+00L, -6.458374852125562576e-01L), C(1.763474912943654736e+00L, -6.458374852125562576e-01L)}},
            {"LLR", 0, {C(1.500000000000000000e+00L, 1.322875655532295295e+00L), C(1.632241882311900196e+00L, -4.052327261871812948e-01L), C(1.632241882311900196e+00L, -4.052327261871812948e-01L)}},
            {"LLR", 1, {C(1.531029514051922868e+00L, 1.265362461179961460e+00L), C(-1.652092229255093930e+00L, 4.067235478511456275e-01L), C(-1.652092229255093930e+00L, 4.067235478511456275e-01L)}},
            {"LLR", 2, {C(1.619298906621570624e+00L, 1.109593627849537530e+00L), C(1.711232389745241426e+00L, -4.015644658431123545e-01L), C(1.711232389745241426e+00L, -4.015644658431123545e-01L)}},
            {"LLR", 3, {C(1.742513147389102065e+00L, 8.918810040988747510e-01L), C(-1.799468636603193617e+00L, 3.680160676209880816e-01L), C(-1.799468636603193617e+00L, 3.680160676209880816e-01L)}}};
        for (auto& p : old) {
            R d = vdist(geometric(p.w, p.th).best.v, p.v);
            check("C1", d < 1e-17L, fmt("%s th=%d selected point unchanged, dist=%.2Le", p.w.c_str(), p.th, d));
        }
    });
    const vector<string> long_words{"LLLR", "LLRR", "LLLRR", "LRLRR", "LLLLR"};
    vector<string> all_words{"LR", "LLR"};
    all_words.insert(all_words.end(), long_words.begin(), long_words.end());
    int d3_nontrivial = 0, d3_cases = 0;
    for (auto& w : all_words) for (int th = 0; th <= 3; th++) guard("D1/D3", [&] {
        // D1: signed sum independent of sgn, of the root t vs 1/t, of p1/p2; odd under conjugation
        // D3: the top-Vol points with the calibrated sign of the signed sum agree in mu and b
        Geo g = geometric(w, th);
        Cand s = select_S(w, th);
        Twist ts = twist_at(w, s.v);
        R dev = 0, odd = 0, dmu = 0, db = 0;
        int nsame = 0;
        for (auto& c : g.all) {
            if (!nonreal(c.v) || c.vol < g.best.vol - 1e-12L) continue;
            R S0 = volume_signed(w, c.v, th, 1);
            for (int sg : {1, -1}) for (bool tf : {false, true}) for (int pm : {0, 1, 2}) dev = max(dev, fabsl(volume_signed(w, c.v, th, sg, nullptr, tf, pm) - S0));
            odd = max(odd, fabsl(S0 + volume_signed(w, V{conj(c.v[0]), conj(c.v[1]), conj(c.v[2])}, th, 1)));
            if ((S0 > 0) == (volume_signed(w, s.v, th, 1) > 0)) {
                nsame++;
                Twist t = twist_at(w, c.v);
                dmu = max(dmu, abs(t.mu - ts.mu)); db = max(db, abs(t.b - ts.b));
            }
        }
        check("D1", dev < 1e-12L && odd < 1e-12L, fmt("%s th=%d branch dev=%.2Le |S(v)+S(conj v)|=%.2Le", w.c_str(), th, dev, odd));
        d3_cases++;
        if (nsame < 2) printf("INFO D3 %s th=%d only %d point with the calibrated sign: nothing to compare\n", w.c_str(), th, nsame);
        else {
            d3_nontrivial++;
            check("D3", dmu < 1e-9L && db < 1e-9L, fmt("%s th=%d %d points with the calibrated sign: max|dmu|=%.2Le max|db|=%.2Le", w.c_str(), th, nsame, dmu, db));
        }
    });
    printf("INFO D3 non-trivial cases (>= 2 points with the calibrated sign): %d of %d\n", d3_nontrivial, d3_cases);
    for (auto& w : words) for (int th = 0; th <= 3; th++) guard("D2", [&] {
        // rule S (sign calibrated on LR theta 0) selects the present Im x > 0 point
        R d = vdist(select_S(w, th).v, geometric(w, th).best.v);
        check("D2", d < 1e-17L, fmt("%s th=%d%s rule S vs Im x > 0 point dist=%.2Le", w.c_str(), th, w == "LR" && th == 0 ? " (calibration case)" : "", d));
    });
    for (auto& w : long_words) guard("D4", [&] {
        // rule S: Im(b) > 0 at theta = 0..3
        string s; int npos = 0;
        for (int th = 0; th <= 3; th++) { R ib = twist_at(w, select_S(w, th).v).b.imag(); npos += ib > 0; s += fmt(" %+.6Lf", ib); }
        check("D4", npos == 4, w + " Im b at theta 0..3:" + s);
    });
    guard("NF1", [&] {
        // G = A o s2 o s1 o N_H o s1^-1 o s2^-1 o A^-1, N_H = (u E(uv), v / E(uv)), E = exp(lam + beta s + cc s^2)
        S u = var(0), v = var(1), s = u * v, one((R)1);
        C lam(0.9L, 0.4L), beta(0.31L, -0.52L), cc(-0.27L, 0.18L), em = exp(lam);
        C a1(0.3L, 0.1L), a2(-0.2L, 0.35L), e1(0.4L, -0.2L), e2(0.15L, 0.25L);
        C A00(1.2L, 0.3L), A01(0.4L, -0.1L), A10(-0.3L, 0.2L), A11 = ((R)1 + A01 * A10) / A00;
        Map NH{u * (em * (one + beta * s + (cc + beta * beta / (R)2) * (s * s))),
               v * (((R)1 / em) * (one - beta * s + (beta * beta / (R)2 - cc) * (s * s)))};
        Map s1{u, v + a1 * (u * u) + a2 * (u * u * u)}, s1i{u, v - a1 * (u * u) - a2 * (u * u * u)};
        Map s2{u + e1 * (v * v) + e2 * (v * v * v), v}, s2i{u - e1 * (v * v) - e2 * (v * v * v), v};
        Map Am{A00 * u + A01 * v, A10 * u + A11 * v}, Ai{A11 * u - A01 * v, -A10 * u + A00 * v};
        Map G = compose(Am, compose(s2, compose(s1, compose(NH, compose(s1i, compose(s2i, Ai))))));
        Twist t = normal_form(G);
        R dmu = abs(t.mu - em), db = abs(t.b - beta), dc = abs(t.c - cc);
        check("NF1", dmu < 1e-10L && db < 1e-10L && dc < 1e-10L,
              fmt("synthetic |mu-e^lam|=%.2Le |b-beta|=%.2Le |c-cc|=%.2Le j1=", dmu, db, dc) + cs(t.j1));
        C cnoj = t.c2p / t.mu - t.b * t.b / (R)2;
        printf("INFO NF1 j1=%s c without j1 correction=%s |c_noj-cc|=%.3Le\n", cs(t.j1).c_str(), cs(cnoj).c_str(), abs(cnoj - cc));
    });
    for (int th : {1, 0}) guard("NF2", [&] {
        // LRL = L^-1 o LLR o L: its fixed point L^-1(p) lies on 2z - xy = 0, where twist_at's (x, y) chart is singular
        V p = select_S("LLR", th).v, q{p[1] * p[0] - p[2], p[1], p[0]};
        printf("INFO NF2 LRL th=%d skipped: fixed_points finds %zu points; q = L^-1(p_LLR): |LRL(q)-q|=%.2Le |2z-xy|=%.2Le\n", th,
               fixed_points("LRL", th).size(), vdist(word_apply("LRL", q), q), abs((R)2 * q[2] - q[0] * q[1]));
    });
    for (auto& [ws, ths] : vector<pair<vector<string>, vector<int>>>{{{"LLR", "RLL"}, {1, 0}}, {{"LR", "RL"}, {0, 1, 2}}})
        for (int th : ths) guard("NF2", [&] {
            // cyclic rotations are conjugate maps: mu, b, c agree (rule S point of each word)
            Twist t0 = twist_at(ws[0], select_S(ws[0], th).v);
            for (size_t k = 1; k < ws.size(); k++) {
                Twist t = twist_at(ws[k], select_S(ws[k], th).v);
                R dmu = abs(t.mu - t0.mu), db = abs(t.b - t0.b), dc = abs(t.c - t0.c);
                check("NF2", dmu < 1e-9L && db < 1e-9L && dc < 1e-9L,
                      fmt("%s vs %s th=%d |dmu|=%.2Le |db|=%.2Le |dc|=%.2Le c=", ws[k].c_str(), ws[0].c_str(), th, dmu, db, dc) + cs(t.c));
            }
        });
    {
        // b from the degree-3 code (before the degree-5 normal form)
        vector<tuple<string, int, C>> oldb{
            {"LR", 0, C(-1.934462610796869539485e-17L, 4.398859193825720770624e-01L)},
            {"LR", 1, C(-5.743584406148464875397e-17L, 3.539855571509915051600e-01L)},
            {"LR", 2, C(-1.463754879340892400002e-17L, 1.975345037207917582936e-01L)},
            {"LR", 3, C(1.377409671725628954852e-17L, 8.359377166241821505860e-02L)},
            {"LLR", 0, C(3.097463284379172179458e-01L, 5.227919813585572863493e-01L)},
            {"LLR", 1, C(2.194844787814619090825e-01L, 5.160062356284023911629e-01L)},
            {"LLR", 2, C(2.353652741714133148734e-02L, 4.212373892439574900282e-01L)},
            {"LLR", 3, C(-8.907318466061656150710e-02L, 2.290632748285385580754e-01L)}};
        for (auto& [w, th, b0] : oldb) guard("NF3", [&] {
            Twist t = twist(w, th);
            check("NF3", abs(t.b - b0) < 1e-14L, fmt("%s th=%d |b-b_old|=%.2Le", w.c_str(), th, abs(t.b - b0)));
            check("NF5", t.consist < 1e-9L && t.jmax < 1e-9L, fmt("%s th=%d u^2 v^3 of N[1] rel=%.2Le density terms max=%.2Le", w.c_str(), th, t.consist, t.jmax));
            printf("INFO NF4 %s th=%d c=%s j1=%s\n", w.c_str(), th, cs(t.c).c_str(), cs(t.j1).c_str());
            if (w == "LR" && th == 0)
                printf("INFO NF4 LR th=0 mu-(5+sqrt21)/2=%s b-16i/(21 sqrt3)=%s\n", cs(t.mu - (5 + sqrtl(21)) / 2).c_str(),
                       cs(t.b - I1 * (16 / (21 * sqrtl(3)))).c_str());
        });
    }
    guard("O2S", [&] {
        // rows from L_n = r Vol/(4 pi) - m log|mu| + (alpha + gamma m^2)/r + D(m)/r^2 + e_n/r^3, D = d0 + d1 m + d2 m^2 + d3 m^3
        const R vol = 2.03L, lmu = 1.567L, A = 0.446L, G = 1.382L, d[4] = {0.13L, -0.21L, 0.37L, 0.52L}, d3p = 0.45L;
        vector<O2Row> rows;
        for (int r : {2001, 3001, 5001, 10001}) {
            O2Row o{"syn", 0, r, 0, G, d3p, {}, {}};
            for (int n = 0; n < 8; n++) {
                R m = n + 0.5L, D = d[0] + m * (d[1] + m * (d[2] + m * d[3]));
                o.L.push_back(r * vol / (4 * PI) - m * lmu + (A + G * m * m) / r + D / ((R)r * r) + (0.8L + 0.1L * n - 0.05L * n * n) / ((R)r * r * r));
            }
            o.res = resid_vec(r, o.L, vol, lmu);
            rows.push_back(o);
        }
        R e3 = 0, e2 = 0;
        for (int n = 0; n < 5; n++) {
            vector<R> z;
            for (auto& o : rows) z.push_back(o2_T(o)[n] - 6 * o.d3p);
            e3 = max(e3, fabsl((fit_r(rows, z, true).p[0] + 6 * d3p) / (6 * d[3]) - 1));
        }
        for (int n = 0; n < 6; n++) {
            vector<R> z;
            for (auto& o : rows) z.push_back(o2_W(o)[n] - 6 * o.d3p * (n + 1.5L));
            e2 = max(e2, fabsl(fit_r(rows, z, true).p[0] - 2 * d[2] - 6 * (d[3] - d3p) * (n + 1.5L)));
        }
        auto bf = b4_fits(rows);
        R ea = 0;
        vector<R> Dn;
        for (int n = 0; n < 8; n++) { ea = max(ea, fabsl(bf[1][n].p[0] - A)); Dn.push_back(bf[1][n].p[1]); }
        auto pf = poly_m(Dn, {0, 1, 2, 3});
        R ed = 0;
        for (int k = 0; k < 4; k++) ed = max(ed, fabsl(pf.p[k] - d[k]));
        check("O2S", e3 < 1e-6L && ea < 1e-8L, fmt("synthetic max rel err 6d3 (B1, 3-term)=%.2Le max|alpha err| (B4, 3-term)=%.2Le", e3, ea));
        check("O2S", e2 < 1e-6L && ed < 1e-6L, fmt("synthetic max|2d2 err| (B3, 3-term)=%.2Le max|d_k err| (B4 cubic in m)=%.2Le", e2, ed));
    });
    printf("failures: %d\n", fails);
    return fails;
}

int cmd_mu(const string& w, R th) {
    Geo g = geometric(w, th);
    printf("word=%s theta=%.12Lf\nfixed point: %s %s %s\nkappa: %s\nmu: %s\nlog|mu|=%.15Lf arg(mu)=%.15Lf\nVol=%.15Lf max-Vol is max-|mu|: %s\n",
           w.c_str(), th, cs(g.best.v[0]).c_str(), cs(g.best.v[1]).c_str(), cs(g.best.v[2]).c_str(), cs(kappa(g.best.v)).c_str(),
           cs(g.best.mu).c_str(), logl(abs(g.best.mu)), arg(g.best.mu), g.best.vol, g.maxmu ? "yes" : "no");
    printf("candidates (x y z | Vol |mu|):\n");
    for (auto& c : g.all)
        printf("  %s %s %s | %.12Lf %.12Lf\n", cs(c.v[0]).c_str(), cs(c.v[1]).c_str(), cs(c.v[2]).c_str(), c.vol, abs(c.mu));
    return 0;
}
int cmd_vol(const string& w, R th) {
    Geo g = geometric(w, th);
    vector<Shape> sh;
    R v = volume(w, g.best.v, th, 1, &sh);
    printf("word=%s theta=%.12Lf Vol=%.15Lf\n", w.c_str(), th, v);
    for (size_t i = 0; i < sh.size(); i++) printf("  w_%zu=%s D=%.15Lf\n", i, cs(sh[i].w).c_str(), sh[i].D);
    return 0;
}
int cmd_resid(const vector<string>& files) {
    printf("word r c theta Vol log|mu| resid_0..resid_7\n");
    for (auto& row : read_raw(files)) {
        Resid o = resid_row(row);
        printf("%s %d %d %.10Lf %.12Lf %.12Lf", row.word.c_str(), row.r, row.c, o.theta, o.vol, o.logmu);
        for (R x : o.res) printf(" %.6Lf", x);
        printf("\n");
    }
    return 0;
}
int cmd_fit(const vector<string>& files) {
    printf("word,theta_nom,alpha,gamma,b_0,b_1,b_2,b_3,b_4,b_5,b_6,b_7,rms,nrows\n");
    for (auto& f : fit(files)) {
        printf("%s,%d,%.6Lf,%.6Lf", f.word.c_str(), f.th, f.alpha, f.gamma);
        for (R x : f.b) printf(",%.6Lf", x);
        printf(",%.3Le,%d\n", f.rms, f.nrows);
    }
    return 0;
}
void fitC(const vector<C>& b, const vector<R>& g, const string& label) {
    vector<vector<R>> X; for (auto& z : b) X.push_back({z.real(), -z.imag()});
    auto p = lsq(X, g);
    C Cc(p[0], p[1]);
    R ss = 0, sg = 0;
    printf("%s: C=%s\n", label.c_str(), cs(Cc).c_str());
    for (size_t i = 0; i < b.size(); i++) {
        R f = (Cc * b[i]).real();
        printf("  gamma=%.5Lf fitted=%.5Lf resid=%+.5Lf\n", g[i], f, g[i] - f);
        ss += (g[i] - f) * (g[i] - f); sg += g[i] * g[i];
    }
    printf("  relative rms=%.4Le\n", sqrtl(ss / sg));
}
int cmd_twist() {
    auto refs = read_fits();
    vector<C> b[3]; vector<R> g[3];
    for (string w : {"LR", "LLR"}) for (int th = 0; th <= 3; th++) {
        R gam = NAN;
        for (auto& rf : refs) if (rf.word == w && rf.th == th) gam = rf.gamma;
        Twist t = twist(w, th);
        printf("%s th=%d mu=%s b=%s gamma=%.5Lf gamma/Re(b)=%.6Lf gamma/Im(b)=%.6Lf gamma/|b|=%.6Lf\n", w.c_str(), th,
               cs(t.mu).c_str(), cs(t.b).c_str(), gam, gam / t.b.real(), gam / t.b.imag(), gam / abs(t.b));
        for (int k : {0, w == "LR" ? 1 : 2}) { b[k].push_back(t.b); g[k].push_back(gam); }
    }
    fitC(b[0], g[0], "all 8");
    fitC(b[1], g[1], "LR only");
    fitC(b[2], g[2], "LLR only");
    return 0;
}
int cmd_robust(const vector<string>& files) {
    auto gs = make_groups(files);
    printf("word,theta_nom,variant,nrows,alpha,gamma,lin,Im_b,gamma/Im_b,gamma/(pi*Im_b)-1,note\n");
    for (auto& g : gs) {
        R imb = twist(g.word, g.th).b.imag();
        auto vs = robust_variants(g);
        for (auto& v : vs) printf("%s\n", robust_line(g, v, imb).c_str());
        for (int strict = 0; strict < 2; strict++) {
            R gmin = INFINITY, gmax = -INFINITY, qmin = INFINITY, qmax = -INFINITY, dmax = 0;
            int k = 0;
            for (auto& v : vs) if (v.ok && !v.singular && !(strict && v.exact)) {
                R q = v.s.gamma / imb;
                gmin = min(gmin, v.s.gamma); gmax = max(gmax, v.s.gamma); qmin = min(qmin, q); qmax = max(qmax, q);
                dmax = max(dmax, fabsl(q / PI - 1)); k++;
            }
            printf("%s,%d,summary(%s),%d variants", g.word.c_str(), g.th, strict ? "without *" : "all", k);
            if (k) printf(",gamma=[%.10Lf %.10Lf],gamma/Im_b=[%.8Lf %.8Lf],max|gamma/(pi*Im_b)-1|=%.3Le", gmin, gmax, qmin, qmax, dmax);
            printf("\n");
        }
    }
    printf("\ndirect: per raw row, second stage over n of res[n] alone; z_row = gamma_row - pi*Im_b(theta_row), b at the row's exact theta\n");
    printf("word,theta_nom,r,c,theta,stage2,gamma_row,pi*Im_b(theta_row),z_row,z_row/gamma_row\n");
    vector<vector<DirectFit>> fits;
    for (auto& g : gs) {
        auto dr = direct_rows(g);
        for (auto& d : dr) for (size_t s = 0; s < STAGES.size(); s++)
            printf("%s,%d,%d,%d,%.10Lf,%s,%.10Lf,%.10Lf,%+.4Le,%+.4Le\n", d.word.c_str(), d.th, d.r, d.c, d.theta, get<0>(STAGES[s]).c_str(),
                   d.gam[s], d.pimb, d.gam[s] - d.pimb, (d.gam[s] - d.pimb) / d.gam[s]);
        fits.push_back(direct_fits(g, dr));
    }
    printf("\ndirect: z_row fitted over the rows of each case; claim z0 = 0\n");
    printf("word,theta_nom,stage2,regressors,nrows,z0,z0/(pi*Im_b(theta_nom)),rms,note\n");
    for (auto& fs : fits) for (auto& f : fs) {
        printf("%s,%d,%s,%s,%d,", f.word.c_str(), f.th, f.stage.c_str(), f.rsq ? "{1,1/r,1/r^2}" : "{1,1/r}", f.nrows);
        if (f.singular) printf("n/a,n/a,n/a,singular\n");
        else printf("%+.4Le,%+.4Le,%.3Le,\n", f.z0, f.rel, f.rms);
    }
    return 0;
}
int cmd_gaps(const string& file) {
    auto rows = read_gaps(file);
    auto refs = read_fits();
    for (int rsq = 0; rsq < 2; rsq++) {
        printf("%sregressors {1,1/r%s} (+ theta-theta_nom for theta_nom != 0)\n", rsq ? "\n" : "", rsq ? ",1/r^2" : "");
        printf("word,theta_nom,nrows,gamma_0,gamma_1,gamma_2,gamma_3,gamma_4,mean,spread,gamma_stored,Im_b,mean/Im_b,rmax*5e-6/2,rounding_bound_mean,note\n");
        for (auto& f : gap_fits(rows, rsq, false)) {
            R st = NAN;
            for (auto& rf : refs) if (rf.word == f.word && rf.th == f.th) st = rf.gamma;
            R imb = twist(f.word, f.th).b.imag();
            printf("%s,%d,%d", f.word.c_str(), f.th, f.nrows);
            if (f.singular) { printf(",n/a,n/a,n/a,n/a,n/a,n/a,n/a,%.5Lf,%.10Lf,n/a,%.2Le,n/a,singular\n", st, imb, f.rmax * 5e-6L / 2); continue; }
            for (R x : f.g) printf(",%.6Lf", x);
            printf(",%.6Lf,%.2Le,%.5Lf,%.10Lf,%.6Lf,%.2Le,%.2Le,\n", f.mean, f.spread, st, imb, f.mean / imb, f.rmax * 5e-6L / 2, f.bmean);
        }
    }
    for (int rsq = 0; rsq < 2; rsq++) {
        printf("\ndirect: w_n = y_n - pi*Im_b(theta_row), regressors {1,1/r%s}, no theta regressor; claim intercept = 0\n", rsq ? ",1/r^2" : "");
        printf("word,theta_nom,nrows,w0_0,w0_1,w0_2,w0_3,w0_4,mean,spread,rounding_bound_mean,pi*Im_b(theta_nom),mean/(pi*Im_b(theta_nom)),note\n");
        for (auto& f : gap_fits(rows, rsq, true)) {
            R pnom = PI * twist(f.word, f.th).b.imag();
            printf("%s,%d,%d", f.word.c_str(), f.th, f.nrows);
            if (f.singular) { printf(",n/a,n/a,n/a,n/a,n/a,n/a,n/a,n/a,%.10Lf,n/a,singular\n", pnom); continue; }
            for (R x : f.g) printf(",%+.3Le", x);
            printf(",%+.3Le,%.2Le,%.2Le,%.10Lf,%+.3Le,\n", f.mean, f.spread, f.bmean, pnom, f.mean / pnom);
        }
    }
    return 0;
}
int cmd_stress(const string& file) {
    auto ss = read_stress(file);
    printf("word,theta,Vol_file,Vol_calc,Vol_calc-file,logmu_file,logmu_calc,logmu_calc-file,maxmu,"
           "S-Vol_calc,S_err,|S-Vol|/S_err,d0-logmu_calc,d0_err,|d0-logmu|/d0_err,d1-logmu_calc,d1_err,|d1-logmu|/d1_err\n");
    for (auto& s : ss) {
        R vol = s.g.best.vol, lm = logl(abs(s.g.best.mu));
        printf("%s,%.0Lf,%.13Lf,%.13Lf,%+.2Le,%.8Lf,%.10Lf,%+.2Le,%s,%+.3Le,%.2Le,%.3Lf,%+.3Le,%.2Le,%.3Lf,%+.3Le,%.2Le,%.3Lf\n",
               s.word.c_str(), s.theta, s.vol_f, vol, vol - s.vol_f, s.logmu_f, lm, lm - s.logmu_f, s.g.maxmu ? "yes" : "no",
               s.S - vol, s.S_err, fabsl(s.S - vol) / s.S_err, s.d0 - lm, s.d0_err, fabsl(s.d0 - lm) / s.d0_err,
               s.d1 - lm, s.d1_err, fabsl(s.d1 - lm) / s.d1_err);
    }
    auto refs = read_fits();
    printf("\nprediction gamma_pred = pi*Im(b); rule S: point of max Vol whose signed sum S has the sign of LR theta 0 (S > 0 there)\n");
    printf("Im x>0: geometric()'s default tie-break, for comparison; same = both rules pick the same point\n");
    printf("word,theta,Vol,S(rule S),b(rule S),gamma_pred(rule S),S(Im x>0),b(Im x>0),gamma_pred(Im x>0),same,gamma_stored,gamma_stored/gamma_pred(rule S)\n");
    vector<string> words{"LR", "LLR"};
    for (auto& s : ss) if (find(words.begin(), words.end(), s.word) == words.end()) words.push_back(s.word);
    for (auto& w : words) for (int th = 0; th <= 3; th++) {
        R st = NAN;
        for (auto& rf : refs) if (rf.word == w && rf.th == th) st = rf.gamma;
        printf("%s,%d,", w.c_str(), th);
        try {
            Geo g = geometric(w, th);
            Cand c = select_S(w, th);
            Twist tS = twist_at(w, c.v), tI = twist_at(w, g.best.v);
            R gS = PI * tS.b.imag(), gI = PI * tI.b.imag();
            printf("%.13Lf,%+.13Lf,%s,%.8Lf,%+.13Lf,%s,%.8Lf,%s,", c.vol, volume_signed(w, c.v, th, 1), cs(tS.b).c_str(), gS,
                   volume_signed(w, g.best.v, th, 1), cs(tI.b).c_str(), gI, vdist(c.v, g.best.v) == 0 ? "yes" : "no");
            if (isnan(st)) printf(",\n"); else printf("%.5Lf,%.6Lf\n", st, st / gS);
        } catch (exception& e) { printf("error: %s\n", e.what()); }
    }
    return 0;
}

// cells "regressors,nrows," then either the values or n/a with the note; * = rows equal regressors (exact)
string rfit_head(const RFit& f, bool rsq, size_t nrows) { return fmt("%s%s,%zu,", rsq ? "{1,1/r,1/r^2}" : "{1,1/r}", f.exact ? "*" : "", nrows); }
string rfit_na(const RFit& f, int ncells) {
    string s;
    for (int i = 0; i < ncells; i++) s += "n/a,";
    return s + (f.ok ? "singular" : "rows<regressors");
}
int cmd_order2(const vector<string>& files) {
    auto gs = make_groups(files);
    vector<vector<O2Row>> all;
    for (auto& g : gs) all.push_back(o2_rows(g));
    printf("B1: T_n = r^2 Delta^3 L_n, claim T_n -> 6 d3_pred, d3_pred = (4 pi^2/3) Re c, c at the point of the direct test\n");
    printf("word,theta_nom,r,theta_row,6d3_pred(theta_row),T_0,T_1,T_2,T_3,T_4\n");
    for (auto& rows : all) for (auto& o : rows) {
        printf("%s,%d,%d,%.10Lf,%.8Lf", o.word.c_str(), o.th, o.r, o.theta, 6 * o.d3p);
        for (R x : o2_T(o)) printf(",%.6Lf", x);
        printf("\n");
    }
    printf("\nB1 fits: z = T_n - 6 d3_pred(theta_row) over the rows; claim z0 = 0\n");
    printf("word,theta_nom,n,regressors,nrows,z0,z0/(6d3_pred(theta_nom)),rms,note\n");
    for (size_t i = 0; i < gs.size(); i++) {
        R p6 = 6 * d3_pred(twist(gs[i].word, gs[i].th));
        for (bool rsq : {false, true}) for (int n = 0; n < 5; n++) {
            vector<R> z;
            for (auto& o : all[i]) z.push_back(o2_T(o)[n] - 6 * o.d3p);
            RFit f = fit_r(all[i], z, rsq);
            printf("%s,%d,%d,%s", gs[i].word.c_str(), gs[i].th, n, rfit_head(f, rsq, all[i].size()).c_str());
            if (!f.ok || f.singular) printf("%s\n", rfit_na(f, 3).c_str());
            else printf("%+.4Le,%+.4Le,%.3Le,\n", f.p[0], f.p[0] / p6, f.rms);
        }
    }
    printf("\nB2: fourth differences, r^2 Delta^4 L_n (claim -> 0) and r^3 Delta^4 L_n\n");
    printf("word,theta_nom,r,r^2D4_0,r^2D4_1,r^2D4_2,r^2D4_3,r^3D4_0,r^3D4_1,r^3D4_2,r^3D4_3\n");
    for (auto& rows : all) for (auto& o : rows) {
        auto d = fdiff(o.L, 4);
        printf("%s,%d,%d", o.word.c_str(), o.th, o.r);
        for (R x : d) printf(",%+.6Le", (R)o.r * o.r * x);
        for (R x : d) printf(",%+.6Le", (R)o.r * o.r * o.r * x);
        printf("\n");
    }
    printf("\nB3: W_n = r (r Delta^2 L_n - 2 pi Im b(theta_row)), claim W_n -> 2 d2 + 6 d3 (n + 3/2)\n");
    printf("word,theta_nom,r,theta_row,W_0,W_1,W_2,W_3,W_4,W_5\n");
    for (auto& rows : all) for (auto& o : rows) {
        printf("%s,%d,%d,%.10Lf", o.word.c_str(), o.th, o.r, o.theta);
        for (R x : o2_W(o)) printf(",%.6Lf", x);
        printf("\n");
    }
    printf("\nB3 fits: V_n = W_n - 6 d3_pred(theta_row)(n + 3/2) over the rows; intercept = estimate of 2 d2\n");
    printf("word,theta_nom,n,regressors,nrows,2d2_est,6d3_pred(theta_nom)(n+3/2),rms,note\n");
    for (size_t i = 0; i < gs.size(); i++) {
        R p6 = 6 * d3_pred(twist(gs[i].word, gs[i].th));
        for (bool rsq : {false, true}) for (int n = 0; n < 6; n++) {
            vector<R> z;
            for (auto& o : all[i]) z.push_back(o2_W(o)[n] - 6 * o.d3p * (n + 1.5L));
            RFit f = fit_r(all[i], z, rsq);
            printf("%s,%d,%d,%s", gs[i].word.c_str(), gs[i].th, n, rfit_head(f, rsq, all[i].size()).c_str());
            if (!f.ok || f.singular) printf("n/a,%.6Lf,n/a,%s\n", p6 * (n + 1.5L), f.ok ? "singular" : "rows<regressors");
            else printf("%+.4Le,%.6Lf,%.3Le,\n", f.p[0], p6 * (n + 1.5L), f.rms);
        }
    }
    printf("\nB4 (theta_nom = 0): y_n = res[n] - pi Im b(theta_row) m^2 over the rows; intercept = alpha, 1/r slope = D(n)\n");
    for (size_t i = 0; i < gs.size(); i++) {
        if (gs[i].th != 0) continue;
        auto bf = b4_fits(all[i]);
        R d3p0 = d3_pred(twist(gs[i].word, 0));
        printf("%s theta_nom=0 nrows=%zu\n", gs[i].word.c_str(), all[i].size());
        printf("n,m,alpha{1,1/r},D{1,1/r},alpha{1,1/r,1/r^2},D{1,1/r,1/r^2}\n");
        for (int n = 0; n < 8; n++) {
            printf("%d,%.1Lf", n, n + 0.5L);
            for (int rsq = 0; rsq < 2; rsq++) {
                if (!bf[rsq][n].ok || bf[rsq][n].singular) printf(",n/a,n/a");
                else printf(",%.8Lf,%.6Lf", bf[rsq][n].p[0], bf[rsq][n].p[1]);
            }
            printf("\n");
        }
        for (int rsq = 0; rsq < 2; rsq++) {
            vector<R> a;
            for (auto& f : bf[rsq]) a.push_back(f.p[0]);
            Spread sp = spread_of(a);
            printf("alpha %s: mean=%.8Lf spread(max-min)=%.3Le\n", rsq ? "{1,1/r,1/r^2}" : "{1,1/r}", sp.mean, sp.spread);
        }
        vector<R> D;
        for (auto& f : bf[1]) D.push_back(f.p[1]);
        auto c3 = poly_m(D, {0, 1, 2, 3}), od = poly_m(D, {1, 3});
        printf("D(n) {1,1/r,1/r^2} cubic in m: d0=%+.6Lf d1=%+.6Lf d2=%+.6Lf d3=%+.6Lf rms=%.3Le d3_pred(0)=%+.6Lf\n",
               c3.p[0], c3.p[1], c3.p[2], c3.p[3], c3.rms, d3p0);
        printf("D(n) {1,1/r,1/r^2} odd powers: d1=%+.6Lf d3=%+.6Lf rms=%.3Le\n", od.p[0], od.p[1], od.rms);
    }
    return 0;
}
// continued-fraction convergents p/q of x with q <= qmax
vector<pair<long long, long long>> convergents(R x, long long qmax) {
    vector<pair<long long, long long>> out;
    long long p1 = 1, p2 = 0, q1 = 0, q2 = 1;
    R y = x;
    for (int it = 0; it < 40; it++) {
        R a = floorl(y);
        long long ai = (long long)a, p = ai * p1 + p2, q = ai * q1 + q2;
        if (q > qmax) break;
        out.push_back({p, q});
        p2 = p1; p1 = p; q2 = q1; q1 = q;
        if (y - a < 1e-15L) break;
        y = 1 / (y - a);
    }
    return out;
}
int cmd_alpha(const vector<string>& files) {
    auto gs = make_groups(files);
    printf("theta_nom = 0, alpha from B4 {1,1/r,1/r^2} (mean over n = 0..7 with spread max-min, and n = 0 alone)\n");
    printf("gamma = pi Im b(0), alpha0 = alpha + gamma/4, S2 = alpha0 + 2 pi Im(b/(mu + 1/mu - 2)), mu and b at theta = 0\n");
    for (auto& g : gs) {
        if (g.th != 0) continue;
        auto rows = o2_rows(g);
        auto bf = b4_fits(rows);
        vector<R> a;
        for (auto& f : bf[1]) a.push_back(f.p[0]);
        Spread sp = spread_of(a);
        Twist t = twist(g.word, 0);
        R gam = PI * t.b.imag(), sh = 2 * PI * (t.b / (t.mu + (R)1 / t.mu - (R)2)).imag();
        printf("\n%s: gamma=%.10Lf 2 pi Im(b/(mu+1/mu-2))=%.10Lf spread=%.3Le\n", g.word.c_str(), gam, sh, sp.spread);
        vector<R> ds = g.word == "LR" ? vector<R>{3, 21} : vector<R>{7, 3, 21};
        for (auto [ename, al] : vector<pair<string, R>>{{"mean", sp.mean}, {"n=0", a[0]}}) {
            R q[3] = {al, al + gam / 4, al + gam / 4 + sh};
            const char* qn[3] = {"alpha", "alpha0", "S2"};
            printf("%s estimate=%s: alpha=%.10Lf alpha0=%.10Lf S2=%.10Lf (spread %.3Le)\n", g.word.c_str(), ename.c_str(), q[0], q[1], q[2], sp.spread);
            if (g.word == "LR") {
                R e1 = q[0] - 31 * PI / (126 * sqrtl(3)), e2 = q[2] - 17 * PI / (18 * sqrtl(3));
                printf("  alpha - 31 pi/(126 sqrt3)=%+.4Le (%+.3Lf spreads)  S2 - 17 pi/(18 sqrt3)=%+.4Le (%+.3Lf spreads)\n",
                       e1, e1 / sp.spread, e2, e2 / sp.spread);
            }
            for (int k = 0; k < 3; k++) for (R d : ds) {
                R x = q[k] * sqrtl(d) / PI, err = sp.spread * sqrtl(d) / PI;
                printf("  %s %s: x=%s*sqrt(%.0Lf)/pi=%.10Lf err=%.3Le convergents (p/q, x-p/q, (x-p/q)/err):", g.word.c_str(), ename.c_str(), qn[k], d, x, err);
                for (auto [p, qq] : convergents(x, 2000)) {
                    R e = x - (R)p / qq;
                    printf(" %lld/%lld %+.3Le %+.3Lf;", p, qq, e, e / err);
                }
                printf("\n");
            }
        }
    }
    printf("\nall cases (direct-test point): word,theta_nom,mu,b,c,6d3_pred,2pi*Im(b/(mu+1/mu-2))\n");
    for (string w : {"LR", "LLR"}) for (int th = 0; th <= 3; th++) {
        Twist t = twist(w, th);
        printf("%s,%d,%s,%s,%s,%+.10Lf,%+.10Lf\n", w.c_str(), th, cs(t.mu).c_str(), cs(t.b).c_str(), cs(t.c).c_str(), 6 * d3_pred(t),
               2 * PI * (t.b / (t.mu + (R)1 / t.mu - (R)2)).imag());
    }
    return 0;
}

int main(int argc, char** argv) {
    vector<string> a(argv + 1, argv + argc);
    try {
        if (a.size() == 1 && a[0] == "selftest") return selftest();
        if (a.size() == 3 && a[0] == "mu") return cmd_mu(a[1], strtold(a[2].c_str(), nullptr));
        if (a.size() == 3 && a[0] == "vol") return cmd_vol(a[1], strtold(a[2].c_str(), nullptr));
        if (a.size() >= 2 && a[0] == "resid") return cmd_resid(vector<string>(a.begin() + 1, a.end()));
        if (a.size() >= 2 && a[0] == "fit") return cmd_fit(vector<string>(a.begin() + 1, a.end()));
        if (a.size() == 1 && a[0] == "twist") return cmd_twist();
        if (a.size() >= 2 && a[0] == "robust") return cmd_robust(vector<string>(a.begin() + 1, a.end()));
        if (a.size() == 2 && a[0] == "gaps") return cmd_gaps(a[1]);
        if (a.size() == 2 && a[0] == "stress") return cmd_stress(a[1]);
        if (a.size() >= 2 && a[0] == "order2") return cmd_order2(vector<string>(a.begin() + 1, a.end()));
        if (a.size() >= 2 && a[0] == "alpha") return cmd_alpha(vector<string>(a.begin() + 1, a.end()));
    } catch (exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    fprintf(stderr, "usage: ladder_tools selftest | mu WORD THETA | vol WORD THETA | resid FILE... | fit FILE... | twist"
                    " | robust FILE... | gaps SUMMARY | stress CSV | order2 FILE... | alpha FILE...\n");
    return 2;
}
