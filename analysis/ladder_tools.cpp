// Trace-map fixed points, multipliers, volumes, twist coefficients and ladder fits.
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
    vector<V> sols;
    for (int start = 0; start < 400; start++) {
        V v;
        for (auto& c : v) { R re = rnd(); c = C(re, rnd()); }
        for (int it = 0; it < 200; it++) {
            V fv = word_apply(w, v);
            M3 J = word_jac(w, v), A;
            A[0] = J[0]; A[0][0] -= (R)1;
            A[1] = J[1]; A[1][1] -= (R)1;
            A[2] = {(R)2 * v[0] - v[1] * v[2], (R)2 * v[1] - v[0] * v[2], (R)2 * v[2] - v[0] * v[1]};
            V F{fv[0] - v[0], fv[1] - v[1], kappa(v) - k0}, dv;
            if (!solve3(A, F, dv)) break;
            for (int i = 0; i < 3; i++) v[i] -= dv[i];
            R n = max({abs(v[0]), abs(v[1]), abs(v[2])});
            if (!(n < 1e6L)) break;
            if (max({abs(dv[0]), abs(dv[1]), abs(dv[2])}) < 1e-18L * (1 + n)) break;
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
R volume(const string& w, V v, R theta, int sgn, vector<Shape>* shapes = nullptr) {
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
        C t = (z + sqrt(z * z - (R)4)) / (R)2;
        M2 A{{{x, -1}, {1, 0}}}, B{{{0, t}, {-(R)1 / t, y}}};
        M2 Ai{{{0, 1}, {-1, x}}}, Bi{{{y, -t}, {(R)1 / t, 0}}};
        M2 K = mm(mm(Ai, Bi), mm(A, B));
        array<C, 2> p1{K[0][1], lam0 - K[0][0]}, p2{lam0 - K[1][1], K[1][0]};
        array<C, 2> p = norm(p1[0]) + norm(p1[1]) >= norm(p2[0]) + norm(p2[1]) ? p1 : p2;
        array<C, 2> q0 = p, q1 = mv(A, p), q2 = mv(mm(A, B), p), q3 = mv(B, p);
        C sh = det(q0, q2) * det(q1, q3) / (det(q0, q3) * det(q1, q2));
        R D = bloch_wigner(sh);
        if (shapes) shapes->push_back({sh, D});
        sum += D;
        v = letter(w[i], v);
    }
    return fabsl(sum);
}

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
struct Twist { C mu, c1, c2, b; R detlin, quadmax; };
// bivariate series truncated at total degree 3: c[i][j] multiplies s^i t^j
struct S {
    C c[4][4] = {};
    S() {}
    S(C a) { c[0][0] = a; }
};
const int DEG = 3;
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

Twist twist(const string& w, R theta) {
    Geo g = geometric(w, theta);
    V p0 = g.best.v;
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
    return out;
}

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
        if (line.empty()) continue;
        vector<string> cells; string cell; istringstream ss(line);
        while (getline(ss, cell, ',')) cells.push_back(cell);
        out.push_back(cells);
    }
    return out;
}
struct Ref { string word; int th; R alpha, gamma, logmu, vol; };
vector<Ref> read_fits() {
    vector<Ref> out;
    for (auto& c : read_csv("data/ladder_subleading_fits.csv"))
        out.push_back({c[0], atoi(c[1].c_str()), strtold(c[2].c_str(), nullptr), strtold(c[3].c_str(), nullptr),
                       strtold(c[12].c_str(), nullptr), strtold(c[13].c_str(), nullptr)});
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

struct Resid { Row row; R theta, vol, logmu; vector<R> res; bool maxmu; };
Resid resid_row(const Row& row) {
    Resid o; o.row = row; o.theta = row_theta(row.r, row.c);
    Geo g = geometric(row.word, o.theta);
    o.vol = g.best.vol; o.logmu = logl(abs(g.best.mu)); o.maxmu = g.maxmu;
    for (int n = 0; n < 8; n++)
        o.res.push_back(row.r * (row.L[n] - row.r * o.vol / (4 * PI) + (n + 0.5L) * o.logmu));
    return o;
}

struct FitOut { string word; int th; R alpha, gamma, rms; vector<R> b; int nrows; };
vector<FitOut> fit(const vector<string>& files) {
    map<pair<string, int>, vector<Resid>> groups;
    for (auto& row : read_raw(files)) {
        Resid o = resid_row(row);
        groups[{row.word, (int)lroundl(o.theta)}].push_back(o);
    }
    vector<FitOut> out;
    for (auto& [key, rs] : groups) {
        FitOut f; f.word = key.first; f.th = key.second; f.nrows = rs.size();
        for (int n = 0; n < 8; n++) {
            vector<vector<R>> X; vector<R> y;
            for (auto& o : rs) {
                vector<R> x{1, 1 / (R)o.row.r};
                if (f.th != 0) x.push_back(o.theta - f.th);
                X.push_back(x); y.push_back(o.res[n]);
            }
            f.b.push_back(lsq(X, y)[0]);
        }
        vector<vector<R>> X; vector<R> y;
        for (int n = 0; n < 8; n++) { X.push_back({1, (n + 0.5L) * (n + 0.5L)}); y.push_back(f.b[n]); }
        auto ag = lsq(X, y);
        f.alpha = ag[0]; f.gamma = ag[1];
        R ss = 0;
        for (int n = 0; n < 8; n++) { R e = f.b[n] - f.alpha - f.gamma * X[n][1]; ss += e * e; }
        f.rms = sqrtl(ss / 8);
        out.push_back(f);
    }
    return out;
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

int main(int argc, char** argv) {
    vector<string> a(argv + 1, argv + argc);
    try {
        if (a.size() == 1 && a[0] == "selftest") return selftest();
        if (a.size() == 3 && a[0] == "mu") return cmd_mu(a[1], strtold(a[2].c_str(), nullptr));
        if (a.size() == 3 && a[0] == "vol") return cmd_vol(a[1], strtold(a[2].c_str(), nullptr));
        if (a.size() >= 2 && a[0] == "resid") return cmd_resid(vector<string>(a.begin() + 1, a.end()));
        if (a.size() >= 2 && a[0] == "fit") return cmd_fit(vector<string>(a.begin() + 1, a.end()));
        if (a.size() == 1 && a[0] == "twist") return cmd_twist();
    } catch (exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    fprintf(stderr, "usage: ladder_tools selftest | mu WORD THETA | vol WORD THETA | resid FILE... | fit FILE... | twist\n");
    return 2;
}
