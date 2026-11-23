# Plan: from data repo to paper

Status: draft, not accepted. Nothing below is implemented. Untracked file, not for commit until decided.

## 1. Goal

A paper that states the spectral ladder as a precise conjecture, explains where it comes from, and backs it
with reproducible numerics. This repo becomes the companion (code + data) the paper cites.

## 2. Facts the plan rests on

Read from the repo:

| # | Fact | Source |
|---|---|---|
| F1 | Claim: abs(lambda_j) = exp(r Vol/4pi) abs(mu)^-(j+1/2) (1 + O(1/r)), j from 1 | README.md:5-7 |
| F2 | Same claim reindexed n = j-1 with exponent (n+1/2) and a 1/r term alpha + gamma (n+1/2)^2 | README.md:14-16 |
| F3 | F1 and F2 disagree on the exponent: j+1/2 = n+3/2, not n+1/2. The data follow F2 (check C1 below). | README.md:7 vs :16 |
| F4 | theta = abs(2pi - 4pi(c+1)/r), printed by the solver with c+1 | README.md:10, solver/ladder3.c:1170 |
| F5 | Cone axis is described as "meridian dSigma" | README.md:10 |
| F6 | alpha, gamma "closed forms not yet identified" | README.md:18 |
| F7 | n=0 statement is tied to the generalized BWY conjecture in one sentence, no bridge given | README.md:20 |
| F8 | High-r data exist for LR and LLR only, theta near 0, 1, 2, 3 | README.md:24, data/r10001_raw.txt:1-8 |
| F9 | Other words stop at r = 101 and are extrapolated | README.md:25 |
| F10 | That table has unusable entries (d4 = -54.8 +- 64.6; d4 = 210 +- 209) | data/stress_7words_r101.csv:5, :18 |
| F11 | Fit model for alpha, gamma: E_j r = b + c/r + d/r^2 + s (theta - theta_nom), r >= 500. No script in the repo produces this file. | data/ladder_subleading_fits.csv:1 |
| F12 | r = 201 to 1001 is present only as a 5 to 6 digit summary, no raw rows | data/r201-1001_summary.txt:1-2 |
| F13 | Solver accepts only the words LR and LLR | solver/ladder3.c:937 |
| F14 | Precision laws are per word (constants for pw = 1 or 2) | solver/ladder3.c:1037, :1040, :1085 |
| F15 | Eigenvalues come from Rayleigh-Ritz on the block; the complex Ritz values exist | solver/ladder3.c:1110-1120 |
| F16 | The sort swaps only the logs, so phase is lost; only logs are printed | solver/ladder3.c:1121, :1172 |
| F17 | Stop rule is a step-size test at 1e-22, not a residual | solver/ladder3.c:1123-1126 |
| F18 | Non-convergence warns on stderr only; stdout line is unchanged | solver/ladder3.c:1168, :1171-1173 |
| F19 | Return value of the eigen-solve is ignored | solver/ladder3.c:1119 |
| F20 | ladder3 uses pread, pwrite, posix_fadvise, fdatasync, getrusage (not available on native Windows) | solver/ladder3.c:49, :365, :510, :514, :529 |
| F21 | Multi-thread runs of ladder3 were not tested | solver/README.md:187 |
| F22 | solver/README.md cites files absent from the repo: ladder.c, r2001.txt, b1-b3.txt, rbig.txt, dev/, big_all.txt, dev/ref/r1001.txt; the sources cite README_ladder2.md | solver/README.md:3, :87, :135, :136, :185, :188; solver/ladder3.c:1 |
| F23 | solver/README.md reports r = 7001 runs; no r = 7001 data file is in data/ | solver/README.md:88, :100 |
| F24 | The Python reference builds S, T_a, T_b from 6j symbols at q = e^(2 pi i/r) (s = 2), even colours | analysis/reference_mprep.py:4-12, :41-53 |
| F25 | The Python reference stores full complex spectra for LR and LLR, so phases are available at small r | analysis/reference_run.py:11 |
| F26 | Trace map, kappa, and fixed-point search; the Jacobian is a finite difference with h = 1e-15; line 28 is dead code | analysis/trace_map.py:3-4, :10, :28-33 |
| F27 | analyze_rows.py execs a text slice of trace_map.py | analysis/analyze_rows.py:3-4 |
| F28 | Vol(theta) is SnapPy b++word with dehn_fill((2pi/theta, 0)) | analysis/analyze_rows.py:5-8 |
| F29 | plot4.py reads analysis_out.txt from the working directory and hardcodes ticks and labels up to r = 2001 | analysis/plot4.py:3, :25-26, :31 |
| F30 | No script produces figures/one_over_r_2001.png, data/extrap_r121.csv or data/stress_7words_r101.csv | analysis/ (five files, none writes them) |
| F31 | No tests, no CI | repo listing |

Hand checks done on the data (arithmetic only, nothing executed):

| # | Check | Result |
|---|---|---|
| C1 | r = 10001 LR theta~0: 4pi L1/r - Vol = -9.7e-4; the (n+1/2) form with n = 0 predicts -2pi log(mu)/r = -9.84e-4 | consistent with F2, not F1 (data/r10001_raw.txt:1) |
| C2 | Same row, L1 - L2 = 1.5665228; prediction log(mu) - 2 gamma/r = 1.566523 | agrees to 1e-6 (data/ladder_subleading_fits.csv:6) |
| C3 | LR theta = 0: log(mu) = 1.5668 gives mu + 1/mu = 5.000, so (1-mu)(1-1/mu) = -3 | data/ladder_subleading_fits.csv:6 |
| C4 | LR theta = 0: d2Vol = -0.14435 = -(2 sqrt 3)/24, the Schlafli value when the cone axis is the fiber boundary (knot longitude) | data/ladder_subleading_fits.csv:6 |

Literature facts. These come from memory or search snippets, so each must be read at the source before any
sentence in the paper depends on it:

| # | Statement to verify | Where |
|---|---|---|
| L1 | TV_r of a mapping torus is a sum over colours of abs(trace of rho_{r,c}(phi))^2 | Detcherry-Kalfagianni, arXiv 1711.03251 |
| L2 | Cone-angle convention theta = abs(2pi - 4pi a/r) for relative RT invariants | Wong-Yang, arXiv 2008.05045 |
| L3 | Asymptotic expansion conjecture with adjoint twisted torsion as 1-loop term | Wong-Yang, arXiv 2103.15056 |
| L4 | Generalized BWY conjecture with cone structures; proved for LR | arXiv 2402.04483 |
| L5 | BWY invariant equals the 1-loop invariant for once-punctured torus bundles | Garoufalidis-Yu, arXiv 2501.00250 |
| L6 | Adjoint torsion of the figure-eight for the longitude has modulus 3 at the complete structure | Porti / Dubois |
| L7 | Stretch factor as a limit of eigenvalues of TQFT matrices for M(0,4) | Andersen-Masbaum-Ueno, arXiv math/0503414 |
| L8 | 1-loop equals torsion for fibered 3-manifolds | arXiv 2304.00469 |

## 3. Work packages

Order: WP0, WP1, then WP2 and WP3 in parallel, then WP4. Each WP is one PR.

### WP0. Make the existing claims reproducible (no new results)

Changes:
1. README.md: fix F3 (one indexing, exponent n+1/2 with n from 0), fix F5 wording (cone axis = fiber boundary,
   the longitude of the fibered knot), soften F7 until the bridge in WP3-T0 exists.
2. solver/README.md: remove or satisfy every reference in F22 and F23. Preferred: add the missing reference
   outputs under solver/dev/ref/ if they still exist (decision D3), else delete the references.
3. No new Python. One dependency-free C++17 file, analysis/ladder_tools.cpp (std::complex<long double>, builds
   with plain g++), with subcommands. It reads the existing raw rows; it needs no solver run.
   - mu: fixed point of the trace map by Newton with the exact chain-rule Jacobian, and the multiplier mu
     (replaces analysis/trace_map.py, F26).
   - vol: Vol(theta) as a dilogarithm sum over the monodromy triangulation, shapes taken from the fixed
     point. The formula is read from the source (Gueritaud, canonical triangulations of once-punctured torus
     bundles) before it is coded. This removes the SnapPy dependency (F28).
   - resid: r (L_n - r Vol/4pi + (n+1/2) log mu) per row and n; should be flat at alpha + gamma (n+1/2)^2.
   - fit: regenerates data/ladder_subleading_fits.csv with the model in F11.
   Double-extended precision is enough: the residuals are about 1e-4 on values near 1600 (data/r10001_raw.txt:1).
4. The existing analysis/*.py stay untouched until the C++ tool reproduces their numbers; they are then
   removed in the same PR, with README.md:36-39 and requirements.txt updated. The extrapolation tables
   (F30) are not regenerated: WP2 replaces them with direct runs, and the junk entries (F10) go with them.
5. Figures are drawn by pgfplots in the paper from the CSVs the tool writes. figures/*.png and
   analysis/plot4.py are removed in the same PR (F29, F30).
6. Add raw rows for r = 201 to 1001 and r = 7001 if they exist (F12, F23; decision D3).

Tests, written first (fixtures are numbers already committed in data/):
- mu, LR, theta = 0: mu + 1/mu = 5 to 15 digits (C3); kappa at the fixed point = -2.
- mu reproduces all 8 logmu values in data/ladder_subleading_fits.csv:2-9 to the printed digits.
- vol reproduces the Vol column of data/extrap_r121.csv:2-17 to 1e-10 (16 values, two words, theta 0 to 3).
- fit reproduces every alpha, gamma in the committed CSV to the printed digits. If it does not, the CSV is
  wrong or the model in F11 is incomplete: stop and report.
- The theta column recomputed from (r, c) matches data/r2001_raw.txt:1-8.

### WP1. Solver changes

Changes in solver/ladder3.c (and ladder2.c where the code is shared):
1. Print arg(lambda_k) next to log abs(lambda_k): sort the Ritz values together with the logs (F15, F16). New
   output columns go at the end of the line so old parsers keep working.
2. Make failure visible on stdout: a status field and a nonzero exit code on non-convergence or a failed
   stage check; check the eigen-solve return (F18, F19).
3. Report an error estimate per eigenvalue from the Ritz residual, and print only the digits it supports (F17).
4. Words L^p R^q: keep K = T_b^-p, apply Theta^q, extend the word parser and the checkpoint signatures (F13).
   The precision laws (F14) are empirical per word, so each new word is validated by a rerun with all
   precisions raised, as solver/README.md:88 did for LR.
5. Reject even r and odd c with a message.

Not in this WP: multi-syllable words such as LRLRR (needs one K store per syllable; a separate plan change if
wanted), Windows or macOS ports (F20).

Tests, written first:
- New solver/ladder_ref.c: a dense small-r reference in C with FLINT (6j symbols, S, T_a, T_b, full
  acb_mat eigen-solve), a direct translation of analysis/reference_mprep.py:41-53 (F24). ladder3 against
  ladder_ref at r = 21, 31, 51, all admissible c, words LR, LLR, LLRR, LLLR: top 8 moduli and phases agree to
  25 digits. This is the test the repo lacks today.
- A run with -maxit 1 exits nonzero and marks the line as not converged.
- ladder3 with -t 4 equals -t 1 to all printed digits at r = 1001 (F21).

### WP2. New data

Runs (who executes is decision D2):
1. Phases for the existing 8 cases at r = 201 to 10001.
2. LLRR, LLLR, LLLRR, LLLLR at theta near 0, 1, 2, 3 and r = 501, 1001, 2001. This replaces the r <= 101
   extrapolations (F9, F10) for four of the five extra words.
3. Raw rows committed under data/ in the existing format plus the new columns; fits regenerated by WP0 scripts.

Acceptance: for every new (word, theta), the n = 0..7 residual from WP0 item 3 (resid) is flat in 1/r. A case
that is not flat is reported as a counterexample candidate, not dropped.

### WP3. Theory checks (cheap, C++ subcommands of ladder_tools, seconds each)

T0. Bridge. Write down exactly which known conjecture the n = 0 statement is. Working hypothesis: trace of
    rho_{r,c}(phi^k) is the relative RT invariant of the k-fold cyclic cover with the puncture orbit coloured
    c (L1, L2), so the top eigenvalue is controlled by Chen-Yang / Wong-Yang, and BWY (L4, L5) is a parallel
    statement for a different representation. Outcome decides the wording of F7.
T1. Trace formula. sum over n of mu^(-k(n+1/2)) = 1/(mu^(k/2) - mu^(-k/2)), whose modulus is
    abs(det(1 - dphi^k))^(-1/2). So the ladder is equivalent to: leading term plus torsion 1-loop term for
    every cyclic cover M_(phi^k), uniformly in k (L3, L8). This is an identity plus a citation; no run.
T2. Twist coefficient. Hypothesis: gamma is the first Birkhoff normal-form (twist) coefficient b of the trace
    map on the level surface kappa = const at the geometric fixed point, times a constant that depends only
    on the normalisation of Planck's constant. Test: compute b for the 8 (word, theta) cases and check that
    gamma / Re(b) (or the analogous real combination if mu is complex) is one constant across all 8
    (data/ladder_subleading_fits.csv:2-9). Gate: if yes, gamma gets a closed form in the conjecture and F6 is
    half resolved; if no, gamma stays a fitted constant and the paper says so.
    Result (ladder_tools twist, analysis/ladder_tools.cpp:225-293, 524-539): passed. With area form
    dx ^ dy / (2z - xy) and normal form (mu u (1 + b uv), v/mu (1 - b uv)), gamma = pi Im(b) in all 8 cases,
    gamma / Im(b) between 3.141577 and 3.141625, which is the precision of the committed gamma. This is
    Bohr-Sommerfeld for H(I) = I log mu + (b/2) I^2 with hbar = 2 pi i / r. b is purely imaginary for LR
    (mu real) and complex for LLR. Independent review: no blockers; b is unchanged to 15 digits under a
    different chart (y,z) and arbitrary eigenvector scaling. Only LLR separates Im(b) from abs(b).
T3. Phase test, needs WP2 item 1: arg(lambda_n) against r CS/4pi - (n+1/2) arg(mu), mod the root-of-unity
    ambiguity of the representation. T2 adds a prediction: an (n+1/2)^2 phase term -pi Re(b) / r, zero for
    LR and nonzero for LLR (sign fixed by the choice of conjugate fixed point).
T4. Torsion check: abs((1-mu)(1-1/mu)) against published adjoint torsion values for LR and LLR (C3, L6).
T5. Stretch goal: alpha from the 2-loop invariant. Only attempted if T2 succeeds.

Tests: T2 and T4 are subcommands of ladder_tools with an assertion and a committed output table. T2 uses
only the 8 committed gamma values; it needs no solver run. T3 needs WP2 item 1.

### WP4. The paper

Sections:
1. Introduction and Conjecture (ladder, in the F2 form; with the gamma closed form if T2 passes).
2. Setup: the SO(3) representation, colour and root conventions (F4, F24).
3. Classical side: cone structures, trace map, mu, torsion (F28, C3, C4).
4. Heuristic: trace formula and normal form (T1, T2).
5. Numerical method (solver/README.md:35-52; WP1 error estimates).
6. Evidence (WP0 residual plot, WP2 tables, T3 phases).
7. Relation to Chen-Yang, Wong-Yang, BWY, AMU (L1 to L8); open problems: proof, higher genus, alpha.

Rule: a claims table is written first. Each sentence that asserts a number or a prior result points to a
data file plus script, or to a literature item L1 to L8 after it was read at the source.

## 4. Decisions needed

| # | Decision | Default if no answer |
|---|---|---|
| D1 | Scope: numerical conjecture paper with a heuristic section, or hold for a proof in the LR case | Conjecture paper |
| D2 | Runs come in two kinds. (a) ladder_tools: seconds, existing data only, any machine with g++. (b) Solver jobs for WP1 tests and WP2: Linux with FLINT, since ladder3 does not build on native Windows (F20). Who runs each. | No default; (a) blocks WP0 and T2, (b) blocks WP1 tests and WP2 |
| D3 | Do ladder.c, dev/ref/, the r = 201 to 1001 raw rows, the r = 7001 rows and the original fit and extrapolation scripts still exist somewhere | Assume no: strip references, rewrite scripts |
| D4 | Where the LaTeX lives: this public repo, or a separate private repo that cites this one | Separate private repo |

## 5. Self-review

- F3 is the one substantive error found in the current statement; C1 and C2 both support the n+1/2 reading.
  If the README's j formula was intended with a different mu normalisation, WP0 item 1 is wrong; needs a
  yes from the author before the edit.
- C3 and C4 use rounded CSV values. They become real tests only in WP0.
- T0 and T2 are hypotheses of this plan, not findings. T2 is designed so the unknown Planck-constant
  normalisation cancels (a ratio across 8 cases), but it can still fail if the 1/r term mixes in
  theta-dependence through c+1 (F4, F11).
- The dilogarithm volume formula for the cone structures is not yet read from a source. If it does not
  reproduce the 16 committed Vol values, the fallback is a one-time Vol table from SnapPy, which is the only
  Python left.
- What the paper can already say from committed data: the conjecture in the F2 form for LR and LLR, moduli
  only, r up to 10001. Runs are needed for three things only: the fits have no generating code (F11), the
  gamma closed-form test (T2), and evidence beyond two manifolds and beyond moduli (WP2).
- Whether mu is real or complex for LLR is not established here; T2 and T3 must handle both.
- WP1 item 4 rests on the review's reading that L^p R^q needs only theta^-p and Theta^q; confirmed for the
  parser (F13) and the precision constants (F14), not for every checkpoint path.
- Novelty: a shallow search found no eigenvalue-level ladder in print. That is not a literature review; a
  proper one is part of WP4 before the introduction is written.
- Nothing in this plan was executed. No commit, push or PR is implied.
