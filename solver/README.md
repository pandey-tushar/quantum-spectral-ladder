# ladder2

Faster, lower-memory replacement for `ladder.c`. It computes the same quantity: the top `ntop` values of
log|lambda| for A = T_b^-1 T_a (LR) or T_b^-2 T_a (LLR), using the SO(3) once-punctured-torus quantum
representation at q = e^{2 pi i/r} with puncture color c. The stdout line has the same format as ladder.c.
`ladder.c` is unchanged.

## Build

    gcc -O3 -march=native -fopenmp ladder2.c -o ladder2 -lflint -lmpfr -lgmp -lm

Tested with FLINT 3.0.1 and gcc on Linux. On macOS with clang, use `-Xpreprocessor -fopenmp -lomp`, or
build with Homebrew gcc.

## Run

    ./ladder2 r c LR|LLR [options]
    OMP_NUM_THREADS=32 ./ladder2 10001 5000 LR > out.txt 2> log.txt

| option | meaning |
|---|---|
| `-t T` | threads (default: `OMP_NUM_THREADS`, or all cores) |
| `-kernel mm` | default. K is stored as residues mod 59-bit primes, and products use parallel `nmod_mat_mul` plus CRT |
| `-kernel acb` | K is stored as `acb_mat`: about half the memory, with a 5x slower iteration stage |
| `-kernel mmsym` | experimental: stores only the lower half of the symmetric K J. About half the memory, but ~1.8x slower because it needs a more accurate K (higher d1, d2) |
| `-ck FILE` | checkpoint. X (and the final d3) are written after every iteration. A rerun with the same FILE resumes from it |
| `-m 12 -ntop 8 -tol 22 -maxit 100` | as in ladder.c, plus an iteration cap |
| `-d1/-d2/-d3 D` | override the automatic precisions (decimal digits) of the three stages |
| `-q` | suppress per-iteration lines on stderr |

stderr reports each stage's precision, its accuracy self-check, the timings and the peak RSS. Lines that
start with `# WARNING` mean a self-check failed after the automatic retries. Treat those results as
suspect.

## Method

1. Y (the tridiagonal loop operator) has real alpha. Its beta are real except for a few purely imaginary
   entries. A diagonal similarity S = diag(i^e_i) makes it real: Y' = S^-1 Y S, with J = S^2 = diag(+-1).
2. **Stage 1 (precision d1, O(n^2), parallel over columns).** Each eigenvector w_b of Y' is computed by
   twisted factorization in real arithmetic. The eigenvectors are never stored. Only rows 0 and n-1 of
   K = S^-1 T_b^-p S = sum_b theta_b^-p w_b w_b^T J / (w_b^T J w_b) are accumulated.
   - Self-check: the same sum with theta = 1 must give rows 0 and n-1 of the identity (the completeness
     residual).
3. **Stage 2 (precision d2, O(n^2)).** K commutes with Y', so every row of K follows from rows 0 and n-1
   by a three-term recurrence. It runs from both ends and meets at the largest |alpha|.
   - Self-check: the mismatch of the two sweeps at the meeting row.
4. **Stage 3 (precision d3).** Block subspace iteration with A' = K Theta, which is similar to A. Each
   iteration does one n x n x m product, also for LLR. The full-mm kernel stores K as per-row fixed point,
   reduced mod about 0.115 d3 primes of 59 bits each. A product is then one `nmod_mat_mul` per prime,
   run in parallel, followed by CRT. log|lambda_1| is read off after 3 iterations. If d3 is too low for
   it, K is rebuilt at a higher d3 and the iteration continues.

### Why the old code needed 0.6 r digits

Computing V's columns stably does not help. The twisted factorization was only about 3 digits better than
the forward recurrence at r = 201, 1001 and 2001.

The real causes:
- V is J-orthogonal (V = S W, W^T J W = I), and its entries reach 10^maxV, with maxV ~ 0.126 r at
  c ~ r/2.
- Computing V from Y loses about 2 maxV digits, because the eigenvectors are badly conditioned.
- Applying V D V^T costs about 2 maxV - log10|lambda| more digits.

ladder2 never multiplies by V. It forms T_b^-p once, through the commuting recurrence, and only that
matrix enters the iteration. Its entries are at most 10^maxV (LR) or about 10^(1.4 maxV) (LLR).

### Precision laws (digits, measured)

| stage | cost | law used (c ~ r/2) | r=2001 | 3001 | 5001 | 7001 |
|---|---|---|---|---|---|---|
| old ladder.c (everything) | | 0.6 r + 20 | 1220 | | | |
| d1 eigenvectors | O(n^2) | ~3 maxV - 0.02 r + 40, about 0.37 r + 40 | 782 | 1153 | 1895 | 2637 |
| d2 recurrence | O(n^2) | 0.06 r + 30 + (1.77 - 1.43 c/r) maxV (+ LLR term) | 434 | 627 | 1014 | 1401 |
| d3 iteration and stored K | O(n^2 m its) | log10 max\|K\| - log10\|lambda_1\| + a r + 30, with a = 0.068 (LR) or 0.098 (LLR) | 290 (LLR 392) | 419 | 664 | 928 |

- For LR at theta ~ 0, d3 is about 0.126 r + 35. That is roughly 5 times fewer digits than the old 0.6 r.
- For LLR, d3 is about 0.19 r.
- All three precisions adapt at run time. maxV is measured in stage 1. A stage is rerun with more
  digits if its self-check fails.

## Validation (agreement with reference values, significant digits of log|lambda_1..8|)

| run | result |
|---|---|
| r = 201, c = 100, 60, 140, LR and LLR (references from ladder.c) | >= 28.6 digits |
| r = 1001, c = 500, 300, LR and LLR (references from ladder.c at 840 digits) | >= 28.9 digits |
| r = 2001, c = 1000, 840, 1318, 522, LR and LLR (r2001.txt, b1-b3.txt) | >= 28.3 digits |
| r = 3001, 5001, 7001 (LR, c = (r-1)/2), rerun with all precisions raised by 60-200 digits | identical to 30 digits |

The r = 1001, 2001 and 7001 references print 30 digits, so agreement of about 28-30 digits is the
maximum this comparison can show.

## Timings on this machine (2 cores, Xeon 2.1 GHz; LR, c = (r-1)/2; default mm kernel)

| r | n | d1/d2/d3 | stage 1 | stage 2 | stage 3 (its) | total | peak RSS |
|---|---|---|---|---|---|---|---|
| 2001 | 500 | 782/434/290 | 2.4 s | 1.2 s | 3.5 s (12) | 6.5 s | 0.16 GB |
| 3001 | 750 | 1153/627/419 | 8.2 s | 4.4 s | 9.6 s (13) | 22 s | 0.46 GB |
| 5001 | 1250 | 1895/1014/664 | 45 s | 20 s | 34 s (13) | 99 s | 1.9 GB |
| 7001 | 1750 | 2637/1401/928 | 181 s | 101 s (+75 s rebuild) | 98 s (13) | 456 s | 5.1 GB |

- At r = 2001, the 8 reference (word, c) cases take 5-18 s each, against 4-8 min for ladder.c. One thread
  takes 13 s for c = 1000 LR, against 6.5 s on two.
- The stage 2 rebuild in the r = 7001 run no longer happens: the max|K| estimate was raised from 1.00 to
  1.02 maxV.

## Projection for r = 10001 (per (word, c))

- Exponents fitted from 5001 to 7001: stage 1 ~ r^4, stage 2 ~ r^4.5, stage 3 ~ r^3.
- LR at theta ~ 0 (c = 5000, n = 2500): about 3100 core-seconds, or roughly 52 core-minutes.
  - Expected wall time is ~4 min on 16 cores, ~2 min on 32 and ~1-1.5 min on 64 (assuming ~80%
    parallel efficiency; stages 1 and 3 are embarrassingly parallel).
  - LLR is about 1.5x that.
  - Colors with theta ~ 3 (c ~ 0.26 r, n ~ 0.37 r) are about 2x.
- Memory is dominated by the stored K: about 2 x 0.115 d3 x 8 n^2 bytes (mm kernel).

| case at r = 10001 | n | d3 | mm (default) | mmsym / acb |
|---|---|---|---|---|
| LR, theta ~ 0 | 2500 | ~1320 | ~15 GB | ~8 GB |
| LLR, theta ~ 0 | 2500 | ~1850 | ~21 GB | ~11-12 GB |
| LR, theta ~ 3 | 3700 | ~1390 | ~35 GB | ~18 GB |
| LLR, theta ~ 3 | 3700 | ~1870 | ~47 GB | ~23-24 GB |

stderr prints the K storage size before allocating it. If it does not fit in RAM, use `-kernel acb`
(slower iteration) or `-kernel mmsym`.

ladder.c would need 6000+ digits and a dense V (tens of GB) at r = 10001, and an estimated days of
single-core time.

## Files

| file | contents |
|---|---|
| `ladder2.c` | the program |
| `rbig.txt` | ladder2 results for r = 3001, 5001, 7001 (LR, theta closest to 0) |
| `dev/` | development tests: `cmp.py` (digit comparison against reference files), `ref/` (reference outputs), validation logs `val2001_final.txt` and `r*_LR.err` |

## ladder3: bounded-memory version

`ladder3.c` computes the same numbers with the same stdout line; `ladder2.c` is unchanged.

    gcc -O3 -march=native -fopenmp ladder3.c -o ladder3 -lflint -lmpfr -lgmp -lm -lpthread
    ./ladder3 10001 5000 LR -t 32 -scratch /fast/disk/run1 -mem 1.5 > out.txt 2> log.txt

| option | meaning |
|---|---|
| `-scratch DIR` | K goes to `DIR/K.bin` and is streamed every iteration (I/O thread prefetches the next chunk). Also checkpoints stage 1 (`stage1.state`), stage 2 (`build.state`, resumes mid-recurrence), the finished K (`K.meta`) and the iteration (`iter.ck`). Rerun the same command to resume. Use one DIR per (r, c, word). |
| `-mem GB` | RAM for K: resident part plus two stream buffers (default 1.5 with `-scratch`; without `-scratch` everything is resident) |
| `-nocache` | drop K.bin pages from the page cache after reading |
| `-ckint SEC` | minimum time between stage-2 checkpoints (600) |
| `-levels NL -rb RB -cb CB -segmb MB` | storage layout (6, 32, 32, 64) |

Other options as ladder2 (`-kernel` is not available).

Changes:
1. K is rounded at one absolute tolerance 2^(log2 max\|K\| - BK) instead of per row (validated: same results). K
   decays by hundreds of digits away from its central peak, so each 32x32 tile is stored only with the
   bits it needs: residues mod np(b) primes with np(b) ~ 2b/58, NL = 6 levels, Y truncated to the same level,
   tiles below the tolerance dropped, thin levels merged by a cost model. One CRT per level and row block.
   Storage relative to ladder2: 0.74 (LR, theta~0), 0.84 (LLR, theta~0), 0.43/0.59 (theta~2, c<r/2),
   0.48/0.56 (theta~3, c<r/2, r=2001), ~1.0 for c>r/2 (small n). Product work drops by the same factor.
2. Out-of-core K (above). Residues are kept (not compact integers) because converting compact integers to
   residues costs ~3x the product itself. K J is symmetric, but storing half of K needs a much more accurate K
   (ladder2's mmsym note), so it is not used.
3. Stage 1 accumulates into shared, chunk-locked sums: 2 vectors of n per thread instead of 9. The H = X^H Z
   step uses no per-thread copies.

Peak RSS is about `-mem` + 0.3-0.5 GB at r = 10001 (32 threads; model of the arrays in the code, matching the
measured 0.03-0.07 GB non-K part at r = 3001).

### Measured (this machine, 1 thread, shared with another 2-thread job; wall s / peak RSS GB)

| r, c, word | ladder2 | ladder3 (RAM) | ladder3 -scratch -mem 0.1 -nocache |
|---|---|---|---|
| 2001 1000 LR | 12.2 / 0.155 | 11.5 / 0.129 | 11.9 / 0.112 |
| 2001 1000 LLR | 17.7 / 0.201 | 17.1 / 0.193 | 16.7 / 0.132 |
| 2001 1318 LR (theta~2) | 9.3 / 0.076 | 8.5 / 0.075 | 9.5 / 0.075 |
| 2001 1318 LLR | 10.3 / 0.099 | 10.4 / 0.112 | 10.6 / 0.130 |
| 3001 1500 LR | 45.7 / 0.463 | 41.6 / 0.353 | 42.5 / 0.133 |
| 3001 1500 LLR | 66.5 / 0.608 | 56.9 / 0.541 | 61.3 / 0.157 |
| 3001 1022 LR (theta~2) | 81.2 / 0.754 | 70.9 / 0.352 | 69.7 / 0.154 |
| 3001 1022 LLR | 96.3 / 1.030 | 87.2 / 0.635 | 90.8 / 0.167 |

Validation: r = 201, 1001, 2001 (all 18 reference cases, RAM mode; r = 2001 also streamed) agree with the
references to >= 28.3 digits, as ladder2. r = 3001, all 8 cases of `big_all.txt`, streamed: identical to
ladder2 to 30 digits. Interrupting and rerunning (in stage 2 and in the iteration) gives identical output.
Multi-thread runs were not tested here (1 thread only); check e.g. `./ladder3 1001 500 LR -t 4` against
`dev/ref/r1001.txt` first.

### Projection for r = 10001

| case | ladder2 RAM | ladder3 K on disk | ladder3 peak RSS (-mem 1.5) | read per run |
|---|---|---|---|---|
| LR, theta~0 | ~15 GB | ~11 GB | ~2 GB | ~130 GB (13 its) |
| LLR, theta~0 | ~21 GB | ~18 GB | ~2 GB | ~170 GB (10 its) |
| LR, theta~3 | ~35 GB | ~17 GB | ~2 GB | ~200 GB |
| LLR, theta~3 | ~47 GB | ~26 GB | ~2 GB | ~250 GB |

- Compute is 0.85-0.95 of ladder2's core-seconds (~47 core-min for LR theta~0).
- Streaming sets a floor on each iteration of (K - mem) / disk bandwidth: about 5 s (LR theta~0) to 12 s (LLR theta~3)
  at 2 GB/s. Below ~6-10 cores that is hidden behind compute. On 32 cores the total is about 1.3-1.5x ladder2's
  wall time at 2 GB/s, and more on a slower disk.
- A larger `-mem`, or free page cache, cuts the reads proportionally.
