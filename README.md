# quantum-spectral-ladder

Numerical study of the spectrum of SO(3) quantum representations of the once-punctured torus mapping class group at the Chen–Yang root q = e^{2πi/r}.

For a pseudo-Anosov word φ in L, R and puncture color c, let λ_1, λ_2, ... be the eigenvalues of ρ_{r,c}(φ) ordered by modulus. The data support

    |λ_j| = exp(r Vol(M_φ(θ)) / 4π) · |μ(θ)|^{-(j+1/2)} · (1 + O(1/r)),

where
- θ = |2π − 4π(c+1)/r| is a cone angle along the core of the cusp (meridian ∂Σ),
- Vol(M_φ(θ)) is the volume of the hyperbolic cone manifold,
- μ(θ) is the expanding eigenvalue of dφ_* on the relative SL2(C) character variety at the geometric fixed point.

With n = j − 1 = 0, 1, 2, ..., the data at r = 201 to 10001 fit the refined form

    log|λ_n| = r Vol/4π − (n + 1/2) log|μ| + [α + γ (n + 1/2)^2] / r + O(1/r^2),

with α, γ given in `data/ladder_subleading_fits.csv` (closed forms not yet identified).

The n = 0 statement contains a spectral form of the generalized Bonahon–Wong–Yang volume conjecture; the spacing log|μ| is the linearization that governs the adjoint torsion (1-loop invariant).

## Evidence in this repo

- `data/r201-1001_summary.txt`, `data/r2001_raw.txt`, `data/r3001-5001_raw.txt`, `data/r10001_raw.txt`: LR and L²R (word LLR), θ ≈ 0, 1, 2, 3, r up to 10001.
- `data/stress_7words_r101.csv`: seven words, r ≤ 101, extrapolated.
- `data/ladder_subleading_fits.csv`: fitted α, γ per case.
- `data/extrap_r121.csv`: LR and LLR on a θ grid, r ≤ 121, extrapolated.
- `figures/`: error versus r; local exponents tend to −1.

Raw file format (one line per run): `r c n theta word iterations log|λ_1| ... log|λ_8|` (30 digits).

## Code

- `solver/ladder3.c`: production solver, C + FLINT 3 (arb/acb) + OpenMP; low memory with disk streaming and checkpoint/resume. See `solver/README.md`.
- `solver/ladder2.c`: in-memory version.
- `analysis/analyze_rows.py raw_files...`: computes Vol(θ) (SnapPy) and log|μ(θ)| (trace map) at the exact θ of each row and prints errors. Writes `analysis_out.txt`.
- `analysis/trace_map.py`: fixed points and linearization of the trace map.
- `analysis/plot4.py`: convergence figure.
- `analysis/reference_*.py`: slow Python/mpmath reference implementation (full spectrum, small r).

Build (needs FLINT 3, MPFR, GMP, OpenMP) and run:

    make            # builds ladder2 and ladder3
    ./ladder3 10001 5000 LR -t 8 -scratch /path/to/scratch -mem 1.5

Python dependencies: `pip install -r requirements.txt`.

## Status

All statements are numerical. Not yet established: a proof, higher genus fibers, and the closed form of the subleading 1/r constants of the ladder gaps.

## License

MIT, see `LICENSE`.
