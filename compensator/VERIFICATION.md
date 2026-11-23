# Formal verification ledger

The analytic proof in `paper.tex` and its PDF covers the full recorded problem,
parts (a)-(d), including maturity. **The full compensator theorem is not yet
formalized in Lean.** The compiler checks precisely the statements included in
the source; their scope is itemized below.

## Pinned environment

- Lean 4.36.0.
- Mathlib v4.36.0, commit `fc4f59e577feb1a4956671bc04ac18aaaf8cd842`.
- 39 theorem statements across three modules, also provided as one combined file.
- Actual compilation with warnings promoted to errors; exact axiom reports in `logs/`.
- No added axioms, `sorry`, `admit`, or unsafe proof shortcuts in the proof source.
- Axioms reported by Lean are the usual foundational `propext`, `Classical.choice`,
  and `Quot.sound` (a subset may suffice for particular theorems).

## Compiled coverage

| Source | Count | Exact scope |
|---|---:|---|
| `lean/CompensatorCore.lean` | 21 | Actual original-payoff atoms, measurability, growth, finite double-sum growth, equivalence of the original indicator, actual spatial primitive derivatives and integration by parts under explicit smoothness hypotheses, generator coefficient, pointwise energy inequality, evaluated time-integral majorant, Brownian finite moments, original terminal payoff in Lp, conditional-payoff martingale. |
| `lean/HeatKernel.lean` | 8 | Actual first/second derivatives and forward/backward heat equations of the logarithmic Gaussian kernel. |
| `lean/QuadraticTest.lean` | 10 | Gaussian second moment and convolution of the quadratic test payoff, actual derivatives/backward heat equation for that test payoff, and supporting moments/bounds. This test is explicitly distinct from the recorded nonsmooth payoff. |

`compensated_path_identity` is a scalar algebra lemma conditional on `P = K + I`;
it does not establish the stochastic decomposition. Likewise,
`terminal_payoff_condExp_martingale` proves a martingale of actual conditional
expectations, not that the proposed `P-K` equals that process.

## Open obligations for the full Lean theorem

1. Prove original-payoff Gaussian convolution smoothing and derivative interchange
   with all domination hypotheses discharged; identify its explicit formula and
   its conditional-expectation version.
2. Construct the predictable L2 Ito integral and prove the isometry and martingale
   theorem for it in Lean.
3. Prove/apply Ito's formula to this time-space primitive, including localization;
   identify `P-K` with the constructed stochastic integral.
4. Prove the original integrand's integrated L2 estimate, compensator
   predictability/finite variation, terminal L2 limit, and the final closed-interval
   `Martingale` theorem, together with uniqueness in the relevant class.

The repository sources inspected for this run contain Brownian laws and
martingale infrastructure. They do not supply a ready-made Ito integral/Ito
formula for this functional. The package adds no assumptions asserting these
missing conclusions. A complete mathematical proof is provided independently in
the paper.

## Numerical verification

`numerics.py` uses the recorded payoff with explicitly chosen illustrative
parameters, not a quadratic replacement: T=1, n=3, a=(-0.75,0,0.85).
It uses 8192 paths, 2048 fine time steps, seed 20261123 and nested coarser grids.
The deterministic checks use independent quadrature. The Monte Carlo diagnostics
include sampling uncertainty and have no role in certifying a general theorem.

See `data/metrics.json`, the CSV tables, the six figure pairs and the run log for
the actual numerical results. Source integrity hashes are in `SHA256SUMS.txt`.
