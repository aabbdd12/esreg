# Brute-force check of the influence functions (1 October 2026)

`ifexact.do <fiml|twostep> [vce(robust)]` computes, on the first 1,000
observations of `data/esr_sim.dta`, the exact influence function of ATT, ATU,
ATE and the mean kappa: the direct sampling term of each average plus the
Jacobian times d(theta)/d(w_i), with d(theta)/d(w_i) = V_oim s_i (FIML) or
-Jac^{-1} g_i (two-step, stacked moments). It prints the current standard
errors (delta + sampling, added) next to sqrt(sum D_i^2) and saves the D_i.

`bf.do <method> <first> <last>` re-runs the whole esreg estimation with the
weight of one observation moved to 1.5 and to 0.5 (iweights): the central
difference is the influence of that observation, by brute force.

Result (observations 1-100): brute force = exact influence function, relative
gap 2-4e-5 (two-step) and 7e-5 to 6e-4 (FIML, the convergence tolerance of ml),
correlation 1.0000000. Ratio current / exact standard error (n = 1,000):

| | ATT | ATU | ATE | kappa |
|---|---|---|---|---|
| two-step | 0.990 | 0.976 | 1.029 | 1.000 |
| FIML, vce(robust) | 0.989 | 0.974 | 1.020 | 1.000 |

The current formula adds the delta-method variance and the sampling variance
as if independent; their covariance moves the standard errors of the effects
by 1-3 percent. Outputs in `out/` (reference values, kept).

## Extended brute force (2 October 2026)

`bf2.do <case> <first> <last>` moves the weight of one observation (1.5, 0.5;
iweights), re-runs the whole estimation and `esrmte if x > 0` (an analysis
subsample, so that esrmte's own sampling part is tested), and records the
central differences of ATT, ATU, ATE, the mean kappa, esrmte's m and kappa
and, case `kx`, the slope of kappa(x). `ifcode.do <case>` stores the influence
functions that the code itself computes (`_esreg_effects`: sampling part U +
parameter part P; `_esreg_mte`; `_esr_psi` for the slope) with the reported
standard errors; `compare2.do <case>` compares the two. The cases
(`bf2_setup.do`), first 1,000 observations:

| case | data | estimation |
|---|---|---|
| `ts` | `esr_sim` | two-step |
| `fiml` | `esr_sim` | FIML |
| `kx` | `esr_kx` | two-step, `kappa(x)` |
| `het` | `esr_kx` | FIML, `hetsigma(x) hetrho(x)` |
| `hm` | `esr_nonlin` | two-step, `hermite(3)` (with the Hermite terms dh2, dh3 of esrmte's curve) |
| `hr` | `esr_nonlin` | two-step, `hermite(3 0)`: the Hermite terms in the treated regime only |

Result (observations 1-100): slope of BF on IF − 1 at most 5e-5 (two-step) and
3.4e-4 (FIML: the convergence tolerance of ml, the BF being a difference of
two maximizations); largest gap relative to the largest influence 2e-5 to
1.9e-4 (two-step), up to 1.9e-3 (FIML, heteroskedastic); correlations
≥ 0.9999998. The reported standard errors equal √ΣD² (two-step) and
√ΣD² × 1.0005 under `vce(robust)` (FIML: ml's n/(n−1) on the parameter part).
Outputs in `out/` (`bf2_*`, `ifcode_*`); `tests/test_effects_se.do` section 8
recomputes the influence functions of the current code and checks them
against the stored brute force.

Case `hr` (3 October 2026, the Hermite terms by regime): slope of BF on IF − 1
at most 4.9e-5, largest gap 7.5e-5 of the largest influence, correlations
1.0000000, reported standard errors = √ΣIF² exactly, for ATT, ATU, ATE, kappa,
esrmte's m and kappa, dh2 and dh3.
