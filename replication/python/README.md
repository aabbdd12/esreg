# esreg — Python reference implementation

Endogenous switching regression with two routes and the treatment effects, written
as the specification of the future Stata command `esreg`. Every number the ado
will produce is first produced here and checked against `movestay` / `msat`.

## Layout

| File | Content |
|---|---|
| `esreg/likelihood.py` | Gaussian ESR log-likelihood per observation and **analytic score**, with `ln sigma_j = Ws a_j`, `atanh rho_j = Wr c_j` (constant by default); Hessian by differences of the score. Parameter order `[b1, b0, g, a1, a0, c1, c0]`. |
| `esreg/probit.py` | Weighted probit, analytic score and Hessian, Newton. |
| `esreg/fiml.py` | FIML: starting values (probit + OLS by regime), BFGS then Newton on the analytic score; `vce` = oim / opg / robust; LR test of rho1 = rho0 = 0. |
| `esreg/twostep.py` | Two-step: probit, OLS by regime on X and the Mills ratio (times `Wk` for kappa(x)); **exact variance of the whole procedure** from the stacked moment conditions (sandwich), replacing the bootstrap. |
| `esreg/effects.py` | Individual effects, ATT / ATU / ATE, kappa; delta-method SEs on the parameter covariance + sampling component; common support of P(Z). |
| `esreg/api.py` | `esreg(df, y, x, d, z, method, weights, hetsigma, hetrho, kappa, vce)` and `Results.summary()`. |
| `validate_esreg.py` | Runs the four datasets of the replication package and prints the Stata reference values next to the output (`validate_esreg.txt`). |
| `check_score.py` | Analytic score against finite differences (max abs. diff 3e-7 on a heterogeneous design). |
| `smoke_het.py` | `kappa(x)` in the two-step route and `hetsigma()/hetrho()` in the likelihood on a design with kappa(x) = 0.9 + 0.4x (runs; not yet validated against Stata). |

## Use

From this folder:

```python
import pandas as pd, sys
sys.path.insert(0, ".")
from esreg import esreg
from paths import DATA     # the data folder (../data/)
df = pd.read_stata(DATA + "esr_sim.dta")
r = esreg(df, "y", ["x"], "d", ["x", "z"], method="fiml")      # or method="twostep"
r.summary()
r.table()      # coefficients as a DataFrame
r.effects()    # ATT, ATU, ATE, kappa
```

Requirements: numpy, scipy, pandas.

## Validation

On `esr_sim` and `esr_annex2` the FIML reproduces `movestay` (SJ build 2.0.0) to every
printed digit: coefficients, standard errors, log-likelihood, LR statistic; and `msat` 1.3.2
for ATT / ATU / ATE / kappa, their standard errors, the two variance components and the
support. The two-step reproduces the `msat, twostep` point estimates on `esr_sim`,
`esr_skew`, `esr_nonlin`. On `esr_skew` the stacked-moment standard errors of the two-step
(0.0517, 0.0723, 0.0459, 0.0675) sit next to the bootstrap-100 values of the note
(0.0539, 0.0747, 0.0446, 0.0742), within bootstrap noise.

Checked later against the Stata commands: `kappa()`, `hetsigma()`/`hetrho()`, the tests and
the pseudo-DiD (`../tests/esrtest_test.do`), the weights and the robust variance
(`../tests/esreg_wt_test.do`, `export_wt.txt`).
