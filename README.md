# esreg — endogenous switching regression with heterogeneous selection on gains

`esreg` is a Stata command (Stata 16 or later) that fits the endogenous
switching regression (Roy model) by full-information maximum likelihood or by a
two-step estimator with the exact variance of the whole procedure, and reports
the treatment effects ATT, ATU and ATE together with the selection-on-gains
parameter κ = ρ₁σ₁ − ρ₀σ₀, constant or heterogeneous in the covariates.
Version **1.0.0** (3 October 2026). It comes with a family of post-estimation
commands:

| command | what it does |
|---|---|
| `esreg` | estimation (FIML or two-step), effects, κ, κ(x), common support; the augmented two-step with Hermite terms when the conditional mean is not linear (`hermite()`, in both regimes or in one); survey weights, `hsize()`, `vce(cluster)`, the `svy:` prefix (FIML) and `vce(svy)` (two-step) |
| `predict` | predictions (effect, potential outcomes, score, Mills ratios, κ_i, ρ_jσ_j) and equation-level scores |
| `esrdiag` | strength, variation and support of the selection equation, with warnings |
| `esrtest` | specification tests: index (`spec`), sufficiency of the instruments (`suff`), constant κ (`kappa`), normality by regime with Hermite controls by regime (`normal`), pseudo-DiD on the index (`pdid`) |
| `esrcurve` | expected effect by quantile group of the score |
| `esrmte` | marginal treatment effect: parametric line (or curve, after `hermite()`) and semiparametric curve |
| `esrreport` | a reading of the results with rule-based notes and the route to report |

The standard errors of the effects come from the influence function of the
whole procedure (the estimation of the parameters, the averaging over the
units, and their covariance), aggregated as the estimation was: by
observation, by cluster or through the survey design. Every printed number is
validated against a Python reference implementation and against Stata's
`etregress, poutcomes`, `movestay` and `msat`, and every standard error
against a brute force. The methods are described in

> Araar, A. (2026). *Endogenous Switching Regression with Heterogeneous
> Selection on Gains: Assumptions, Tests, and Estimation.* Zenodo.
> <https://doi.org/10.5281/zenodo.22717029>

## Installation

From Stata:

```stata
net install esreg, from("https://raw.githubusercontent.com/aabbdd12/esreg/main/") replace
net get esreg, from("https://raw.githubusercontent.com/aabbdd12/esreg/main/") replace
help esreg
```

`net get` copies the example data (`esr_kx.dta`, `esr_nonlin.dta`,
`esr_wt.dta`) into the current folder. The examples of `help esreg` run from
their links (in the command window, in the dialog box filled in, or as a
do-file); without the copy, they read the data from GitHub. The package will
also be submitted to the SSC archive (`ssc install esreg, all`).

## Quick start

```stata
webuse union3
esreg ln_wage age grade smsa black tenure, select(union = south black tenure)
esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep)
esrdiag
esrtest, spec
esrtest, normal
esrtest, pdid nq(4)
predict double p, pr
esrcurve, rank(p) nq(5) graph
esrmte, semipar graph
esrreport
```

Heterogeneous selection on gains: `use esr_kx`, then
`esreg y x, select(d = x z) method(twostep) kappa(x)` and `esrtest, kappa`.
A conditional mean not linear in u: `use esr_nonlin`, then `esrtest, normal`
after the two-step and `esreg y x, select(d = x z) method(twostep) hermite(3 0)`.
Survey design: `use esr_wt`, `svyset [pweight = wt], strata(region)`, then
`esreg income educ, select(treatment = educ i.region inst) method(twostep) vce(svy)`
or `svy: esreg income educ, select(treatment = educ i.region inst)`.

## Files

`src/`: `esreg.ado` (estimation), `esreg_lf1.ado` (likelihood evaluator),
`esreg_engine.ado` (loader), `esreg_mata.ado` (the Mata engine), `_esreg_getest.ado`,
`_esreg_data.ado`, `_esreg_effects_post.ado`, `_esreg_esample.ado`,
`_esreg_ifcov.ado` (internal helpers), `esreg_p.ado` (predict), `esrdiag.ado`,
`esrtest.ado`, `esrcurve.ado`, `esrmte.ado`, `_esreg_pwr.ado` (percentile-weights
regression engine), `esrreport.ado`, `esreg_examples.ado` (the examples of the
help); `esreg.dlg` (the dialog box, opened by `db esreg`); help files `*.sthlp`;
`esreg_returns.txt` (the layout of `e()`). `examples/`: the example data.
`paper/`: the technical note. `replication/`: the replication material of the
paper (not copied by `net install`). `esreg.pkg`, `stata.toc` at the root.

## Replication

The folder [`replication/`](replication/) reproduces every table and figure of
the paper, in Python and in Stata, from the same data: the Python reference
library, the Monte Carlo scripts with their outputs, the simulated samples, the
do-files and their logs (see its README).

The first version of the paper (September 2026) was produced with esreg 0.4.1,
which stays installable from the tag `v0.4.1`:

```stata
net install esreg, from("https://raw.githubusercontent.com/aabbdd12/esreg/v0.4.1/") replace
```

and its replication package is <https://github.com/aabbdd12/esregrepl>.

## License

MIT. Author: Abdelkrim Araar, Université Laval and Partnership for Economic
Policy (PEP).
