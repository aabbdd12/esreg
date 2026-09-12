# esreg — endogenous switching regression with heterogeneous selection on gains

`esreg` is a Stata command (Stata 16 or later) that fits the endogenous
switching regression (Roy model) by full-information maximum likelihood or by a
two-step estimator with the exact variance of the whole procedure, and reports
the treatment effects ATT, ATU and ATE together with the selection-on-gains
parameter κ = ρ₁σ₁ − ρ₀σ₀, constant or heterogeneous in the covariates.
It comes with a family of post-estimation commands:

| command | what it does |
|---|---|
| `esreg` | estimation (FIML or two-step), effects, κ, κ(x), common support; `svy:` prefix (FIML), survey weights, `hsize()` |
| `predict` | predictions (effect, potential outcomes, score, Mills ratios, κ_i, ρ_jσ_j) and equation-level scores |
| `esrdiag` | strength, variation and support of the selection equation, with warnings |
| `esrtest` | specification tests: index (`spec`), sufficiency of the instruments (`suff`), constant κ (`kappa`), normality by regime (`normal`), pseudo-DiD on the index (`pdid`) |
| `esrcurve` | expected effect by quantile group of the score |
| `esrmte` | marginal treatment effect: parametric line and semiparametric curve |
| `esrreport` | a reading of the results with rule-based notes and the route to report |

Every printed number is validated against a Python reference implementation
and against Stata's `etregress, poutcomes`, `movestay` and `msat`. The
methods are described in

> Araar, A. (2026). *Endogenous Switching Regression with Heterogeneous
> Selection on Gains: Assumptions, Tests, and Estimation.* Working paper,
> Université Laval and PEP.

## Installation

From Stata:

```stata
net install esreg, from("https://raw.githubusercontent.com/aabbdd12/esreg/main/") replace
help esreg
```

(The package will also be submitted to the SSC archive: `ssc install esreg`.)

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

Heterogeneous selection on gains: `esreg y x, select(d = x z) method(twostep) kappa(x)`
then `esrtest, kappa`. Survey design: `svyset psu [pw = wt], strata(s)` then
`svy: esreg y x, select(d = x z)`.

## Files

`esreg.ado` (estimation), `esreg_lf1.ado` (likelihood evaluator),
`esreg_engine.ado` (loader), `esreg_mata.ado` (the Mata engine), `_esreg_getest.ado`,
`_esreg_data.ado`, `_esreg_effects_post.ado` (internal helpers), `esreg_p.ado` (predict),
`esrdiag.ado`, `esrtest.ado`, `esrcurve.ado`, `esrmte.ado`, `_esreg_pwr.ado`
(percentile-weights regression engine), `esrreport.ado`; help files `*.sthlp`;
`esreg_returns.txt` (the layout of `e()`); `esreg.pkg`, `stata.toc`.

The replication package of the paper (Python reference library, Monte Carlo
scripts, simulated samples, do-files and logs of every table and figure) is a
separate repository.

## License

MIT. Author: Abdelkrim Araar, Université Laval and Partnership for Economic
Policy (PEP).
