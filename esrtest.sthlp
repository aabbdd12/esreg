{smcl}
{* *! version 1.0.0  12sep2026}{...}
{vieweralsosee "esreg" "help esreg"}{...}
{vieweralsosee "esrdiag" "help esrdiag"}{...}
{vieweralsosee "esrreport" "help esrreport"}{...}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{bf:esrtest} {hline 2}}Specification tests after esreg{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:esrtest} [{cmd:,} {it:subtest} {opt est(name)} {opt l:evel(#)}]

{synoptset 28}{...}
{synopthdr:subtest}
{synoptline}
{synopt:{opt spec} [{opt ord:er(#)} {opt add(varlist)}]}specification of the selection index: link test and gamma contrast{p_end}
{synopt:{opt suff} [{opt keep(varname)}]}index sufficiency of the excluded instruments{p_end}
{synopt:{opt kap:pa}}constant selection on gains, H0: kappa(x) = kappa (the default subtest){p_end}
{synopt:{opt norm:al}}normality by regime: Hausman contrast FIML vs two-step, and Hermite controls{p_end}
{synopt:{opt pdid} [{opt nq(#)}]}pseudo-difference-in-differences on strata of the selection index{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:esrtest} runs one of the specification tests of the paper on a stored
{helpb esreg} estimation. Every subtest works from the probit index and the two
regime regressions; those that re-estimate do so on the sample of the stored
estimation and restore the user's {cmd:e()} on exit.

{pstd}
{cmd:spec} (a) refits the probit of D on the fitted index v = Zg and its powers
v^2 ... v^{it:order} (default 2; the link test of Pregibon): under a correctly
specified index the powers have zero coefficients; Wald chi2(order - 1). With
{cmd:add(}{it:varlist}{cmd:)} the listed terms (polynomials, interactions) are
added to the probit instead and tested by Wald and LR. (b) contrasts the probit
and the FIML estimates of the selection coefficients gamma with the
influence-function covariance of the difference: this is the gamma block of the
Hausman contrast of {cmd:normal}, what trivariate normality imposes on the
selection equation. A rejection of (b) with a rejection of (a) points to the
index; with a clean (a), to the joint law.

{pstd}
{cmd:suff} tests that the excluded instruments carry no information on the regime
errors beyond the index, E[w_j | X, Z, u] = E[w_j | X, u]: the excluded
instruments other than {cmd:keep()} (default: the strongest in the probit) are
added to the two regime regressions beside the Mills ratio; under H0 their
coefficients are zero. Wald on the stacked variance, joint and by regime. A
rejection says that an instrument has a direct effect on the outcome (the
exclusion is wrong) or that participation has a second index. With a single
excluded instrument the test is not available (the instrument would test itself).

{pstd}
{cmd:kappa} tests that selection on gains is constant in the covariates. After
{cmd:method(twostep) kappa(W)}, H0 is that the non-constant entries of
theta_1 - theta_0 are zero (Wald on the stacked variance; the table gives each
slope d kappa/d w). After {cmd:method(fiml)} with {cmd:hetsigma()} or
{cmd:hetrho()}, H0 is that all non-constant coefficients of ln sigma_j and
atanh rho_j are zero, a sufficient condition for a constant kappa.

{pstd}
{cmd:normal} has two pieces. (a) The Hausman contrast between the FIML and the
two-step estimates of q = (gamma, b_1, rho_1 sigma_1, b_0, rho_0 sigma_0), with
the covariance of the difference built from the influence functions of the two
estimators (positive semi-definite, valid under the alternative), chi2 on the rank,
and the contrast on kappa with its standard error: it rejects when the regime
errors are not normal (the likelihood is then inconsistent, the two-step is not).
(b) The Hermite controls E[u^2 - 1 | D, Z] and E[u^3 - 3u | D, Z], added beside
the Mills ratio in each regime regression, with a Wald chi2(4) joint and chi2(2)
by regime: they reject when E[w_j | u] is not linear in u, in which case both
routes fit the wrong curve and the semiparametric MTE of {helpb esrmte} is the
estimate to report. Skewed regime errors with a linear conditional mean leave (b)
clean and make (a) reject.

{pstd}
{cmd:pdid} cuts the sample into {cmd:nq()} strata of the probit index (default 2),
computes the X-adjusted treated-untreated gap in each stratum, and the two
within-group contrasts between the top and bottom strata: C1 among the treated,
whose expectation is rho_1 sigma_1 times the contrast of mean lambda_1, and C0
among the untreated, whose expectation is -rho_0 sigma_0 times the contrast of mean
lambda_0. Their signs read the signs of rho_1 sigma_1 and rho_0 sigma_0 with the
linear conditional mean alone (no normality); the pattern C1 < 0, C0 > 0 is the
signature of selection on gains, the same sign a common selection on the level of
the outcome. kappa_dd = C1/dlambda_1 + C0/dlambda_0 is a consistent estimator of
kappa that never uses the joint law. Never cut the strata on the outcome.


{title:Options}

{phang}{cmd:est(}{it:name}{cmd:)} uses the estimation stored under {it:name}.{p_end}
{phang}{cmd:order(#)} is the highest power of the index in the link test (default 2).{p_end}
{phang}{cmd:add(}{it:varlist}{cmd:)} tests the listed terms added to the probit instead of the powers.{p_end}
{phang}{cmd:keep(}{it:varname}{cmd:)} is the excluded instrument left out of the sufficiency test.{p_end}
{phang}{cmd:nq(#)} is the number of strata of the index in the pseudo-DiD (default 2).{p_end}
{phang}{cmd:level(#)} sets the confidence level.{p_end}


{title:Stored results}

{pstd}All subtests store {cmd:r(chi2)}, {cmd:r(df)}, {cmd:r(p)} and {cmd:r(test)}.
In addition: {cmd:spec} stores {cmd:r(link)} (est se z p of the tested terms),
{cmd:r(gamma)} (probit, fiml, se_diff by coefficient), {cmd:r(chi2_g)},
{cmd:r(df_g)}, {cmd:r(p_g)}, {cmd:r(lr)}, {cmd:r(p_lr)}, {cmd:r(tested)};
{cmd:suff} stores {cmd:r(chi2_1)}, {cmd:r(df_1)}, {cmd:r(p_1)}, {cmd:r(chi2_0)},
{cmd:r(df_0)}, {cmd:r(p_0)}, {cmd:r(coefs)}, {cmd:r(keep)}, {cmd:r(tested)};
{cmd:kappa} stores {cmd:r(slopes)} (two-step); {cmd:normal} stores {cmd:r(q)},
{cmd:r(hermite)}, {cmd:r(kappa_2s)}, {cmd:r(kappa_fiml)}, {cmd:r(kappa_diff)},
{cmd:r(kappa_se)}, {cmd:r(chi2_h)}, {cmd:r(df_h)}, {cmd:r(p_h)}, {cmd:r(chi2_h1)},
{cmd:r(chi2_h0)}; {cmd:pdid} stores {cmd:r(strata)}, {cmd:r(C1)}, {cmd:r(se_C1)},
{cmd:r(C0)}, {cmd:r(se_C0)}, {cmd:r(dlam1)}, {cmd:r(dlam0)}, {cmd:r(rhosig1)},
{cmd:r(rhosig0)}, {cmd:r(kappa_dd)}, {cmd:r(se_kappa_dd)}, {cmd:r(dd)},
{cmd:r(se_dd)}, {cmd:r(nq)}.


{title:Examples}

{phang2}{cmd:. esreg y x, select(d = x z1 z2) method(twostep)}{p_end}
{phang2}{cmd:. esrtest, spec}{p_end}
{phang2}{cmd:. esrtest, spec add(c.z1#c.z1 c.x#c.z1)}{p_end}
{phang2}{cmd:. esrtest, suff keep(z1)}{p_end}
{phang2}{cmd:. esrtest, normal}{p_end}
{phang2}{cmd:. esrtest, pdid nq(4)}{p_end}
{phang2}{cmd:. esreg y x, select(d = x z1 z2) method(twostep) kappa(x)}{p_end}
{phang2}{cmd:. esrtest, kappa}{p_end}
