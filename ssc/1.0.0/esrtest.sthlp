{smcl}
{* *! version 1.0.0  03oct2026}{...}
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
{synopt:{opt pdid} [{opt nq(#)} {opt reps(#)}]}pseudo-difference-in-differences on strata of the selection index{p_end}
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
routes fit the wrong curve. The augmented two-step
({cmd:esreg, method(twostep) hermite(}{it:#}{cmd:)}) fits it, of order 3 if the
cubic terms (u^3 - 3u, both regimes, Wald chi2(2), reported on a line of their own)
reject, else 2; the semiparametric MTE of {helpb esrmte} checks it. Skewed regime
errors with a linear conditional mean leave (b) clean and make (a) reject.

{pstd}
{cmd:pdid} cuts the sample into {cmd:nq()} strata of the probit index v = Zg
(weighted quantile groups, default 2) and computes the X-adjusted
treated-untreated gap in each stratum. In the two extreme strata, it regresses,
within each group j, the outcome on X and the top-stratum indicator: the
coefficient of the indicator is the contrast C1 among the treated and C0 among the
untreated. The same regressions with the Mills ratio lambda_j(v) in place of the
outcome give pi_1 and pi_0, the X-adjusted contrasts of the Mills ratios (pi_1 < 0,
pi_0 > 0). When the regime errors have a conditional mean linear in u,
C1 = rho_1 sigma_1 pi_1 and C0 = -rho_0 sigma_0 pi_0, whatever the distribution of
X across the strata: rho_1 sigma_1 = C1/pi_1, rho_0 sigma_0 = -C0/pi_0, and
kappa_dd = C1/pi_1 + C0/pi_0 estimates kappa without the joint law of the errors.
The signs of C1 and C0 need even less, a conditional mean increasing in u: the
pattern C1 < 0, C0 > 0 is the signature of selection on gains (kappa > 0); C1 and C0
of the same sign point to a common selection on the level of the outcome. Never
cut the strata on the outcome.

{pstd}
The standard errors of {cmd:pdid} come from the influence function of the whole
procedure: the probit, the boundaries of the strata (which move with the
estimated index) and the regressions, aggregated as in the estimation
(independent observations, by cluster, or through the survey design). With
{cmd:reps(#)} they come from a bootstrap of the whole procedure instead. Two notes
flag where the standard errors are approximate: a boundary between strata on a
block of tied values of the index that a resample can move (a discrete index:
the gaps and the double difference are the most affected, kappa_dd the least),
and a Mills contrast pi_j with a coefficient of variation above 10%, for which the
ratios C_j/pi_j and kappa_dd are far from normal while the sign of C_j stays
informative.


{title:Options}

{phang}{cmd:est(}{it:name}{cmd:)} uses the estimation stored under {it:name}.{p_end}
{phang}{cmd:order(#)} is the highest power of the index in the link test (default 2).{p_end}
{phang}{cmd:add(}{it:varlist}{cmd:)} tests the listed terms added to the probit instead of the powers.{p_end}
{phang}{cmd:keep(}{it:varname}{cmd:)} is the excluded instrument left out of the sufficiency test.{p_end}
{phang}{cmd:nq(#)} is the number of strata of the index in the pseudo-DiD (default 2).
With tied values of the index fewer strata may form; a note says so.{p_end}
{phang}{cmd:reps(#)} replaces the standard errors of the pseudo-DiD by a bootstrap of
the whole procedure with # replications (by cluster when the estimation was; use
{helpb set seed} for reproducible results). Not available after {cmd:svy}.{p_end}
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
{cmd:r(chi2_h0)}, {cmd:r(chi2_h3)}, {cmd:r(p_h3)} (the cubic terms); {cmd:pdid} stores {cmd:r(strata)} (n1 n0 v_mean p_mean gap
se_gap lambda1 lambda0 by stratum), {cmd:r(C1)}, {cmd:r(se_C1)}, {cmd:r(C0)},
{cmd:r(se_C0)}, {cmd:r(pi1)}, {cmd:r(se_pi1)}, {cmd:r(pi0)}, {cmd:r(se_pi0)},
{cmd:r(rhosig1)}, {cmd:r(se_rhosig1)}, {cmd:r(rhosig0)}, {cmd:r(se_rhosig0)},
{cmd:r(kappa_dd)}, {cmd:r(se_kappa_dd)}, {cmd:r(dd)}, {cmd:r(se_dd)},
{cmd:r(dlam1)} and {cmd:r(dlam0)} (the raw contrasts of the mean Mills ratios,
top - bottom, descriptive), {cmd:r(nq)} (the strata formed), {cmd:r(boundaries)}
(by boundary: the share below, its target, the shares of the tied values on
each side, and a flag when a resample can move them), {cmd:r(ties)},
{cmd:r(cv_pi1)}, {cmd:r(cv_pi0)}, {cmd:r(reps)} and {cmd:r(setype)}.


{title:Examples}

{phang2}{cmd:. esreg y x, select(d = x z1 z2) method(twostep)}{p_end}
{phang2}{cmd:. esrtest, spec}{p_end}
{phang2}{cmd:. esrtest, spec add(c.z1#c.z1 c.x#c.z1)}{p_end}
{phang2}{cmd:. esrtest, suff keep(z1)}{p_end}
{phang2}{cmd:. esrtest, normal}{p_end}
{phang2}{cmd:. esrtest, pdid nq(4)}{p_end}
{phang2}{cmd:. set seed 1}{p_end}
{phang2}{cmd:. esrtest, pdid nq(4) reps(500)}{p_end}
{phang2}{cmd:. esreg y x, select(d = x z1 z2) method(twostep) kappa(x)}{p_end}
{phang2}{cmd:. esrtest, kappa}{p_end}
