{smcl}
{* *! version 1.0.0  03oct2026}{...}
{vieweralsosee "esreg postestimation" "help esreg_postestimation"}{...}
{vieweralsosee "esrdiag" "help esrdiag"}{...}
{vieweralsosee "esrtest" "help esrtest"}{...}
{vieweralsosee "esrcurve" "help esrcurve"}{...}
{vieweralsosee "esrmte" "help esrmte"}{...}
{vieweralsosee "esrreport" "help esrreport"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[CAUSAL] etregress" "help etregress"}{...}
{vieweralsosee "[SVY] svy estimation" "help svy estimation"}{...}
{viewerjumpto "Syntax" "esreg##syntax"}{...}
{viewerjumpto "Description" "esreg##description"}{...}
{viewerjumpto "Options" "esreg##options"}{...}
{viewerjumpto "Weights and survey design" "esreg##svy"}{...}
{viewerjumpto "Stored results" "esreg##results"}{...}
{viewerjumpto "Examples" "esreg##examples"}{...}
{viewerjumpto "References" "esreg##references"}{...}
{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{bf:esreg} {hline 2}}Endogenous switching regression with heterogeneous selection on gains{p_end}
{p2colreset}{...}
{p 4 4 2}{txt}Package {cmd:esreg}, version {res}1.0.0{txt} (03/10/2026) {c |} Stata {res}16{txt} or later {c |} first release {res}0.4.1{txt} (12/09/2026){p_end}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
[{cmd:by} {varlist}{cmd::}] {cmd:esreg} {depvar} [{indepvars}] {ifin} [{it:{help esreg##weight:weight}}]{cmd:,}
{cmdab:sel:ect(}{it:treatvar} {cmd:=} {varlist}{cmd:)} [{it:options}]

{p 8 16 2}
{cmd:svy:} {cmd:esreg} {depvar} [{indepvars}]{cmd:,} {cmd:select(}{it:treatvar} {cmd:=} {varlist}{cmd:)} [{it:options}]   (FIML; the two-step has {cmd:vce(svy)})

{p 8 16 2}
{cmd:esreg} [{cmd:,} {opt eff:ects} {opt l:evel(#)}]   (replay; {cmd:effects} recomputes the effects from the current {cmd:e(V)})

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {cmdab:sel:ect(}{it:treatvar} {cmd:=} {varlist}{cmd:)}}selection equation: the 0/1 treatment and its regressors{p_end}
{synopt:{cmdab:meth:od(fiml)}}full-information maximum likelihood (the default){p_end}
{synopt:{cmdab:meth:od(twostep)}}probit, then OLS by regime with the Mills ratio; exact stacked-moment variance{p_end}
{synopt:{cmdab:hets:igma(}{varlist}{cmd:)}}FIML: ln sigma_j linear in {it:varlist} (heterogeneous regime variances){p_end}
{synopt:{cmdab:hetr:ho(}{varlist}{cmd:)}}FIML: atanh rho_j linear in {it:varlist} (heterogeneous correlations){p_end}
{synopt:{cmdab:kap:pa(}{varlist}{cmd:)}}two-step: rho_j sigma_j(x) = W theta_j, hence kappa(x) = W(theta_1 - theta_0){p_end}
{synopt:{cmdab:her:mite}[{cmd:(}{it:#} [{it:#}]{cmd:)}]}two-step, the augmented model: Hermite terms of order 2 or 3 in E[w_j | u],
in both regimes or by regime (treated, untreated; 0 = none); {cmd:hermite} alone is {cmd:hermite(3)}{p_end}

{syntab:Weights}
{synopt:{opth hs:ize(varname)}}multiply the weight by the household size (effects per individual){p_end}
{synopt:{opt nosvy:set}}do not take the {cmd:svyset} weight when no weight is given{p_end}

{syntab:SE}
{synopt:{opth vce(vcetype)}}FIML: {opt oim} (default), {opt r:obust}, {opt cl:uster} {it:clustvar}, {opt opg};
two-step: {opt r:obust} (default, the stacked moments), {opt cl:uster} {it:clustvar}, {opt svy}{p_end}

{syntab:Reporting}
{synopt:{opt noeff:ects}}do not compute the effects and kappa{p_end}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt nolog}}suppress the iteration log{p_end}
{synopt:{opt sto:re(name)}}also store the estimation under {it:name} ({helpb estimates store}){p_end}

{syntab:Maximization}
{synopt:{opt iter:ate(#)}, {opt dif:ficult}}passed to {helpb ml}{p_end}
{synoptline}
{p 4 6 2}* {cmd:select()} is required.{p_end}
{p 4 6 2}{it:indepvars} and the variables of {cmd:select()}, {cmd:hetsigma()}, {cmd:hetrho()} and {cmd:kappa()} may contain factor variables; see {help fvvarlist}.{p_end}
{marker weight}{...}
{p 4 6 2}{cmd:pweight}s and {cmd:iweight}s are allowed; see {help weight}. {cmd:fweight}s are not: a unit of a survey
stands for its sampling weight, it is not a replicated record (replicated records can be {helpb expand}ed first).
{cmd:svy} (FIML), {cmd:bootstrap} and {cmd:jackknife} are allowed as prefixes; see {help prefix} and
{help esreg##svy:Weights and survey design} below.{p_end}
{p 4 6 2}{cmd:by} is allowed; see {help by}.{p_end}
{p 4 6 2}A dialog box is available: type {cmd:db esreg}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:esreg} fits the endogenous switching regression (Roy model)

{p 12 12 2}Y_1 = X b_1 + w_1,   Y_0 = X b_0 + w_0,   D = 1{Z g + u > 0},   Y = D Y_1 + (1 - D) Y_0,{p_end}

{pstd}
with (w_1, w_0, u) jointly normal (FIML) or with only the linear conditional means
E[w_j | u] = rho_j sigma_j u (two-step), and reports, next to the coefficients, the
treatment effects ATT, ATU and ATE and the parameter of selection on gains

{p 12 12 2}kappa = Cov(w_1 - w_0, u) = rho_1 sigma_1 - rho_0 sigma_0,{p_end}

{pstd}
with standard errors from the influence function of the whole procedure: the
estimation of the parameters (delta method on {cmd:e(V)}), the averaging over the
units, and the covariance of the two, aggregated as the estimation was (by
observation, by cluster, or by the survey design).
ATT - ATU = kappa (mean lambda_1 + mean lambda_0): kappa > 0 means that those who
select in gain more than those who stay out. The command also reports the common
support of the score P(Z) = Phi(Z g).

{pstd}
With the option {cmd:poutcomes} and the treatment interacted with every covariate,
Stata's {helpb etregress} fits the same likelihood; {cmd:esreg} adds the two-step
route with its exact variance, kappa and kappa(x), heterogeneous regime laws, and
the family of post-estimation commands: {helpb esrdiag} (strength, variation and
support of the selection equation), {helpb esrtest} (specification tests: index,
sufficiency of the instruments, constant kappa, normality by regime, pseudo-DiD),
{helpb esrcurve} (effect by quantile of the score), {helpb esrmte} (marginal
treatment effect, parametric line and semiparametric curve), {helpb esrreport}
(a reading of the results for the practitioner) and {helpb esreg postestimation}
({cmd:predict}).

{pstd}
The estimation is stored automatically as {cmd:_esreg} (the last {cmd:esreg} of
the session); the post-estimation commands look for {cmd:est(}{it:name}{cmd:)},
then the current {cmd:e()}, then {cmd:_esreg}, and leave the user's {cmd:e()}
untouched.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmd:select(}{it:treatvar} {cmd:=} {varlist}{cmd:)} names the 0/1 treatment variable
and the regressors of the selection equation. The variables of {it:varlist} that are
not in {it:indepvars} are the excluded instruments; at least one is needed for
identification beyond functional form ({helpb esrdiag} says how much rests on the
form).

{phang}
{cmd:method(fiml)}, the default, maximizes the log-likelihood of (Y, D) given Z with
{helpb ml} (method lf1, analytic score), in the parameterization ln sigma_j and
atanh rho_j. The LR test of independent equations (rho_1 = rho_0 = 0) is reported.

{phang}
{cmd:method(twostep)} estimates the probit of D on Z, then regresses Y on X and
lambda_1 = phi(Zg)/Phi(Zg) among the treated and on X and -lambda_0 =
-phi(Zg)/(1 - Phi(Zg)) among the untreated. The coefficients of the Mills terms are
rho_1 sigma_1 and rho_0 sigma_0. The variance is the sandwich of the stacked moment
conditions of the whole procedure (probit score, two sets of normal equations),
which carries the first-step correction and the heteroskedasticity induced by the
truncation; it replaces the bootstrap.

{phang}
{cmd:hetsigma(}{varlist}{cmd:)} and {cmd:hetrho(}{varlist}{cmd:)} (FIML) make
ln sigma_j and atanh rho_j linear in the listed variables, so that
rho_j sigma_j and kappa vary across units; {helpb esrtest}{cmd:, kappa} tests the
homogeneity.

{phang}
{cmd:kappa(}{varlist}{cmd:)} (two-step) interacts the Mills ratio of each regime
with the listed variables W, so that rho_j sigma_j(x) = W theta_j and
kappa(x) = W(theta_1 - theta_0); {helpb esrtest}{cmd:, kappa} tests that the
non-constant entries of theta_1 - theta_0 are zero.

{phang}
{cmd:hermite}[{cmd:(}{it:#} [{it:#}]{cmd:)}] (two-step) fits the augmented model, in which
the conditional means of the regime errors are no longer linear in u:
E[w_j | u] = rho_j sigma_j u + h_j2 (u^2 - 1) + h_j3 (u^3 - 3u), with the second
and third Hermite polynomials. Each regime regression adds the truncated moments
E[u^2 - 1 | D, Z] and, at order 3, E[u^3 - 3u | D, Z] beside the Mills ratio;
{it:#} is 2 (the quadratic term only) or 3 (both; {cmd:hermite} alone).
{cmd:hermite(}{it:#1 #0}{cmd:)} sets the order by regime, {it:#1} for the treated
and {it:#0} for the untreated, 0 for no Hermite terms in that regime: for example
{cmd:hermite(3 0)} when only the treated regime departs from the line
({helpb esrtest}{cmd:, normal} tests each regime). ATT depends on the untreated
regime's terms only and ATU on the treated regime's (the treated mean of the
fitted Y_1 is the mean outcome of the treated, whatever the regressors), so that
terms added where a regime is linear cost precision for nothing:
{cmd:hermite(3 0)} has the ATT of the line, and the ATU of {cmd:hermite(3)}. ATT and
ATU add the Hermite terms of the truncated moments, with dh_k = h_1k - h_0k, and
their standard errors come from the same influence function of the whole
procedure; ATE = mean X(b_1 - b_0), as the Hermite terms have mean zero; kappa =
Cov(w_1 - w_0, u) is the linear part. {helpb esrmte} evaluates the curve
MTE(u) = m + kappa v + dh_2 (v^2 - 1) + dh_3 (v^3 - 3v), v = invnormal(1 - u),
with its band. The order is guided by {helpb esrtest}{cmd:, normal}: 3 if its
cubic terms reject, else 2; {helpb esrreport} proposes this route when the Hermite
controls reject with strong instruments. Choosing the augmented model after these
tests is a choice to state when reporting, and the polynomial extrapolates poorly
beyond the support of the score: read the effects and the curve on the common
support, with {cmd:esrmte, semipar} as the check.

{dlgtab:Weights}

{phang}
{cmd:hsize(}{it:varname}{cmd:)} multiplies the weight (the one given, or the
{cmd:svyset} one, or 1) by the household size, so that the effects are averages
over individuals when the observations are households. Not allowed under the
{cmd:svy} prefix.

{phang}
{cmd:nosvyset} declines the {cmd:svyset} pweight, which is otherwise taken by
default (with a note) when no weight is given in the command. An explicit weight is
always used alone: it is never multiplied by the {cmd:svyset} one.

{dlgtab:SE}

{phang}
{cmd:vce(}{it:vcetype}{cmd:)} with {cmd:method(fiml)}: {cmd:oim} (default),
{cmd:robust} (the sandwich, also used automatically with pweights),
{cmd:cluster} {it:clustvar}, {cmd:opg}; for the survey design use the {cmd:svy}
prefix. With {cmd:method(twostep)}: {cmd:robust} (default; the sandwich of the
stacked moments of the whole procedure), {cmd:cluster} {it:clustvar} (the
influence functions summed within the clusters, with the factor G/(G - 1) of
{helpb ml}), {cmd:svy} (the design-based variance; see
{help esreg##svy:Weights and survey design}). The standard errors of the effects
follow the same aggregation. For {cmd:bootstrap} and {cmd:jackknife} use the
prefixes.

{dlgtab:Reporting}

{phang}
{cmd:noeffects} skips the effects, kappa, the support and the ancillary matrices
(faster; the post-estimation commands still work).

{phang}
{cmd:store(}{it:name}{cmd:)} stores the estimation under {it:name} in addition to
{cmd:_esreg}; the post-estimation commands accept {cmd:est(}{it:name}{cmd:)}.

{phang}
{cmd:effects} (replay only) recomputes the effects, kappa and their standard errors
from the current {cmd:e(b)} and {cmd:e(V)}; after {cmd:svy: esreg} this is done
automatically at the first replay, from the linearized variance.


{marker svy}{...}
{title:Weights and survey design}

{pstd}
Estimation with pweights weights the probit, the regime regressions and the
likelihood, uses the robust (FIML) or stacked (two-step) variance, and averages
the effects with the weights. When the data are {cmd:svyset} and no weight is
given, the design weight is used as a pweight and a note is printed.

{pstd}
The FIML route supports the {cmd:svy} prefix ({cmd:svy linearized},
{cmd:svy bootstrap}, {cmd:svy jackknife}, {cmd:svy brr}): {cmd:esreg} accepts the iweights
the prefix passes and {cmd:predict}{cmd:, scores} returns the seven equation-level
scores of the likelihood. After {cmd:svy: esreg}, {cmd:svy} has replaced
{cmd:e(V)} by the design-based variance; the effects and kappa are recomputed from
it the first time the results are replayed (type {cmd:esreg}), or with
{cmd:esreg, effects}; {cmd:e(eff_vce)} says which variance they rest on. The LR
test of independent equations is not available under {cmd:svy}.

{pstd}
The two-step route has no likelihood score for the {cmd:svy} prefix to linearize
(as {cmd:etregress, twostep}); its option {cmd:vce(svy)} gives the design-based
variance of the whole procedure directly. The influence functions of the stacked
moments (probit score, the two sets of normal equations) are totaled through the
design of {cmd:svyset} (strata, primary sampling units, finite-population
correction; {cmd:svy linearized: total}), with the {cmd:svyset} pweight as the
weight; the coefficient table then uses the design degrees of freedom. With PSUs
and no strata, {cmd:vce(svy)} equals {cmd:vce(cluster} {it:psu}{cmd:)}.
Poststratification and calibration are not supported.

{pstd}
The standard errors of the effects and kappa add the variance of the parameters
(delta method on {cmd:e(V)}), the variance of the average over the units, and
their covariance, all three from the influence function of the whole procedure;
under {cmd:svy} or {cmd:vce(svy)} the last two are design-based (linearized) and
the effects table uses t with the design degrees of freedom. After
{cmd:svy, subpop():} {cmd:esreg}, the effects are those of the subpopulation and
their variance is that of a domain.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:esreg} stores the following in {cmd:e()} (see also {cmd:esreg_returns.txt}
in the package):

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}, {cmd:e(N_treated)}, {cmd:e(N_untreated)}, {cmd:e(sum_w)}}observations, by group, and sum of weights{p_end}
{synopt:{cmd:e(ll)}, {cmd:e(ll_indep)}, {cmd:e(lr_indep)}, {cmd:e(lr_df)}, {cmd:e(p_indep)}}log-likelihoods and LR test of rho_1 = rho_0 = 0 (FIML){p_end}
{synopt:{cmd:e(att)}, {cmd:e(atu)}, {cmd:e(ate)}, {cmd:e(kappa)}}effects and kappa{p_end}
{synopt:{cmd:e(se_att)}, {cmd:e(se_atu)}, {cmd:e(se_ate)}, {cmd:e(se_kappa)}}their standard errors{p_end}
{synopt:{cmd:e(sigma1)}, {cmd:e(sigma0)}, {cmd:e(rho1)}, {cmd:e(rho0)}, {cmd:e(rhosig1)}, {cmd:e(rhosig0)}}ancillary parameters (means when heterogeneous){p_end}
{synopt:{cmd:e(supp_lo)}, {cmd:e(supp_hi)}, {cmd:e(p_min1)}, {cmd:e(p_max1)}, {cmd:e(p_min0)}, {cmd:e(p_max0)}}common support and ranges of P(Z){p_end}
{synopt:{cmd:e(ml1)}, {cmd:e(ml0)}}mean lambda_1 (treated), mean lambda_0 (untreated){p_end}
{synopt:{cmd:e(k_x)}, {cmd:e(k_z)}, {cmd:e(k_hs)}, {cmd:e(k_hr)}, {cmd:e(k_kap)}, {cmd:e(k_h)}}columns of the design matrices ({cmd:e(k_h)}: Hermite controls per regime){p_end}
{synopt:{cmd:e(hermite)}, {cmd:e(hermite1)}, {cmd:e(hermite0)}}order of the augmented two-step, 2 or 3 (0: the linear model): the larger one, and by regime{p_end}
{synopt:{cmd:e(k_h1)}, {cmd:e(k_h0)}}Hermite controls of the treated and of the untreated regime{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters ({cmd:vce(cluster)}){p_end}
{synopt:{cmd:e(N_strata)}, {cmd:e(N_psu)}, {cmd:e(df_r)}}strata, PSUs and design degrees of freedom (two-step, {cmd:vce(svy)}){p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}, {cmd:e(cmdline)}, {cmd:e(version)}, {cmd:e(method)}}{cmd:esreg}; {cmd:fiml} or {cmd:twostep}{p_end}
{synopt:{cmd:e(depvar)}, {cmd:e(treat)}}outcome and treatment{p_end}
{synopt:{cmd:e(xvars)}, {cmd:e(zvars)}, {cmd:e(hetsigma)}, {cmd:e(hetrho)}, {cmd:e(kappavars)}}expanded variable lists{p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}, {cmd:e(hsize)}, {cmd:e(vce)}, {cmd:e(vcetype)}, {cmd:e(clustvar)}}weights and variance; two-step {cmd:e(vce)}: {cmd:stacked}, {cmd:cluster} or {cmd:svy}{p_end}
{synopt:{cmd:e(eff_vce)}}aggregation of the standard errors of the effects: {cmd:iid}, {cmd:cluster}, {cmd:svy} ({cmd:none} with {cmd:noeffects}){p_end}
{synopt:{cmd:e(eqnames)}, {cmd:e(predict)}, {cmd:e(store)}, {cmd:e(title)}}{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients and variance (FIML: equations y_1, y_0, treat, lnsigma_1, lnsigma_0, atanhrho_1, atanhrho_0){p_end}
{synopt:{cmd:e(effects)}}4 x 5: ATT, ATU, ATE, kappa by est, se, var_param, var_samp, cov_ps (twice the covariance of the parameter and sampling parts); se^2 = var_param + var_samp + cov_ps{p_end}
{synopt:{cmd:e(b_sel)}, {cmd:e(b_1)}, {cmd:e(b_0)}}selection and regime coefficient blocks{p_end}
{synopt:{cmd:e(anc)}, {cmd:e(support)}, {cmd:e(lambda)}}(sigma1 sigma0 rho1 rho0), (p_min1 p_max1 p_min0 p_max0), (mean lambda1, mean lambda0){p_end}

{p2col 5 22 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
Every example runs from its links: in the command window, in the dialog box
filled in, or as a do-file in the Do-file Editor; the data in memory are given
back at the end. Examples 1 to 3 use the manual's {cmd:union3} data, read with
{cmd:webuse} (an internet connection is needed); examples 4 to 6 use simulated
samples shipped with the package ({stata "ssc install esreg, all replace"} or
{cmd:net get esreg} copies them into the current folder; otherwise they are
read from the SSC archive or from GitHub). The syntax is
{cmd:esreg_examples} {it:#} [{cmd:, db} | {cmd:do}].

{pstd}{bf:Example 1.} The manual's example of {cmd:etregress}, union membership
and wages: the switching regression by full-information maximum likelihood{p_end}
{phang2}{cmd:. webuse union3}{p_end}
{phang2}{cmd:. esreg ln_wage age grade smsa black tenure, select(union = south black tenure)}{p_end}
{p 8 8 2}{txt}({stata "esreg_examples 1":example 1: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 1, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 1, do":open as a do-file}){p_end}

{pstd}{bf:Example 2.} The two-step, the diagnostics of the selection equation,
the tests and the reading{p_end}
{phang2}{cmd:. webuse union3}{p_end}
{phang2}{cmd:. esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep)}{p_end}
{phang2}{cmd:. esrdiag}{p_end}
{phang2}{cmd:. esrtest, normal}{p_end}
{phang2}{cmd:. esrtest, pdid}{p_end}
{phang2}{cmd:. esrreport}{p_end}
{p 8 8 2}{txt}({stata "esreg_examples 2":example 2: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 2, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 2, do":open as a do-file}){p_end}

{pstd}On {cmd:union3} the two routes differ widely (ATT 0.39 by maximum likelihood,
1.10 by the two-step, log wage): {cmd:esrtest, normal} and {cmd:esrreport} say which
one the tests support, the two-step, and both rest on {cmd:south}, the only excluded
instrument, whose exclusion no test can check here (Section 7 of the technical note, References).{p_end}

{pstd}{bf:Example 3.} The effect along the score and the marginal treatment
effect{p_end}
{phang2}{cmd:. webuse union3}{p_end}
{phang2}{cmd:. esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep)}{p_end}
{phang2}{cmd:. predict double P, pr}{p_end}
{phang2}{cmd:. esrcurve, rank(P) nq(5) graph}{p_end}
{phang2}{cmd:. esrmte, semipar graph}{p_end}
{p 8 8 2}{txt}({stata "esreg_examples 3":example 3: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 3, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 3, do":open as a do-file}){p_end}

{pstd}{bf:Example 4.} Selection on gains that varies with x, kappa(x), and its
test (simulated sample {cmd:esr_kx}, kappa(x) = 1.08 + 0.4x){p_end}
{phang2}{cmd:. use esr_kx}{p_end}
{phang2}{cmd:. esreg y x, select(d = x z) method(twostep) kappa(x)}{p_end}
{phang2}{cmd:. esrtest, kappa}{p_end}
{p 8 8 2}{txt}({stata "esreg_examples 4":example 4: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 4, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 4, do":open as a do-file}){p_end}

{pstd}{bf:Example 5.} A conditional mean not linear in u: the Hermite test,
then the augmented two-step with the Hermite terms in the treated regime
(simulated sample {cmd:esr_nonlin}, whose gain is cubic in u){p_end}
{phang2}{cmd:. use esr_nonlin}{p_end}
{phang2}{cmd:. esreg y x, select(d = x z) method(twostep)}{p_end}
{phang2}{cmd:. esrtest, normal}{p_end}
{phang2}{cmd:. esreg y x, select(d = x z) method(twostep) hermite(3 0)}{p_end}
{phang2}{cmd:. esrmte, semipar graph}{p_end}
{p 8 8 2}{txt}({stata "esreg_examples 5":example 5: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 5, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 5, do":open as a do-file}){p_end}

{pstd}{bf:Example 6.} A survey design: the two-step with {cmd:vce(svy)}, the
likelihood with the {cmd:svy} prefix (simulated sample {cmd:esr_wt}, with
sampling weights){p_end}
{phang2}{cmd:. use esr_wt}{p_end}
{phang2}{cmd:. svyset [pweight = wt], strata(region)}{p_end}
{phang2}{cmd:. esreg income educ, select(treatment = educ i.region inst) method(twostep) vce(svy)}{p_end}
{phang2}{cmd:. svy: esreg income educ, select(treatment = educ i.region inst)}{p_end}
{phang2}{cmd:. esreg}{p_end}
{p 8 8 2}{txt}({stata "esreg_examples 6":example 6: click to run in command window}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 6, db":click to run in dialog box}){p_end}
{p 8 8 2}{txt}({stata "esreg_examples 6, do":open as a do-file}){p_end}


{marker references}{...}
{title:References}

{phang}Araar, A. (2026). Endogenous switching regression with heterogeneous
selection on gains: assumptions, tests, and estimation. Zenodo,
doi:10.5281/zenodo.22717029.{p_end}

{phang}Araar, A. (2026). Estimating treatment effects under selection on gains:
models, assumptions, and policy implications. Zenodo, doi:10.5281/zenodo.22672713.{p_end}

{phang}Heckman, J. J., and E. Vytlacil (2005). Structural equations, treatment
effects, and econometric policy evaluation. {it:Econometrica} 73(3), 669-738.{p_end}

{phang}Lokshin, M., and Z. Sajaia (2004). Maximum likelihood estimation of
endogenous switching regression models. {it:Stata Journal} 4(3), 282-289.{p_end}

{phang}Maddala, G. S. (1983).
{it:Limited-Dependent and Qualitative Variables in Econometrics}.
Cambridge University Press.{p_end}


{title:Author}

{pstd}Abdelkrim Araar, Universite Laval and Partnership for Economic Policy (PEP).
Bug reports and suggestions: see the package page on GitHub.{p_end}
{pstd}Version 1.0.0. Requires Stata 16 or later. License: MIT.{p_end}
