{smcl}
{* *! version 1.0.0  12sep2026}{...}
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


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:esreg} {depvar} [{indepvars}] {ifin} [{it:{help esreg##weight:weight}}]{cmd:,}
{cmdab:sel:ect(}{it:treatvar} {cmd:=} {varlist}{cmd:)} [{it:options}]

{p 8 16 2}
{cmd:svy:} {cmd:esreg} {depvar} [{indepvars}]{cmd:,} {cmd:select(}{it:treatvar} {cmd:=} {varlist}{cmd:)} [{it:options}]   (FIML only)

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
{synopt:{opt her:mite}}two-step: add the Hermite controls E[u^2-1|D,Z], E[u^3-3u|D,Z] beside the Mills ratio (used by {helpb esrtest})){p_end}

{syntab:Weights}
{synopt:{opth hs:ize(varname)}}multiply the weight by the household size (effects per individual){p_end}
{synopt:{opt nosvy:set}}do not take the {cmd:svyset} weight when no weight is given{p_end}

{syntab:SE}
{synopt:{opth vce(vcetype)}}FIML: {opt oim} (default), {opt r:obust}, {opt cl:uster} {it:clustvar}, {opt opg}; the two-step always reports the stacked-moment variance{p_end}

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
{p 4 6 2}{cmd:pweight}s, {cmd:fweight}s and {cmd:iweight}s are allowed; see {help weight}. {cmd:svy}, {cmd:bootstrap} and {cmd:jackknife} are allowed as prefixes; see {help prefix} and {help esreg##svy:Weights and survey design} below.{p_end}


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
with standard errors by the delta method on {cmd:e(V)} plus a sampling component.
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
{cmd:vce(}{it:vcetype}{cmd:)} applies to FIML: {cmd:oim} (default), {cmd:robust}
(the sandwich, also used automatically with pweights), {cmd:cluster} {it:clustvar},
{cmd:opg}. For {cmd:bootstrap} and {cmd:jackknife} use the prefixes.

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
The FIML route supports the {cmd:svy} prefix ({cmd:svy linearized}, {cmd:svy
bootstrap}, {cmd:svy jackknife}, {cmd:svy brr}): {cmd:esreg} accepts the iweights
the prefix passes and {cmd:predict}{cmd:, scores} returns the seven equation-level
scores of the likelihood. After {cmd:svy: esreg}, {cmd:svy} has replaced
{cmd:e(V)} by the design-based variance; the effects and kappa are recomputed from
it the first time the results are replayed (type {cmd:esreg}), or with
{cmd:esreg, effects}; {cmd:e(eff_vce)} says which variance they rest on. The LR
test of independent equations is not available under {cmd:svy}.

{pstd}
The two-step route is not {cmd:svy}-able (as {cmd:etregress, twostep} is not): a
sequential estimator has no single likelihood score to linearize. Use pweights
(the stacked variance is design-consistent for weights) or the {cmd:bootstrap}
prefix with {cmd:cluster()} for the primary sampling units.


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
{synopt:{cmd:e(k_x)}, {cmd:e(k_z)}, {cmd:e(k_hs)}, {cmd:e(k_hr)}, {cmd:e(k_kap)}, {cmd:e(k_h)}}columns of the design matrices{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}, {cmd:e(cmdline)}, {cmd:e(version)}, {cmd:e(method)}}{cmd:esreg}; {cmd:fiml} or {cmd:twostep}{p_end}
{synopt:{cmd:e(depvar)}, {cmd:e(treat)}}outcome and treatment{p_end}
{synopt:{cmd:e(xvars)}, {cmd:e(zvars)}, {cmd:e(hetsigma)}, {cmd:e(hetrho)}, {cmd:e(kappavars)}}expanded variable lists{p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}, {cmd:e(hsize)}, {cmd:e(vce)}, {cmd:e(vcetype)}, {cmd:e(eff_vce)}}weights and variances{p_end}
{synopt:{cmd:e(eqnames)}, {cmd:e(predict)}, {cmd:e(store)}, {cmd:e(title)}}{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients and variance (FIML: equations y_1, y_0, treat, lnsigma_1, lnsigma_0, atanhrho_1, atanhrho_0){p_end}
{synopt:{cmd:e(effects)}}4 x 4: ATT, ATU, ATE, kappa by est, se, var_param, var_samp{p_end}
{synopt:{cmd:e(b_sel)}, {cmd:e(b_1)}, {cmd:e(b_0)}}selection and regime coefficient blocks{p_end}
{synopt:{cmd:e(anc)}, {cmd:e(support)}, {cmd:e(lambda)}}(sigma1 sigma0 rho1 rho0), (p_min1 p_max1 p_min0 p_max0), (mean lambda1, mean lambda0){p_end}

{p2col 5 22 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}{p_end}


{marker examples}{...}
{title:Examples}

{pstd}The manual example of {cmd:etregress}: union membership and wages{p_end}
{phang2}{cmd:. webuse union3}{p_end}
{phang2}{cmd:. esreg ln_wage age grade smsa black tenure, select(union = south black tenure)}{p_end}
{phang2}{cmd:. esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep)}{p_end}
{phang2}{cmd:. esrdiag}{p_end}
{phang2}{cmd:. esrtest, normal}{p_end}
{phang2}{cmd:. esrreport}{p_end}

{pstd}Heterogeneous selection on gains{p_end}
{phang2}{cmd:. esreg y x, select(d = x z) method(twostep) kappa(x)}{p_end}
{phang2}{cmd:. esrtest, kappa}{p_end}
{phang2}{cmd:. esreg y x, select(d = x z) hetsigma(x) hetrho(x)}{p_end}

{pstd}Survey design{p_end}
{phang2}{cmd:. svyset psu [pw = wt], strata(strata)}{p_end}
{phang2}{cmd:. svy: esreg y x, select(d = x z)}{p_end}
{phang2}{cmd:. esreg}{p_end}
{phang2}{cmd:. esreg y x, select(d = x z) method(twostep) hsize(hhsize)}{p_end}


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

{phang}Maddala, G. S. (1983). {it:Limited-Dependent and Qualitative Variables in
Econometrics}. Cambridge University Press.{p_end}


{title:Author}

{pstd}Abdelkrim Araar, Universite Laval and Partnership for Economic Policy (PEP).
Bug reports and suggestions: see the package page on GitHub.{p_end}
