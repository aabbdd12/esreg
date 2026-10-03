{smcl}
{* *! version 1.0.0  12sep2026}{...}
{vieweralsosee "esreg" "help esreg"}{...}
{vieweralsosee "esrtest" "help esrtest"}{...}
{vieweralsosee "esrreport" "help esrreport"}{...}
{title:Title}

{p2colset 5 16 18 2}{...}
{p2col:{bf:esrdiag} {hline 2}}Diagnostics of the selection equation after esreg{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:esrdiag} [{cmd:,} {opt est(name)}]


{title:Description}

{pstd}
{cmd:esrdiag} reports what the selection equation of an {helpb esreg} estimation
can bear, in four blocks.

{pstd}
{it:Fit}: the log likelihood and the pseudo-R2 of the selection equation on its
own, that is of the probit of the treatment on the whole of Z, together with the
chi2 against the intercept-only model. This block is always reported, whether or
not the selection equation holds a variable excluded from the outcome equations.

{pstd}
{it:Strength}: the LR test of the probit with and without the excluded
instruments, the increment of the pseudo-R2, and a Wald chi2(1) per instrument.

{pstd}
{it:Variation}: the share of the variance of the score P(Z) not explained by X,
Var(P | X)/Var(P), and the variance-inflation factors of the Mills ratios
lambda_1 (among the treated) and lambda_0 (among the untreated) on X. A large VIF
says that the Mills ratio is nearly a linear function of X, so that
identification rests on the curvature of lambda alone (identification by form).

{pstd}
{it:Support}: the range of P(Z) in each group, the common support, the weighted
share of the treated (untreated) whose counterfactual is extrapolated outside the
other group's range, and the ATT and ATU restricted to the common support next to
the full ones.

{pstd}
Warnings are printed when the LR statistic is below 10, when a VIF exceeds 10, and
when more than 10 percent of a group is extrapolated. The score and the Mills
ratios are those of the stored estimation (the probit gamma after
{cmd:method(twostep)}, the FIML gamma after {cmd:method(fiml)}).


{title:Options}

{phang}
{cmd:est(}{it:name}{cmd:)} uses the estimation stored under {it:name}; the default
is the current {cmd:e()} if it is an {cmd:esreg} estimation, else {cmd:_esreg}.


{title:Stored results}

{pstd}{cmd:esrdiag} stores in {cmd:r()}: {cmd:r(lr_excl)}, {cmd:r(df_excl)},
{cmd:r(p_excl)}, {cmd:r(r2_full)}, {cmd:r(r2_restr)}, {cmd:r(r2_incr)},
{cmd:r(ll_sel)}, {cmd:r(ll0_sel)}, {cmd:r(df_sel)},
{cmd:r(varP_share)}, {cmd:r(vif1)}, {cmd:r(vif0)}, {cmd:r(supp_lo)},
{cmd:r(supp_hi)}, {cmd:r(p_min1)}, {cmd:r(p_max1)}, {cmd:r(p_min0)},
{cmd:r(p_max0)}, {cmd:r(share_att)}, {cmd:r(share_atu)}, {cmd:r(att)},
{cmd:r(atu)}, {cmd:r(att_cs)}, {cmd:r(atu_cs)}, the flags {cmd:r(weak)} and
{cmd:r(form)}, the matrix {cmd:r(instr)} (coef, se, chi2 per instrument) and the
macros {cmd:r(excluded)}, {cmd:r(method)}.


{title:Example}

{phang2}{cmd:. webuse union3}{p_end}
{phang2}{cmd:. esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep)}{p_end}
{phang2}{cmd:. esrdiag}{p_end}
