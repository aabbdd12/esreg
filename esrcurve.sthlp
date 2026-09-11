{smcl}
{* *! version 1.0.0  12sep2026}{...}
{vieweralsosee "esreg" "help esreg"}{...}
{vieweralsosee "esrmte" "help esrmte"}{...}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col:{bf:esrcurve} {hline 2}}Expected treatment effect by quantile group after esreg{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:esrcurve} {ifin}{cmd:,} {opth rank(varname)} [{opt nq(#)} {opt est(name)} {opt gr:aph} {opt nose} {opt l:evel(#)} {opt ti:tle(string)} {opt na:me(string)} {opt sav:ing(filename)}]


{title:Description}

{pstd}
{cmd:esrcurve} averages the expected individual effect E[Y_1 - Y_0 | X, Z, D] of
a stored {helpb esreg} estimation within quantile groups of {cmd:rank()}
(typically the score, {cmd:predict}{cmd:, pr}), separately for the treated, the
untreated and all, with standard errors by the delta method on {cmd:e(V)} plus
the within-cell sampling component. Groups are weighted quantiles of {cmd:rank()}.
Along the score, the profile of the effect is the parametric marginal treatment
effect read through the groups; see {helpb esrmte} for the line and the
semiparametric curve.


{title:Options}

{phang}{cmd:rank(}{it:varname}{cmd:)} is the ranking variable (required).{p_end}
{phang}{cmd:nq(#)} is the number of groups (default 10).{p_end}
{phang}{cmd:est(}{it:name}{cmd:)} uses the estimation stored under {it:name}.{p_end}
{phang}{cmd:graph} draws the three profiles; {cmd:title()}, {cmd:name()} and {cmd:saving()} apply to the graph.{p_end}
{phang}{cmd:nose} skips the standard errors (faster).{p_end}


{title:Stored results}

{pstd}{cmd:r(table)}: one row per group with columns rank_mean, treated, se_treated,
untreated, se_untreated, all, se_all, n; {cmd:r(nq)}, {cmd:r(rank)}, {cmd:r(method)}.


{title:Example}

{phang2}{cmd:. esreg y x, select(d = x z) method(twostep)}{p_end}
{phang2}{cmd:. predict double p, pr}{p_end}
{phang2}{cmd:. esrcurve, rank(p) nq(5) graph}{p_end}
