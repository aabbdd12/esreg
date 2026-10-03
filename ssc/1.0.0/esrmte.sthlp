{smcl}
{* *! version 1.0.0  03oct2026}{...}
{vieweralsosee "esreg" "help esreg"}{...}
{vieweralsosee "esrcurve" "help esrcurve"}{...}
{title:Title}

{p2colset 5 15 17 2}{...}
{p2col:{bf:esrmte} {hline 2}}Marginal treatment effect after esreg: parametric line and semiparametric curve{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:esrmte} {ifin} [{cmd:,} {opt at(numlist)} {opt semi:par} {opt est(name)} {opt gr:aph} {opt l:evel(#)} {opt na:me(string)} {opt sav:ing(filename)}]


{title:Description}

{pstd}
{cmd:esrmte} evaluates the marginal treatment effect of a stored {helpb esreg}
estimation at percentiles u of the participation unobservable (u = 0 the most
eager): the parametric line MTE(u) = m + kappa invnormal(1 - u), with
m = mean(X)(b_1 - b_0), its standard error from the influence function of the
whole procedure (the parameters, the averages m and kappa over the units, and
their covariance; by cluster or by the survey design when the estimation was), and
a flag saying
whether u lies inside the common support of P(Z) (outside it the line is an
extrapolation).

{pstd}
After the augmented two-step ({cmd:esreg, method(twostep) hermite(}{it:#}{cmd:)}),
the parametric MTE is the curve MTE(u) = m + kappa v + dh_2 (v^2 - 1) +
dh_3 (v^3 - 3v), v = invnormal(1 - u), where dh_k = h_1k - h_0k are the differences
of the Hermite coefficients of the two regimes (a term absent from a regime,
{cmd:hermite(}{it:#1 #0}{cmd:)}, counting as zero), with its band from the covariance
of all its coefficients (the same influence function). Outside the common support
the polynomial extrapolates its Hermite terms, poorly beyond the data;
{cmd:semipar} is the check.

{pstd}
With {cmd:semipar}, the marginal treatment effect is also estimated without the
parametric line, as the derivative of E[Y | X, P] in P by percentile-weights
regression on the probit score at the quantiles of P (a weighted regression of Y on
X centred, X x P, P and P^2 at each quantile, bandwidth by a derivative-aware rule);
the curve is reported next to the line, with the effective sample and the width of
each window. The line is trusted where the curve follows it; the curve needs
variation in the score, and with a single indicator instrument it is uninformative.


{title:Options}

{phang}{cmd:at(}{it:numlist}{cmd:)} are the percentiles u (default .05 .1 .25 .5 .75 .9 .95).{p_end}
{phang}{cmd:semipar} adds the semiparametric curve.{p_end}
{phang}{cmd:est(}{it:name}{cmd:)} uses the estimation stored under {it:name}.{p_end}
{phang}{cmd:graph} draws the line with its band, the common support, the ATE, and the curve if requested.{p_end}


{title:Stored results}

{pstd}{cmd:r(mte)}: rows u with columns u, mte, se, support; with {cmd:semipar},
{cmd:r(mte_sp)}: columns tau, p_tau, parametric, pwr, se, h, neff, width_p;
{cmd:r(curve)}: the coefficients of the parametric MTE (m, kappa, and dh2, dh3
after {cmd:hermite}), {cmd:r(V_curve)} their covariance; {cmd:r(hermite)} the order
of the augmented two-step (0 for the line); {cmd:r(m)}, {cmd:r(kappa)},
{cmd:r(se_m)}, {cmd:r(se_kappa)}.


{title:Example}

{phang2}{cmd:. esreg y x, select(d = x z) method(twostep)}{p_end}
{phang2}{cmd:. esrmte, semipar graph}{p_end}
