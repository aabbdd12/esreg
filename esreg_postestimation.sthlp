{smcl}
{* *! version 1.0.0  12sep2026}{...}
{vieweralsosee "esreg" "help esreg"}{...}
{vieweralsosee "esrdiag" "help esrdiag"}{...}
{vieweralsosee "esrtest" "help esrtest"}{...}
{vieweralsosee "esrcurve" "help esrcurve"}{...}
{vieweralsosee "esrmte" "help esrmte"}{...}
{vieweralsosee "esrreport" "help esrreport"}{...}
{title:Title}

{p2colset 5 30 32 2}{...}
{p2col:{bf:esreg postestimation} {hline 2}}Postestimation tools for esreg{p_end}
{p2colreset}{...}


{title:Postestimation commands}

{pstd}The following commands are available after {helpb esreg}:

{synoptset 16}{...}
{synopt:{helpb esrdiag}}strength, variation and support of the selection equation{p_end}
{synopt:{helpb esrtest}}specification tests: index, sufficiency, constant kappa, normality by regime, pseudo-DiD{p_end}
{synopt:{helpb esrcurve}}expected effect by quantile group of a ranking variable (the score){p_end}
{synopt:{helpb esrmte}}marginal treatment effect: parametric line and semiparametric curve{p_end}
{synopt:{helpb esrreport}}a reading of the results with rule-based notes{p_end}
{synopt:{cmd:predict}}predictions and scores (below){p_end}
{synopt:{helpb estimates}, {helpb test}, {helpb lincom}, {helpb nlcom}}standard tools on {cmd:e(b)}, {cmd:e(V)}{p_end}


{title:Syntax for predict}

{p 8 16 2}
{cmd:predict} [{it:type}] {newvar} {ifin} [{cmd:,} {it:statistic}]

{p 8 16 2}
{cmd:predict} [{it:type}] {it:stub}{cmd:*} {ifin}{cmd:,} {opt sc:ores}   (FIML only)

{synoptset 12}{...}
{synopthdr:statistic}
{synoptline}
{synopt:{opt eff:ect}}expected individual effect: E[Y_1 - Y_0 | X, Z, D] (the default){p_end}
{synopt:{opt xb1}, {opt xb0}}linear predictions X b_1, X b_0{p_end}
{synopt:{opt xbs:el}}selection index Z g{p_end}
{synopt:{opt pr}}P(D = 1 | Z) = Phi(Z g){p_end}
{synopt:{opt lam:bda1}, {opt lam:bda0}}phi(Zg)/Phi(Zg), phi(Zg)/(1 - Phi(Zg)){p_end}
{synopt:{opt yc11}, {opt yc10}}E[Y_1 | D = 1, X, Z], E[Y_1 | D = 0, X, Z]{p_end}
{synopt:{opt yc01}, {opt yc00}}E[Y_0 | D = 1, X, Z], E[Y_0 | D = 0, X, Z]{p_end}
{synopt:{opt kap:pa}}kappa_i = rho_1 sigma_1(x) - rho_0 sigma_0(x) (constant unless heterogeneous){p_end}
{synopt:{opt rhos:ig1}, {opt rhos:ig0}}rho_j sigma_j(x){p_end}
{synopt:{opt sig:ma1}, {opt sig:ma0}, {opt rho1}, {opt rho0}}sigma_j(x), rho_j(x) (FIML only){p_end}
{synoptline}

{pstd}
{cmd:scores} creates the seven equation-level scores of the FIML, in the order of
{cmd:e(b)}: y_1, y_0, the selection equation, lnsigma_1, lnsigma_0, atanhrho_1,
atanhrho_0. They are what the {cmd:svy} prefix uses for its linearized variance.

{pstd}
Correspondence with {cmd:mspredict} after {cmd:movestay}: xb1/xb2 = xb1/xb0,
psel = pr, mills1/mills2 = lambda1/lambda0, yc1_1 = yc11, yc2_1 = yc01,
yc1_2 = yc10, yc2_2 = yc00.


{title:Examples}

{phang2}{cmd:. esreg y x, select(d = x z)}{p_end}
{phang2}{cmd:. predict double eff, effect}{p_end}
{phang2}{cmd:. predict double p, pr}{p_end}
{phang2}{cmd:. esrcurve, rank(p) nq(10) graph}{p_end}
{phang2}{cmd:. predict double s*, scores}{p_end}
