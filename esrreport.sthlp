{smcl}
{* *! version 1.0.0  12sep2026}{...}
{vieweralsosee "esreg" "help esreg"}{...}
{vieweralsosee "esrdiag" "help esrdiag"}{...}
{vieweralsosee "esrtest" "help esrtest"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{bf:esrreport} {hline 2}}A reading of an esreg estimation for the practitioner{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:esrreport} [{cmd:,} {opt est(name)} {opt nq(#)} {opt lr:min(#)} {opt vif:max(#)} {opt extrap:max(#)} {opt alpha(#)} {opt exp:ort(filename)} {opt replace} {opt notests}]


{title:Description}

{pstd}
{cmd:esrreport} runs the diagnostics and the tests of the family on a stored
{helpb esreg} estimation and prints an ordered reading in four steps: (1) the
selection equation (strength, form, support, link test, sufficiency); (2)
selection on gains (kappa and its sign, the pseudo-DiD pattern, kappa(x)); (3) the
joint law and the route to report (Hausman contrast, gamma contrast, Hermite
controls); (4) the effects (full and on the common support, profile by quantile
of the score). Numbered notes are printed when a threshold is crossed; each note
names the test and the section of the paper it rests on. The route is decided
in the order of the assumptions: a rejected link test comes first ("index":
respecify the selection equation), then a rejected sufficiency ("instruments"),
then rejected Hermite controls ("semipar"), then a rejected Hausman contrast
("twostep"), else "fiml".

{pstd}
The notes are statements and questions, not verdicts: the thresholds are
conventions set by the options. {cmd:export()} writes the figures and the notes to
a plain-text sheet, which an assistant, human or not, can read together with the
context of the study.


{title:Options}

{phang}{cmd:est(}{it:name}{cmd:)} uses the estimation stored under {it:name}.{p_end}
{phang}{cmd:nq(#)} groups of the score for the profile of effects (default 5).{p_end}
{phang}{cmd:lrmin(#)} LR of the excluded instruments below which they are called weak (default 10).{p_end}
{phang}{cmd:vifmax(#)} VIF of the Mills ratio above which identification is called "by form" (default 10).{p_end}
{phang}{cmd:extrapmax(#)} share of a group outside the other's support above which the effects on the common support are recommended (default 0.10).{p_end}
{phang}{cmd:alpha(#)} size of the tests in the reading (default 0.05).{p_end}
{phang}{cmd:export(}{it:filename}{cmd:)} writes the reading sheet; {cmd:replace} overwrites.{p_end}
{phang}{cmd:notests} prints the diagnostics and the effects only (no re-estimation).{p_end}


{title:Stored results}

{pstd}{cmd:r(route)} ({cmd:fiml}, {cmd:twostep}, {cmd:semipar}, {cmd:index} or
{cmd:instruments}), {cmd:r(n_notes)}, {cmd:r(note1)}, {cmd:r(note2)}, ..., and the
figures used: {cmd:r(kappa)}, {cmd:r(se_kappa)}, {cmd:r(kappa_dd)}, {cmd:r(p_link)},
{cmd:r(p_gamma)}, {cmd:r(p_suff)}, {cmd:r(p_hausman)}, {cmd:r(p_hermite)},
{cmd:r(lr_excl)}, {cmd:r(vif)}, {cmd:r(share_extrap)}, {cmd:r(att)}, {cmd:r(atu)},
{cmd:r(ate)}, {cmd:r(att_cs)}, {cmd:r(atu_cs)}, {cmd:r(eff_first)}, {cmd:r(eff_last)}.


{title:Example}

{phang2}{cmd:. webuse union3}{p_end}
{phang2}{cmd:. esreg ln_wage age grade smsa black tenure, select(union = south black tenure)}{p_end}
{phang2}{cmd:. esrreport, export(reading.txt) replace}{p_end}
