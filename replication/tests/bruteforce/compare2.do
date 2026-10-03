* compare2.do <case> -- the brute force (bf2.do) against the influence functions
* of the code (ifcode.do): for each statistic, the largest gap between the two on
* the brute-forced observations relative to the largest influence, their
* correlation, and the standard error sqrt(sum IF^2) on the 1,000 observations
* against the reported one.  Cases: see bf2_setup.do.
args case
version 16
clear all
set more off
set linesize 120
local W = subinstr("`c(pwd)'", "\", "/", .) + "/tests/bruteforce/out"
local bf : dir "`W'" files "bf2_`case'_*.dta"
local first 1
foreach f of local bf {
    if (`first') use "`W'/`f'", clear
    else         append using "`W'/`f'"
    local first 0
}
merge 1:1 id using "`W'/ifcode_`case'.dta", nogen
qui count if B_att < .
local nbf = r(N)
di as txt _n "{hline 104}" _n "case `case': brute force (" as res `nbf' ///
   as txt " observations) against the influence functions of the code" _n "{hline 104}"
di as txt "  statistic  {c |} max|BF-IF|/max|IF|  corr(BF, IF)  slope-1 (BF on IF) {c |} se(IF, code)  se reported    ratio"
di as txt "{hline 13}{c +}{hline 52}{c +}{hline 38}"
foreach s in att atu ate kap m km slope dh2 dh3 {
    cap confirm variable IF_`s' B_`s'
    if (_rc) continue
    qui count if IF_`s' < .
    if (r(N) == 0) continue
    tempvar g a
    qui gen double `g' = abs(B_`s' - IF_`s') if B_`s' < .
    qui gen double `a' = abs(IF_`s') if B_`s' < .
    qui summarize `g', meanonly
    local gmax = r(max)
    qui summarize `a', meanonly
    local rel = `gmax' / r(max)
    qui correlate B_`s' IF_`s'
    local rho = r(rho)
    * a systematic error (a missing or wrong term) moves the slope; the
    * convergence tolerance of ml only adds noise around it
    qui regress B_`s' IF_`s', noconstant
    local sl = _b[IF_`s'] - 1
    tempvar q
    qui gen double `q' = IF_`s'^2
    qui summarize `q', meanonly
    local seif = sqrt(r(sum))
    qui summarize rep_`s', meanonly
    local rep = r(mean)
    di as txt %11s "`s'" "  {c |}" as res %17.2e `rel' %15.8f `rho' %17.1e `sl' as txt "   {c |}" ///
       as res %12.6f `seif' %13.6f `rep' %10.5f `rep'/`seif'
}
di as txt "{hline 104}"
di as txt "BF: the whole estimation re-run with the weight of one observation at 1.5 and 0.5 (`nbf' observations)."
di as txt "FIML cases: the reported standard errors use vce(robust), whose parameter part carries n/(n-1)."
