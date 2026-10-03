* mc_summary.do <design> [routes] -- combines the blocks of mc_se.do (mc_aug.do)
* and reports, for each statistic: the true value, the mean estimate (bias), the
* standard deviation of the estimates over the replications, the root mean square
* of the analytic standard errors, their ratio (with its Monte Carlo precision
* 1/sqrt(2R)), the coverage of the 95% interval (with its precision), and the root
* mean square error of the estimates.  routes (designs AN, AC): other routes of
* mc_aug.do computed later on the same samples (files mc_<design>-<routes>_*.dta),
* merged by replication with the stored ones.
* True values (closed forms, checked by a simulation of 4 x 10^7 draws):
*   N: ATT 1.1778739, ATU -0.3546754, ATE 0.5, kappa 1.08, kappa_DD 1.08,
*      MTE(u) = 0.5 + 1.08 invnormal(1 - u): 1.884076, 0.5, -0.884076
*   K: ATT 1.1713871, ATU -0.3464967, ATE 0.5, mean kappa 1.08, slope 0.4, MTE as N
*   AN (mc_aug.do, design N): as N, and dh2 = dh3 = 0
*   AC (mc_aug.do, design C1): ATT 1.0349929, ATU -0.1745285, ATE 0.5, kappa 1.08,
*      dh2 0, dh3 0.54, MTE(u) = 0.5 + 1.08 v + 0.54 (v^3 - 3v), v = invnormal(1 - u):
*      0.9445474, 0.5, 0.0554526
args design routes
version 16
clear all
set more off
set linesize 120
local W = subinstr("`c(pwd)'", "\", "/", .) + "/tests/mc/out"
local files : dir "`W'" files "mc_`design'_*.dta"
local first 1
foreach f of local files {
    if (`first') use "`W'/`f'", clear
    else         append using "`W'/`f'"
    local first 0
}
if ("`routes'" != "") {
    local rtag = subinstr("`routes'", " ", "", .)
    tempfile base xr
    qui save `base'
    local xf : dir "`W'" files "mc_`design'-`rtag'_*.dta"
    local first 1
    foreach f of local xf {
        if (`first') use "`W'/`f'", clear
        else         append using "`W'/`f'"
        local first 0
    }
    qui save `xr'
    use `base', clear
    qui merge 1:1 rep using `xr', keep(match) nogen
}
duplicates report rep
local R = _N
local t_att  = cond(inlist("`design'", "N", "AN"), 1.1778739, cond("`design'" == "AC", 1.0349929, 1.1713871))
local t_atu  = cond(inlist("`design'", "N", "AN"), -0.3546754, cond("`design'" == "AC", -0.1745285, -0.3464967))
local t_ate  = 0.5
local t_kap  = 1.08
local t_kdd  = 1.08
local t_slope = 0.4
local t_m10  = 0.5 + 1.08*invnormal(.9)
local t_m50  = 0.5
local t_m90  = 0.5 + 1.08*invnormal(.1)
local t_dh2  = 0
local t_dh3  = cond("`design'" == "AC", 0.54, 0)
if ("`design'" == "AC") {
    foreach p in 10 50 90 {
        local v = invnormal(1 - `p'/100)
        local t_m`p' = 0.5 + 1.08*`v' + 0.54*(`v'^3 - 3*`v')
    }
}
if ("`design'" == "N")     local stats "ts_att ts_atu ts_ate ts_kap ts_m10 ts_m50 ts_m90 kdd fi_att fi_atu fi_ate fi_kap fi_m10 fi_m50 fi_m90"
else if ("`design'" == "K") local stats "kx_att kx_atu kx_ate kx_kap kx_slope kx_m10 kx_m50 kx_m90"
else local stats "ln_att ln_atu ln_ate ln_kap ln_m10 ln_m50 ln_m90 h3_att h3_atu h3_ate h3_kap h3_m10 h3_m50 h3_m90 h3_dh2 h3_dh3"
foreach route of local routes {
    local stats "`stats' `route'_att `route'_atu `route'_ate `route'_kap `route'_m10 `route'_m50 `route'_m90 `route'_dh2 `route'_dh3"
}
tempname S
mat `S' = J(`: word count `stats'', 10, .)
di as txt _n "{hline 118}"
di as txt "Monte Carlo, design `design': " as res `R' as txt " replications; ratio = rms(se) / sd(estimates), precision about +-" ///
   as res %4.1f 100/sqrt(2*`R') as txt "%; coverage precision about +-" as res %4.1f 100*sqrt(.95*.05/`R') as txt "%"
di as txt "{hline 118}"
di as txt "  statistic {c |}     true      mean      bias        sd   rms(se)  mean(se)     ratio  coverage95   valid      rmse"
di as txt "{hline 11}{c +}{hline 106}"
local i = 0
foreach s of local stats {
    local ++i
    local k = substr("`s'", strpos("`s'", "_") + 1, .)
    if ("`s'" == "kdd") local k "kdd"
    local tv = `t_`k''
    qui count if `s' < . & se_`s' < .
    local nv = r(N)
    qui summarize `s' if se_`s' < .
    local mn = r(mean)
    local sd = r(sd)
    tempvar q c
    qui gen double `q' = se_`s'^2 if `s' < .
    qui summarize `q', meanonly
    local rms = sqrt(r(mean))
    qui summarize se_`s' if `s' < ., meanonly
    local mse = r(mean)
    qui gen byte `c' = abs(`s' - `tv') <= invnormal(.975)*se_`s' if `s' < . & se_`s' < .
    qui summarize `c', meanonly
    local cov = r(mean)
    di as txt %10s "`s'" " {c |}" as res %9.4f `tv' %10.4f `mn' %10.4f `mn' - `tv' %10.4f `sd' %10.4f `rms' ///
       %10.4f `mse' %10.4f `rms'/`sd' %12.4f `cov' %8.0f `nv' %10.4f sqrt((`mn' - `tv')^2 + `sd'^2)
    mat `S'[`i', 1] = `tv'
    mat `S'[`i', 2] = `mn'
    mat `S'[`i', 3] = `sd'
    mat `S'[`i', 4] = `rms'
    mat `S'[`i', 5] = `mse'
    mat `S'[`i', 6] = `rms'/`sd'
    mat `S'[`i', 7] = `cov'
    mat `S'[`i', 8] = `nv'
    mat `S'[`i', 9] = `R'
    mat `S'[`i', 10] = sqrt((`mn' - `tv')^2 + `sd'^2)
}
di as txt "{hline 118}"
mat rownames `S' = `stats'
mat colnames `S' = true mean sd rms_se mean_se ratio coverage valid R rmse
mat list `S', format(%9.4f)
* identities of the routes by regime, replication by replication: ATT depends on
* the untreated regime's Hermite terms only, ATU on the treated regime's
foreach route of local routes {
    if ("`route'" == "r1") local pairs "r1_att ln_att r1_atu h3_atu"
    if ("`route'" == "r0") local pairs "r0_att h3_att r0_atu ln_atu"
    if (!inlist("`route'", "r1", "r0")) continue
    while ("`pairs'" != "") {
        gettoken a pairs : pairs
        gettoken b pairs : pairs
        foreach v in "" "se_" {
            tempvar g
            qui gen double `g' = abs(`v'`a' - `v'`b')
            qui summarize `g', meanonly
            di as txt "identity `v'`a' = `v'`b': largest |difference| over the replications " as res %9.2e r(max)
        }
    }
}
