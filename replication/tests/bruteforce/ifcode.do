* ifcode.do <case> -- the influence functions that esreg's own code computes (the
* functions behind the reported standard errors), at the estimate, for ATT, ATU,
* ATE and the mean kappa (_esreg_effects: sampling part U + parameter part P),
* esrmte's m and kappa on the analysis sample x > 0 (_esreg_mte) and, case kx, the slope of kappa(x) (_esr_psi
* times the contrast of the two lambda_x coefficients); with the reported
* standard errors (FIML with vce(robust)).  Cases: see bf2_setup.do.
* Writes out/ifcode_<case>.dta, compared with the brute force by compare2.do.
args case
version 16
clear all
set more off
* run from the folder that holds data/ and tests/: the root of the development
* project, or replication/ of the public repository (the package in ../src/)
local P = subinstr("`c(pwd)'", "\", "/", .)
local S "`P'/src"
capture confirm file "`S'/esreg.ado"
if (_rc) local S "`P'/../src"
local W "`P'/tests/bruteforce/out"
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrmte {
    run "`S'/`f'.ado"
}
mata:
// the influence function of the kappa(x) slope (two-step, kappa(x)): rows of
// -Jac^{-1} g_i (the stacked moments), contrast of the two lambda_x coefficients
void _bf2_slope(real scalar j1, real scalar j0)
{
    real colvector y, d, th
    real matrix X, Z, W1, Psi
    real scalar n
    y = st_data(., "y", "smp"); d = st_data(., "d", "smp"); n = rows(y)
    X  = st_data(., "x", "smp"), J(n, 1, 1)
    Z  = st_data(., ("x", "z"), "smp"), J(n, 1, 1)
    W1 = st_data(., "x", "smp"), J(n, 1, 1)
    th = st_matrix("e(b)")'
    Psi = _esr_psi(2, th, y, X, Z, d, W1, J(n, 0, .), J(n, 1, 1))
    st_store(., "IF_slope", "smp", Psi[., j1] - Psi[., j0])
}
end

do "`P'/tests/bruteforce/bf2_setup.do" `case'
local vce = cond(inlist("`case'", "fiml", "het"), "vce(robust)", "")
esreg y x, select(d = x z) $BF2OPT `vce' nolog
matrix E = e(effects)
matrix list E
qui esrmte if x > 0, at(.5)
local se_m  = r(se_m)
local se_km = r(se_kappa)
tempname VC
mat `VC' = r(V_curve)
local q = colsof(`VC')
local se_dh2 = cond(`q' >= 3, sqrt(`VC'[min(3, `q'), min(3, `q')]), .)
local se_dh3 = cond(`q' >= 4, sqrt(`VC'[min(4, `q'), min(4, `q')]), .)
di as txt "esrmte: m = " r(m) " (" `se_m' "), kappa = " r(kappa) " (" `se_km' ")"
local se_slope = .
if ("`case'" == "kx") {
    lincom [y_1]lambda_x - [y_0]lambda_x
    local se_slope = r(se)
}
gen byte smp = e(sample)
gen double one = 1
_esreg_data if smp
local xl "`r(x)'"
local zl "`r(z)'"
local hs "`r(hs)'"
local hr "`r(hr)'"
local kl "`r(kap)'"
mata: _esreg_effects("y", "`xl'", "`zl'", "d", "`hs'", "`hr'", "`kl'", "one", "smp", "E2", ///
                     "U_att U_atu U_ate U_kap P_att P_atu P_ate P_kap")
gen byte av = smp & x > 0
local un "U_m U_km"
local pn "P_m P_km"
if (`q' >= 3) {
    local un "`un' U_dh2"
    local pn "`pn' P_dh2"
}
if (`q' >= 4) {
    local un "`un' U_dh3"
    local pn "`pn' P_dh3"
}
mata: _esreg_mte("y", "`xl'", "`zl'", "d", "`hs'", "`hr'", "`kl'", "one", "smp", "av", "MK", ///
                 "`un' `pn'", "V2")
gen double IF_slope = .
if ("`case'" == "kx") {
    local cn : colfullnames e(b)
    local j1 : list posof "y_1:lambda_x" in cn
    local j0 : list posof "y_0:lambda_x" in cn
    mata: _bf2_slope(`j1', `j0')
}
foreach s in att atu ate kap m km dh2 dh3 {
    cap confirm variable U_`s'
    if (_rc == 0) gen double IF_`s' = U_`s' + P_`s'
}
local i = 0
foreach s in att atu ate kap {
    local ++i
    gen double rep_`s' = E[`i', 2]
}
gen double rep_m = `se_m'
gen double rep_km = `se_km'
gen double rep_slope = `se_slope'
gen double rep_dh2 = `se_dh2'
gen double rep_dh3 = `se_dh3'
keep id IF_* rep_*
save "`W'/ifcode_`case'.dta", replace

