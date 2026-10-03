* bf2.do <case> <first> <last> -- extended brute force: the influence of
* observations first..last on ATT, ATU, ATE, the mean kappa, esrmte's m and kappa
* (on the analysis sample x > 0, to test its own sampling part)
* and (case kx) the slope of kappa(x), (cases hm, hr) the Hermite coefficients dh2 dh3
* of esrmte's curve: the whole estimation (esreg, then esrmte)
* re-run with the weight of one observation moved to 1.5 and to 0.5 (iweights),
* central difference.  Cases: see bf2_setup.do.
* Writes out/bf2_<case>_<first>.dta (id B_att B_atu B_ate B_kap B_m B_km B_slope
* B_dh2 B_dh3).
args case first last
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
do "`P'/tests/bruteforce/bf2_setup.do" `case'
tempname R
matrix `R' = J(`last' - `first' + 1, 10, .)
local r 0
forvalues i = `first'/`last' {
    local ++r
    foreach s in 1.5 0.5 {
        qui replace wt = 1
        qui replace wt = `s' in `i'
        qui esreg y x [iw = wt], select(d = x z) $BF2OPT nolog
        local t = 10*`s'
        local a`t' = e(att)
        local u`t' = e(atu)
        local e`t' = e(ate)
        local k`t' = e(kappa)
        local s`t' = .
        if ("`case'" == "kx") local s`t' = _b[y_1:lambda_x] - _b[y_0:lambda_x]
        qui esrmte if x > 0, at(.5)
        local m`t' = r(m)
        local q`t' = r(kappa)
        local h2`t' = .
        local h3`t' = .
        if (inlist("`case'", "hm", "hr")) {
            tempname CV
            mat `CV' = r(curve)
            local h2`t' = `CV'[1, 3]
            local h3`t' = `CV'[1, 4]
        }
    }
    matrix `R'[`r', 1] = `i'
    matrix `R'[`r', 2] = `a15' - `a5'
    matrix `R'[`r', 3] = `u15' - `u5'
    matrix `R'[`r', 4] = `e15' - `e5'
    matrix `R'[`r', 5] = `k15' - `k5'
    matrix `R'[`r', 6] = `m15' - `m5'
    matrix `R'[`r', 7] = `q15' - `q5'
    matrix `R'[`r', 8] = `s15' - `s5'
    matrix `R'[`r', 9] = `h215' - `h25'
    matrix `R'[`r', 10] = `h315' - `h35'
    if (mod(`r', 10) == 0) di as txt "bf2 `case': observation `i' done"
}
clear
svmat double `R'
rename (`R'1 `R'2 `R'3 `R'4 `R'5 `R'6 `R'7 `R'8 `R'9 `R'10) (id B_att B_atu B_ate B_kap B_m B_km B_slope B_dh2 B_dh3)
save "`W'/bf2_`case'_`first'.dta", replace
