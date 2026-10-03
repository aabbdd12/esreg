* bf.do -- brute-force influence of observations first..last on ATT, ATU, ATE
* and mean kappa: the whole esreg estimation re-run with the weight of one
* observation moved to 1.5 and to 0.5 (iweights), central difference.
* Writes bf_<method>_<first>.dta (id, B_att, B_atu, B_ate, B_kap).
args method first last
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
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p {
    run "`S'/`f'.ado"
}
use "`P'/data/esr_sim.dta", clear
keep in 1/1000
gen long id = _n
gen double wt = 1
tempname R
matrix `R' = J(`last' - `first' + 1, 5, .)
local r 0
forvalues i = `first'/`last' {
    local ++r
    foreach s in 1.5 0.5 {
        qui replace wt = 1
        qui replace wt = `s' in `i'
        qui esreg y x [iw = wt], select(d = x z) method(`method') nolog
        local a`=10*`s'' = e(att)
        local u`=10*`s'' = e(atu)
        local e`=10*`s'' = e(ate)
        local k`=10*`s'' = e(kappa)
    }
    matrix `R'[`r', 1] = `i'
    matrix `R'[`r', 2] = (`a15' - `a5')
    matrix `R'[`r', 3] = (`u15' - `u5')
    matrix `R'[`r', 4] = (`e15' - `e5')
    matrix `R'[`r', 5] = (`k15' - `k5')
}
clear
svmat double `R'
rename (`R'1 `R'2 `R'3 `R'4 `R'5) (id B_att B_atu B_ate B_kap)
save "`W'/bf_`method'_`first'.dta", replace
