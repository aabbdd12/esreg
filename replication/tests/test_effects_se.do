* test_effects_se.do -- the standard errors of esreg (effects, two-step
* clusters and survey design) against oracles, with tolerances and a pass/fail
* summary.  Run from the project root:  tools\run_stata.bat tests\test_effects_se.do
*
*  1. iid: the standard errors of ATT, ATU, ATE, kappa equal sqrt(sum D_i^2), D_i
*     the exact influence function of tests/bruteforce/out (itself checked against
*     the brute force: the whole estimation re-run with the weight of one
*     observation moved); two-step, and FIML vce(robust) (ml's robust e(V) carries
*     n/(n-1), removed from the parameter part for the comparison).
*  2. clusters: the influence function of theta aggregated by cluster equals ml's
*     vce(cluster) e(V) (FIML); the effects by cluster equal sqrt(G/(G-1) sum_g
*     D_g^2) for both routes (two-step: the new vce(cluster)).
*  3. survey design: the influence function of theta through svy: total equals
*     the svy prefix's e(V) (FIML, strata and PSUs); two identities of the
*     two-step vce(svy): with PSUs and no strata it equals vce(cluster psu); with
*     one unit per PSU and no strata it equals n/(n-1) times the default.
*  4. subpopulation: after svy, subpop() the effects, esrdiag and esrtest are
*     those of the subpopulation (= esreg [pw] if subpop).
*  5. refusals: fweights, vce(svy) with fiml, vce(oim) with twostep, svy prefix
*     with twostep.
*  6. post-commands: esrmte se(kappa) = e(se_kappa) under constant kappa; esrcurve
*     and esrmte run by cluster and by design.
*  7. pseudo-DiD (esrtest, pdid): estimates and influence-function standard
*     errors = the Python reference; by cluster and by design with one unit per
*     cluster or PSU = sqrt(n/(n-1)) times the default; the influence function
*     against the bootstrap of the whole procedure (reps()).
*  8. extended brute force (tests/bruteforce/bf2.do, stored): the influence
*     functions of the code for two-step kappa(), FIML hetsigma()/hetrho(), the
*     augmented two-step hermite(3) and esrmte on a subsample (with the Hermite
*     terms of its curve), observation by observation.
version 16
clear all
set more off
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrcurve esrmte esrdiag esrtest {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}
global npass 0
global nfail 0
cap program drop _chk
program define _chk
    args label got want tol
    local rd = reldif(`got', `want')
    if (`rd' <= `tol' & `got' < . & `want' < .) {
        di as txt "  pass  `label'" _col(58) "got " as res %12.8g `got' as txt "  want " as res %12.8g `want' as txt "  reldif " as res %8.1e `rd'
        global npass = $npass + 1
    }
    else {
        di as err "  FAIL  `label'" _col(58) "got " %12.8g `got' "  want " %12.8g `want' "  reldif " %8.1e `rd'
        global nfail = $nfail + 1
    }
end
cap program drop _chkrc
program define _chkrc
    args label rc want
    if (`rc' == `want') {
        di as txt "  pass  `label' (rc = `rc')"
        global npass = $npass + 1
    }
    else {
        di as err "  FAIL  `label': rc = `rc', expected `want'"
        global nfail = $nfail + 1
    }
end

* ============================================================================
di as txt _n "{hline 78}" _n "1. iid: standard errors = sqrt(sum D^2) of the exact influence function" _n "{hline 78}"
use "data/esr_sim.dta", clear
keep in 1/1000
gen long id = _n
gen long cl = ceil(id / 10)
merge 1:1 id using "tests/bruteforce/out/ifexact_twostep.dta", nogen
rename D_* T_*
merge 1:1 id using "tests/bruteforce/out/ifexact_fiml_robust.dta", nogen
rename D_* F_*
local rows "att atu ate kap"
* exact variances, iid and by cluster (G/(G-1) sum over the clusters)
foreach r in T F {
    local i = 0
    foreach v of local rows {
        local ++i
        tempvar d2 cs c2
        qui gen double `d2' = `r'_`v'^2
        qui summarize `d2', meanonly
        local `r'iid`i' = r(sum)
        qui egen double `cs' = total(`r'_`v'), by(cl)
        qui bysort cl: gen double `c2' = `cs'^2 if _n == 1
        qui summarize `c2', meanonly
        local `r'cl`i' = r(sum) * 100 / 99
        sort id
    }
}
esreg y x, select(d = x z) method(twostep) nolog
mat E = e(effects)
local i = 0
foreach v of local rows {
    local ++i
    _chk "two-step se(`v')" E[`i',2] sqrt(`Tiid`i'') 1e-6
}
esreg y x, select(d = x z) vce(robust) nolog
mat E = e(effects)
local n = e(N)
local i = 0
foreach v of local rows {
    local ++i
    _chk "FIML robust se(`v') (param x (n-1)/n)" "sqrt(E[`i',3]*(`n'-1)/`n' + E[`i',4] + E[`i',5])" "sqrt(`Fiid`i'')" 1e-5
}

* ============================================================================
di as txt _n "{hline 78}" _n "2. clusters (100 clusters of 10)" _n "{hline 78}"
esreg y x, select(d = x z) vce(cluster cl) nolog
mat Vml = e(V)
* the influence function of theta aggregated by cluster, against ml's vce(cluster)
mata:
    y = st_data(., "y"); d = st_data(., "d"); n = rows(y)
    X = st_data(., "x"), J(n, 1, 1); Z = st_data(., ("x", "z")), J(n, 1, 1)
    W1 = J(n, 1, 1)
    th = st_matrix("e(b)")'
    Psi = _esr_psi(1, th, y, X, Z, d, W1, W1, J(n, 1, 1))
    Us = _esr_clsum(Psi, st_data(., "cl"))
    G = rows(Us)
    Vc = cross(Us, Us) * G / (G - 1)
    st_matrix("Vc", Vc)
end
local p = colsof(Vml)
forvalues j = 1/`p' {
    _chk "FIML cluster V[`j',`j']: Psi by cluster = ml" Vc[`j',`j'] Vml[`j',`j'] 1e-5
}
mat E = e(effects)
local i = 0
foreach v of local rows {
    local ++i
    _chk "FIML cluster se(`v')" E[`i',2] sqrt(`Fcl`i'') 1e-5
}
esreg y x, select(d = x z) method(twostep) vce(cluster cl) nolog
mat E = e(effects)
_chk "two-step cluster: e(N_clust)" e(N_clust) 100 0
local i = 0
foreach v of local rows {
    local ++i
    _chk "two-step cluster se(`v')" E[`i',2] sqrt(`Tcl`i'') 1e-6
}

* ============================================================================
di as txt _n "{hline 78}" _n "3. survey design (data/esr_wt.dta, 250 PSUs in 5 strata)" _n "{hline 78}"
use "data/esr_wt.dta", clear
gen int psu = ceil(_n / 20)
gen byte str = mod(psu, 5) + 1
set seed 12345
gen byte sub = runiform() < 0.7
svyset psu [pw = wt], strata(str)
svy: esreg income educ, select(treatment = educ inst)
mat Vsvy = e(V)
mata:
    y = st_data(., "income"); d = st_data(., "treatment"); n = rows(y)
    X = st_data(., "educ"), J(n, 1, 1); Z = st_data(., ("educ", "inst")), J(n, 1, 1)
    W1 = J(n, 1, 1)
    th = st_matrix("e(b)")'
    Psi = _esr_psi(1, th, y, X, Z, d, W1, W1, st_data(., "wt"))
    (void) st_addvar("double", "psi" :+ strofreal(1..cols(Psi)))
    st_store(., "psi" :+ strofreal(1..cols(Psi)), Psi)
end
_esreg_ifcov psi*, svy
mat Vm = r(V)
local p = colsof(Vsvy)
forvalues j = 1/`p' {
    _chk "FIML svy V[`j',`j']: svy total of Psi/w = svy prefix" Vm[`j',`j'] Vsvy[`j',`j'] 1e-5
}
drop psi*
_chk "FIML svy: effects use t (e(df_r))" e(df_r) 245 0
* the two-step vce(svy): PSUs without strata = vce(cluster psu)
svyset psu [pw = wt]
esreg income educ, select(treatment = educ inst) method(twostep) vce(svy) nolog
mat Vs = e(V)
mat Es = e(effects)
esreg income educ [pw = wt], select(treatment = educ inst) method(twostep) vce(cluster psu) nolog
mat Vc = e(V)
mat Ec = e(effects)
local p = colsof(Vs)
forvalues j = 1/`p' {
    _chk "two-step vce(svy) = vce(cluster psu): V[`j',`j']" Vs[`j',`j'] Vc[`j',`j'] 1e-6
}
forvalues i = 1/4 {
    _chk "two-step vce(svy) = vce(cluster psu): effect `i' se" Es[`i',2] Ec[`i',2] 1e-6
}
* one unit per PSU, no strata: n/(n-1) times the default
svyset _n [pw = wt]
esreg income educ, select(treatment = educ inst) method(twostep) vce(svy) nolog
mat Vs = e(V)
mat Es = e(effects)
local n = e(N)
esreg income educ [pw = wt], select(treatment = educ inst) method(twostep) nolog
mat Vd = e(V)
mat Ed = e(effects)
local p = colsof(Vs)
forvalues j = 1/`p' {
    _chk "two-step vce(svy), SRS = n/(n-1) default: V[`j',`j']" Vs[`j',`j'] Vd[`j',`j']*`n'/(`n'-1) 1e-6
}
forvalues i = 1/4 {
    _chk "two-step vce(svy), SRS = n/(n-1) default: var effect `i'" Es[`i',2]^2 Ed[`i',2]^2*`n'/(`n'-1) 1e-6
}

* ============================================================================
di as txt _n "{hline 78}" _n "4. subpopulation" _n "{hline 78}"
svyset psu [pw = wt], strata(str)
svy, subpop(sub): esreg income educ, select(treatment = educ inst)
local a1 = e(att)
local k1 = e(kappa)
local s1 = e(se_att)
esreg income educ [pw = wt] if sub, select(treatment = educ inst) nolog
_chk "svy subpop: ATT = esreg [pw] if sub" `a1' e(att) 1e-6
_chk "svy subpop: kappa = esreg [pw] if sub" `k1' e(kappa) 1e-6
_chk "svy subpop: se(ATT) defined" "(`s1' < .)" 1 0
* the post-commands after svy, subpop() work on the subpopulation too
qui svy, subpop(sub): esreg income educ, select(treatment = educ inst)
qui esrdiag
local lr1 = r(lr_excl)
qui esrtest, spec
local c1 = r(chi2)
qui esrtest, pdid
local k1 = r(kappa_dd)
qui esreg income educ [pw = wt] if sub, select(treatment = educ inst) nolog
qui esrdiag
_chk "svy subpop: esrdiag LR = esreg [pw] if sub" `lr1' r(lr_excl) 1e-6
qui esrtest, spec
_chk "svy subpop: esrtest spec chi2 = esreg [pw] if sub" `c1' r(chi2) 1e-6
qui esrtest, pdid
_chk "svy subpop: esrtest pdid kappa_dd = esreg [pw] if sub" `k1' r(kappa_dd) 1e-6

* ============================================================================
di as txt _n "{hline 78}" _n "5. refusals" _n "{hline 78}"
gen int fw = 1 + mod(_n, 3)
cap esreg income educ [fw = fw], select(treatment = educ inst) nolog
local rc = _rc
_chkrc "fweight refused" `rc' 101
cap esreg income educ, select(treatment = educ inst) vce(svy) nolog
local rc = _rc
_chkrc "vce(svy) with method(fiml) refused" `rc' 198
cap esreg income educ, select(treatment = educ inst) method(twostep) vce(oim) nolog
local rc = _rc
_chkrc "vce(oim) with method(twostep) refused" `rc' 198
cap svy: esreg income educ, select(treatment = educ inst) method(twostep)
local rc = (_rc != 0)
_chkrc "svy prefix with method(twostep) refused" `rc' 1

* ============================================================================
di as txt _n "{hline 78}" _n "6. post-commands" _n "{hline 78}"
esreg income educ [pw = wt], select(treatment = educ inst) method(twostep) nolog
local sk = e(se_kappa)
esrmte, at(.5)
_chk "esrmte se(kappa) = e(se_kappa) (constant kappa)" r(se_kappa) `sk' 1e-8
esreg income educ [pw = wt], select(treatment = educ inst) method(twostep) vce(cluster psu) nolog
local sk = e(se_kappa)
esrmte, at(.5)
_chk "esrmte se(kappa) = e(se_kappa), by cluster" r(se_kappa) `sk' 1e-8
esrcurve, rank(educ) nq(4)
mat R = r(table)
_chk "esrcurve by cluster: se defined" "(R[1,3] < . & R[4,7] < .)" 1 0
svyset psu [pw = wt], strata(str)
esreg income educ, select(treatment = educ inst) method(twostep) vce(svy) nolog
local sk = e(se_kappa)
esrmte, at(.5)
_chk "esrmte se(kappa) = e(se_kappa), survey design" r(se_kappa) `sk' 1e-8
esrcurve, rank(educ) nq(4)
mat R = r(table)
_chk "esrcurve survey design: se defined" "(R[1,3] < . & R[4,7] < .)" 1 0

* ============================================================================
di as txt _n "{hline 78}" _n "7. pseudo-DiD (esrtest, pdid)" _n "{hline 78}"
* the estimates and the influence-function standard errors against the Python
* reference (replication/python/esreg/tests.py: pseudo_did, an independent
* implementation of the same estimating equations), esr_sim, nq = 2 and 4
use "data/esr_sim.dta", clear
qui esreg y x, select(d = x z) method(twostep) nolog
qui esrtest, pdid
foreach s in "C1 -0.4846735411 0.02472017532" "C0 0.2452658675 0.0248215396" ///
             "pi1 -0.6596334444 0.006812713775" "pi0 0.7766444383 0.008829046088" ///
             "rhosig1 0.7347619276 0.03787168079" "rhosig0 -0.3158020008 0.03200808242" ///
             "kappa_dd 1.050563928 0.05013674959" "dd -0.4698060433 0.03763929233" {
    tokenize `s'
    _chk "pdid nq=2: `1' = Python" r(`1') `2' 1e-6
    _chk "pdid nq=2: se(`1') = Python" r(se_`1') `3' 1e-6
}
qui esrtest, pdid nq(4)
foreach s in "C1 -0.8655469378 0.04253851768" "C0 0.4990154893 0.04656259321" ///
             "kappa_dd 1.101150477 0.0508871514" "dd -0.8116029707 0.07141332312" {
    tokenize `s'
    _chk "pdid nq=4: `1' = Python" r(`1') `2' 1e-6
    _chk "pdid nq=4: se(`1') = Python" r(se_`1') `3' 1e-6
}
* by cluster and by design, one observation per cluster or PSU: sqrt(n/(n-1))
* times the default
local n = _N
gen long id = _n
qui esrtest, pdid
local sk = r(se_kappa_dd)
local sd = r(se_dd)
qui esreg y x, select(d = x z) method(twostep) vce(cluster id) nolog
qui esrtest, pdid
_chk "pdid by cluster (1 obs/cluster): se(kappa_dd)" r(se_kappa_dd) `sk'*sqrt(`n'/(`n'-1)) 1e-6
_chk "pdid by cluster (1 obs/cluster): se(dd)" r(se_dd) `sd'*sqrt(`n'/(`n'-1)) 1e-6
svyset id
qui esreg y x, select(d = x z) method(twostep) vce(svy) nolog
qui esrtest, pdid
_chk "pdid survey design (1 obs/PSU): se(kappa_dd)" r(se_kappa_dd) `sk'*sqrt(`n'/(`n'-1)) 1e-6
_chk "pdid survey design (1 obs/PSU): se(dd)" r(se_dd) `sd'*sqrt(`n'/(`n'-1)) 1e-6
* the influence function against the bootstrap of the whole procedure (reps(),
* B = 200, fixed seed; the bootstrap noise on a standard error is about 5%):
* ratio within 15% of 1 (reldif(ratio, 1) = |ratio - 1|/2 <= .075)
use "data/esr_sim.dta", clear
keep in 1/5000
qui esreg y x, select(d = x z) method(twostep) nolog
qui esrtest, pdid
foreach k in C1 C0 kappa_dd {
    local if_`k' = r(se_`k')
}
set seed 20261002
qui esrtest, pdid reps(200)
foreach k in C1 C0 kappa_dd {
    _chk "pdid se(`k'): influence function / bootstrap" `if_`k''/r(se_`k') 1 .075
}
local ok = ("`e(cmd)'" == "esreg")
_chk "pdid reps(): e() restored" `ok' 1 0

* ============================================================================
di as txt _n "{hline 78}" _n "8. extended brute force: kappa(), hetsigma/hetrho, hermite(), esrmte" _n "{hline 78}"
* the influence functions of the code (the effects, esrmte's m and kappa on the
* analysis sample x > 0, the kappa(x) slope, the Hermite terms dh2 dh3 of the
* augmented two-step's MTE curve) against the brute force stored in
* tests/bruteforce/out/bf2_<case>_*.dta (the whole estimation re-run with the
* weight of one observation moved, observations 1-100; tests/bruteforce/bf2.do):
* slope of BF on IF = 1 and largest gap relative to the largest influence, within
* the convergence noise of the estimation (two-step 1e-3, FIML 5e-3)
cap mata: mata drop _t8_slope()
mata:
void _t8_slope(real scalar j1, real scalar j0)
{
    real colvector y, d
    real matrix X, Z, W1, Psi
    real scalar n
    y = st_data(., "y", "smp"); d = st_data(., "d", "smp"); n = rows(y)
    X  = st_data(., "x", "smp"), J(n, 1, 1)
    Z  = st_data(., ("x", "z"), "smp"), J(n, 1, 1)
    W1 = st_data(., "x", "smp"), J(n, 1, 1)
    Psi = _esr_psi(2, st_matrix("e(b)")', y, X, Z, d, W1, J(n, 0, .), J(n, 1, 1))
    st_store(., "IF_slope", "smp", Psi[., j1] - Psi[., j0])
}
end
foreach c in ts fiml kx het hm hr {
    qui do "tests/bruteforce/bf2_setup.do" `c'
    qui esreg y x, select(d = x z) $BF2OPT nolog
    qui gen byte smp = e(sample)
    qui gen byte av = smp & x > 0
    qui gen double one = 1
    qui _esreg_data if smp
    local xl "`r(x)'"
    local zl "`r(z)'"
    local hs "`r(hs)'"
    local hr "`r(hr)'"
    local kl "`r(kap)'"
    mata: _esreg_effects("y", "`xl'", "`zl'", "d", "`hs'", "`hr'", "`kl'", "one", "smp", "E8", ///
                         "U_att U_atu U_ate U_kap P_att P_atu P_ate P_kap")
    local un "U_m U_km"
    local pn "P_m P_km"
    local stats "att atu ate kap m km"
    if (inlist("`c'", "hm", "hr")) {
        local un "`un' U_dh2 U_dh3"
        local pn "`pn' P_dh2 P_dh3"
        local stats "`stats' dh2 dh3"
    }
    mata: _esreg_mte("y", "`xl'", "`zl'", "d", "`hs'", "`hr'", "`kl'", "one", "smp", "av", "MK8", ///
                     "`un' `pn'", "V8")
    if ("`c'" == "kx") {
        qui gen double IF_slope = .
        local cn : colfullnames e(b)
        local j1 : list posof "y_1:lambda_x" in cn
        local j0 : list posof "y_0:lambda_x" in cn
        mata: _t8_slope(`j1', `j0')
        local stats "`stats' slope"
    }
    foreach s in att atu ate kap m km dh2 dh3 {
        cap confirm variable U_`s'
        if (_rc == 0) qui gen double IF_`s' = U_`s' + P_`s'
    }
    local bf : dir "tests/bruteforce/out" files "bf2_`c'_*.dta"
    foreach f of local bf {
        qui merge 1:1 id using "tests/bruteforce/out/`f'", nogen update
    }
    local tol = cond(inlist("`c'", "fiml", "het"), 5e-3, 1e-3)
    foreach s of local stats {
        tempvar g a
        qui gen double `g' = abs(B_`s' - IF_`s') if B_`s' < .
        qui gen double `a' = abs(IF_`s') if B_`s' < .
        qui summarize `g', meanonly
        local gmax = r(max)
        qui summarize `a', meanonly
        local rel = `gmax' / r(max)
        qui regress B_`s' IF_`s', noconstant
        _chk "brute force `c': slope of BF on IF(`s') = 1" _b[IF_`s'] 1 `=`tol'/10'
        _chk "brute force `c': max|BF - IF(`s')| / max|IF|" `rel' 0 `tol'
    }
}
* a bare hermite is hermite(3), and e(cmdline) keeps the command as typed
qui esreg y x, select(d = x z) method(twostep) hermite nolog
local ok = (`"`e(cmdline)'"' == "esreg y x, select(d = x z) method(twostep) hermite nolog")
_chk "e(cmdline) as typed (bare hermite)" `ok' 1 0
_chk "bare hermite = hermite(3)" e(hermite) 3 0
* the Hermite terms by regime, hermite(#1 #0): ATT depends on the untreated
* regime's terms only and ATU on the treated regime's (the treated mean of the
* fitted Y_1 is the treated mean of y whatever the treated regime's regressors):
* hermite(3 0) has the line's ATT and hermite(3)'s ATU, hermite(0 3) the line's
* ATU and hermite(3)'s ATT, estimates and standard errors
qui esreg y x, select(d = x z) method(twostep) nolog
scalar t8_ln_att = e(att)
scalar t8_ln_sat = e(se_att)
scalar t8_ln_atu = e(atu)
scalar t8_ln_sau = e(se_atu)
qui esreg y x, select(d = x z) method(twostep) hermite(3) nolog
scalar t8_hb_att = e(att)
scalar t8_hb_sat = e(se_att)
scalar t8_hb_atu = e(atu)
scalar t8_hb_sau = e(se_atu)
qui esreg y x, select(d = x z) method(twostep) hermite(3 0) nolog
_chk "hermite(3 0): e(hermite1), e(hermite0) = 3, 0" `=10*e(hermite1) + e(hermite0)' 30 0
_chk "hermite(3 0): ATT = the line's" e(att) scalar(t8_ln_att) 1e-10
_chk "hermite(3 0): se(ATT) = the line's" e(se_att) scalar(t8_ln_sat) 1e-7
_chk "hermite(3 0): ATU = hermite(3)'s" e(atu) scalar(t8_hb_atu) 1e-10
_chk "hermite(3 0): se(ATU) = hermite(3)'s" e(se_atu) scalar(t8_hb_sau) 1e-7
qui esreg y x, select(d = x z) method(twostep) hermite(0 3) nolog
_chk "hermite(0 3): ATT = hermite(3)'s" e(att) scalar(t8_hb_att) 1e-10
_chk "hermite(0 3): se(ATT) = hermite(3)'s" e(se_att) scalar(t8_hb_sat) 1e-7
_chk "hermite(0 3): ATU = the line's" e(atu) scalar(t8_ln_atu) 1e-10
_chk "hermite(0 3): se(ATU) = the line's" e(se_atu) scalar(t8_ln_sau) 1e-7

di as txt _n "{hline 78}"
di as txt "test_effects_se: " as res $npass as txt " passed, " as res $nfail as txt " failed"
if ($nfail > 0) exit 9
