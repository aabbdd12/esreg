* esreg_test.do -- step 2 of the progressive build: esreg.ado against the
* Python reference (esreg v0.1) and against movestay / msat.
* Run from the root of the replication folder (sogrepl): datasets in data/.
* Put esreg.ado in the current folder or on the adopath, then:  do esreg_test.do

version 16
clear all
set more off
* reload the three ado files from disk on every run (they drop their own previous definitions)
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}

* ---- 1. esr_sim, FIML -----------------------------------------------------
use "data/esr_sim.dta", clear
esreg y x, select(d = x z)
di as txt _n "  expected (Python / movestay): ll -39940.2911, LR 787.01"
di as txt   "  y_1: x 1.205903 (.0121553) _cons .9831363 (.0189707); y_0: x .3865108 _cons .4803906"
di as txt   "  d: x .292636 z .9081515 _cons .2042513; lnsigma_1 .2549284 lnsigma_0 -.0094671"
di as txt   "  atanhrho_1 .699128 (.0276441) atanhrho_0 -.323829 (.0263308)"
di as txt   "  ATT 1.17919 (.0342351) ATU -.343525 (.0380631) ATE .508129 (.0278976) kappa 1.08916 (.0373076)"
di as txt   "  var comp ATT .00111 / 5.84e-05 ; support [0.0213, 0.9910]"

* ---- 2. esr_sim, two-step ---------------------------------------------------
esreg y x, select(d = x z) method(twostep)
di as txt _n "  expected (Python): d: x .292311 (.0105781) z .913279 (.0129668) _cons .203296 (.0102485)"
di as txt   "  y_1: x 1.20116 (.0121475) lambda .733444 (.0311533) _cons 1.00788 (.020892)"
di as txt   "  y_0: x .386429 (.0109357) lambda -.31011 (.025222) _cons .480755 (.0197401)"
di as txt   "  ATT 1.17991 (.0345843) ATU -.287578 (.0424281) ATE .53319 (.0293374) kappa 1.04355 (.0406198)"
di as txt   "  sigma_1 1.28158 sigma_0 .990469 rho_1 .572297 rho_0 -.313094 ; support [0.0207, 0.9912]"

* ---- 3. esr_annex2, FIML with a factor variable -----------------------------
use "data/esr_annex2.dta", clear
esreg income educ, select(treatment = educ i.region inst)
di as txt _n "  expected: ll -4229.3298, LR 172.38; income_1: educ .5039501 _cons 62.8872; income_0: educ .503302 _cons 60.9076"
di as txt   "  treatment: educ -.0034487 2.region .2872393 3.region .3312872 inst .7008038 _cons .3195647"
di as txt   "  ATT 1.98115 (.0291947) ATU 1.98233 (.0267532) ATE 1.98155 (.0230826) kappa -.000887548 (.0263384)"
di as txt   "  support [0.0421, 0.9754]"

* ---- 4. esr_skew: FIML then two-step -----------------------------------------
use "data/esr_skew.dta", clear
esreg y x, select(d = x z)
di as txt _n "  expected: ll -20670.0850, LR 732.70; ATT 2.33027 (.0544608) ATU -1.63347 (.0622513) ATE .57433 (.0375635) kappa 2.82488 (.0511241)"
esreg y x, select(d = x z) method(twostep)
di as txt _n "  expected: ATT 1.25885 (.0522976) ATU -.344916 (.0733844) ATE .548381 (.0450786) kappa 1.16692 (.0674883)"
di as txt   "  (bootstrap-100 of msat: .0539 .0747 .0446 .0742)"

* ---- 5. esr_nonlin, two-step ---------------------------------------------------
use "data/esr_nonlin.dta", clear
esreg y x, select(d = x z) method(twostep)
di as txt _n "  expected: ATT 1.08803 (.0561481) ATU -.210066 (.117275) ATE .517385 (.0592842) kappa .87816 (.109231)"

* ---- 6. replay, stored results, automatic and named storage --------------------
esreg
ereturn list
estimates dir
use "data/esr_sim.dta", clear
esreg y x, select(d = x z) method(twostep) store(ts) nolog
regress y x                       // some other estimation in between
estimates restore _esreg         // the last esreg comes back
esreg
estimates restore ts
di as txt "ATT from the stored two-step: " as res e(att) as txt "   (expected 1.17991)"
estimates dir
