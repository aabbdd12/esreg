* esrdiag_test.do -- diagnostics of the selection equation (Section 4.4) against the
* Python reference (esreg/diag.py).  Run from the project root:
*   tools\run_stata.bat tests\esrdiag_test.do
version 16
clear all
set more off
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrdiag {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}

* ---- esr_sim, FIML ---------------------------------------------------------------------
use "data/esr_sim.dta", clear
esreg y x, select(d = x z) nolog
esrdiag
di as txt _n "expected (Python): LR chi2(1) = 6816.98; pseudo-R2 .26832 / .01993 (incr .2484); z: chi2 4931.49"
di as txt   "   Var(P|X)/Var(P) .91537; VIF 1.06695 / 1.05808; shares extrapolated .02414 / .02417"
di as txt   "   ATT on common support 1.18076 (full 1.17919); ATU on common support -.34765 (full -.34353)"
return list

* ---- esr_suff, two-step, two instruments -------------------------------------------------
use "data/esr_suff.dta", clear
esreg y x, select(d = x z1 z2) method(twostep) nolog
esrdiag
di as txt _n "expected (Python): LR chi2(2) = 1775.85; pseudo-R2 .27754 / .01872; z1 chi2 924.39, z2 chi2 692.28"
di as txt   "   Var(P|X)/Var(P) .92373; VIF 1.04563 / 1.06642; shares .02468 / .02813; ATU on support -.20274 (full -.2016)"

* ---- esr_annex2 with a factor variable in the selection equation ------------------------------
use "data/esr_annex2.dta", clear
esreg income educ, select(treatment = educ i.region inst) nolog
esrdiag
di as txt _n "expected (Python, region dummies): LR chi2(3) = 1126.64; pseudo-R2 .17688 / .00001;"
di as txt   "   2.region chi2 40.64, 3.region 57.05, inst 887.55; Var(P|X)/Var(P) .99979; VIF 1.00027 / 1.00001;"
di as txt   "   shares .03333 / .0012; ATT on support 1.98114"

* ---- a weak instrument: the warnings ------------------------------------------------------------
use "data/esr_sim.dta", clear
set seed 2026
gen double noise = rnormal()
esreg y x, select(d = x noise) method(twostep) nolog
esrdiag
di as txt _n "expected: LR chi2(1) close to 0, VIF of lambda in the tens or more, both warnings printed"
di as txt "e(cmd) after esrdiag: `e(cmd)'   e(att) = " e(att) "   (the user's estimation is kept)"
