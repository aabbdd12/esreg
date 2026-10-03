* esreg_wt_test.do -- weights: pweights, svyset auto-detection, hsize(), the svy prefix,
* predict scores; against the Python references of export_wt.py (data/esr_wt.dta)
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrdiag esrtest {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}
use "data/esr_wt.dta", clear
* ---- 1. explicit pweights, FIML (robust) and two-step ----------------------------------
esreg income educ [pw = wt], select(treatment = educ i.region inst) nolog
di as txt _n "expected (Python, robust): educ_1 .506362 (.005588)  cons_1 62.878491 (.020202)  inst .698443 (.025182)"
di as txt   "   atanhrho_1 .595651 (.070358)  atanhrho_0 .554330 (.066949); ATT 1.97979 (.03102) ATU 1.97070 (.02977) kappa .00675 (.02884)"
esreg income educ [pw = wt], select(treatment = educ i.region inst) method(twostep) nolog
di as txt _n "expected (Python): inst .692478 (.025678)  lambda1 .191507 (.024015)  -lambda0 .183995 (.023050);"
di as txt   "   ATT 1.97998 (.03326) ATU 1.96986 (.03434) kappa .00751 (.03220)"
* ---- 2. svyset auto-detection and hsize() -----------------------------------------------
svyset [pw = wt]
esreg income educ, select(treatment = educ i.region inst) nolog
di as txt _n "expected: same as the explicit [pw = wt] estimation above (note printed)"
esreg income educ, select(treatment = educ i.region inst) hsize(hsize) nolog
di as txt _n "expected (Python, weight wt*hsize): educ_1 .507171 (.006108)  inst .721402 (.028689)  atanhrho_1 .553533 (.077230);"
di as txt   "   ATT 1.97828 (.03209) ATU 1.98353 (.03246) kappa -.00408 (.03101)"
esreg income educ, select(treatment = educ i.region inst) hsize(hsize) method(twostep) nolog
di as txt _n "expected (Python): inst .717161 (.028904); ATT 1.97093 (.03499) ATU 1.98633 (.03500) kappa -.01176 (.03375)"
esreg income educ, select(treatment = educ i.region inst) nosvyset nolog
di as txt _n "expected: unweighted estimation (ll -4229.3298, kappa -.000888)"
* ---- 3a. the equation-level scores, called directly (the error, if any, is shown) ----------
esreg income educ, select(treatment = educ i.region inst) nolog
predict double sc*, scores
summarize sc* [aw = wt]
di as txt "expected: seven scores, weighted means ~ 0 (first-order conditions)"
drop sc*
* ---- 3b. the svy prefix (FIML): linearized = sandwich * sqrt(n/(n-1)) under pweights only
svy: esreg income educ, select(treatment = educ i.region inst)
di as txt _n "expected: coefficients as in 1., s.e. = Python robust * 1.0001 (educ_1 .005589, inst .025185)"
di as txt "e(prefix) = `e(prefix)'  e(vce) = `e(vce)'  e(eff_vce) = `e(eff_vce)'"
esreg
di as txt _n "expected: effects recomputed from the linearized e(V): ATT 1.97979 (~.03102) kappa .00675 (~.02884); e(eff_vce) = `e(eff_vce)'"
esreg, effects
estimates dir
* ---- 4. the post-commands after a weighted estimation ----------------------------------
esreg income educ [pw = wt], select(treatment = educ i.region inst) nolog
esrdiag
esrtest, spec
esrtest, pdid
* ---- 5. svy with the two-step: the prefix is refused (no scores), vce(svy) instead -------
cap noi svy: esreg income educ, select(treatment = educ i.region inst) method(twostep)
di as txt "expected: an error (scores are available after method(fiml) only; use vce(svy))"
esreg income educ, select(treatment = educ i.region inst) method(twostep) vce(svy)
di as txt "expected: pweight-only design: s.e. = the two-step ones of 1. times sqrt(n/(n-1)) = 1.0001"
