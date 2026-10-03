* esrtest_test.do -- esrtest, kappa (constant selection on gains) against the
* Python reference on data/esr_kx.dta (design K: kappa(x) = 1.08 + 0.4 x, n = 5000).
* Run from the project root: tools\run_stata.bat tests\esrtest_test.do
version 16
clear all
set more off
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrtest {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}

use "data/esr_kx.dta", clear
summarize kappa_true

* ---- two-step with kappa(x) --------------------------------------------------------
esreg y x, select(d = x z) method(twostep) kappa(x)
di as txt _n "expected (Python): y_1: x 1.16936 lambda_x .368271 (.0700732) lambda .775691 (.0645897) _cons 1.00769"
di as txt   "                   y_0: x .384036 lambda_x -.0639327 (.0495426) lambda -.266896 (.0503117) _cons .517259"
di as txt   "                   ATT 1.13062 (.0724412) ATU -.324953 (.0856042) kappa 1.0394 (.0821178)"
esrtest, kappa
di as txt _n "expected (Python): slope .43220 (.08598), chi2(1) = 25.27, p = 5.0e-07"
predict kx, kappa
summarize kx kappa_true
correlate kx kappa_true

* ---- FIML with hetsigma(x) hetrho(x) ---------------------------------------------------
esreg y x, select(d = x z) hetsigma(x) hetrho(x) nolog
di as txt _n "expected (Python): ll -10151.0278; lnsigma_1: x .161898 (.0141098); lnsigma_0: x .0131259 (.0169904)"
di as txt   "                   atanhrho_1: x .262743 (.0494046); atanhrho_0: x -.0870503 (.0585534)"
di as txt   "                   ATT 1.14599 (.0734285) ATU -.390395 (.0767194) kappa 1.0979 (.0760317)"
esrtest, kappa
di as txt _n "expected (Python): chi2(4) = 133.94"
predict kx2, kappa
correlate kx kx2 kappa_true

* ---- constant-kappa estimations refuse the test with a clear message ----------------------
esreg y x, select(d = x z) method(twostep)
cap noisily esrtest, kappa
esreg y x, select(d = x z) nolog
cap noisily esrtest, kappa

* =====================================================================================
* esrtest, suff  (index sufficiency, Section 4.2) against the Python reference
* =====================================================================================
foreach f in esrtest {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}
* ---- sufficiency holds (two instruments, z1 kept, z2 tested) -------------------------
use "data/esr_suff.dta", clear
esreg y x, select(d = x z1 z2) method(twostep) nolog
esrtest, suff
di as txt _n "expected (Python): keep z1, tested z2; coef .00189 (.02982) regime 1, .02837 (.02644) regime 0;"
di as txt   "                   joint chi2(2) = 1.16 (p .560), regime 1 chi2(1) = 0.00, regime 0 chi2(1) = 1.15"
return list
* the choice of the kept instrument can be forced
esrtest, suff keep(z2)
* ---- double hurdle: a second index in participation ---------------------------------------
use "data/esr_dh.dta", clear
esreg y x, select(d = x z1 z2) method(twostep) nolog
esrtest, suff
di as txt _n "expected (Python): coef .28318 (.04084) regime 1, .05318 (.02137) regime 0;"
di as txt   "                   joint chi2(2) = 54.98, regime 1 chi2(1) = 48.09, regime 0 chi2(1) = 6.19"
* ---- single instrument (esr_sim): the instrument tests itself, warning printed -------------
use "data/esr_sim.dta", clear
esreg y x, select(d = x z) nolog
esrtest, suff
di as txt _n "expected (Python): coef -.04668 (.04301) regime 1, -.00414 (.04509) regime 0; joint chi2(2) = 1.19 (p .552)"
di as txt "e(cmd) after esrtest: `e(cmd)'   e(method) = `e(method)'   (expected esreg, fiml: the user's estimation is kept)"
estimates dir

* =====================================================================================
* esrtest, normal  (normality by regime, Section 4.5) against the Python reference
* =====================================================================================
foreach f in esreg esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esrtest {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}
use "data/esr_sim.dta", clear
esreg y x, select(d = x z) nolog
esrtest, normal
di as txt _n "expected (Python): Hausman chi2(9) = 11.00 (p .276); kappa 1.044 vs 1.089, diff -.046 (s.e. .018)"
di as txt   "                   Hermite: coefs .030 -.082 (regime 1) .011 -.013 (regime 0); joint chi2(4) = 9.10, regime 1 7.79, regime 0 1.31"
return list
use "data/esr_skew.dta", clear
esreg y x, select(d = x z) method(twostep) nolog
esrtest, normal
di as txt _n "expected (Python): Hausman chi2(9) = 397.70; kappa 1.167 vs 2.825, diff -1.658 (s.e. .107); Hermite chi2(4) = 3.33 (p .50)"
use "data/esr_nonlin.dta", clear
esreg y x, select(d = x z) method(twostep) nolog
esrtest, normal
di as txt _n "expected (Python): Hausman chi2(9) = 8.30 (p .50); kappa .878 vs 1.233, diff -.354 (s.e. .619)"
di as txt   "                   Hermite: coefs -.45 .881 (regime 1) .174 .044 (regime 0); joint chi2(4) = 28.57, regime 1 24.99, regime 0 4.32"
di as txt "e(cmd) after esrtest: `e(cmd)'  e(method) = `e(method)'  e(att) = " e(att) "  (expected esreg, twostep, 1.08803)"
estimates dir

* =====================================================================================
* esrtest, pdid  (pseudo-DiD on the selection index, Section 4.6) against Python
* =====================================================================================
foreach f in esrtest {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}
use "data/esr_sim.dta", clear
esreg y x, select(d = x z) method(twostep) nolog
esrtest, pdid
di as txt _n "expected (Python, nq=2): gaps .9527 (.0253) and .4829 (.0283); C1 -.4847 (.0247) C0 .2453 (.0248);"
di as txt   "   pi1 -.6596 pi0 .7766; rho1sigma1 .7348 rho0sigma0 -.3158; kappa_dd 1.0506 (.0501); DD -.4698 (.0376)"
esrtest, pdid nq(4)
di as txt _n "expected (Python, nq=4): gaps 1.1356 .8667 .6288 .3240; C1 -.8655 (.0425) C0 .4990 (.0466);"
di as txt   "   pi1 -1.1815 pi0 1.3540; kappa_dd 1.1012 (.0509); DD -.8116 (.0714)"
use "data/esr_annex2.dta", clear
esreg income educ, select(treatment = educ i.region inst) nolog
esrtest, pdid
di as txt _n "expected (Python, nq=2): gaps 2.2578 2.2983; C1 -.0809 (.0122) C0 -.1212 (.0195); pi1 -.4867 pi0 .7304;"
di as txt   "   rho1sigma1 .1662 rho0sigma0 .1659; kappa_dd .0003 (.0361); DD .0405 (.0229)"
use "data/esr_skew.dta", clear
esreg y x, select(d = x z) method(twostep) nolog
esrtest, pdid nq(4)
di as txt _n "expected (Python, nq=4): C1 -.9235 (.0744) C0 .5246 (.0710); pi1 -1.1261 pi0 1.3166; kappa_dd 1.2185 (.0864)"
return list

* =====================================================================================
* esrtest, spec  (specification of the selection equation, Section 4.1) against Python
* =====================================================================================
foreach f in esreg esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esrtest {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}
* ---- omitted quadratic in the index: D = 1{.2+.3x+.9z+.5(z^2-1)+u>0}, working probit on (x z)
use "data/esr_link.dta", clear
esreg y x, select(d = x z) nolog
esrtest, spec
di as txt _n "expected (Python): link v^2 .84348 (.04807), Wald chi2(1) = 307.93;"
di as txt   "                   gamma probit .24969 .6053 .05555 vs FIML .24864 .59472 .05145, se diff .0021 .00273 .00104, chi2(3) = 31.20 (p .000)"
return list
esrtest, spec order(3)
di as txt _n "expected (Python): v^2 .84487 (.06508) v^3 .00182 (.05749), Wald chi2(2) = 307.79"
esrtest, spec add(z2)
di as txt _n "expected (Python): z2 .49394 (.02288), Wald chi2(1) = 466.18, LR chi2(1) = 593.70"
esrtest, spec add(c.z#c.z)
di as txt _n "expected: same as add(z2)"
* ---- the right index: the quadratic in the selection equation
esreg y x, select(d = x z z2) nolog
esrtest, spec
di as txt _n "expected (Python): link v^2 -.03028 (.03103), Wald chi2(1) = 0.95 (p .329);"
di as txt   "                   gamma probit .28854 .90247 .49394 -.29689 vs FIML .2891 .90496 .50293 -.30105, chi2(4) = 2.64 (p .620)"
di as txt   "                   kappa: select(x z) .9779 (fiml) .9322 (twostep); select(x z z2) 1.1567 / 1.1301  (true 1.08)"
esreg y x, select(d = x z z2) method(twostep) nolog
esrtest, spec
di as txt "e(cmd) after esrtest: `e(cmd)'  e(method) = `e(method)'  (expected esreg, twostep: the user's estimation is kept)"
* ---- esr_sim: index correctly specified, joint normal
use "data/esr_sim.dta", clear
esreg y x, select(d = x z) nolog
esrtest, spec
di as txt _n "expected (Python): link v^2 -.00297 (.01334), chi2(1) = 0.05 (p .824); gamma chi2(3) = 6.66 (p .084)"
esrtest, spec order(3)
di as txt _n "expected (Python): v^2 -.00287 (.01385) v^3 -.00028 (.0106), chi2(2) = 0.05 (p .975)"
* ---- factor terms in the selection equation (esr_annex2): spec, suff and normal
* (the internal re-estimations now run on the plain variables behind 2.region 3.region)
use "data/esr_annex2.dta", clear
esreg income educ, select(treatment = educ i.region inst) nolog
esrtest, spec
di as txt _n "expected (Python): link v^2 -.01928 (.03843), chi2(1) = 0.25 (p .616);"
di as txt   "                   gamma probit -.0057 .3258 .361 .6972 .3007 vs FIML -.0034 .2872 .3313 .7008 .3196,"
di as txt   "                   se diff .0026 .0209 .0191 .004 .0148, chi2(5) = 7.17 (p .208)"
return list
esrtest, suff
di as txt _n "expected (Python, keep inst, tested 2.region 3.region): coefs .01577 (.01602) .01192 (.01489) regime 1,"
di as txt   "                   .05245 (.0215) .03352 (.0205) regime 0; joint chi2(4) = 6.94, regime 1 1.06, regime 0 6.25"
esrtest, normal
di as txt _n "expected (Python): Hausman chi2(11) = 11.36 (p .413); kappa -.0055 vs -.0009, diff -.0046 (s.e. .0130);"
di as txt   "                   Hermite coefs .1022 -.014 (regime 1) .1312 .0447 (regime 0); joint chi2(4) = 8.02, regime 1 6.21, regime 0 1.81"
di as txt "e(cmd) after esrtest: `e(cmd)'  e(method) = `e(method)'  (expected esreg, fiml)"
estimates dir
