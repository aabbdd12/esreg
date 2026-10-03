* etregress_check.do -- esreg against Stata's own etregress (Stata 19, [CAUSAL] etregress)
* With poutcomes and the treatment interacted with every outcome covariate, etregress fits
* exactly the unconstrained ESR (separate beta_j, sigma_j, rho_j): the log likelihoods and
* the coefficients must agree to every digit, and margins, predict(cte) subpop() must give
* esreg's ATT.  cfunction (one-step GMM on stacked moments) is a second reference for the
* stacked-moment variance of the two-step route.
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrdiag esrtest {
    cap confirm file "`f'.ado"
    if (_rc == 0) run "`f'.ado"
    else {
        cap confirm file "src/`f'.ado"
        if (_rc == 0) run "src/`f'.ado"
        else          run "../src/`f'.ado"
    }
}
* ---------------------------------------------------------------- drugexp (manual, example 2 and 4)
use "data/drugexp.dta", clear
etregress lndrug i.ins#c.(chron age lninc) chron age lninc, treat(ins = age married lninc work) poutcomes nolog
estimates store et_ml
* ATT (cte on the treated) and ATE; vce(unconditional) needs vce(robust), so the delta-method
* variance here is the "parameter" component of esreg's, the sampling component apart
margins, predict(cte) subpop(ins)
margins r.ins, contrast(nowald)
di as txt "expected from esreg / Python: ll = -12465.1123; ATT -.71550 ATU -.90863 ATE -.76198 kappa .11797"
esreg lndrug chron age lninc, select(ins = age married lninc work) nolog
di as txt "esreg FIML: ll = " e(ll) "   ATT = " e(att) "  ATU = " e(atu) "  kappa = " e(kappa)
etregress lndrug i.ins#c.(chron age lninc) chron age lninc, treat(ins = age married lninc work) poutcomes cfunction nolog
margins, predict(cte) subpop(ins) vce(unconditional)
esreg lndrug chron age lninc, select(ins = age married lninc work) method(twostep) nolog
di as txt "expected (Python two-step): probit .021142 .084631 .102303 .288418 -.622993 (= etregress's ins equation);"
di as txt "   ATT -.63055 (.53236) ATU -1.2645 (.45150) kappa .37954 (.39468); sigma1 1.21642 sigma0 1.35122"
* ---------------------------------------------------------------- union3 (manual, examples 1, 3, 5)
use "data/union3.dta", clear
etregress wage i.union#c.(age grade smsa black tenure) age grade smsa black tenure, ///
    treat(union = south black tenure) poutcomes nolog
margins, predict(cte) subpop(union)
di as txt "expected from esreg / Python: ll = -3008.4210; ATT 4.1863 ATU 2.3262 kappa 1.1431; rho0 = -.9039"
esreg wage age grade smsa black tenure, select(union = south black tenure) nolog
esrdiag
esrtest, spec
esrtest, normal
di as txt "expected (Python): Hausman chi2(18) = 55.84; gamma contrast chi2(4) = 22.79; link 1.26; Hermite 3.66"
esreg ln_wage age grade smsa black tenure, select(union = south black tenure) nolog
esrtest, normal
di as txt "expected (Python, ln_wage): Hausman 29.72 (p .040); kappa 2S .2102 vs FIML -.0571; Hermite 4.86"
esreg ln_wage age grade smsa black tenure, select(union = south black tenure) method(twostep) nolog
esrtest, pdid nq(4)
di as txt "expected (Python, xtile strata): kappa_dd .7820 (.2989); C1 .0839 (.0632) C0 .2052 (.0442); pi1 -.4033 pi0 .2073"
