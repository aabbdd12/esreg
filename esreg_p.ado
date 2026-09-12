*! esreg_p 0.4.1  11sep2026  A. Araar
*! predict after esreg.  One statistic per call (default: effect).
*!   xb1 xb0        linear predictions X b_1, X b_0
*!   xbsel          selection index Z g            pr         P(D = 1 | Z)
*!   lambda1        phi(Zg)/Phi(Zg)                lambda0    phi(Zg)/(1 - Phi(Zg))
*!   yc11 yc10      E[Y_1 | D = 1, X], E[Y_1 | D = 0, X]
*!   yc01 yc00      E[Y_0 | D = 1, X], E[Y_0 | D = 0, X]
*!   effect         expected individual effect: yc11 - yc01 if treated, yc10 - yc00 if not
*!   kappa          kappa_i = rho1 sigma1(x) - rho0 sigma0(x)   (constant unless heterogeneous)
*!   rhosig1 rhosig0  rho_j sigma_j(x)
*!   sigma1 sigma0 rho1 rho0   (fiml only) sigma_j(x), rho_j(x)
*! Correspondence with mspredict after movestay: xb1/xb2 = xb1/xb0, psel = pr,
*! mills1/mills2 = lambda1/lambda0, yc1_1 = yc11, yc2_1 = yc01 (treated only),
*! yc1_2 = yc10, yc2_2 = yc00 (untreated only); predict computes all four for everyone.
* all program drops first: an autoloaded ado-file defines its programs only while
* it reads program definitions, so nothing else may stand between them
cap program drop esreg_p
cap program drop _esreg_p_lin
cap program drop _esreg_p_rs
cap program drop _esreg_p_scores
program define esreg_p
    version 16
    if ("`e(cmd)'" != "esreg") {
        di as err "esreg estimates not found"
        exit 301
    }
    * ---- equation-level scores (svy prefix): predict stub* [if] [in], scores ----------
    local 0orig `0'
    syntax [anything] [if] [in] [, SCores *]
    if ("`scores'" != "") {
        if ("`e(method)'" != "fiml") {
            di as err "scores are available after method(fiml) only (the svy prefix needs them)"
            exit 198
        }
        _score_spec `0orig'
        local vl "`s(varlist)'"
        local vt : word 1 of `s(typlist)'
        marksample touse, novarlist
        _esreg_p_scores `vl' if `touse', vtyp(`vt')
        exit
    }
    local 0 `0orig'
    syntax newvarname [if] [in] [, XB1 XB0 XBSel PR LAMbda1 LAMbda0 YC11 YC10 YC01 YC00 ///
                                   EFFect KAPpa RHOSig1 RHOSig0 SIGma1 SIGma0 RHO1 RHO0 ]
    local stat "`xb1'`xb0'`xbsel'`pr'`lambda1'`lambda0'`yc11'`yc10'`yc01'`yc00'`effect'`kappa'`rhosig1'`rhosig0'`sigma1'`sigma0'`rho1'`rho0'"
    local nst : word count `xb1' `xb0' `xbsel' `pr' `lambda1' `lambda0' `yc11' `yc10' `yc01' `yc00' `effect' `kappa' `rhosig1' `rhosig0' `sigma1' `sigma0' `rho1' `rho0'
    if (`nst' > 1) {
        di as err "only one statistic may be specified"
        exit 198
    }
    if ("`stat'" == "") {
        local stat effect
        di as txt "(option effect assumed)"
    }
    marksample touse, novarlist
    local y   "`e(depvar)'"
    local dv  "`e(treat)'"
    local method "`e(method)'"
    if (inlist("`stat'", "sigma1", "sigma0", "rho1", "rho0") & "`method'" != "fiml") {
        di as err "`stat' is available after method(fiml) only; use rhosig1/rhosig0 after method(twostep)"
        exit 198
    }

    * ---- building blocks ------------------------------------------------------
    tempvar xb1v xb0v zg l1 l0 rs1 rs0
    _esreg_p_lin `xb1v' if `touse', eq(`y'_1)
    _esreg_p_lin `xb0v' if `touse', eq(`y'_0)
    _esreg_p_lin `zg'   if `touse', eq(`dv')
    qui gen double `l1' = exp(lnnormalden(`zg') - lnnormal(`zg'))  if `touse'
    qui gen double `l0' = exp(lnnormalden(`zg') - lnnormal(-`zg')) if `touse'
    if ("`method'" == "fiml") {
        tempvar ls1 ls0 ar1 ar0
        _esreg_p_lin `ls1' if `touse', eq(lnsigma_1)
        _esreg_p_lin `ls0' if `touse', eq(lnsigma_0)
        _esreg_p_lin `ar1' if `touse', eq(atanhrho_1)
        _esreg_p_lin `ar0' if `touse', eq(atanhrho_0)
        qui gen double `rs1' = tanh(`ar1') * exp(`ls1') if `touse'
        qui gen double `rs0' = tanh(`ar0') * exp(`ls0') if `touse'
    }
    else {
        _esreg_p_rs `rs1' if `touse', eq(`y'_1)
        _esreg_p_rs `rs0' if `touse', eq(`y'_0)
    }

    * ---- the requested statistic ------------------------------------------------
    if ("`stat'" == "xb1")     qui gen `typlist' `varlist' = `xb1v' if `touse'
    if ("`stat'" == "xb0")     qui gen `typlist' `varlist' = `xb0v' if `touse'
    if ("`stat'" == "xbsel")   qui gen `typlist' `varlist' = `zg'   if `touse'
    if ("`stat'" == "pr")      qui gen `typlist' `varlist' = normal(`zg') if `touse'
    if ("`stat'" == "lambda1") qui gen `typlist' `varlist' = `l1' if `touse'
    if ("`stat'" == "lambda0") qui gen `typlist' `varlist' = `l0' if `touse'
    if ("`stat'" == "yc11")    qui gen `typlist' `varlist' = `xb1v' + `rs1' * `l1' if `touse'
    if ("`stat'" == "yc10")    qui gen `typlist' `varlist' = `xb1v' - `rs1' * `l0' if `touse'
    if ("`stat'" == "yc01")    qui gen `typlist' `varlist' = `xb0v' + `rs0' * `l1' if `touse'
    if ("`stat'" == "yc00")    qui gen `typlist' `varlist' = `xb0v' - `rs0' * `l0' if `touse'
    if ("`stat'" == "effect") {
        qui gen `typlist' `varlist' = cond(`dv' == 1, ///
                                   (`xb1v' - `xb0v') + (`rs1' - `rs0') * `l1', ///
                                   (`xb1v' - `xb0v') - (`rs1' - `rs0') * `l0') if `touse'
    }
    if ("`stat'" == "kappa")   qui gen `typlist' `varlist' = `rs1' - `rs0' if `touse'
    if ("`stat'" == "rhosig1") qui gen `typlist' `varlist' = `rs1' if `touse'
    if ("`stat'" == "rhosig0") qui gen `typlist' `varlist' = `rs0' if `touse'
    if ("`stat'" == "sigma1")  qui gen `typlist' `varlist' = exp(`ls1')  if `touse'
    if ("`stat'" == "sigma0")  qui gen `typlist' `varlist' = exp(`ls0')  if `touse'
    if ("`stat'" == "rho1")    qui gen `typlist' `varlist' = tanh(`ar1') if `touse'
    if ("`stat'" == "rho0")    qui gen `typlist' `varlist' = tanh(`ar0') if `touse'
    label variable `varlist' "esreg `stat'"
end

* linear index of one equation of e(b), skipping the lambda columns of the two-step
program define _esreg_p_lin
    version 16
    syntax newvarname [if], eq(string)
    marksample touse, novarlist
    tempname b sub
    mat `b' = e(b)
    local cn : colfullnames `b'
    local first = 1
    foreach c of local cn {
        gettoken ceq cname : c, parse(":")
        if ("`ceq'" != "`eq'") continue
        local cname = substr("`cname'", 2, .)
        if (substr("`cname'", 1, 6) == "lambda") continue
        if (`first') {
            mat `sub' = `b'[1, "`c'"]
            local first = 0
        }
        else mat `sub' = `sub', `b'[1, "`c'"]
    }
    qui matrix score double `varlist' = `sub' if `touse', equation(`eq')
end

* rho_j sigma_j (x) of the two-step route: coefficients of the lambda columns of
* equation eq, applied to e(kappavars) and the constant
program define _esreg_p_rs
    version 16
    syntax newvarname [if], eq(string)
    marksample touse, novarlist
    tempname b sub
    mat `b' = e(b)
    local cn : colfullnames `b'
    local kv "`e(kappavars)' _cons"
    local first = 1
    local j = 0
    foreach c of local cn {
        gettoken ceq cname : c, parse(":")
        if ("`ceq'" != "`eq'") continue
        local cname = substr("`cname'", 2, .)
        if (substr("`cname'", 1, 6) != "lambda") continue
        local ++j
        local newname : word `j' of `kv'
        tempname one
        mat `one' = `b'[1, "`c'"]
        mat colnames `one' = `newname'
        if (`first') {
            mat `sub' = `one'
            local first = 0
        }
        else mat `sub' = `sub', `one'
    }
    qui matrix score double `varlist' = `sub' if `touse'
end

* ---------------------------------------------------------------------------
* equation-level scores of the FIML, d lnL_i / d(x_i b_eq), in the order of e(b):
* y_1, y_0, <d>, lnsigma_1, lnsigma_0, atanhrho_1, atanhrho_0 (unweighted; the svy
* prefix applies the design)
program define _esreg_p_scores
    version 16
    syntax newvarlist(min=7 max=7) [if] [, vtyp(string)]
    marksample touse, novarlist
    local y  "`e(depvar)'"
    local dv "`e(treat)'"
    tempvar xb1 xb0 zg ls1 ls0 ar1 ar0 e1 e0 r1 r0 q1 q0 m1 m0
    _esreg_p_lin `xb1' if `touse', eq(`y'_1)
    _esreg_p_lin `xb0' if `touse', eq(`y'_0)
    _esreg_p_lin `zg'  if `touse', eq(`dv')
    _esreg_p_lin `ls1' if `touse', eq(lnsigma_1)
    _esreg_p_lin `ls0' if `touse', eq(lnsigma_0)
    _esreg_p_lin `ar1' if `touse', eq(atanhrho_1)
    _esreg_p_lin `ar0' if `touse', eq(atanhrho_0)
    qui gen double `e1' = (`y' - `xb1') / exp(`ls1') if `touse'
    qui gen double `e0' = (`y' - `xb0') / exp(`ls0') if `touse'
    qui gen double `r1' = tanh(`ar1') if `touse'
    qui gen double `r0' = tanh(`ar0') if `touse'
    qui gen double `q1' = sqrt(1 - `r1'^2) if `touse'
    qui gen double `q0' = sqrt(1 - `r0'^2) if `touse'
    qui gen double `m1' = exp(lnnormalden((`zg' + `r1'*`e1')/`q1') - lnnormal((`zg' + `r1'*`e1')/`q1'))   if `touse'
    qui gen double `m0' = exp(lnnormalden((`zg' + `r0'*`e0')/`q0') - lnnormal(-(`zg' + `r0'*`e0')/`q0')) if `touse'
    if ("`vtyp'" == "") local vtyp double
    tokenize `varlist'
    qui gen `vtyp' `1' = cond(`dv' == 1, (`e1' - `m1'*`r1'/`q1') / exp(`ls1'), 0) if `touse'
    qui gen `vtyp' `2' = cond(`dv' == 0, (`e0' + `m0'*`r0'/`q0') / exp(`ls0'), 0) if `touse'
    qui gen `vtyp' `3' = cond(`dv' == 1, `m1'/`q1', -`m0'/`q0') if `touse'
    qui gen `vtyp' `4' = cond(`dv' == 1, `e1'^2 - 1 - `m1'*`r1'*`e1'/`q1', 0) if `touse'
    qui gen `vtyp' `5' = cond(`dv' == 0, `e0'^2 - 1 + `m0'*`r0'*`e0'/`q0', 0) if `touse'
    qui gen `vtyp' `6' = cond(`dv' == 1, `m1'*(`e1' + `r1'*`zg')/`q1', 0) if `touse'
    qui gen `vtyp' `7' = cond(`dv' == 0, -`m0'*(`e0' + `r0'*`zg')/`q0', 0) if `touse'
    local i = 0
    foreach eq in `y'_1 `y'_0 `dv' lnsigma_1 lnsigma_0 atanhrho_1 atanhrho_0 {
        local ++i
        label variable ``i'' "equation-level score from esreg: `eq'"
    }
end
