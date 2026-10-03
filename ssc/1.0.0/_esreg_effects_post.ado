*! _esreg_effects_post.ado 1.0.0  03oct2026  A. Araar
*! Internal to the esreg family: computes the effects, kappa, the support and the
*! ancillary results from the current e(b) and e(V) and posts them in e(): at the
*! end of esreg, after the svy prefix has replaced e(V) by the linearized
*! variance, and on request (esreg, effects).
*! The standard error of each effect adds three terms of its influence function:
*! the parameter part (delta method on e(V), so the vce of the estimation), the
*! sampling part (the average over the units) and twice their covariance; the
*! last two are aggregated as the estimation was: independent observations, by
*! cluster (e(clustvar)), or by the survey design (svy prefix, or vce(svy) of the
*! two-step route).  After svy, subpop(), the effects are those of the
*! subpopulation and the design variance is that of a domain.
*! Options (esreg itself): plain xl() zl() hsl() hrl() kl() wvar() give the plain
*! variables and the weight of the estimation, so that _esreg_data (which drops
*! and recreates __esr_fv*, possibly in use by the caller) is not called.

cap program drop _esreg_effects_post
program define _esreg_effects_post, eclass
    version 16
    syntax [, PLAIN XL(string) ZL(string) HSL(string) HRL(string) KL(string) WVar(varname) ]
    tempvar smp dom w
    qui gen byte `smp' = e(sample)
    * the fitted observations (the subpopulation after svy, subpop()) and how the
    * estimation aggregated them
    _esreg_esample `dom'
    local mode "`r(mode)'"
    local subpop = cond(r(domain), "domain(`dom')", "")
    if ("`wvar'" != "") qui gen double `w' = `wvar' if `smp'
    else if ("`e(wtype)'" != "") {
        cap qui gen double `w' `e(wexp)' if `smp'
        if (_rc) {
            * under svy the weight passed to esreg was a temporary variable: use the design weight
            cap qui svyset
            if ("`r(wvar)'" != "") qui gen double `w' = `r(wvar)' if `smp'
            else                   qui gen double `w' = 1 if `smp'
        }
    }
    else qui gen double `w' = 1 if `smp'
    if ("`plain'" == "") {
        _esreg_data if `dom'
        local xl  "`r(x)'"
        local zl  "`r(z)'"
        local hsl "`r(hs)'"
        local hrl "`r(hr)'"
        local kl  "`r(kap)'"
    }
    tempname E A S L U C
    * the influence functions of the four estimates: sampling part (4), parameter part (4)
    local ifv ""
    forvalues j = 1/8 {
        tempvar f`j'
        local ifv "`ifv' `f`j''"
    }
    mata: _esreg_effects("`e(depvar)'", "`xl'", "`zl'", "`e(treat)'", "`hsl'", "`hrl'", "`kl'", "`w'", "`dom'", "`E'", "`ifv'")
    if ("`plain'" == "") cap drop __esr_fv*
    * keep the r() of Mata before any other rclass call
    mat `A' = (r(sigma1), r(sigma0), r(rho1), r(rho0))
    mat `S' = (r(p_min1), r(p_max1), r(p_min0), r(p_max0))
    mat `L' = (r(ml1), r(ml0))
    mat `U' = (r(supp_lo), r(supp_hi))
    if ("`mode'" == "cluster") _esreg_ifcov `ifv' if `smp', cluster(`e(clustvar)')
    if ("`mode'" == "svy")     _esreg_ifcov `ifv' if `smp', svy `subpop'
    if ("`mode'" != "iid") {
        mat `C' = r(V)
        forvalues j = 1/4 {
            mat `E'[`j', 4] = `C'[`j', `j']
            mat `E'[`j', 5] = 2 * `C'[`j', `j' + 4]
            mat `E'[`j', 2] = sqrt(`E'[`j', 3] + `E'[`j', 4] + `E'[`j', 5])
        }
    }
    mat rownames `E' = ATT ATU ATE kappa
    mat colnames `E' = est se var_param var_samp cov_ps
    ereturn matrix effects = `E', copy
    local j = 0
    foreach s in att atu ate kappa {
        local ++j
        ereturn scalar `s'    = `E'[`j', 1]
        ereturn scalar se_`s' = `E'[`j', 2]
    }
    ereturn scalar sigma1  = `A'[1,1]
    ereturn scalar sigma0  = `A'[1,2]
    ereturn scalar rho1    = `A'[1,3]
    ereturn scalar rho0    = `A'[1,4]
    ereturn scalar rhosig1 = `A'[1,3] * `A'[1,1]
    ereturn scalar rhosig0 = `A'[1,4] * `A'[1,2]
    ereturn scalar supp_lo = `U'[1,1]
    ereturn scalar supp_hi = `U'[1,2]
    ereturn scalar ml1     = `L'[1,1]
    ereturn scalar ml0     = `L'[1,2]
    ereturn scalar p_min1  = `S'[1,1]
    ereturn scalar p_max1  = `S'[1,2]
    ereturn scalar p_min0  = `S'[1,3]
    ereturn scalar p_max0  = `S'[1,4]
    mat colnames `A' = sigma1 sigma0 rho1 rho0
    ereturn matrix anc = `A'
    mat colnames `S' = p_min1 p_max1 p_min0 p_max0
    ereturn matrix support = `S'
    mat colnames `L' = lambda1_treated lambda0_untreated
    ereturn matrix lambda = `L'
    ereturn local eff_vce "`mode'"
end
