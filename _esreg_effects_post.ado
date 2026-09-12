*! _esreg_effects_post.ado 0.4.4  12sep2026  A. Araar
*! Internal to the esreg family: recomputes the effects, kappa, the support and the
*! ancillary results from the current e(b) and e(V) and reposts them (used after the
*! svy prefix has replaced e(V) by the linearized variance, and by -esreg, effects-).

cap program drop _esreg_effects_post
* ---------------------------------------------------------------------------
* recompute the effects, kappa, support and ancillary results from the current
* e(b) and e(V): after the svy prefix has replaced e(V) by its linearized
* variance, or on request (esreg, effects); reposts them in e()
program define _esreg_effects_post, eclass
    version 16
    tempvar smp w
    qui gen byte `smp' = e(sample)
    if ("`e(wtype)'" != "") {
        cap qui gen double `w' `e(wexp)' if `smp'
        if (_rc) {
            * under svy the weight passed to esreg was a temporary variable: use the design weight
            cap qui svyset
            if ("`r(wvar)'" != "") qui gen double `w' = `r(wvar)' if `smp'
            else                   qui gen double `w' = 1 if `smp'
        }
    }
    else qui gen double `w' = 1 if `smp'
    _esreg_data if `smp'
    local xl  "`r(x)'"
    local zl  "`r(z)'"
    local hsl "`r(hs)'"
    local hrl "`r(hr)'"
    local kl  "`r(kap)'"
    tempname E A S L
    mata: _esreg_effects("`e(depvar)'", "`xl'", "`zl'", "`e(treat)'", "`hsl'", "`hrl'", "`kl'", "`w'", "`smp'", "`E'")
    cap drop __esr_fv*
    mat rownames `E' = ATT ATU ATE kappa
    mat colnames `E' = est se var_param var_samp
    ereturn matrix effects = `E'
    foreach s in att atu ate kappa se_att se_atu se_ate se_kappa ///
                 sigma1 sigma0 rho1 rho0 supp_lo supp_hi ml1 ml0 ///
                 p_min1 p_max1 p_min0 p_max0 {
        ereturn scalar `s' = r(`s')
    }
    ereturn scalar rhosig1 = r(rho1) * r(sigma1)
    ereturn scalar rhosig0 = r(rho0) * r(sigma0)
    mat `A' = (r(sigma1), r(sigma0), r(rho1), r(rho0))
    mat colnames `A' = sigma1 sigma0 rho1 rho0
    ereturn matrix anc = `A'
    mat `S' = (r(p_min1), r(p_max1), r(p_min0), r(p_max0))
    mat colnames `S' = p_min1 p_max1 p_min0 p_max0
    ereturn matrix support = `S'
    mat `L' = (r(ml1), r(ml0))
    mat colnames `L' = lambda1_treated lambda0_untreated
    ereturn matrix lambda = `L'
    local ev = cond("`e(prefix)'" == "svy", "svy", "model")
    ereturn local eff_vce "`ev'"
end
