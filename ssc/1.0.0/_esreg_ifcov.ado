*! _esreg_ifcov.ado 1.0.0  03oct2026  A. Araar
*! Internal to the esreg family: the covariance of the sums of influence functions
*! held in Stata variables (one variable per estimate, the weight of the
*! observation included), for independent observations, by cluster, or by the
*! survey design (svyset).
*!
*!   _esreg_ifcov varlist [if] [, CLuster(varname) SVY DOMain(varname) ]
*!
*!   default : sum_i u_i u_i'
*!   cluster : sum_g u_g u_g' * G/(G-1), u_g the sum within cluster g (as ml);
*!             r(N_clust)
*!   svy     : the linearized design variance of the totals of u_i / w_i, w the
*!             svyset sampling weight (svy linearized: total), on the observations
*!             of the if condition; domain() is the subpopulation (svy, subpop()),
*!             the influence functions being zero outside it; r(N_strata),
*!             r(N_psu), r(df_r)
*! r(V) is k x k, in the order of varlist.  A missing value of an influence
*! function counts as zero (it is zero outside the estimation sample).

cap program drop _esreg_ifcov
program define _esreg_ifcov, rclass
    version 16
    syntax varlist(numeric) [if] [, CLuster(varname) SVY DOMain(varname numeric) ]
    marksample smp, novarlist
    tempname C
    if ("`svy'" == "") {
        local cid ""
        if ("`cluster'" != "") {
            tempvar cid
            qui egen long `cid' = group(`cluster') if `smp'
            qui replace `smp' = 0 if missing(`cid')
        }
        local zl ""
        local k = 0
        foreach v of local varlist {
            local ++k
            tempvar z`k'
            qui gen double `z`k'' = cond(missing(`v'), 0, `v') if `smp'
            local zl "`zl' `z`k''"
        }
        mata: _esr_ifcov_m("`zl'", "`smp'", "`cid'", "`C'")
        mat rownames `C' = `varlist'
        mat colnames `C' = `varlist'
        if ("`cluster'" != "") return scalar N_clust = `ncl'
        return matrix V = `C'
        exit
    }
    * ---- survey design ----------------------------------------------------------
    qui svyset
    if ("`r(settings)'" == ", clear" | "`r(settings)'" == "") {
        di as err "the data are not svyset"
        exit 119
    }
    if (`"`r(poststrata)'`r(rake)'`r(regress)'"' != "") {
        di as err "the survey-design variance of esreg does not handle poststratification or calibration"
        exit 198
    }
    local wv "`r(wvar)'"
    local zl ""
    local k = 0
    foreach v of local varlist {
        local ++k
        tempvar z`k'
        if ("`wv'" != "") qui gen double `z`k'' = cond(missing(`v') | `wv' <= 0, 0, `v' / `wv') if `smp'
        else              qui gen double `z`k'' = cond(missing(`v'), 0, `v') if `smp'
        local zl "`zl' `z`k''"
    }
    tempname hold
    _estimates hold `hold', restore nullok
    if ("`domain'" != "") qui svy linearized, subpop(`domain'): total `zl' if `smp'
    else                  qui svy linearized: total `zl' if `smp'
    mat `C' = e(V)
    local ns  = e(N_strata)
    local np  = e(N_psu)
    local dfr = e(df_r)
    _estimates unhold `hold'
    mat rownames `C' = `varlist'
    mat colnames `C' = `varlist'
    return scalar N_strata = `ns'
    return scalar N_psu    = `np'
    return scalar df_r     = `dfr'
    return matrix V = `C'
end
