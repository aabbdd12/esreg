*! _esreg_data.ado 0.4.4  12sep2026  A. Araar
*! Internal to the esreg family: plain variables (__esr_fv#) behind the factor terms
*! of e(xvars), e(zvars), e(hetsigma), e(hetrho), e(kappavars), so that probit,
*! regress and the internal esreg calls do not re-merge 2.region 3.region.

cap program drop _esreg_data
* the plain variable lists of the current esreg estimation, for Mata
* returns r(x) r(z) r(hs) r(hr) r(kap); variables __esr_fv* are created as needed
program define _esreg_data, rclass
    version 16
    syntax [if]
    cap drop __esr_fv*
    local j = 0
    local made_terms ""
    local made_vars ""
    foreach o in xvars zvars hetsigma hetrho kappavars {
        local out ""
        foreach v in `e(`o')' {
            cap confirm variable `v'
            if (_rc == 0) {
                local out "`out' `v'"
                continue
            }
            * the same factor term gets the same plain variable in every list
            local pos : list posof "`v'" in made_terms
            if (`pos' > 0) {
                local out "`out' `: word `pos' of `made_vars''"
                continue
            }
            local ++j
            local made_terms "`made_terms' `v'"
            local made_vars  "`made_vars' __esr_fv`j'"
            fvrevar `v' `if'
            local r "`r(varlist)'"
            qui summarize `r' `if', meanonly
            if (r(min) == r(max)) {
                local lev = substr("`v'", 1, strpos("`v'", ".") - 1)
                local var = substr("`v'", strpos("`v'", ".") + 1, .)
                qui gen double __esr_fv`j' = (`var' == `lev') `if'
            }
            else qui gen double __esr_fv`j' = `r' `if'
            local out "`out' __esr_fv`j'"
        }
        local `o' "`out'"
    }
    return local x   "`xvars'"
    return local z   "`zvars'"
    return local hs  "`hetsigma'"
    return local hr  "`hetrho'"
    return local kap "`kappavars'"
end
