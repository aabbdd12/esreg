*! _esreg_esample.ado 1.0.0  03oct2026  A. Araar
*! Internal to the esreg family: the observations the current esreg estimation
*! was fitted on, and how its variance aggregates them.
*!
*!   _esreg_esample newvarname
*!
*! creates newvarname (byte): 1 on e(sample), within the subpopulation after
*! svy, subpop() (svy then reports the whole design sample in e(sample)).
*! Returns r(mode) = iid | cluster | svy, r(clustvar), and r(domain) = 1 after
*! svy, subpop() (the design variance is then that of a domain).

cap program drop _esreg_esample
program define _esreg_esample, rclass
    version 16
    syntax newvarname
    qui gen byte `varlist' = e(sample)
    local mode "iid"
    if ("`e(prefix)'" == "svy" | "`e(vce)'" == "svy") local mode "svy"
    else if ("`e(clustvar)'" != "") local mode "cluster"
    local domain = 0
    if ("`e(prefix)'" == "svy" & `"`e(subpop)'"' != "") {
        * e(subpop) is "varname", "varname if exp" or "if exp"
        local domain = 1
        local spec `"`e(subpop)'"'
        gettoken first rest : spec
        if ("`first'" == "if") {
            qui replace `varlist' = 0 if !(`rest')
        }
        else {
            qui replace `varlist' = 0 if `first' == 0 | missing(`first')
            gettoken iff rest2 : rest
            if ("`iff'" == "if") qui replace `varlist' = 0 if !(`rest2')
        }
    }
    return local mode "`mode'"
    return local clustvar "`e(clustvar)'"
    return scalar domain = `domain'
end
