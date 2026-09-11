*! esrmte 0.3.5  11sep2026  A. Araar
*! Marginal treatment effect after esreg.
*!   Parametric line  MTE(u) = m + kappa * invnormal(1 - u), m = E[X(b1 - b0)],
*!   at chosen percentiles u of the participation unobservable (0 = most eager),
*!   with delta-method standard errors and common-support flags;
*!   with semipar, the semiparametric MTE by percentile-weights regression on
*!   the probit score (engine _esreg_pwr, derivative-aware bandwidth) at the
*!   same percentiles, next to the line.
*!
*! esrmte [if] [in] [, AT(numlist) SEMIpar EST(name) GRaph Level(#)
*!                     NAme(string) SAVing(string) ]
*! Default at(): .05 .1 .25 .5 .75 .9 .95.  r(mte) = (u mte se support),
*! r(mte_sp) = (tau p_tau parametric pwr se h neff width_p), r(m), r(kappa).
cap program drop esrmte
program define esrmte, rclass
    version 16
    syntax [if] [in] [, AT(numlist >0 <1 sort) SEMIpar EST(name) GRaph Level(cilevel) ///
                        NAme(string) SAVing(string asis) ]
    esreg_engine
    * keep the user's e() whatever happens below (probit and _esreg_pwr post their own)
    tempname hold
    _estimates hold `hold', restore copy
    _esreg_getest `est'
    if ("`at'" == "") local at ".05 .1 .25 .5 .75 .9 .95"

    marksample touse, novarlist
    qui replace `touse' = 0 if !e(sample)
    local y      "`e(depvar)'"
    local dv     "`e(treat)'"
    local zv     "`e(zvars)'"
    local method "`e(method)'"
    local wtype  "`e(wtype)'"
    local wexp0  "`e(wexp)'"
    local ate    = e(ate)
    local plo    = e(supp_lo)
    local phi    = e(supp_hi)
    local z = invnormal(1 - (1 - `level'/100)/2)

    tempvar w
    if ("`wtype'" != "") qui gen double `w' `wexp0' if `touse'
    else                 qui gen double `w' = 1 if `touse'

    * ---- m, kappa and their covariance --------------------------------------------
    _esreg_data if `touse'
    local xl  "`r(x)'"
    local zl  "`r(z)'"
    local hs  "`r(hs)'"
    local hr  "`r(hr)'"
    local kap "`r(kap)'"
    tempname MK
    mata: _esreg_mte("`xl'", "`zl'", "`dv'", "`hs'", "`hr'", "`kap'", "`w'", "`touse'", "`MK'")
    local m   = `MK'[1,1]
    local kp  = `MK'[1,2]
    local vm  = `MK'[1,3] + `MK'[1,6]
    local vk  = `MK'[1,4]
    local vmk = `MK'[1,5]

    * ---- parametric line at the percentiles -------------------------------------------
    local nm : word count `at'
    tempname MT
    mat `MT' = J(`nm', 4, .)
    di
    di as txt "Marginal treatment effect at percentiles of the participation unobservable"
    di as txt "(" as res "`method'" as txt " estimation; MTE(u) = m + kappa*invnormal(1-u), m = " ///
       as res %8.5g `m' as txt ", kappa = " as res %8.5g `kp' as txt ")"
    di as txt "{hline 8}{c TT}{hline 24}{c TT}{hline 12}"
    di as txt "   u    {c |}   MTE(u)     Std. err.{c |}  support"
    di as txt "{hline 8}{c +}{hline 24}{c +}{hline 12}"
    local i = 0
    foreach u of numlist `at' {
        local ++i
        local cu  = invnormal(1 - `u')
        local est = `m' + `kp'*`cu'
        local se  = sqrt(`vm' + (`cu')^2*`vk' + 2*(`cu')*`vmk')
        local ins = (`u' >= `plo' & `u' <= `phi')
        mat `MT'[`i',1] = `u'
        mat `MT'[`i',2] = `est'
        mat `MT'[`i',3] = `se'
        mat `MT'[`i',4] = `ins'
        if (`ins') local tag "  observed"
        else       local tag "  extrapol."
        di as txt %7.3f `u' " {c |}" as res %11.4f `est' %12.4f `se' as txt "{c |}`tag'"
    }
    di as txt "{hline 8}{c BT}{hline 24}{c BT}{hline 12}"
    di as txt "support: u inside the common support of P(Z) [" %6.4f `plo' ", " %6.4f `phi' "]"
    mat colnames `MT' = u mte se support

    * ---- semiparametric MTE on the probit score ------------------------------------------
    local nsp = 0
    if ("`semipar'" != "") {
        * is the engine available (as a file on the adopath or already in memory)?
        cap _esreg_pwr
        if (_rc == 199) {
            di as txt _n "semipar: _esreg_pwr.ado not found (it ships with esreg) -- semiparametric MTE skipped"
        }
        else {
            local wexp ""
            if ("`wtype'" != "") local wexp "[`wtype' `wexp0']"
            local wexpa = subinstr("`wexp'", "pweight", "aweight", 1)
            tempvar pzs pz2
            qui probit `dv' `zv' `wexp' if `touse'
            qui predict double `pzs' if `touse', pr
            qui gen double `pz2' = `pzs'^2 if `touse'
            local wsp ""
            foreach v of local xl {
                tempvar c i
                qui summarize `v' `wexpa' if `touse', meanonly
                qui gen double `c' = `v' - r(mean) if `touse'
                qui gen double `i' = `c' * `pzs' if `touse'
                local wsp "`wsp' `c' `i'"
            }
            local nsp : word count `at'
            tempname SP
            mat `SP' = J(`nsp', 8, .)
            local i = 0
            foreach u of numlist `at' {
                local ++i
                cap qui _esreg_pwr `y' `wsp' `pzs' `pz2' `wexp' if `touse', ///
                    per(`u') rankvar(`pzs') target(`pzs' `pz2')
                if (_rc) {
                    di as txt "semipar: percentile-weights regression failed at tau = `u' (rc = " _rc ")"
                    continue
                }
                local qt = e(q_tau)
                local hh = e(h)
                local ne = e(N_eff)
                local lo = max(`u' - `hh'*sqrt(2), 0.005)
                local hi = min(`u' + `hh'*sqrt(2), 0.995)
                _pctile `pzs' if `touse', p(`=100*`lo'' `=100*`hi'')
                local wid = r(r2) - r(r1)
                qui lincom _b[`pzs'] + 2*`qt'*_b[`pz2']
                mat `SP'[`i',1] = `u'
                mat `SP'[`i',2] = `qt'
                mat `SP'[`i',3] = `m' + `kp'*invnormal(1 - `qt')
                mat `SP'[`i',4] = r(estimate)
                mat `SP'[`i',5] = r(se)
                mat `SP'[`i',6] = `hh'
                mat `SP'[`i',7] = `ne'
                mat `SP'[`i',8] = `wid'
            }
            di as txt _n "Semiparametric MTE (percentile-weights regression on the probit score, derivative-aware bandwidth)"
            di as txt "{hline 7}{c TT}{hline 8}{c TT}{hline 11}{c TT}{hline 22}{c TT}{hline 8}{c TT}{hline 8}{c TT}{hline 9}"
            di as txt "  tau  {c |} p_tau  {c |} parametric{c |}     PWR     Std. err.{c |}    h   {c |}  N_eff {c |} width p"
            di as txt "{hline 7}{c +}{hline 8}{c +}{hline 11}{c +}{hline 22}{c +}{hline 8}{c +}{hline 8}{c +}{hline 9}"
            forvalues j = 1/`nsp' {
                di as txt %6.2f `SP'[`j',1] " {c |}" as res %7.3f `SP'[`j',2] as txt " {c |}" ///
                   as res %10.4f `SP'[`j',3] as txt " {c |}" ///
                   as res %10.4f `SP'[`j',4] %11.4f `SP'[`j',5] as txt " {c |}" ///
                   as res %7.3f `SP'[`j',6] as txt " {c |}" as res %7.0f `SP'[`j',7] as txt " {c |}" ///
                   as res %8.3f `SP'[`j',8]
            }
            di as txt "{hline 7}{c BT}{hline 8}{c BT}{hline 11}{c BT}{hline 22}{c BT}{hline 8}{c BT}{hline 8}{c BT}{hline 9}"
            di as txt "tau: quantile of the probit score; p_tau: the score at that quantile; parametric: m + kappa*invnormal(1-p_tau)."
            di as txt "PWR std. err. exclude the estimation of the score."
            mat colnames `SP' = tau p_tau parametric pwr se h neff width_p
        }
    }

    * ---- graph -----------------------------------------------------------------------------
    if ("`graph'" != "") {
        preserve
        qui drop _all
        qui set obs 99
        qui gen double u   = _n/100
        qui gen double cu  = invnormal(1-u)
        qui gen double mte = `m' + `kp'*cu
        qui gen double se  = sqrt(`vm' + cu^2*`vk' + 2*cu*`vmk')
        qui gen double lb  = mte - `z'*se
        qui gen double ub  = mte + `z'*se
        local spplot ""
        local splegend ""
        if (`nsp' > 0) {
            qui set obs `=99+`nsp''
            qui gen double qsp  = .
            qui gen double msp  = .
            qui gen double lbsp = .
            qui gen double ubsp = .
            forvalues j = 1/`nsp' {
                qui replace qsp  = `SP'[`j',2] in `=99+`j''
                qui replace msp  = `SP'[`j',4] in `=99+`j''
                qui replace lbsp = `SP'[`j',4] - `z'*`SP'[`j',5] in `=99+`j''
                qui replace ubsp = `SP'[`j',4] + `z'*`SP'[`j',5] in `=99+`j''
            }
            local spplot (rcap lbsp ubsp qsp, lcolor(gs6)) (scatter msp qsp, mcolor(black) msymbol(O))
            local splegend 4 "PWR on the probit score (semiparametric)"
        }
        if ("`name'" != "") local name "name(`name', replace)"
        else                 local name "name(esrmte, replace)"
        if (`"`saving'"' != "") local saving `"saving(`saving')"'
        twoway (rarea lb ub u, color(gs13) lwidth(none))                  ///
               (line mte u, lcolor(black))                                ///
               `spplot'                                                   ///
               , xline(`plo' `phi', lpattern(shortdash) lcolor(gs8))      ///
                 yline(`ate', lpattern(dot) lcolor(gs8))                  ///
                 legend(order(2 "MTE(u), parametric" `splegend') rows(1) size(small)) ///
                 ytitle("Marginal treatment effect", margin(r=3)) ylabel(, angle(horizontal) format(%5.2f)) ///
                 xtitle("u : percentile of the participation unobservable (0 = most eager)", margin(t=2)) ///
                 title("Marginal treatment effect", size(medsmall))       ///
                 note("Dashed verticals: common support of P(Z). Dotted: ATE. `level'% band: delta method.", size(vsmall)) ///
                 graphregion(color(white)) plotregion(margin(small)) `name' `saving'
        restore
    }

    * ---- returns -----------------------------------------------------------------------------
    return matrix mte = `MT'
    if (`nsp' > 0) return matrix mte_sp = `SP'
    return scalar m     = `m'
    return scalar kappa = `kp'
    return scalar se_m  = sqrt(`vm')
    return scalar se_kappa = sqrt(`vk')
    return scalar supp_lo = `plo'
    return scalar supp_hi = `phi'
    return local method "`method'"
end
