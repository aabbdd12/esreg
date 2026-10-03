*! esrmte 1.0.0  03oct2026  A. Araar
*! Corrected 02oct2026: after esreg, method(twostep) hermite(#), the MTE curve
*! with the Hermite terms, MTE(u) = m + kappa v + sum_k dh_k H_k(v), v =
*! invnormal(1 - u), and its band from the covariance of all its coefficients.
*! Corrected 03oct2026: with Hermite terms in one regime only (hermite(#1 #0)), the
*! absent terms count as zero in dh_k = h_1k - h_0k.
*! Marginal treatment effect after esreg.
*!   Parametric line  MTE(u) = m + kappa * invnormal(1 - u), m = E[X(b1 - b0)],
*!   at chosen percentiles u of the participation unobservable (0 = most eager),
*!   with standard errors from the influence function of the whole procedure
*!   (parameters, sampling of m and kappa, their covariance; by cluster or by the
*!   survey design as the estimation) and common-support flags;
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
    * the fitted observations (within the subpopulation after svy, subpop()), the
    * design sample, and how the estimation aggregated them
    tempvar esmp smp
    qui gen byte `smp' = e(sample)
    _esreg_esample `esmp'
    local mode "`r(mode)'"
    local domopt = cond(r(domain), "domain(`esmp')", "")
    qui replace `touse' = 0 if !`esmp'
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
    if ("`mode'" == "svy" & e(df_r) < .) local z = invttail(e(df_r), (1 - `level'/100)/2)

    tempvar w
    if ("`wtype'" != "") qui gen double `w' `wexp0' if `esmp'
    else                 qui gen double `w' = 1 if `esmp'

    * ---- m, kappa and their covariance --------------------------------------------
    * (the influence function of theta on the fitted observations, the means of m
    * and kappa over the analysis sample)
    _esreg_data if `esmp'
    local xl  "`r(x)'"
    local zl  "`r(z)'"
    local hs  "`r(hs)'"
    local hr  "`r(hr)'"
    local kap "`r(kap)'"
    * the coefficients of the curve: b = (m, kappa) for the line; with hermite(#)
    * (two-step) also dh_k = h_1k - h_0k, the differences of the coefficients of
    * the Hermite controls: MTE(v) = m + kappa v + dh2 (v^2 - 1) + dh3 (v^3 - 3v)
    local ph = 0
    if ("`method'" == "twostep" & e(k_h) < .) local ph = e(k_h)
    local q = 2 + `ph'
    tempname MK V2 C CF B VT
    local ifv ""
    forvalues j = 1/`=2*`q'' {
        tempvar f`j'
        local ifv "`ifv' `f`j''"
    }
    mata: _esreg_mte("`y'", "`xl'", "`zl'", "`dv'", "`hs'", "`hr'", "`kap'", "`w'", "`esmp'", "`touse'", "`MK'", "`ifv'", "`V2'", "`CF'")
    cap drop __esr_fv*
    mat `B'  = `CF'[1..., 1]'
    mat `VT' = `CF'[1..., 2...]
    if ("`mode'" != "iid") {
        * the sampling and covariance terms by cluster or by the survey design
        * (order of the influence functions: U_1..U_q, P_1..P_q)
        if ("`mode'" == "cluster") _esreg_ifcov `ifv' if `smp', cluster(`e(clustvar)')
        else                       _esreg_ifcov `ifv' if `smp', svy `domopt'
        mat `C' = r(V)
        local q1 = `q' + 1
        local q2 = 2*`q'
        mat `VT' = `V2' + `C'[1..`q', 1..`q'] + `C'[1..`q', `q1'..`q2'] + `C'[`q1'..`q2', 1..`q']
    }
    local m   = `B'[1,1]
    local kp  = `B'[1,2]
    local vm  = `VT'[1,1]
    local vk  = `VT'[2,2]
    local vmk = `VT'[1,2]
    local dh2 = cond(`ph' >= 1, `B'[1, min(3, `q')], 0)
    local dh3 = cond(`ph' >= 2, `B'[1, min(4, `q')], 0)

    * ---- parametric line at the percentiles -------------------------------------------
    local nm : word count `at'
    tempname MT GV GE GS
    mat `MT' = J(`nm', 4, .)
    di
    di as txt "Marginal treatment effect at percentiles of the participation unobservable"
    if (`ph' == 0) {
        di as txt "(" as res "`method'" as txt " estimation; MTE(u) = m + kappa*invnormal(1-u), m = " ///
           as res %8.5g `m' as txt ", kappa = " as res %8.5g `kp' as txt ")"
    }
    else {
        local dh3t ""
        if (`ph' >= 2) local dh3t = ", dh3 = " + strofreal(`dh3', "%8.5g")
        local hspec "`=e(hermite)'"
        if (e(hermite1) < . & e(hermite0) < . & e(hermite1) != e(hermite0)) local hspec "`=e(hermite1)' `=e(hermite0)'"
        di as txt "(two-step, hermite(`hspec'); MTE(u) = m + kappa*v + dh2*(v^2-1)" ///
           cond(`ph' >= 2, " + dh3*(v^3-3v)", "") ", v = invnormal(1-u);"
        di as txt " m = " as res %8.5g `m' as txt ", kappa = " as res %8.5g `kp' as txt ///
           ", dh2 = " as res %8.5g `dh2' as txt "`dh3t')"
    }
    di as txt "{hline 8}{c TT}{hline 24}{c TT}{hline 12}"
    di as txt "   u    {c |}   MTE(u)     Std. err.{c |}  support"
    di as txt "{hline 8}{c +}{hline 24}{c +}{hline 12}"
    local i = 0
    foreach u of numlist `at' {
        local ++i
        local cu  = invnormal(1 - `u')
        * the gradient of MTE(u) in the coefficients of the curve
        mat `GV' = (1, `cu')
        if (`ph' >= 1) mat `GV' = `GV', ((`cu')^2 - 1)
        if (`ph' >= 2) mat `GV' = `GV', ((`cu')^3 - 3*(`cu'))
        mat `GE' = `GV' * `B''
        mat `GS' = `GV' * `VT' * `GV''
        local est = `GE'[1,1]
        local se  = sqrt(`GS'[1,1])
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
    if (`ph' > 0) di as txt "outside the support the polynomial curve extrapolates its Hermite terms (poorly beyond the data)"
    mat colnames `MT' = u mte se support

    * ---- semiparametric MTE on the probit score ------------------------------------------
    local nsp = 0
    if ("`semipar'" != "") {
        * is the engine available (as a file on the adopath or already in memory)?
        cap _esreg_pwr
        local rcp = _rc
        if (`rcp' != 199) {
            * the ado is there; make sure its Mata functions are compiled (an autoloaded
            * ado does not execute its mata: block): compile the file explicitly if not
            mata: st_local("epwr_ok", strofreal(findexternal("_epwr_main()") != NULL & findexternal("_epwr_revCumSumT()") != NULL))
            if ("`epwr_ok'" != "1") {
                cap findfile _esreg_pwr.ado
                if (_rc == 0) cap noisily version `c(stata_version)': run "`r(fn)'"
            }
        }
        if (`rcp' == 199) {
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
                local cq = invnormal(1 - `qt')
                mat `SP'[`i',3] = `m' + `kp'*`cq' + `dh2'*(`cq'^2 - 1) + `dh3'*(`cq'^3 - 3*`cq')
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
            di as txt "tau: quantile of the probit score; p_tau: the score at that quantile; parametric: the curve above at u = p_tau."
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
        qui gen double mte = `m' + `kp'*cu + `dh2'*(cu^2 - 1) + `dh3'*(cu^3 - 3*cu)
        * the variance of the curve: g' VT g, g = (1, v, v^2 - 1, v^3 - 3v)[1..q]
        qui gen double g1 = 1
        qui gen double g2 = cu
        qui gen double g3 = cu^2 - 1
        qui gen double g4 = cu^3 - 3*cu
        qui gen double se = 0
        forvalues i = 1/`q' {
            forvalues j = 1/`q' {
                qui replace se = se + g`i'*g`j'*`VT'[`i', `j']
            }
        }
        qui replace se = sqrt(se)
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
                 note("Dashed verticals: common support of P(Z). Dotted: ATE. `level'% band: influence function.", size(vsmall)) ///
                 graphregion(color(white)) plotregion(margin(small)) `name' `saving'
        restore
    }

    * ---- returns -----------------------------------------------------------------------------
    return matrix mte = `MT'
    if (`nsp' > 0) return matrix mte_sp = `SP'
    local cn "m kappa"
    if (`ph' >= 1) local cn "`cn' dh2"
    if (`ph' >= 2) local cn "`cn' dh3"
    mat colnames `B' = `cn'
    mat rownames `VT' = `cn'
    mat colnames `VT' = `cn'
    return matrix curve = `B'
    return matrix V_curve = `VT'
    return scalar hermite = cond(`ph' > 0, `ph' + 1, 0)
    return scalar m     = `m'
    return scalar kappa = `kp'
    return scalar se_m  = sqrt(`vm')
    return scalar se_kappa = sqrt(`vk')
    return scalar supp_lo = `plo'
    return scalar supp_hi = `phi'
    return local method "`method'"
end
