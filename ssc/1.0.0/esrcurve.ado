*! esrcurve 1.0.0  03oct2026  A. Araar
*! Expected treatment effect by quantile group of a ranking variable, after esreg:
*! treated, untreated and all, with standard errors from the influence function
*! of the whole procedure (parameters, within-cell sampling, their covariance; by
*! cluster or by the survey design as the estimation), table, r(table) and
*! optional graph.
*!
*! esrcurve [if] [in], RANK(varname) [ NQ(#) EST(name) GRaph noSE Level(#)
*!                                     TItle(string) NAme(string) SAVing(string) ]
*!
*! Works from est(name) if given, else from the current e() if it is an esreg
*! estimation, else from _esreg (the last esreg of the session).
cap program drop esrcurve
program define esrcurve, rclass
    version 16
    syntax [if] [in], RANK(varname numeric) [ NQ(integer 10) EST(name) GRaph noSE ///
                                             Level(cilevel) TItle(string) NAme(string) SAVing(string asis) ]
    esreg_engine
    tempname hold
    _estimates hold `hold', restore copy
    _esreg_getest `est'
    if (`nq' < 2) {
        di as err "nq() must be at least 2"
        exit 198
    }

    marksample touse, novarlist
    * the fitted observations (within the subpopulation after svy, subpop()), the
    * design sample, and how the estimation aggregated them
    tempvar esmp smp
    qui gen byte `smp' = e(sample)
    _esreg_esample `esmp'
    local mode "`r(mode)'"
    local domopt = cond(r(domain), "domain(`esmp')", "")
    qui replace `touse' = 0 if !`esmp'
    markout `touse' `rank'
    local y  "`e(depvar)'"
    local dv "`e(treat)'"

    * weights of the estimation
    tempvar w
    if ("`e(wtype)'" != "") qui gen double `w' `e(wexp)' if `esmp'
    else                    qui gen double `w' = 1 if `esmp'

    * quantile groups of the ranking variable on the sample
    tempvar g
    if ("`e(wtype)'" != "") qui xtile `g' = `rank' [aw = `w'] if `touse', nq(`nq')
    else                    qui xtile `g' = `rank' if `touse', nq(`nq')
    qui summarize `g' if `touse', meanonly
    local G = r(max)

    * plain variables and the Mata call: the influence function of theta on the
    * fitted observations, the cell means over the analysis sample (g is missing
    * outside it)
    _esreg_data if `esmp'
    local xl  "`r(x)'"
    local zl  "`r(z)'"
    local hs  "`r(hs)'"
    local hr  "`r(hr)'"
    local kap "`r(kap)'"
    tempname R VP C
    local ifv ""
    local nif = 6*`G'
    forvalues j = 1/`nif' {
        tempvar f`j'
        local ifv "`ifv' `f`j''"
    }
    mata: _esreg_cells("`y'", "`xl'", "`zl'", "`dv'", "`hs'", "`hr'", "`kap'", "`w'", "`g'", `G', "`esmp'", "`R'", "`ifv'", "`VP'")
    cap drop __esr_fv*
    if ("`mode'" != "iid") {
        * the sampling and covariance terms by cluster or by the survey design
        * (influence functions: U then P, each treated 1..G, untreated, all)
        if ("`mode'" == "cluster") _esreg_ifcov `ifv' if `smp', cluster(`e(clustvar)')
        else                       _esreg_ifcov `ifv' if `smp', svy `domopt'
        mat `C' = r(V)
        forvalues q = 1/`G' {
            forvalues c = 0/2 {
                local k = `c'*`G' + `q'
                mat `R'[`q', 2 + 2*`c'] = sqrt(`VP'[1,`k'] + `C'[`k',`k'] + 2*`C'[`k', 3*`G' + `k'])
            }
        }
    }

    * mean of the ranking variable by group
    tempname M
    mat `M' = J(`G', 1, .)
    forvalues q = 1/`G' {
        qui summarize `rank' [aw = `w'] if `touse' & `g' == `q', meanonly
        mat `M'[`q', 1] = r(mean)
    }
    mat `R' = `M', `R'
    mat colnames `R' = rank_mean treated se_treated untreated se_untreated all se_all n
    local rn ""
    forvalues q = 1/`G' {
        local rn "`rn' q`q'"
    }
    mat rownames `R' = `rn'

    * ---- display ----------------------------------------------------------------
    di
    di as txt "Expected treatment effect by quantile of " as res "`rank'" as txt " (nq = `G'), " ///
       as res "`e(method)'" as txt " estimation of " as res "`y'" as txt " on treatment " as res "`dv'"
    if ("`se'" == "") {
        di as txt "{hline 84}"
        di as txt "  q  {c |}  mean(`rank') {c |}   Treated   Std. err.{c |}  Untreated  Std. err.{c |}     All     Std. err."
        di as txt "{hline 5}{c +}{hline 14}{c +}{hline 22}{c +}{hline 22}{c +}{hline 22}"
        forvalues q = 1/`G' {
            di as txt %4.0f `q' " {c |}" as res %13.4f `R'[`q',1] " " as txt "{c |}" ///
               as res %10.4f `R'[`q',2] %11.4f `R'[`q',3] as txt " {c |}" ///
               as res %10.4f `R'[`q',4] %11.4f `R'[`q',5] as txt " {c |}" ///
               as res %10.4f `R'[`q',6] %11.4f `R'[`q',7]
        }
        di as txt "{hline 84}"
    }
    else {
        di as txt "{hline 52}"
        di as txt "  q  {c |}  mean(`rank') {c |}   Treated  Untreated        All"
        di as txt "{hline 5}{c +}{hline 14}{c +}{hline 31}"
        forvalues q = 1/`G' {
            di as txt %4.0f `q' " {c |}" as res %13.4f `R'[`q',1] " " as txt "{c |}" ///
               as res %10.4f `R'[`q',2] %11.4f `R'[`q',4] %11.4f `R'[`q',6]
        }
        di as txt "{hline 52}"
    }
    di as txt "Std. err.: influence function (parameters + within-cell sampling + covariance); groups by weighted quantiles of `rank'"

    * ---- graph ------------------------------------------------------------------
    if ("`graph'" != "") {
        preserve
        qui clear
        qui svmat double `R', names(col)
        qui gen q = _n
        if ("`title'" == "") local title "Expected treatment effect by quantile of `rank'"
        if ("`name'" != "") local name "name(`name', replace)"
        if (`"`saving'"' != "") local saving `"saving(`saving')"'
        twoway (connected treated q, lcolor(navy) mcolor(navy)) ///
               (connected untreated q, lcolor(maroon) mcolor(maroon) msymbol(square)) ///
               (connected all q, lcolor(black) mcolor(black) msymbol(triangle) lpattern(dash)), ///
               legend(order(1 "Treated" 2 "Untreated" 3 "All") rows(1)) ///
               xtitle("Quantile group of `rank'", margin(t=2)) ///
               ytitle("Expected effect", margin(r=3)) ylabel(, angle(horizontal) format(%5.2f)) ///
               title(`"`title'"', size(medsmall)) xlabel(1(1)`G') ///
               graphregion(color(white)) plotregion(margin(small)) `name' `saving'
        restore
    }

    return matrix table = `R'
    return scalar nq = `G'
    return local rank "`rank'"
    return local method "`e(method)'"
end
