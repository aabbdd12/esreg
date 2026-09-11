*! esreg 0.4.0  11sep2026  A. Araar (Universite Laval / PEP)
*! Endogenous switching regression: FIML (ml lf1, analytic score) or two-step
*! (probit + OLS with Mills ratios, exact stacked-moment variance), with the
*! treatment effects ATT / ATU / ATE and kappa = rho1*sigma1 - rho0*sigma0.
*! Files: esreg.ado (this), esreg_lf1.ado (ml evaluator), esreg_engine.ado
*! (shared Mata engine).  Mirrors the Python reference esreg v0.1.
*!
*! esreg depvar [indepvars] [if] [in] [pw fw iw], SELect(treatvar = varlist)
*!        [ METHod(fiml|twostep) HETSigma(varlist) HETRho(varlist) KAPpa(varlist)
*!          vce(oim|opg|robust|cluster clustvar) noEFFects Level(#) noLOg
*!          ITERate(#) DIFficult STORe(name) ]
*!
*! Stored results (see esreg_returns.txt for the full list).  The estimation is
*! always stored as -estimates store _esreg- (last esreg of the session); a
*! named copy is added with store().  Post-commands look for est(name), then the
*! current e() if e(cmd) == "esreg", then _esreg.
*!   scalars : N N_treated N_untreated sum_w k_x k_z k_hs k_hr k_kap
*!             ll ll_indep lr_indep lr_df p_indep (fiml)  converged
*!             att atu ate kappa  se_att se_atu se_ate se_kappa
*!             sigma1 sigma0 rho1 rho0 rhosig1 rhosig0 (means when heterogeneous)
*!             supp_lo supp_hi p_min1 p_max1 p_min0 p_max0 ml1 ml0
*!   macros  : cmd cmdline version method depvar treat xvars zvars hetsigma hetrho
*!             kappavars wtype wexp vce vcetype eqnames predict title store
*!   matrices: b V effects(4x4: est se var_param var_samp) b_sel b_1 b_0 anc(1x4)
*!             support(1x4) lambda(1x2)
*!   function: sample


* repeated-test convenience: drop previous definitions before redefining
cap program drop esreg
cap program drop _esreg_display
cap program drop _esreg_llindep
cap program drop _esreg_fvx
cap program drop _esreg_fvone

program define esreg, eclass properties(svyb svyj svyr)
    version 16
    if replay() {
        if ("`e(cmd)'" != "esreg") error 301
        syntax [, Level(cilevel) EFFects]
        * after svy: the linearized e(V) replaces the model-based one; the effects,
        * kappa and their standard errors are recomputed from it once (or on request)
        if ("`effects'" != "" | ("`e(prefix)'" == "svy" & "`e(eff_vce)'" != "svy")) {
            esreg_engine
            _esreg_effects_post
        }
        _esreg_display, level(`level')
        exit
    }
    syntax varlist(min=1 numeric fv) [if] [in] [pw fw iw] , SELect(string) ///
        [ METHod(string) HETSigma(varlist numeric fv) HETRho(varlist numeric fv) ///
          KAPpa(varlist numeric fv) vce(passthru) noEFFects Level(cilevel) ///
          noLOg ITERate(passthru) DIFficult STORe(name) noSTORE2 HERmite ///
          HSize(varname numeric) noSVYset ]

    * shared engine (esreg_engine.ado): Mata functions and helpers, loaded on demand
    esreg_engine
    * ---- parse -----------------------------------------------------------
    gettoken y xvars : varlist
    gettoken dvar zvars : select, parse("=")
    gettoken eq zvars : zvars, parse("=")
    if ("`eq'" != "=") {
        di as err "select() must be select(treatvar = varlist)"
        exit 198
    }
    local dvar = trim("`dvar'")
    confirm variable `dvar'
    local zvars = trim("`zvars'")
    if ("`method'" == "")            local method fiml
    if ("`method'" == "ml")          local method fiml
    if ("`method'" == "2s" | "`method'" == "twostep" | "`method'" == "2step") local method twostep
    if (!inlist("`method'", "fiml", "twostep")) {
        di as err "method() must be fiml or twostep"
        exit 198
    }
    if ("`method'" == "twostep" & ("`hetsigma'`hetrho'" != "")) {
        di as err "hetsigma() and hetrho() are options of method(fiml); use kappa() with method(twostep)"
        exit 198
    }
    if ("`method'" == "fiml" & "`kappa'" != "") {
        di as err "kappa() is an option of method(twostep); use hetsigma()/hetrho() with method(fiml)"
        exit 198
    }
    if ("`hermite'" != "" & "`method'" != "twostep") {
        di as err "hermite is an option of method(twostep)"
        exit 198
    }
    local ph = 0
    if ("`hermite'" != "") {
        local ph = 2
        local effects noeffects
    }
    marksample touse
    markout `touse' `dvar'
    foreach o in zvars hetsigma hetrho kappa {
        if ("``o''" != "") {
            fvrevar ``o'', list
            markout `touse' `r(varlist)'
        }
    }
    qui count if `touse'
    local N = r(N)
    qui count if `touse' & `dvar' == 1
    local N1 = r(N)
    local N0 = `N' - `N1'
    cap assert inlist(`dvar', 0, 1) if `touse'
    if (_rc) {
        di as err "`dvar' must be 0/1"
        exit 450
    }
    * ---- weights: svyset pweight taken by default when none is given (nosvyset
    *      declines it); hsize() multiplies the weight (effects per individual when
    *      the unit is the household); under svy: (iweights) hsize() is not allowed
    if ("`weight'" == "" & "`svyset'" == "") {
        cap qui svyset
        if (_rc == 0 & "`r(wvar)'" != "") {
            local weight "pweight"
            local exp "= `r(wvar)'"
            di as txt "(svyset weight " as res "`r(wvar)'" as txt " used as pweight; design-based standard errors need the svy: prefix)"
        }
    }
    else if ("`weight'" != "" & "`weight'" != "iweight") {
        * an explicit weight is used once: the svyset weight, if any, is not applied on top
        cap qui svyset
        if (_rc == 0 & "`r(wvar)'" != "") {
            local exp0 = trim(subinstr("`exp'", "=", "", 1))
            if ("`exp0'" != "`r(wvar)'") di as txt "(explicit weight [`weight'`exp'] used; the svyset weight " as res "`r(wvar)'" as txt " is not applied on top of it)"
        }
    }
    if ("`hsize'" != "") {
        if ("`weight'" == "iweight") {
            di as err "hsize() is not allowed under the svy: prefix (the design weights are already in use)"
            exit 198
        }
        if ("`weight'" == "") {
            local weight "pweight"
            local exp "= `hsize'"
        }
        else {
            local exp0 = subinstr("`exp'", "=", "", 1)
            local exp "= (`exp0') * `hsize'"
        }
        markout `touse' `hsize'
    }
    local wexp ""
    if ("`weight'" != "") local wexp "[`weight'`exp']"
    tempvar wv
    if ("`weight'" != "") qui gen double `wv' `exp' if `touse'
    else                  qui gen double `wv' = 1     if `touse'

    tempvar tuse
    qui gen byte `tuse' = `touse'
    qui summarize `wv' if `touse', meanonly
    local sumw = r(sum)

    * expand factor variables once, so that Mata sees plain columns
    _esreg_fvx `xvars' if `touse'
    local xlist "`r(varlist)'"
    _esreg_fvx `zvars' if `touse'
    local zlist "`r(varlist)'"
    foreach o in hetsigma hetrho kappa {
        if ("``o''" != "") {
            _esreg_fvx ``o'' if `touse'
            local `o'list "`r(varlist)'"
        }
    }
    * plain (temporary) variables behind factor terms, for Mata and for generate;
    * one term at a time, otherwise 2.region 3.region is re-merged with 2 as base
    foreach o in xlist zlist hetsigmalist hetrholist kappalist {
        local `o'_m ""
        foreach v of local `o' {
            cap confirm variable `v'
            if (_rc == 0) {
                local `o'_m "``o'_m' `v'"
            }
            else {
                tempvar tv
                _esreg_fvone `v' `tv' if `touse'
                local `o'_m "``o'_m' `tv'"
            }
        }
    }

    * ======================================================================
    if ("`method'" == "fiml") {
        * starting values: probit for gamma, OLS by regime for beta and ln sigma
        tempname b0mat
        qui probit `dvar' `zlist_m' `wexp' if `touse'
        tempname bg
        mat `bg' = e(b)
        qui regress `y' `xlist_m' `wexp' if `touse' & `dvar' == 1
        tempname b1
        mat `b1' = e(b)
        local ls1 = ln(e(rmse) * sqrt((e(N) - e(df_m) - 1) / e(N)))
        qui regress `y' `xlist_m' `wexp' if `touse' & `dvar' == 0
        tempname b00
        mat `b00' = e(b)
        local ls0 = ln(e(rmse) * sqrt((e(N) - e(df_m) - 1) / e(N)))

        local nhs : word count `hetsigmalist'
        local nhr : word count `hetrholist'
        tempname z1 z2
        mat `z1' = J(1, `nhs' + 1, 0)
        mat `z2' = J(1, `nhr' + 1, 0)
        mat `z1'[1, `nhs' + 1] = `ls1'
        mat `b0mat' = `b1', `b00', `bg', `z1'
        mat `z1'[1, `nhs' + 1] = `ls0'
        mat `b0mat' = `b0mat', `z1', `z2', `z2'

        * log-likelihood of the independent-equations model (closed form), before ml posts e()
        tempname llind
        _esreg_llindep `y' `dvar' `wexp', touse(`touse') xlist(`xlist_m') zlist(`zlist_m')
        scalar `llind' = r(ll)

        ml model lf1 esreg_lf1 (`y'_1: `y' `dvar' = `xlist_m') (`y'_0: `xlist_m') ///
            (`dvar': `zlist_m') (lnsigma_1: `hetsigmalist_m') (lnsigma_0: `hetsigmalist_m') ///
            (atanhrho_1: `hetrholist_m') (atanhrho_0: `hetrholist_m') ///
            `wexp' if `touse', `vce' `iterate' `difficult' maximize init(`b0mat', copy) ///
            search(off) `log' nooutput title(Endogenous switching regression -- full information ML)
        * put the expanded factor-variable names back on e(b) and e(V)
        local cn ""
        foreach v of local xlist {
            local cn `cn' `y'_1:`v'
        }
        local cn `cn' `y'_1:_cons
        foreach v of local xlist {
            local cn `cn' `y'_0:`v'
        }
        local cn `cn' `y'_0:_cons
        foreach v of local zlist {
            local cn `cn' `dvar':`v'
        }
        local cn `cn' `dvar':_cons
        foreach e in lnsigma_1 lnsigma_0 {
            foreach v of local hetsigmalist {
                local cn `cn' `e':`v'
            }
            local cn `cn' `e':_cons
        }
        foreach e in atanhrho_1 atanhrho_0 {
            foreach v of local hetrholist {
                local cn `cn' `e':`v'
            }
            local cn `cn' `e':_cons
        }
        tempname bb VV
        mat `bb' = e(b)
        mat `VV' = e(V)
        mat colnames `bb' = `cn'
        mat colnames `VV' = `cn'
        mat rownames `VV' = `cn'
        ereturn repost b = `bb' V = `VV', rename
        ereturn scalar ll_indep = `llind'
        ereturn scalar lr_indep = 2 * (e(ll) - `llind')
        ereturn scalar lr_df    = 2 * (`nhr' + 1)
        ereturn scalar p_indep  = chi2tail(e(lr_df), e(lr_indep))
        ereturn local  method   "fiml"
        ereturn local  eqnames  "`y'_1 `y'_0 `dvar' lnsigma_1 lnsigma_0 atanhrho_1 atanhrho_0"
    }
    * ======================================================================
    else {
        * ---- two-step -------------------------------------------------------
        tempvar zg l1 l0
        qui probit `dvar' `zlist_m' `wexp' if `touse'
        tempname bg
        mat `bg' = e(b)
        qui predict double `zg' if `touse', xb
        qui gen double `l1' =  exp(lnnormalden(`zg') - lnnormal(`zg'))  if `touse'
        qui gen double `l0' =  exp(lnnormalden(`zg') - lnnormal(-`zg')) if `touse'
        local lam1 ""
        local lam0 ""
        local knames ""
        local j = 0
        foreach v of local kappalist_m {
            local ++j
            tempvar m1_`j' m0_`j'
            qui gen double `m1_`j'' =  `l1' * `v' if `touse'
            qui gen double `m0_`j'' = -`l0' * `v' if `touse'
            local lam1 "`lam1' `m1_`j''"
            local lam0 "`lam0' `m0_`j''"
            local vn : word `j' of `kappalist'
            local vn = subinstr("`vn'", ".", "_", .)
            local knames "`knames' lambda_`vn'"
        }
        tempvar ml0
        qui gen double `ml0' = -`l0' if `touse'
        local lam1 "`lam1' `l1'"
        local lam0 "`lam0' `ml0'"
        local knames "`knames' lambda"
        * Hermite controls E[u^2 - 1 | D, Z], E[u^3 - 3u | D, Z] (esrtest, normal)
        local hnames ""
        if (`ph' > 0) {
            tempvar h2_1 h3_1 h2_0 h3_0
            qui gen double `h2_1' = -`zg' * `l1' if `touse'
            qui gen double `h3_1' = (`zg'^2 - 1) * `l1' if `touse'
            qui gen double `h2_0' = `zg' * `l0' if `touse'
            qui gen double `h3_0' = -(`zg'^2 - 1) * `l0' if `touse'
            local lam1 "`lam1' `h2_1' `h3_1'"
            local lam0 "`lam0' `h2_0' `h3_0'"
            local hnames "h2 h3"
        }
        qui regress `y' `xlist_m' `lam1' `wexp' if `touse' & `dvar' == 1
        tempname c1
        mat `c1' = e(b)
        qui regress `y' `xlist_m' `lam0' `wexp' if `touse' & `dvar' == 0
        tempname c0
        mat `c0' = e(b)
        tempname b V
        mat `b' = `bg', `c1', `c0'
        * column names: d:z ... ; y_1: x ... lambda ; y_0: x ... lambda
        local cn ""
        foreach v of local zlist {
            local cn `cn' `dvar':`v'
        }
        local cn `cn' `dvar':_cons
        foreach v of local xlist {
            local cn `cn' `y'_1:`v'
        }
        foreach v of local knames {
            local cn `cn' `y'_1:`v'
        }
        foreach v of local hnames {
            local cn `cn' `y'_1:`v'
        }
        local cn `cn' `y'_1:_cons
        foreach v of local xlist {
            local cn `cn' `y'_0:`v'
        }
        foreach v of local knames {
            local cn `cn' `y'_0:`v'
        }
        foreach v of local hnames {
            local cn `cn' `y'_0:`v'
        }
        local cn `cn' `y'_0:_cons
        * reorder c1, c0 so that lambda terms come before _cons? No: keep Stata's
        * regress order (x..., lambda..., _cons), which the Mata code expects.
        mat colnames `b' = `cn'
        * stacked-moment variance in Mata
        mata: _esreg_ts_var("`y'", "`xlist_m'", "`zlist_m'", "`dvar'", "`kappalist_m'", "`wv'", "`touse'", "`b'", "`V'", `ph')
        mat rownames `V' = `cn'
        mat colnames `V' = `cn'
        ereturn post `b' `V' `wexp', depname(`y') obs(`N') esample(`touse')
        ereturn local method  "twostep"
        ereturn local title   "Endogenous switching regression -- two-step"
        ereturn local eqnames "`dvar' `y'_1 `y'_0"
        ereturn local vce     "stacked"
        ereturn local vcetype ""
        ereturn scalar converged = 1
    }

    * ---- common e() ---------------------------------------------------------
    ereturn local cmd       "esreg"
    ereturn local cmdline   "esreg `0'"
    ereturn local version   "0.4.0"
    ereturn local depvar    "`y'"
    ereturn local treat     "`dvar'"
    ereturn local xvars     "`xlist'"
    ereturn local zvars     "`zlist'"
    ereturn local hetsigma  "`hetsigmalist'"
    ereturn local hetrho    "`hetrholist'"
    ereturn local kappavars "`kappalist'"
    ereturn local wtype     "`weight'"
    ereturn local wexp      "`exp'"
    ereturn local hsize     "`hsize'"
    ereturn local eff_vce   "model"
    ereturn local predict   "esreg_p"
    ereturn scalar N_treated   = `N1'
    ereturn scalar N_untreated = `N0'
    ereturn scalar sum_w       = `sumw'
    ereturn scalar k_x   = `: word count `xlist'' + 1
    ereturn scalar k_z   = `: word count `zlist'' + 1
    ereturn scalar k_hs  = `: word count `hetsigmalist'' + 1
    ereturn scalar k_hr  = `: word count `hetrholist'' + 1
    ereturn scalar k_kap = `: word count `kappalist'' + 1
    ereturn scalar k_h   = `ph'
    * coefficient blocks by equation, for the post-commands
    tempname bb bs b1m b0m
    mat `bb' = e(b)
    mat `bs' = `bb'[1, "`dvar':"]
    mat `b1m' = `bb'[1, "`y'_1:"]
    mat `b0m' = `bb'[1, "`y'_0:"]
    ereturn matrix b_sel = `bs'
    ereturn matrix b_1   = `b1m'
    ereturn matrix b_0   = `b0m'

    * ---- effects ------------------------------------------------------------
    if ("`effects'" == "") {
        tempname E A S L
        mata: _esreg_effects("`y'", "`xlist_m'", "`zlist_m'", "`dvar'", "`hetsigmalist_m'", "`hetrholist_m'", "`kappalist_m'", "`wv'", "`tuse'", "`E'")
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
    }
    * ---- keep the estimation in memory: always as _esreg, and under store() ----
    ereturn local store "`store'"
    if ("`store2'" == "") qui estimates store _esreg
    if ("`store'" != "") qui estimates store `store'
    _esreg_display, level(`level')
end

* ---------------------------------------------------------------------------
* ---------------------------------------------------------------------------
* expand factor variables and drop base / omitted levels
program define _esreg_fvx, rclass
    syntax [varlist(default=none fv)] [if]
    local out ""
    if ("`varlist'" != "") {
        fvexpand `varlist' `if'
        foreach v in `r(varlist)' {
            if (strpos("`v'", "b.") == 0 & strpos("`v'", "o.") == 0) local out "`out' `v'"
        }
    }
    return local varlist "`out'"
end

* ---------------------------------------------------------------------------
* fill a caller-owned variable with one factor term (fvrevar creates temporaries
* that vanish when this subprogram ends, hence the copy); falls back to a
* generated indicator for a simple #.varname term if fvrevar returns a constant
program define _esreg_fvone
    syntax anything(name=args) [if]
    gettoken term target : args
    fvrevar `term' `if'
    local v "`r(varlist)'"
    qui summarize `v' `if', meanonly
    if (r(min) == r(max)) {
        local lev = substr("`term'", 1, strpos("`term'", ".") - 1)
        local var = substr("`term'", strpos("`term'", ".") + 1, .)
        cap confirm number `lev'
        if (_rc) {
            di as err "cannot expand factor term `term'"
            exit 198
        }
        qui gen double `target' = (`var' == `lev') `if'
    }
    else qui gen double `target' = `v' `if'
end

* ---------------------------------------------------------------------------
program define _esreg_display
    syntax [, Level(cilevel)]
    if ("`e(method)'" == "fiml" & "`e(prefix)'" == "svy") {
        di
        di as txt "Endogenous switching regression -- full information ML, survey design (svy)"
        di as txt "Number of obs = " as res %9.0f e(N) as txt "   Design df = " as res %5.0f e(df_r) ///
           as txt "   Std. err.: linearized"
        ereturn display, level(`level')
        di as txt "(the LR test of independent equations is not available under svy)"
    }
    else if ("`e(method)'" == "fiml") {
        di
        ml display, level(`level')
        di as txt "LR test of independent equations (rho1 = rho0 = 0): chi2(" as res e(lr_df) ///
           as txt ") = " as res %8.2f e(lr_indep) as txt "   Prob > chi2 = " as res %6.4f chi2tail(e(lr_df), e(lr_indep))
    }
    else {
        di
        di as txt "`e(title)'" _col(49) as txt "Number of obs = " as res %9.0f e(N)
        di as txt "Outcome: " as res "`e(depvar)'" as txt "   Treatment: " as res "`e(treat)'" ///
           _col(49) as txt "Treated       = " as res %9.0f e(N_treated)
        di as txt "Std. err.: exact variance of the two-step procedure (stacked moments)"
        ereturn display, level(`level')
    }
    di as txt "sigma_1 = " as res %8.6g e(sigma1) as txt "   sigma_0 = " as res %8.6g e(sigma0) ///
       as txt "   rho_1 = " as res %8.6g e(rho1) as txt "   rho_0 = " as res %8.6g e(rho0) ///
       cond("`e(hetsigma)'`e(hetrho)'" != "", "   (means; heterogeneous)", "")
    cap confirm matrix e(effects)
    if (_rc) exit
    tempname E
    mat `E' = e(effects)
    local z = invnormal(1 - (1 - `level'/100)/2)
    di
    di as txt "Treatment effects (delta method on e(V) + sampling component)"
    di as txt "{hline 76}"
    di as txt "Index    {c |}   Estimate   Std. err.        z    P>|z|     [`level'% conf. interval]"
    di as txt "{hline 9}{c +}{hline 66}"
    foreach r in ATT ATU ATE kappa {
        local i = rownumb(`E', "`r'")
        local b = `E'[`i', 1]
        local s = `E'[`i', 2]
        di as txt %8s "`r'" " {c |}" as res %11.6g `b' "  " %10.6g `s' "  " %7.2f `b'/`s' ///
           "  " %6.3f 2*normal(-abs(`b'/`s')) "   " %10.6g `b' - `z'*`s' "  " %10.6g `b' + `z'*`s'
    }
    di as txt "{hline 76}"
    di as txt "kappa = rho1*sigma1 - rho0*sigma0;  ATT - ATU = kappa*(mean lambda1 + mean lambda0)"
    di as txt "Variance components (ATT): parameter " as res %8.3g `E'[1,3] as txt "   sampling " as res %8.3g `E'[1,4]
    di as txt "Common support of P(Z): [" as res %6.4f e(supp_lo) as txt ", " as res %6.4f e(supp_hi) as txt "]"
end

* ---------------------------------------------------------------------------
* log-likelihood of the model with rho1 = rho0 = 0: probit + normal OLS by regime
program define _esreg_llindep, rclass
    syntax varlist(min=2 max=2) [pw fw iw], touse(string) xlist(string) zlist(string)
    gettoken y d : varlist
    local wexp ""
    if ("`weight'" != "") local wexp "[`weight'`exp']"
    tempvar w
    if ("`weight'" != "") qui gen double `w' `exp' if `touse'
    else                  qui gen double `w' = 1     if `touse'
    qui probit `d' `zlist' `wexp' if `touse'
    local ll = e(ll)
    foreach j in 1 0 {
        qui regress `y' `xlist' `wexp' if `touse' & `d' == `j'
        tempvar r
        qui predict double `r' if `touse' & `d' == `j', resid
        qui summarize `w' if `touse' & `d' == `j', meanonly
        local nj = r(sum)
        tempvar r2
        qui gen double `r2' = `w' * `r'^2 if `touse' & `d' == `j'
        qui summarize `r2' if `touse' & `d' == `j', meanonly
        local s2 = r(sum) / `nj'
        local ll = `ll' - 0.5 * `nj' * (ln(2*_pi) + ln(`s2') + 1)
    }
    return scalar ll = `ll'
end

