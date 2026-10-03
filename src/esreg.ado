*! esreg 1.0.0  03oct2026  A. Araar (Universite Laval / PEP)
*! Endogenous switching regression: FIML (ml lf1, analytic score) or two-step
*! (probit + OLS with Mills ratios, exact stacked-moment variance), with the
*! treatment effects ATT / ATU / ATE and kappa = rho1*sigma1 - rho0*sigma0.
*! Files: esreg.ado (this), esreg_lf1.ado (ml evaluator), esreg_engine.ado
*! (shared Mata engine).  Mirrors the Python reference esreg v0.1.
*!
*! [by varlist:] esreg depvar [indepvars] [if] [in] [pw iw], SELect(treatvar = varlist)
*!        [ METHod(fiml|twostep) HETSigma(varlist) HETRho(varlist) KAPpa(varlist)
*!          vce(vcetype) noEFFects Level(#) noLOg ITERate(#) DIFficult STORe(name) ]
*!   vce(): fiml   oim (default) | opg | robust | cluster clustvar  (svy: the svy prefix)
*!          twostep  robust (default: the stacked moments) | cluster clustvar | svy
*! Corrected 02oct2026: the standard errors of the effects include the covariance between the
*! parameter and sampling parts (influence function of the whole procedure) and
*! follow the clusters and the survey design; method(twostep) honours
*! vce(cluster) and has vce(svy); the effects after svy, subpop() are those of
*! the subpopulation; fweights are refused.
*! Corrected 03oct2026: the augmented two-step method(twostep) hermite(#) (Hermite
*! terms of order 2 or 3 in E[w_j | u]), in both regimes or by regime with
*! hermite(#1 #0) (treated, untreated; 0 = none): e(hermite1) e(hermite0) e(k_h1) e(k_h0).
*!
*! Stored results (see esreg_returns.txt for the full list).  The estimation is
*! always stored as -estimates store _esreg- (last esreg of the session); a
*! named copy is added with store().  Post-commands look for est(name), then the
*! current e() if e(cmd) == "esreg", then _esreg.
*!   scalars : N N_treated N_untreated sum_w k_x k_z k_hs k_hr k_kap
*!             hermite hermite1 hermite0 k_h k_h1 k_h0 (two-step, hermite())
*!             ll ll_indep lr_indep lr_df p_indep (fiml)  converged
*!             att atu ate kappa  se_att se_atu se_ate se_kappa
*!             sigma1 sigma0 rho1 rho0 rhosig1 rhosig0 (means when heterogeneous)
*!             supp_lo supp_hi p_min1 p_max1 p_min0 p_max0 ml1 ml0
*!             N_clust (cluster)  N_strata N_psu df_r (twostep, vce(svy))
*!   macros  : cmd cmdline version method depvar treat xvars zvars hetsigma hetrho
*!             kappavars wtype wexp vce vcetype clustvar eff_vce eqnames predict
*!             title store
*!   matrices: b V effects(4x5: est se var_param var_samp cov_ps) b_sel b_1 b_0
*!             anc(1x4) support(1x4) lambda(1x2)
*!   function: sample


* repeated-test convenience: drop previous definitions before redefining
cap program drop esreg
cap program drop _esreg_display
cap program drop _esreg_llindep
cap program drop _esreg_fvx
cap program drop _esreg_fvone
cap program drop _esreg_hermopt

program define esreg, eclass byable(recall) properties(svyb svyj svyr)
    version 16
    if replay() {
        if ("`e(cmd)'" != "esreg") error 301
        syntax [, Level(cilevel) EFFects]
        * after svy: the linearized e(V) replaces the model-based one; the effects,
        * kappa and their standard errors are recomputed from it once (or on request)
        if ("`effects'" != "" | ("`e(prefix)'" == "svy" & !inlist("`e(eff_vce)'", "svy", "none"))) {
            esreg_engine
            _esreg_effects_post
        }
        _esreg_display, level(`level')
        exit
    }
    * hermite alone is hermite(3) (syntax has no option with an optional argument);
    * e(cmdline) keeps the command as typed
    local cmdline0 `"`0'"'
    _esreg_hermopt `0'
    local 0 `"`r(cmd)'"'
    syntax varlist(min=1 numeric fv) [if] [in] [pw fw iw] , SELect(string) ///
        [ METHod(string) HETSigma(varlist numeric fv) HETRho(varlist numeric fv) ///
          KAPpa(varlist numeric fv) VCE(string) noEFFects Level(cilevel) ///
          noLOg ITERate(passthru) DIFficult STORe(name) noSTORE2 HERmite(numlist integer min=1 max=2 >=0 <=3) ///
          HSize(varname numeric) noSVYset ]
    * fweights are refused: a unit of the survey stands for its sampling weight,
    * it is not a replicated record (both variances count the units)
    if ("`weight'" == "fweight") {
        di as err "fweights are not allowed: sampling weights go in {bf:[pweight=]};"
        di as err "a frequency weight built from a sampling weight, such as int(pw*10000), is a pweight;"
        di as err "replicated records (true frequencies) can be expanded first: {bf:expand} {it:fwvar}"
        exit 101
    }

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
    * hermite(#): the same order in both regimes; hermite(#1 #0): the order in the
    * treated regime, then in the untreated one, 0 for no Hermite terms there
    local herm1 = 0
    local herm0 = 0
    if ("`hermite'" != "") {
        local herm1 : word 1 of `hermite'
        local herm0 : word 2 of `hermite'
        if ("`herm0'" == "") local herm0 = `herm1'
    }
    if ((`herm1' > 0 | `herm0' > 0) & "`method'" != "twostep") {
        di as err "hermite() is an option of method(twostep)"
        exit 198
    }
    if (!inlist(`herm1', 0, 2, 3) | !inlist(`herm0', 0, 2, 3)) {
        di as err "hermite(#) or hermite(#1 #0): the order of the Hermite expansion of E[w_j | u], 2 or 3,"
        di as err "in both regimes, or in the treated (#1) and the untreated (#0) regime, 0 for none"
        exit 198
    }
    * ph1, ph0 Hermite controls beside the Mills ratio in each regime:
    * E[u^2 - 1 | D, Z] (order 2), and E[u^3 - 3u | D, Z] (order 3)
    local ph1 = cond(`herm1' > 0, `herm1' - 1, 0)
    local ph0 = cond(`herm0' > 0, `herm0' - 1, 0)
    * ---- vce: method(fiml) hands it to ml; method(twostep) aggregates the stacked
    *      moments by observation (default; robust says the same), by cluster, or
    *      by the survey design (vce(svy): the strata, PSUs, fpc and weight of svyset)
    local vce = strtrim(`"`vce'"')
    gettoken vt vrest : vce
    local vt = lower("`vt'")
    local vrest = strtrim(`"`vrest'"')
    local clustvar ""
    local vsvy ""
    if ("`vt'" == "") {
    }
    else if ("`vt'" == substr("cluster", 1, max(2, strlen("`vt'")))) {
        if (`: word count `vrest'' != 1) {
            di as err "vce(cluster clustvar): one cluster variable"
            exit 198
        }
        unab clustvar : `vrest', max(1)
    }
    else if ("`vt'" == "svy" & `"`vrest'"' == "") {
        if ("`method'" == "fiml") {
            di as err "vce(svy) is an option of method(twostep); with method(fiml), use the svy prefix: {bf:svy: esreg ...}"
            exit 198
        }
        local vsvy "svy"
    }
    else if ("`vt'" == substr("robust", 1, max(1, strlen("`vt'"))) & `"`vrest'"' == "") {
    }
    else if (inlist("`vt'", "oim", "opg") & `"`vrest'"' == "") {
        if ("`method'" == "twostep") {
            di as err "vce(`vt') is an option of method(fiml); the two-step variance is the sandwich of the stacked"
            di as err "moments (default), by cluster with vce(cluster clustvar), by the survey design with vce(svy)"
            exit 198
        }
    }
    else if ("`method'" == "twostep") {
        di as err "vce(`vce') not allowed with method(twostep): use robust (default), cluster clustvar or svy"
        exit 198
    }
    local vceml ""
    if ("`method'" == "fiml" & `"`vce'"' != "") local vceml `"vce(`vce')"'
    marksample touse
    markout `touse' `dvar'
    if ("`clustvar'" != "") markout `touse' `clustvar', strok
    foreach o in zvars hetsigma hetrho kappa {
        if ("``o''" != "") {
            fvrevar ``o'', list
            markout `touse' `r(varlist)'
        }
    }
    * ---- vce(svy) (two-step): the design and the weight of svyset ------------------
    if ("`vsvy'" != "") {
        if ("`weight'" != "") {
            di as err "vce(svy): the sampling weight comes from svyset; do not specify weights"
            exit 198
        }
        if ("`svyset'" != "") {
            di as err "vce(svy) uses the svyset design: nosvyset is not allowed with it"
            exit 198
        }
        qui svyset
        if ("`r(settings)'" == ", clear" | "`r(settings)'" == "") {
            di as err "vce(svy): the data are not svyset"
            exit 119
        }
        if ("`r(wvar)'" != "" & "`r(wtype)'" != "pweight") {
            di as err "vce(svy): the svyset weight must be a pweight"
            exit 198
        }
        if (`"`r(poststrata)'`r(rake)'`r(regress)'"' != "") {
            di as err "vce(svy): poststratification and calibration are not supported"
            exit 198
        }
        local svyvars "`r(wvar)' `r(strata1)' `r(fpc1)'"
        if ("`r(su1)'" != "_n") local svyvars "`svyvars' `r(su1)'"
        markout `touse' `svyvars', strok
        if ("`r(wvar)'" != "") {
            local weight "pweight"
            local exp "= `r(wvar)'"
        }
    }
    * ---- weights: svyset pweight taken by default when none is given (nosvyset
    *      declines it); hsize() multiplies the weight (effects per individual when
    *      the unit is the household); under svy: (iweights) hsize() is not allowed
    if ("`weight'" == "" & "`svyset'" == "") {
        cap qui svyset
        if (_rc == 0 & "`r(wvar)'" != "") {
            local weight "pweight"
            local exp "= `r(wvar)'"
            markout `touse' `r(wvar)'
            local hint = cond("`method'" == "fiml", "the svy: prefix", "vce(svy)")
            di as txt "(svyset weight " as res "`r(wvar)'" as txt " used as pweight; design-based standard errors need `hint')"
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
    * counts on the final sample (after the design variables and hsize())
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
            `wexp' if `touse', `vceml' `iterate' `difficult' maximize init(`b0mat', copy) ///
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
        * Hermite controls E[u^2 - 1 | D, Z] (order 2), and E[u^3 - 3u | D, Z] (order 3),
        * in each regime as its order says
        local hn1 ""
        local hn0 ""
        if (`ph1' >= 1) {
            tempvar h2_1
            qui gen double `h2_1' = -`zg' * `l1' if `touse'
            local lam1 "`lam1' `h2_1'"
            local hn1 "h2"
        }
        if (`ph1' >= 2) {
            tempvar h3_1
            qui gen double `h3_1' = (`zg'^2 - 1) * `l1' if `touse'
            local lam1 "`lam1' `h3_1'"
            local hn1 "h2 h3"
        }
        if (`ph0' >= 1) {
            tempvar h2_0
            qui gen double `h2_0' = `zg' * `l0' if `touse'
            local lam0 "`lam0' `h2_0'"
            local hn0 "h2"
        }
        if (`ph0' >= 2) {
            tempvar h3_0
            qui gen double `h3_0' = -(`zg'^2 - 1) * `l0' if `touse'
            local lam0 "`lam0' `h3_0'"
            local hn0 "h2 h3"
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
        foreach v of local hn1 {
            local cn `cn' `y'_1:`v'
        }
        local cn `cn' `y'_1:_cons
        foreach v of local xlist {
            local cn `cn' `y'_0:`v'
        }
        foreach v of local knames {
            local cn `cn' `y'_0:`v'
        }
        foreach v of local hn0 {
            local cn `cn' `y'_0:`v'
        }
        local cn `cn' `y'_0:_cons
        * reorder c1, c0 so that lambda terms come before _cons? No: keep Stata's
        * regress order (x..., lambda..., _cons), which the Mata code expects.
        mat colnames `b' = `cn'
        * stacked-moment variance in Mata: by observation, by cluster (sums of the
        * influence functions within the clusters), or by the survey design (svy
        * linearized total of the influence functions)
        local cid ""
        if ("`clustvar'" != "") {
            tempvar cid
            qui egen long `cid' = group(`clustvar') if `touse'
        }
        local psiv ""
        if ("`vsvy'" != "") {
            local np = colsof(`b')
            forvalues j = 1/`np' {
                tempvar psi`j'
                local psiv "`psiv' `psi`j''"
            }
        }
        mata: _esreg_ts_var("`y'", "`xlist_m'", "`zlist_m'", "`dvar'", "`kappalist_m'", "`wv'", "`touse'", "`b'", "`V'", (`ph1', `ph0'), "`cid'", "`psiv'")
        if ("`vsvy'" != "") {
            _esreg_ifcov `psiv' if `touse', svy
            mat `V' = r(V)
            local svy_df = r(df_r)
            local svy_ns = r(N_strata)
            local svy_np = r(N_psu)
            mata: st_local("vmiss", strofreal(hasmissing(st_matrix("`V'"))))
            if (`vmiss') {
                di as err "vce(svy): the design variance is missing (a stratum with a single PSU in the estimation"
                di as err "sample?); see the singleunit() option of svyset"
                exit 459
            }
        }
        mat rownames `V' = `cn'
        mat colnames `V' = `cn'
        if ("`vsvy'" != "") ereturn post `b' `V' `wexp', depname(`y') obs(`N') esample(`touse') dof(`svy_df')
        else                ereturn post `b' `V' `wexp', depname(`y') obs(`N') esample(`touse')
        ereturn local method  "twostep"
        ereturn local title   "Endogenous switching regression -- two-step"
        ereturn local eqnames "`dvar' `y'_1 `y'_0"
        if ("`clustvar'" != "") {
            ereturn local  vce      "cluster"
            ereturn local  vcetype  "Robust"
            ereturn local  clustvar "`clustvar'"
            ereturn scalar N_clust  = `esr_ncl'
        }
        else if ("`vsvy'" != "") {
            ereturn local  vce      "svy"
            ereturn local  vcetype  "Linearized"
            ereturn scalar N_strata = `svy_ns'
            ereturn scalar N_psu    = `svy_np'
        }
        else {
            ereturn local vce     "stacked"
            ereturn local vcetype ""
        }
        ereturn scalar converged = 1
    }

    * ---- common e() ---------------------------------------------------------
    ereturn local cmd       "esreg"
    ereturn local cmdline   `"esreg `cmdline0'"'
    ereturn local version   "1.0.0"
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
    ereturn local eff_vce   "none"
    ereturn local predict   "esreg_p"
    ereturn scalar N_treated   = `N1'
    ereturn scalar N_untreated = `N0'
    ereturn scalar sum_w       = `sumw'
    ereturn scalar k_x   = `: word count `xlist'' + 1
    ereturn scalar k_z   = `: word count `zlist'' + 1
    ereturn scalar k_hs  = `: word count `hetsigmalist'' + 1
    ereturn scalar k_hr  = `: word count `hetrholist'' + 1
    ereturn scalar k_kap = `: word count `kappalist'' + 1
    ereturn scalar k_h1  = `ph1'
    ereturn scalar k_h0  = `ph0'
    ereturn scalar k_h   = max(`ph1', `ph0')
    ereturn scalar hermite1 = `herm1'
    ereturn scalar hermite0 = `herm0'
    ereturn scalar hermite  = max(`herm1', `herm0')
    * coefficient blocks by equation, for the post-commands
    tempname bb bs b1m b0m
    mat `bb' = e(b)
    mat `bs' = `bb'[1, "`dvar':"]
    mat `b1m' = `bb'[1, "`y'_1:"]
    mat `b0m' = `bb'[1, "`y'_0:"]
    ereturn matrix b_sel = `bs'
    ereturn matrix b_1   = `b1m'
    ereturn matrix b_0   = `b0m'

    * ---- effects (_esreg_effects_post: the influence function of the whole
    *      procedure, aggregated as the estimation was) -----------------------------
    if ("`effects'" == "") {
        _esreg_effects_post, plain xl(`xlist_m') zl(`zlist_m') hsl(`hetsigmalist_m') ///
            hrl(`hetrholist_m') kl(`kappalist_m') wvar(`wv')
    }
    * ---- keep the estimation in memory: always as _esreg, and under store() ----
    ereturn local store "`store'"
    if ("`store2'" == "") qui estimates store _esreg
    if ("`store'" != "") qui estimates store `store'
    _esreg_display, level(`level')
end

* ---------------------------------------------------------------------------
* a bare hermite (or an abbreviation of it) among the options becomes hermite(3)
program define _esreg_hermopt, rclass
    _parse comma lhs rhs : 0
    if (`"`rhs'"' == "") {
        return local cmd `"`0'"'
        exit
    }
    gettoken comma rhs : rhs, parse(",")
    local out ""
    while (`"`rhs'"' != "") {
        gettoken tok rhs : rhs, bind
        if inlist(lower(`"`tok'"'), "her", "herm", "hermi", "hermit", "hermite") local tok "hermite(3)"
        local out `"`out' `tok'"'
    }
    return local cmd `"`lhs', `out'"'
end

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
        if ("`e(vce)'" == "svy") {
            di as txt _col(49) "Strata        = " as res %9.0f e(N_strata)
            di as txt _col(49) "PSUs          = " as res %9.0f e(N_psu)
            di as txt _col(49) "Design df     = " as res %9.0f e(df_r)
            di as txt "Std. err.: exact two-step variance, linearized by the survey design"
        }
        else if ("`e(vce)'" == "cluster") {
            di as txt _col(49) "Clusters      = " as res %9.0f e(N_clust)
            di as txt "Std. err.: exact two-step variance (stacked moments), by cluster of `e(clustvar)'"
        }
        else di as txt "Std. err.: exact variance of the two-step procedure (stacked moments)"
        ereturn display, level(`level')
    }
    if (e(sigma1) < .) {
        di as txt "sigma_1 = " as res %8.6g e(sigma1) as txt "   sigma_0 = " as res %8.6g e(sigma0) ///
           as txt "   rho_1 = " as res %8.6g e(rho1) as txt "   rho_0 = " as res %8.6g e(rho0) ///
           cond("`e(hetsigma)'`e(hetrho)'" != "", "   (means; heterogeneous)", "")
    }
    cap confirm matrix e(effects)
    if (_rc) exit
    tempname E
    mat `E' = e(effects)
    * t with the design degrees of freedom (survey design), z otherwise
    local df = .
    if ("`e(eff_vce)'" == "svy" & e(df_r) < .) local df = e(df_r)
    if (`df' < .) {
        local crit = invttail(`df', (1 - `level'/100)/2)
        local st "t"
    }
    else {
        local crit = invnormal(1 - (1 - `level'/100)/2)
        local st "z"
    }
    di
    di as txt "Treatment effects (std. err.: influence function of the whole procedure)"
    di as txt "{hline 76}"
    di as txt "Index    {c |}   Estimate   Std. err.        `st'    P>|`st'|     [`level'% conf. interval]"
    di as txt "{hline 9}{c +}{hline 66}"
    foreach r in ATT ATU ATE kappa {
        local i = rownumb(`E', "`r'")
        local b = `E'[`i', 1]
        local s = `E'[`i', 2]
        if (`df' < .) local p = 2*ttail(`df', abs(`b'/`s'))
        else          local p = 2*normal(-abs(`b'/`s'))
        di as txt %8s "`r'" " {c |}" as res %11.6g `b' "  " %10.6g `s' "  " %7.2f `b'/`s' ///
           "  " %6.3f `p' "   " %10.6g `b' - `crit'*`s' "  " %10.6g `b' + `crit'*`s'
    }
    di as txt "{hline 76}"
    if (e(k_h) > 0 & e(k_h) < .) {
        local hspec = cond(e(hermite1) == e(hermite0), "hermite(" + strofreal(e(hermite1)) + ")", ///
                           "hermite(" + strofreal(e(hermite1)) + " " + strofreal(e(hermite0)) + ")")
        di as txt "`hspec': E[w1 - w0 | u] = kappa*u + dh2*(u^2-1)" cond(e(k_h) >= 2, " + dh3*(u^3-3u)", "") ///
           ", dh_k = h_1k - h_0k;"
        if (e(hermite1) != e(hermite0)) {
            di as txt "Hermite terms of order " e(hermite1) " in the treated regime, " e(hermite0) ///
               " in the untreated one (h_jk = 0 for a term a regime does not have);"
        }
        di as txt "kappa = Cov(w1 - w0, u), the linear part; the MTE curve: esrmte"
    }
    else di as txt "kappa = rho1*sigma1 - rho0*sigma0;  ATT - ATU = kappa*(mean lambda1 + mean lambda0)"
    if (colsof(`E') >= 5) {
        local agg = cond("`e(eff_vce)'" == "svy", "survey design", cond("`e(eff_vce)'" == "cluster", "by cluster", "by observation"))
        di as txt "Variance of ATT: parameters " as res %8.3g `E'[1,3] as txt " + sampling " as res %8.3g `E'[1,4] ///
           as txt " + covariance " as res %8.3g `E'[1,5] as txt " (`agg')"
    }
    else di as txt "Variance components (ATT): parameter " as res %8.3g `E'[1,3] as txt "   sampling " as res %8.3g `E'[1,4]
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

