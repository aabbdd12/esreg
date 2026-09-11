*! esrtest 0.5.1  11sep2026  A. Araar
*! Specification tests after esreg (Section 4 of the esreg paper): spec (4.1),
*! suff (4.2), kappa (4.3), normal (4.5), pdid (4.6); strength and support are
*! in esrdiag (4.4).
*!
*!   esrtest [, KAPpa SUFF [KEEP(varname)] SPEC [ORDer(#) ADD(varlist)] NORMal PDID [NQ(#)]
*!               EST(name) Level(#)]
*!
*! kappa : constant selection on gains, H0: kappa(x) = kappa.
*!   after method(twostep) with kappa(W):  kappa(x) = W(t1 - t0); H0: the non-
*!       constant entries of t1 - t0 are zero (exact null); Wald on the stacked
*!       two-step variance; the table gives each slope d kappa / d w with its s.e.
*!   after method(fiml) with hetsigma()/hetrho():  H0: the non-constant
*!       coefficients of ln sigma_j and atanh rho_j are all zero (homogeneous
*!       regime laws), a sufficient condition for kappa(x) constant.
*! suff  : index sufficiency, H0: E[omega_j | X, Z, u] = E[omega_j | X, u].  The
*!   excluded instruments other than keep() (default: the strongest in the probit)
*!   are added to the two regime regressions of the two-step route; under H0
*!   their coefficients are zero.  Wald on the stacked variance, joint and by
*!   regime.  With a single excluded instrument the test rests on the curvature
*!   of lambda alone (a warning is printed).
*! spec  : specification of the selection equation (Section 4.1), two pieces.  (a) link
*!   test: probit of D on v = Z g and v^2..v^order (default order 2), Wald on the
*!   powers (H0: the index is linear in Z and the link is the normal cdf); with
*!   add(varlist) instead, a RESET-type test: the terms are added to the probit,
*!   Wald and LR on their coefficients.  (b) gamma contrast: probit against FIML
*!   on the selection coefficients, the gamma block of the Hausman contrast of
*!   normal, covariance from the influence functions (H0: A3; also rejects when
*!   the index is misspecified).  Both re-estimate on the sample of e().
*! normal: normality by regime (Section 4.5), two pieces.  (a) Hausman contrast of
*!   the FIML and two-step estimates of q = (gamma, beta_1, rho1 sigma1, beta_0,
*!   rho0 sigma0), covariance of the difference from the influence functions of
*!   the two estimators (A3: trivariate normality), with the kappa contrast
*!   kappa_2S - kappa_FIML; (b) Hermite controls E[u^2-1 | D,Z], E[u^3-3u | D,Z]
*!   added beside the Mills ratio in each regime, Wald joint and by regime
*!   (A2: E[omega_j | u] linear in u).  Both re-estimate on the sample of e().
*! pdid  : pseudo-DiD on the selection index (Section 4.6).  nq() strata of the
*!   probit index (default 2); within each stratum the X-adjusted treated-untreated
*!   gap; the within-treated contrast C1 (top vs bottom stratum) has expectation
*!   rho1 sigma1 (dlam1 < 0), the within-untreated contrast C0 has expectation
*!   -rho0 sigma0 dlam0 (dlam0 > 0); their signs are model-free under A2; kappa_dd =
*!   C1/dlam1 + C0/dlam0 needs the Mills factors (A1) but no trivariate normality.
*!   First-step variation of gamma ignored in the standard errors.
*! Returns r(chi2) r(df) r(p) r(test) and, after twostep, r(slopes) (est, se, z, p);
*! for spec also r(link) r(gamma) r(chi2_g) r(df_g) r(p_g) r(lr) r(p_lr) r(tested);
*! for suff also r(chi2_1) r(df_1) r(p_1) r(chi2_0) r(df_0) r(p_0) r(keep) r(tested) r(coefs).
cap program drop esrtest
program define esrtest, rclass
    version 16
    syntax [, KAPpa SUFF KEEP(varname) SPEC ORDer(integer 2) ADD(varlist numeric fv) NORMal PDID NQ(integer 2) EST(name) Level(cilevel)]
    esreg_engine
    tempname hold
    _estimates hold `hold', restore copy
    _esreg_getest `est'
    if ("`kappa'`suff'`spec'`normal'`pdid'" == "") local kappa kappa
    local y      "`e(depvar)'"
    local method "`e(method)'"

    * =========================================================================
    if ("`kappa'" != "") {
        di
        di as txt "Test of constant selection on gains:  H0: kappa(x) = kappa"
        if ("`method'" == "twostep") {
            local kv "`e(kappavars)'"
            if ("`kv'" == "") {
                di as err "the estimation has no kappa() variables: re-estimate with esreg ..., method(twostep) kappa(varlist)"
                exit 198
            }
            di as txt "(two-step estimation, kappa(x) = W(t1 - t0); Wald on the non-constant entries of t1 - t0)"
            local nk : word count `kv'
            tempname S
            mat `S' = J(`nk', 4, .)
            local tl ""
            local rn ""
            di as txt "{hline 68}"
            di as txt "  d kappa / d w        {c |}      Coef.   Std. err.        z    P>|z|"
            di as txt "{hline 24}{c +}{hline 43}"
            local i = 0
            foreach v of local kv {
                local ++i
                local vn = subinstr("`v'", ".", "_", .)
                qui lincom [`y'_1]lambda_`vn' - [`y'_0]lambda_`vn'
                mat `S'[`i',1] = r(estimate)
                mat `S'[`i',2] = r(se)
                mat `S'[`i',3] = r(estimate)/r(se)
                mat `S'[`i',4] = 2*normal(-abs(r(estimate)/r(se)))
                di as txt %22s "`v'" "  {c |}" as res %11.5g `S'[`i',1] %12.5g `S'[`i',2] ///
                   %9.2f `S'[`i',3] %9.3f `S'[`i',4]
                local tl "`tl' ([`y'_1]lambda_`vn' = [`y'_0]lambda_`vn')"
                local rn "`rn' `vn'"
            }
            di as txt "{hline 68}"
            qui test `tl'
            mat colnames `S' = est se z p
            mat rownames `S' = `rn'
            return matrix slopes = `S'
        }
        else {
            local hs "`e(hetsigma)'"
            local hr "`e(hetrho)'"
            if ("`hs'`hr'" == "") {
                di as err "the estimation has no hetsigma()/hetrho() variables: re-estimate with esreg ..., hetsigma(varlist) hetrho(varlist)"
                exit 198
            }
            di as txt "(fiml estimation; Wald that the non-constant coefficients of ln sigma_j and atanh rho_j are zero:"
            di as txt " homogeneous regime laws, a sufficient condition for kappa(x) constant)"
            local tl ""
            foreach v of local hs {
                local tl "`tl' [lnsigma_1]`v' [lnsigma_0]`v'"
            }
            foreach v of local hr {
                local tl "`tl' [atanhrho_1]`v' [atanhrho_0]`v'"
            }
            qui test `tl'
        }
        local chi2 = r(chi2)
        local df   = r(df)
        local p    = r(p)
        di as txt "  Wald chi2(" as res `df' as txt ") = " as res %8.2f `chi2' ///
           as txt "     Prob > chi2 = " as res %6.4f `p'
        di as txt "  kappa at the sample mean (esreg): " as res %8.5g e(kappa) as txt " (s.e. " as res %7.4g e(se_kappa) as txt ")"
        return scalar chi2 = `chi2'
        return scalar df   = `df'
        return scalar p    = `p'
        return local  test "kappa"
    }

    * =========================================================================
    if ("`suff'" != "") {
        local dv    "`e(treat)'"
        local xv    "`e(xvars)'"
        local zv    "`e(zvars)'"
        local wtype "`e(wtype)'"
        local wexp0 "`e(wexp)'"
        local wexp ""
        if ("`wtype'" != "") local wexp "[`wtype' `wexp0']"
        tempvar smp
        qui gen byte `smp' = e(sample)
        local excl : list zv - xv
        if ("`excl'" == "") {
            di as err "no excluded instrument: the selection equation has no variable outside the outcome equations"
            exit 198
        }
        di
        di as txt "Test of index sufficiency:  H0: E[omega_j | X, Z, u] = E[omega_j | X, u]"
        * plain variables behind factor terms (2.region 3.region would be re-merged
        * with 2 as base by probit and by the internal esreg call); same order as e()
        _esreg_data if `smp'
        local xv_m "`r(x)'"
        local zv_m "`r(z)'"
        local excl_m : list zv_m - xv_m
        * the instrument kept excluded: keep() or the strongest in the probit
        if ("`keep'" == "") {
            qui probit `dv' `zv_m' `wexp' if `smp'
            local best = -1
            local i = 0
            foreach v of local excl_m {
                local ++i
                cap local st = (_b[`v'] / _se[`v'])^2
                if (_rc) continue
                if (`st' > `best') {
                    local best = `st'
                    local keep : word `i' of `excl'
                }
            }
        }
        else {
            local ok : list keep in excl
            if (!`ok') {
                di as err "keep(`keep') is not an excluded instrument (excluded: `excl')"
                exit 198
            }
        }
        local tested : list excl - keep
        if ("`tested'" == "") {
            di as txt "note: a single excluded instrument (`keep'); it is tested itself, so that"
            di as txt "      identification rests on the curvature of lambda alone -- read with care"
            local tested "`keep'"
            local keep ""
        }
        di as txt "(two-step re-estimation with " as res "`tested'" as txt " added to both regime equations; kept excluded: " as res cond("`keep'" == "", "none", "`keep'") as txt ")"
        local tested_m ""
        foreach v of local tested {
            local pos : list posof "`v'" in zv
            local tested_m "`tested_m' `: word `pos' of `zv_m''"
        }
        qui esreg `y' `xv_m' `tested_m' `wexp' if `smp', select(`dv' = `zv_m') method(twostep) nostore2
        tempname C
        local nt : word count `tested'
        mat `C' = J(2*`nt', 4, .)
        local tl1 ""
        local tl0 ""
        local rn ""
        di as txt "{hline 68}"
        di as txt "  added instrument      {c |}      Coef.   Std. err.        z    P>|z|"
        di as txt "{hline 24}{c +}{hline 43}"
        local i = 0
        foreach j in 1 0 {
            di as txt "  regime `j' (`y'_`j')"
            foreach v of local tested {
                local ++i
                local pos : list posof "`v'" in tested
                local vm : word `pos' of `tested_m'
                mat `C'[`i',1] = _b[`y'_`j':`vm']
                mat `C'[`i',2] = _se[`y'_`j':`vm']
                mat `C'[`i',3] = `C'[`i',1] / `C'[`i',2]
                mat `C'[`i',4] = 2*normal(-abs(`C'[`i',3]))
                di as txt %22s "`v'" "  {c |}" as res %11.5g `C'[`i',1] %12.5g `C'[`i',2] ///
                   %9.2f `C'[`i',3] %9.3f `C'[`i',4]
                local tl`j' "`tl`j'' [`y'_`j']`vm'"
                local vn = subinstr("`v'", ".", "_", .)
                local rn "`rn' r`j'_`vn'"
            }
        }
        di as txt "{hline 68}"
        qui test `tl1' `tl0'
        local chi2 = r(chi2)
        local df   = r(df)
        local p    = r(p)
        qui test `tl1'
        local chi21 = r(chi2)
        local df1   = r(df)
        local p1    = r(p)
        qui test `tl0'
        local chi20 = r(chi2)
        local df0   = r(df)
        local p0    = r(p)
        di as txt "  joint      Wald chi2(" as res `df'  as txt ") = " as res %8.2f `chi2'  as txt "   Prob > chi2 = " as res %6.4f `p'
        di as txt "  regime 1   Wald chi2(" as res `df1' as txt ") = " as res %8.2f `chi21' as txt "   Prob > chi2 = " as res %6.4f `p1'
        di as txt "  regime 0   Wald chi2(" as res `df0' as txt ") = " as res %8.2f `chi20' as txt "   Prob > chi2 = " as res %6.4f `p0'
        di as txt "A rejection says that the added instrument(s) carry information on the regime errors beyond the"
        di as txt "selection index: a direct effect on the outcome (exclusion) or a second index in participation."
        mat colnames `C' = est se z p
        mat rownames `C' = `rn'
        return matrix coefs = `C'
        return scalar chi2   = `chi2'
        return scalar df     = `df'
        return scalar p      = `p'
        return scalar chi2_1 = `chi21'
        return scalar df_1   = `df1'
        return scalar p_1    = `p1'
        return scalar chi2_0 = `chi20'
        return scalar df_0   = `df0'
        return scalar p_0    = `p0'
        cap drop __esr_fv*
        return local  keep   "`keep'"
        return local  tested "`tested'"
        return local  test   "suff"
    }

    * =========================================================================
    if ("`spec'" != "") {
        if (`order' < 2) {
            di as err "order() must be at least 2"
            exit 198
        }
        local dv    "`e(treat)'"
        local xv    "`e(xvars)'"
        local zv    "`e(zvars)'"
        local wtype "`e(wtype)'"
        local wexp0 "`e(wexp)'"
        local wexp ""
        if ("`wtype'" != "") local wexp "[`wtype' `wexp0']"
        tempvar smp w v
        qui gen byte `smp' = e(sample)
        if ("`wtype'" != "") qui gen double `w' `wexp0' if `smp'
        else                 qui gen double `w' = 1 if `smp'
        _esreg_data if `smp'
        local xl "`r(x)'"
        local zl "`r(z)'"
        di
        di as txt "Specification of the selection equation"
        * ---- (a) link test of the probit index, or RESET-type test of added terms ------
        qui probit `dv' `zl' `wexp' if `smp'
        local ll0 = e(ll)
        qui predict double `v' if `smp', xb
        local lr = .
        local plr = .
        if ("`add'" == "") {
            local tested ""
            local tlab ""
            forvalues j = 2/`order' {
                tempvar v`j'
                qui gen double `v`j'' = `v'^`j' if `smp'
                local tested "`tested' `v`j''"
                local tlab "`tlab' v^`j'"
            }
            qui probit `dv' `v' `tested' `wexp' if `smp'
            local vb  = _b[`v']
            local vse = _se[`v']
            di as txt "(a) Link test: probit of `dv' on v = Z g and its powers 2..`order' (H0: the powers are zero)"
        }
        else {
            local tested ""
            local tlab ""
            foreach t of local add {
                cap confirm variable `t'
                if (_rc == 0) {
                    local tested "`tested' `t'"
                    local tlab "`tlab' `t'"
                    continue
                }
                fvexpand `t' if `smp'
                foreach tt in `r(varlist)' {
                    if (strpos("`tt'", "b.") > 0 | strpos("`tt'", "o.") > 0) continue
                    tempvar tv
                    _esreg_fvone `tt' `tv' if `smp'
                    local tested "`tested' `tv'"
                    local tlab "`tlab' `tt'"
                }
            }
            qui probit `dv' `zl' `tested' `wexp' if `smp'
            local lr  = 2 * (e(ll) - `ll0')
            local nadd : word count `tested'
            local plr = chi2tail(`nadd', `lr')
            di as txt "(a) Added terms in the probit of `dv' (H0: their coefficients are zero)"
        }
        tempname L
        local nt : word count `tested'
        mat `L' = J(`nt', 4, .)
        di as txt "{hline 66}"
        di as txt "  term                 {c |}      Coef.   Std. err.        z    P>|z|"
        di as txt "{hline 22}{c +}{hline 43}"
        if ("`add'" == "") di as txt %20s "v" "  {c |}" as res %11.5g `vb' %12.5g `vse' %9.2f `vb'/`vse' %9.3f 2*normal(-abs(`vb'/`vse'))
        local i = 0
        foreach t of local tested {
            local ++i
            mat `L'[`i',1] = _b[`t']
            mat `L'[`i',2] = _se[`t']
            mat `L'[`i',3] = `L'[`i',1] / `L'[`i',2]
            mat `L'[`i',4] = 2*normal(-abs(`L'[`i',3]))
            local lab : word `i' of `tlab'
            di as txt %20s "`lab'" "  {c |}" as res %11.5g `L'[`i',1] %12.5g `L'[`i',2] %9.2f `L'[`i',3] %9.3f `L'[`i',4]
        }
        di as txt "{hline 66}"
        qui test `tested'
        local lchi2 = r(chi2)
        local ldf   = r(df)
        local lp    = r(p)
        di as txt "  Wald chi2(" as res `ldf' as txt ") = " as res %8.2f `lchi2' as txt "   Prob > chi2 = " as res %6.4f `lp'
        if ("`add'" != "") di as txt "  LR   chi2(" as res `ldf' as txt ") = " as res %8.2f `lr' as txt "   Prob > chi2 = " as res %6.4f `plr'
        mat colnames `L' = est se z p
        local rn = subinstr("`tlab'", ".", "_", .)
        local rn = subinstr("`rn'", "^", "", .)
        local rn = subinstr("`rn'", "#", "_", .)
        mat rownames `L' = `rn'
        * ---- (b) gamma: probit against FIML (the gamma block of the Hausman contrast) -----
        tempname bF VF bS R Q VD bg Vg
        qui esreg `y' `xl' `wexp' if `smp', select(`dv' = `zl') nostore2 noeffects nolog
        mat `bF' = e(b)
        mat `VF' = e(V)
        qui esreg `y' `xl' `wexp' if `smp', select(`dv' = `zl') method(twostep) nostore2 noeffects
        mat `bS' = e(b)
        mata: _esreg_hausman("`bF'", "`VF'", "`bS'", "`y'", "`xl'", "`zl'", "`dv'", "`w'", "`smp'", "`R'", "`Q'", "`VD'")
        cap drop __esr_fv*
        local m : word count `zv'
        local mg = `m' + 1
        mat `bg' = (`Q'[1..`mg', 1] - `Q'[1..`mg', 2])'
        mat `Vg' = `VD'[1..`mg', 1..`mg']
        tempname Vgi
        mat `Vgi' = syminv(`Vg')
        local gdf = `mg' - diag0cnt(`Vgi')
        tempname Wg
        mat `Wg' = `bg' * `Vgi' * `bg''
        local gchi2 = `Wg'[1,1]
        local gp    = chi2tail(`gdf', `gchi2')
        di
        di as txt "(b) gamma: probit against FIML (H0: trivariate normality, A3; the gamma block of esrtest, normal)"
        di as txt "{hline 70}"
        di as txt "                          {c |}     probit        FIML   difference   Std. err."
        di as txt "{hline 26}{c +}{hline 43}"
        local qn ""
        foreach vv of local zv {
            local qn "`qn' `vv'"
        }
        local qn "`qn' _cons"
        forvalues i = 1/`mg' {
            local nm : word `i' of `qn'
            di as txt %24s "`nm'" "  {c |}" as res %11.4f `Q'[`i',1] %12.4f `Q'[`i',2] ///
               %12.4f `Q'[`i',1] - `Q'[`i',2] %12.4f `Q'[`i',3]
        }
        di as txt "{hline 70}"
        di as txt "  Wald chi2(" as res `gdf' as txt ") = " as res %8.2f `gchi2' as txt "   Prob > chi2 = " as res %6.4f `gp'
        di as txt "Reading: (a) rejects when the index is not linear in Z (missing polynomials or"
        di as txt "interactions) or the link is not the normal cdf; both routes then use the wrong"
        di as txt "Mills ratios.  (b) rejects when the outcomes, read through trivariate normality,"
        di as txt "pull gamma away from the probit: A3 fails, or the index of (a) is misspecified."
        tempname G
        mat `G' = `Q'[1..`mg', 1..3]
        mat colnames `G' = probit fiml se_diff
        mat rownames `G' = `qn'
        return matrix gamma = `G'
        return matrix link  = `L'
        return scalar chi2     = `lchi2'
        return scalar df       = `ldf'
        return scalar p        = `lp'
        return scalar lr       = `lr'
        return scalar p_lr     = `plr'
        return scalar chi2_g   = `gchi2'
        return scalar df_g     = `gdf'
        return scalar p_g      = `gp'
        return local  tested   "`tlab'"
        return local  test     "spec"
    }

    * =========================================================================
    if ("`normal'" != "") {
        local dv    "`e(treat)'"
        local xv    "`e(xvars)'"
        local zv    "`e(zvars)'"
        local wtype "`e(wtype)'"
        local wexp0 "`e(wexp)'"
        local wexp ""
        if ("`wtype'" != "") local wexp "[`wtype' `wexp0']"
        tempvar smp w
        qui gen byte `smp' = e(sample)
        if ("`wtype'" != "") qui gen double `w' `wexp0' if `smp'
        else                 qui gen double `w' = 1 if `smp'
        di
        di as txt "Normality by regime"
        * ---- (a) Hausman: FIML against two-step on the common parameters --------------
        tempname bF VF bS R Q
        * plain variables behind factor terms, for the internal esreg calls and for Mata
        _esreg_data if `smp'
        local xl "`r(x)'"
        local zl "`r(z)'"
        qui esreg `y' `xl' `wexp' if `smp', select(`dv' = `zl') nostore2 noeffects nolog
        mat `bF' = e(b)
        mat `VF' = e(V)
        qui esreg `y' `xl' `wexp' if `smp', select(`dv' = `zl') method(twostep) nostore2 noeffects
        mat `bS' = e(b)
        mata: _esreg_hausman("`bF'", "`VF'", "`bS'", "`y'", "`xl'", "`zl'", "`dv'", "`w'", "`smp'", "`R'", "`Q'")
        local nq = rowsof(`Q')
        local k : word count `xv'
        local m : word count `zv'
        local qn ""
        foreach v of local zv {
            local qn "`qn' `dv':`v'"
        }
        local qn "`qn' `dv':_cons"
        foreach j in 1 0 {
            foreach v of local xv {
                local qn "`qn' `y'_`j':`v'"
            }
            local qn "`qn' `y'_`j':_cons `y'_`j':rho`j'sigma`j'"
        }
        di as txt "(a) Hausman contrast, FIML vs two-step (H0: trivariate normality, A3)"
        di as txt "{hline 70}"
        di as txt "                          {c |}   two-step        FIML   difference   Std. err."
        di as txt "{hline 26}{c +}{hline 43}"
        forvalues i = 1/`nq' {
            local nm : word `i' of `qn'
            di as txt %24s "`nm'" "  {c |}" as res %11.4f `Q'[`i',1] %12.4f `Q'[`i',2] ///
               %12.4f `Q'[`i',1] - `Q'[`i',2] %12.4f `Q'[`i',3]
        }
        di as txt "{hline 70}"
        di as txt "  Hausman chi2(" as res `R'[1,2] as txt ") = " as res %8.2f `R'[1,1] ///
           as txt "   Prob > chi2 = " as res %6.4f `R'[1,3]
        di as txt "  kappa: two-step " as res %8.5g `R'[1,4] as txt "   FIML " as res %8.5g `R'[1,5] ///
           as txt "   difference " as res %8.5g `R'[1,6] as txt " (s.e. " as res %7.4g `R'[1,7] ///
           as txt ", z = " as res %6.2f `R'[1,6]/`R'[1,7] as txt ")"
        di as txt "  Std. err. of the differences from the influence functions of the two estimators."
        * ---- (b) Hermite controls: linearity of E[omega_j | u] ---------------------------
        qui esreg `y' `xl' `wexp' if `smp', select(`dv' = `zl') method(twostep) hermite nostore2
        cap drop __esr_fv*
        tempname H
        mat `H' = J(4, 4, .)
        local i = 0
        di
        di as txt "(b) Hermite controls beside the Mills ratio (H0: E[omega_j | u] linear in u, A2)"
        di as txt "{hline 68}"
        di as txt "  control                {c |}      Coef.   Std. err.        z    P>|z|"
        di as txt "{hline 24}{c +}{hline 43}"
        foreach j in 1 0 {
            di as txt "  regime `j' (`y'_`j')"
            foreach v in h2 h3 {
                local ++i
                mat `H'[`i',1] = _b[`y'_`j':`v']
                mat `H'[`i',2] = _se[`y'_`j':`v']
                mat `H'[`i',3] = `H'[`i',1] / `H'[`i',2]
                mat `H'[`i',4] = 2*normal(-abs(`H'[`i',3]))
                local lab = cond("`v'" == "h2", "E[u^2-1 | D,Z]", "E[u^3-3u | D,Z]")
                di as txt %22s "`lab'" "  {c |}" as res %11.5g `H'[`i',1] %12.5g `H'[`i',2] ///
                   %9.2f `H'[`i',3] %9.3f `H'[`i',4]
            }
        }
        di as txt "{hline 68}"
        qui test [`y'_1]h2 [`y'_1]h3 [`y'_0]h2 [`y'_0]h3
        local hchi2 = r(chi2)
        local hdf   = r(df)
        local hp    = r(p)
        qui test [`y'_1]h2 [`y'_1]h3
        local h1chi2 = r(chi2)
        local h1p    = r(p)
        qui test [`y'_0]h2 [`y'_0]h3
        local h0chi2 = r(chi2)
        local h0p    = r(p)
        di as txt "  joint      Wald chi2(" as res `hdf' as txt ") = " as res %8.2f `hchi2' as txt "   Prob > chi2 = " as res %6.4f `hp'
        di as txt "  regime 1   Wald chi2(2) = " as res %8.2f `h1chi2' as txt "   Prob > chi2 = " as res %6.4f `h1p'
        di as txt "  regime 0   Wald chi2(2) = " as res %8.2f `h0chi2' as txt "   Prob > chi2 = " as res %6.4f `h0p'
        di as txt "Reading: (a) rejects when the regime laws are not normal (the likelihood is then"
        di as txt "inconsistent, the two-step is not); (b) rejects when E[omega_j | u] bends (both routes"
        di as txt "lose the curve; the semiparametric MTE of esrmte follows it)."
        mat colnames `H' = est se z p
        mat rownames `H' = r1_h2 r1_h3 r0_h2 r0_h3
        mat colnames `Q' = twostep fiml se_diff
        return matrix hermite = `H'
        return matrix q       = `Q'
        return scalar chi2       = `R'[1,1]
        return scalar df         = `R'[1,2]
        return scalar p          = `R'[1,3]
        return scalar kappa_2s   = `R'[1,4]
        return scalar kappa_fiml = `R'[1,5]
        return scalar kappa_diff = `R'[1,6]
        return scalar kappa_se   = `R'[1,7]
        return scalar chi2_h     = `hchi2'
        return scalar df_h       = `hdf'
        return scalar p_h        = `hp'
        return scalar chi2_h1    = `h1chi2'
        return scalar chi2_h0    = `h0chi2'
        return local  test "normal"
    }

    * =========================================================================
    if ("`pdid'" != "") {
        if (`nq' < 2) {
            di as err "nq() must be at least 2"
            exit 198
        }
        local dv    "`e(treat)'"
        local xv    "`e(xvars)'"
        local zv    "`e(zvars)'"
        local wtype "`e(wtype)'"
        local wexp0 "`e(wexp)'"
        local kap2  = e(kappa)
        local skap2 = e(se_kappa)
        local wexp  ""
        local wexpa ""
        if ("`wtype'" != "") {
            local wexp  "[`wtype' `wexp0']"
            local wexpa "[aweight `wexp0']"
        }
        tempvar smp w v l1 l0 s top bot
        qui gen byte `smp' = e(sample)
        if ("`wtype'" != "") qui gen double `w' `wexp0' if `smp'
        else                 qui gen double `w' = 1 if `smp'
        _esreg_data if `smp'
        local xv_m "`r(x)'"
        local zv_m "`r(z)'"
        qui probit `dv' `zv_m' `wexp' if `smp'
        qui predict double `v' if `smp', xb
        qui gen double `l1' = exp(lnnormalden(`v') - lnnormal(`v'))  if `smp'
        qui gen double `l0' = exp(lnnormalden(`v') - lnnormal(-`v')) if `smp'
        qui xtile `s' = `v' [aw = `w'] if `smp', nq(`nq')
        qui summarize `s' if `smp', meanonly
        local G = r(max)
        qui gen byte `top' = (`s' == `G') if `smp'
        qui gen byte `bot' = (`s' == 1)   if `smp'
        tempname T
        mat `T' = J(`G', 8, .)
        forvalues q = 1/`G' {
            qui count if `smp' & `s' == `q' & `dv' == 1
            mat `T'[`q',1] = r(N)
            qui count if `smp' & `s' == `q' & `dv' == 0
            mat `T'[`q',2] = r(N)
            qui summarize `v' [aw = `w'] if `smp' & `s' == `q', meanonly
            mat `T'[`q',3] = r(mean)
            mat `T'[`q',4] = normal(r(mean))
            qui regress `y' `xv_m' `dv' `wexp' if `smp' & `s' == `q'
            mat `T'[`q',5] = _b[`dv']
            mat `T'[`q',6] = _se[`dv']
            qui summarize `l1' [aw = `w'] if `smp' & `s' == `q' & `dv' == 1, meanonly
            mat `T'[`q',7] = r(mean)
            qui summarize `l0' [aw = `w'] if `smp' & `s' == `q' & `dv' == 0, meanonly
            mat `T'[`q',8] = r(mean)
        }
        * within-group contrasts, top against bottom stratum, X-adjusted
        qui regress `y' `xv_m' `top' `wexp' if `smp' & `dv' == 1 & (`top' | `bot')
        local C1  = _b[`top']
        local sC1 = _se[`top']
        qui regress `y' `xv_m' `top' `wexp' if `smp' & `dv' == 0 & (`top' | `bot')
        local C0  = _b[`top']
        local sC0 = _se[`top']
        cap drop __esr_fv*
        local dlam1 = `T'[`G',7] - `T'[1,7]
        local dlam0 = `T'[`G',8] - `T'[1,8]
        local rs1 = `C1' / `dlam1'
        local rs0 = -`C0' / `dlam0'
        local kdd = `rs1' - `rs0'
        local skdd = sqrt((`sC1'/`dlam1')^2 + (`sC0'/`dlam0')^2)
        local dd  = `T'[`G',5] - `T'[1,5]
        local sdd = sqrt(`T'[`G',6]^2 + `T'[1,6]^2)
        * ---- display ----
        di
        di as txt "Pseudo-DiD on the selection index (H0: no selection on gains, kappa = 0)"
        di as txt "(strata of the probit index Z'g; X-adjusted gap treated - untreated within each)"
        di as txt "{hline 78}"
        di as txt "stratum {c |}    n1     n0   mean P {c |}     gap  Std. err. {c |} mean lambda1  mean lambda0"
        di as txt "{hline 8}{c +}{hline 22}{c +}{hline 20}{c +}{hline 26}"
        forvalues q = 1/`G' {
            di as txt %5.0f `q' "   {c |}" as res %6.0f `T'[`q',1] %7.0f `T'[`q',2] %9.3f `T'[`q',4] as txt " {c |}" ///
               as res %8.4f `T'[`q',5] %11.4f `T'[`q',6] as txt " {c |}" as res %13.4f `T'[`q',7] %14.4f `T'[`q',8]
        }
        di as txt "{hline 78}"
        di as txt "  double difference of the gaps, top - bottom:    " as res %8.4f `dd' as txt "  (s.e. " as res %6.4f `sdd' as txt ")"
        di as txt "  within treated,   top - bottom (C1):           " as res %8.4f `C1' as txt "  (s.e. " as res %6.4f `sC1' as txt ", z = " as res %6.2f `C1'/`sC1' as txt ")   sign of rho1*sigma1 = sign(C1/dlam1), dlam1 = " as res %6.3f `dlam1'
        di as txt "  within untreated, top - bottom (C0):           " as res %8.4f `C0' as txt "  (s.e. " as res %6.4f `sC0' as txt ", z = " as res %6.2f `C0'/`sC0' as txt ")   sign of rho0*sigma0 = -sign(C0/dlam0), dlam0 = " as res %6.3f `dlam0'
        di as txt "  rho1*sigma1 = C1/dlam1 = " as res %8.4f `rs1' as txt "    rho0*sigma0 = -C0/dlam0 = " as res %8.4f `rs0'
        di as txt "  kappa_dd = " as res %8.4f `kdd' as txt "  (s.e. " as res %6.4f `skdd' as txt ", z = " as res %6.2f `kdd'/`skdd' ///
           as txt ")     esreg kappa = " as res %8.4f `kap2' as txt " (s.e. " as res %6.4f `skap2' as txt ")"
        di as txt "Reading: with kappa > 0 the gap falls across the strata (the eager treated at low P gain"
        di as txt "more, the reluctant untreated at high P would have gained less); the signs of C1 and C0 need"
        di as txt "only A2; kappa_dd uses the Mills factors of A1. Never split on the outcome itself (collider)."
        mat colnames `T' = n1 n0 v_mean p_mean gap se_gap lambda1 lambda0
        return matrix strata = `T'
        return scalar C1     = `C1'
        return scalar se_C1  = `sC1'
        return scalar C0     = `C0'
        return scalar se_C0  = `sC0'
        return scalar dlam1  = `dlam1'
        return scalar dlam0  = `dlam0'
        return scalar rhosig1 = `rs1'
        return scalar rhosig0 = `rs0'
        return scalar kappa_dd = `kdd'
        return scalar se_kappa_dd = `skdd'
        return scalar dd     = `dd'
        return scalar se_dd  = `sdd'
        return scalar nq     = `G'
        return local  test "pdid"
    }
end
