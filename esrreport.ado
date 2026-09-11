*! esrreport 0.1.1  12sep2026  A. Araar
*! Reading of an esreg estimation for the practitioner (skeleton).  Runs the
*! diagnostics and tests of the family on the stored estimation, then prints an
*! ordered reading with rule-based notes: (1) the selection equation (strength,
*! form, support, index specification, sufficiency); (2) selection on gains
*! (kappa, its sign, the pseudo-DiD pattern, kappa(x)); (3) the joint law and the
*! route to retain (Hausman, gamma contrast, Hermite); (4) the effects to report
*! (full and on the common support, profile by quantile of the score).
*! The notes are statements and questions, not verdicts: the thresholds are
*! conventions, set by options, and every line names the test it rests on.
*!
*!   esrreport [, EST(name) NQ(#) LRmin(#) VIFmax(#) EXTRAPmax(#) ALPHA(#)
*!                EXPort(filename) REPLACE noTESTs]
*!
*!   nq(#)        quantile groups of the score for the profile of effects (default 5)
*!   lrmin(#)     LR of the excluded instruments below which they are called weak (10)
*!   vifmax(#)    VIF of lambda above which identification is called "by form" (10)
*!   extrapmax(#) share of a group outside the other's support above which the
*!                effects on the common support are recommended (0.10)
*!   alpha(#)     size of the tests in the reading (0.05)
*!   export()     write the figures and the notes to a plain-text file (the sheet
*!                an assistant, human or not, can read together with the study's
*!                context); replace overwrites
*!   notests      diagnostics and effects only (no re-estimation)
*! Returns r(route) ("fiml", "twostep", "semipar", or "index" / "instruments" when the
*! selection equation must be fixed first), r(n_notes), r(note1)... and
*! the figures used: r(kappa) r(se_kappa) r(kappa_dd) r(p_link) r(p_gamma) r(p_suff)
*! r(p_hausman) r(p_hermite) r(lr_excl) r(vif) r(share_extrap) r(att) r(atu) r(ate)
*! r(att_cs) r(atu_cs) r(eff_first) r(eff_last).
cap program drop esrreport
program define esrreport, rclass
    version 16
    syntax [, EST(name) NQ(integer 5) LRmin(real 10) VIFmax(real 10) EXTRAPmax(real 0.10) ///
              ALPHA(real 0.05) EXPort(string) REPLACE noTESTs]
    esreg_engine
    tempname hold
    _estimates hold `hold', restore copy
    _esreg_getest `est'
    local y      "`e(depvar)'"
    local dv     "`e(treat)'"
    local xv     "`e(xvars)'"
    local zv     "`e(zvars)'"
    local method "`e(method)'"
    local kv     "`e(kappavars)'"
    local hs     "`e(hetsigma)'"
    local hr     "`e(hetrho)'"
    local N      = e(N)
    local N1     = e(N_treated)
    local att    = e(att)
    local atu    = e(atu)
    local ate    = e(ate)
    local kappa  = e(kappa)
    local skappa = e(se_kappa)
    local rho1   = e(rho1)
    local rho0   = e(rho0)
    local excl : list zv - xv
    local nexcl : word count `excl'
    local nnotes = 0

    * ------------------------------------------------------------------ 1. selection equation
    tempvar P
    qui predict double `P' if e(sample), pr
    qui esrdiag
    local lr     = r(lr_excl)
    local plr    = r(p_excl)
    local vif    = max(r(vif1), r(vif0))
    local varP   = r(varP_share)
    local sh_att = r(share_att)
    local sh_atu = r(share_atu)
    local att_cs = r(att_cs)
    local atu_cs = r(atu_cs)
    local plink = .
    local pgam  = .
    local psuff = .
    local psuff1 = .
    local psuff0 = .
    local keep  ""
    if ("`tests'" == "") {
        qui esrtest, spec
        local plink = r(p)
        local pgam  = r(p_g)
        if (`nexcl' > 1) {
            qui esrtest, suff
            local psuff  = r(p)
            local psuff1 = r(p_1)
            local psuff0 = r(p_0)
            local keep   "`r(keep)'"
        }
    }
    di
    di as txt "{hline 78}"
    di as txt "Reading of the esreg estimation of " as res "`y'" as txt " on " as res "`dv'" ///
       as txt " (" as res "`method'" as txt ", N = " as res `N' as txt ", treated = " as res `N1' as txt ")"
    di as txt "{hline 78}"
    di as txt "1. The selection equation"
    if (`nexcl' == 0) {
        _esrr_note `nnotes' "no excluded instrument: the selection equation has no variable outside the outcome equations; kappa and the effects are identified by the curvature of lambda alone (identification by form, Section 3). Add an instrument before reading further."
    }
    else {
        di as txt "   excluded instrument(s): " as res "`excl'" as txt "   LR chi2(" as res `nexcl' as txt ") = " ///
           as res %6.2f `lr' as txt "   Var(P|X)/Var(P) = " as res %5.3f `varP' as txt "   max VIF(lambda) = " as res %5.2f `vif'
        if (`lr' < `lrmin') {
            _esrr_note `nnotes' "weak instruments (LR " + string(`lr', "%5.1f") + " < " + string(`lrmin') + "): kappa and the effects will be imprecise and the likelihood may return a confident wrong number; look for a stronger instrument (esrdiag)."
        }
        if (`vif' > `vifmax') {
            _esrr_note `nnotes' "the Mills ratio is nearly collinear with X (VIF " + string(`vif', "%5.1f") + "): identification rests on the functional form of lambda; read kappa as a form-dependent quantity (esrdiag)."
        }
        if (`nexcl' == 1) {
            _esrr_note `nnotes' "a single excluded instrument: index sufficiency cannot be tested (the instrument tests itself, Section 4.2) and part of the identification is by form; a second instrument would make the exclusion testable."
        }
        if (`plink' < `alpha') {
            _esrr_note `nnotes' "the link test rejects (p = " + string(`plink', "%5.3f") + "): the selection index is not linear in Z or the link is not the normal cdf; add polynomials or interactions to the selection equation (esrtest, spec add()) before anything else, since both routes use the Mills ratios of this index."
        }
        else if (`plink' < .) {
            di as txt "   link test of the index: p = " as res %5.3f `plink' as txt " (no sign of a misspecified index)"
        }
        if (`psuff' < `alpha') {
            local where = cond(`psuff1' < `alpha' & `psuff0' >= `alpha', "through the treated regime", ///
                          cond(`psuff0' < `alpha' & `psuff1' >= `alpha', "through the untreated regime", "in both regimes"))
            _esrr_note `nnotes' "index sufficiency rejected (p = " + string(`psuff', "%5.3f") + ", `where'): an instrument other than `keep' carries information on the regime errors beyond the index -- a direct effect on the outcome (the exclusion is wrong) or a second index in participation; kappa is biased (Section 3). Reconsider the excluded instruments."
        }
        else if (`psuff' < .) {
            di as txt "   index sufficiency (kept: " as res "`keep'" as txt "): p = " as res %5.3f `psuff' as txt " (no sign of information beyond the index)"
        }
    }
    di as txt "   support of P(Z): " as res "[" %5.3f e(supp_lo) ", " %5.3f e(supp_hi) "]" ///
       as txt "   extrapolated: treated " as res %5.3f `sh_att' as txt ", untreated " as res %5.3f `sh_atu'
    if (max(`sh_att', `sh_atu') > `extrapmax') {
        _esrr_note `nnotes' "more than " + string(100*`extrapmax', "%3.0f") + "% of one group has no counterpart in the other: report the effects on the common support (ATT " + string(`att_cs', "%8.4g") + ", ATU " + string(`atu_cs', "%8.4g") + ") next to the full ones."
    }

    * ------------------------------------------------------------------ 2. selection on gains
    di as txt "2. Selection on gains"
    local zk = `kappa' / `skappa'
    di as txt "   kappa = rho1 sigma1 - rho0 sigma0 = " as res %8.4g `kappa' as txt " (s.e. " as res %7.4g `skappa' ///
       as txt ", z = " as res %5.2f `zk' as txt ")   rho1 = " as res %6.3f `rho1' as txt "   rho0 = " as res %6.3f `rho0'
    local kdd = .
    if ("`tests'" == "") {
        qui esrtest, pdid nq(4)
        local kdd  = r(kappa_dd)
        local skdd = r(se_kappa_dd)
        local C1   = r(C1)
        local sC1  = r(se_C1)
        local C0   = r(C0)
        local sC0  = r(se_C0)
        di as txt "   pseudo-DiD on the index (4 strata): kappa_dd = " as res %8.4g `kdd' as txt " (s.e. " as res %7.4g `skdd' ///
           as txt ")   C1 = " as res %7.4g `C1' as txt " (z " as res %4.1f `C1'/`sC1' as txt ")   C0 = " as res %7.4g `C0' as txt " (z " as res %4.1f `C0'/`sC0' as txt ")"
        if (`C1'/`sC1' < -1.96 & `C0'/`sC0' > 1.96) {
            _esrr_note `nnotes' "the pseudo-DiD pattern C1 < 0, C0 > 0 is the signature of selection on gains (kappa > 0): the eager treated gain more, the reluctant untreated would gain less; the ATT exceeds the ATU and a program extended to the untreated would yield less than the ATT (Section 4.6)."
        }
        else if (`C1'/`sC1' > 1.96 & `C0'/`sC0' < -1.96) {
            _esrr_note `nnotes' "the pattern C1 > 0, C0 < 0 says kappa < 0: those who select in are those who gain least -- the targeting story of a program reaching the disadvantaged; the ATT understates the effect on the untreated."
        }
        else if (abs(`C1'/`sC1') > 1.96 & abs(`C0'/`sC0') > 1.96) {
            _esrr_note `nnotes' "C1 and C0 have the same sign: rho1 and rho0 have the same sign, a common selection on the level of the outcome rather than a selection on gains; kappa is the difference of two like quantities and is imprecise by construction."
        }
        else if (abs(`C1'/`sC1') > 1.96 | abs(`C0'/`sC0') > 1.96) {
            local one = cond(abs(`C1'/`sC1') > 1.96, "C1 (treated): rho1 sigma1 is " + cond(`C1' < 0, "positive", "negative"), "C0 (untreated): rho0 sigma0 is " + cond(`C0' > 0, "negative", "positive"))
            _esrr_note `nnotes' "only one within-group contrast is significant, `one'; the other regime shows no selection on unobservables in the strata, and selection on gains (their difference) is not established without the parametric model."
        }
        else {
            _esrr_note `nnotes' "no pseudo-DiD pattern: neither within-group contrast is significant; no selection on unobservables is visible across the strata of the index."
        }
    }
    if (abs(`zk') < 1.96) {
        _esrr_note `nnotes' "kappa is not significantly different from zero: the ATT, ATU and ATE do not differ beyond sampling error; the model of a common effect is not rejected on this ground, but a wide standard error is not evidence of kappa = 0 (see the strength of the instruments above)."
    }
    else {
        _esrr_note `nnotes' "kappa " + cond(`kappa' > 0, "> 0: selection on gains, ATT > ATU; the ATT is the effect of the program as it is, the ATU what its extension would yield.", "< 0: those who participate gain less than those who do not; the ATT understates the potential effect of the program on the untreated.")
    }
    if ("`kv'" != "" | "`hs'`hr'" != "") {
        qui esrtest, kappa
        local pkx = r(p)
        if (`pkx' < `alpha') _esrr_note `nnotes' "kappa varies with the covariates of kappa()/hetsigma()/hetrho() (p = " + string(`pkx', "%5.3f") + "): report kappa(x) and the effects by these covariates rather than a single kappa (Section 4.3)."
        else di as txt "   kappa(x): p = " as res %5.3f `pkx' as txt " (no evidence of heterogeneity in the listed covariates)"
    }

    * ------------------------------------------------------------------ 3. joint law and route
    local route "`method'"
    local phaus = .
    local pherm = .
    if ("`tests'" == "") {
        di as txt "3. The joint law and the route"
        qui esrtest, normal
        local phaus = r(p)
        local pherm = r(p_h)
        local k2s   = r(kappa_2s)
        local kF    = r(kappa_fiml)
        local zkd   = r(kappa_diff) / r(kappa_se)
        di as txt "   Hausman FIML vs two-step: p = " as res %5.3f `phaus' as txt "   gamma contrast: p = " as res %5.3f `pgam' ///
           as txt "   Hermite controls: p = " as res %5.3f `pherm'
        di as txt "   kappa two-step " as res %8.4g `k2s' as txt "   FIML " as res %8.4g `kF' as txt "   z of the difference " as res %5.2f `zkd'
        if (`plink' < `alpha') {
            local route "index"
            _esrr_note `nnotes' "the index is misspecified (link test above): the Mills ratios, the Hermite controls and the contrasts are all functions of the wrong index, so the tests of this block are not readable yet. Respecify the selection equation first, then run esrreport again."
        }
        else if (`psuff' < `alpha') {
            local route "instruments"
            _esrr_note `nnotes' "index sufficiency fails (note above): both routes are biased whatever the joint law (Section 3), and the contrasts of this block compare two biased estimators. Reconsider the excluded instruments first, then run esrreport again."
        }
        else if (`pherm' < `alpha') {
            local route "semipar"
            _esrr_note `nnotes' "the Hermite controls reject (p = " + string(`pherm', "%5.3f") + "): E[omega | u] is not linear in u; both parametric routes fit the wrong curve and the ATT/ATU depend on where the curve is read. Use the semiparametric MTE curve (esrmte, semipar) and the effects on the common support; treat kappa as a local slope."
        }
        else if (`phaus' < `alpha') {
            local route "twostep"
            _esrr_note `nnotes' "the Hausman contrast rejects (p = " + string(`phaus', "%5.3f") + ") with a clean link test and clean Hermite controls: the regime errors are not normal (A3); the likelihood is inconsistent and the two-step is not. Report the two-step estimates (method(twostep)); a transformation of the outcome may remove the rejection."
        }
        else {
            local route "fiml"
            _esrr_note `nnotes' "no rejection of trivariate normality: the two routes agree and the likelihood is the efficient one; report FIML, with the two-step as a check."
        }
        if (`pgam' < `alpha' & `phaus' >= `alpha') {
            _esrr_note `nnotes' "the gamma contrast rejects (p = " + string(`pgam', "%5.3f") + ") while the full Hausman does not: the disagreement is confined to the selection equation; read it with the link test (index) before A3."
        }
    }

    * ------------------------------------------------------------------ 4. effects and profile
    di as txt "4. The effects"
    di as txt "   ATT = " as res %8.4g `att' as txt "   ATU = " as res %8.4g `atu' as txt "   ATE = " as res %8.4g `ate' ///
       as txt "   (common support: ATT " as res %8.4g `att_cs' as txt ", ATU " as res %8.4g `atu_cs' as txt ")"
    qui esrcurve, rank(`P') nq(`nq') nose
    tempname T
    mat `T' = r(table)
    local G = r(nq)
    local e1 = `T'[1, 6]
    local eG = `T'[`G', 6]
    di as txt "   effect by quantile of P(Z), " as res `G' as txt " groups: first " as res %8.4g `e1' as txt "   last " as res %8.4g `eG'
    if ("`route'" != "`method'" & inlist("`route'", "fiml", "twostep")) {
        di as txt "   (the effects above are those of the current estimation, " as res "`method'" as txt "; the route to report is " as res "`route'" as txt ": re-estimate with method(`route'))"
    }
    else if (inlist("`route'", "index", "instruments")) {
        di as txt "   (the effects above are not to be reported before the selection equation is fixed)"
    }
    else if ("`route'" == "semipar") {
        di as txt "   (report the effects on the common support and the semiparametric curve of esrmte, semipar)"
    }
    if (`kappa' != 0 & abs(`zk') >= 1.96) {
        _esrr_note `nnotes' "the profile of the effect along the score (first to last group: " + string(`e1', "%8.4g") + " to " + string(`eG', "%8.4g") + ") is the marginal treatment effect read through the parametric line; where the semiparametric curve of esrmte leaves the line, trust the curve (Section 5.3)."
    }
    di as txt "{hline 78}"
    di as txt "Route to report: " as res "`route'" as txt "   (" as res `nnotes' as txt " notes above)"
    di as txt "{hline 78}"

    * ------------------------------------------------------------------ export
    if ("`export'" != "") {
        tempname fh
        cap confirm file "`export'"
        if (_rc == 0 & "`replace'" == "") {
            di as err "file `export' already exists; use replace"
            exit 602
        }
        cap erase "`export'"
        file open `fh' using "`export'", write text
        file write `fh' "esreg reading sheet" _n
        file write `fh' "outcome `y'; treatment `dv'; method `method'; N `N'; treated `N1'" _n
        file write `fh' "x: `xv'" _n "z: `zv'" _n "excluded: `excl'" _n
        file write `fh' "lr_excl `lr'; p_excl `plr'; varP_share `varP'; vif `vif'" _n
        file write `fh' "supp_lo " (e(supp_lo)) "; supp_hi " (e(supp_hi)) "; share_att_extrap `sh_att'; share_atu_extrap `sh_atu'" _n
        file write `fh' "p_link `plink'; p_gamma `pgam'; p_suff `psuff'; keep `keep'" _n
        file write `fh' "kappa `kappa'; se_kappa `skappa'; rho1 `rho1'; rho0 `rho0'; kappa_dd `kdd'" _n
        file write `fh' "p_hausman `phaus'; p_hermite `pherm'" _n
        file write `fh' "att `att'; atu `atu'; ate `ate'; att_cs `att_cs'; atu_cs `atu_cs'; eff_first `e1'; eff_last `eG'" _n
        file write `fh' "route `route'" _n
        forvalues i = 1/`nnotes' {
            file write `fh' "note`i': `note`i''" _n
        }
        file close `fh'
        di as txt "(reading sheet written to " as res "`export'" as txt ")"
    }

    * ------------------------------------------------------------------ returns
    return local  route      "`route'"
    return scalar n_notes    = `nnotes'
    forvalues i = 1/`nnotes' {
        return local note`i' "`note`i''"
    }
    return scalar kappa      = `kappa'
    return scalar se_kappa   = `skappa'
    return scalar kappa_dd   = `kdd'
    return scalar p_link     = `plink'
    return scalar p_gamma    = `pgam'
    return scalar p_suff     = `psuff'
    return scalar p_hausman  = `phaus'
    return scalar p_hermite  = `pherm'
    return scalar lr_excl    = `lr'
    return scalar vif        = `vif'
    return scalar share_extrap = max(`sh_att', `sh_atu')
    return scalar att        = `att'
    return scalar atu        = `atu'
    return scalar ate        = `ate'
    return scalar att_cs     = `att_cs'
    return scalar atu_cs     = `atu_cs'
    return scalar eff_first  = `e1'
    return scalar eff_last   = `eG'
end

* ---------------------------------------------------------------------------
* print one numbered note, wrapped, and keep it in the caller's locals
cap program drop _esrr_note
program define _esrr_note
    gettoken n 0 : 0
    local txt = `0'
    local i = `n' + 1
    c_local nnotes `i'
    c_local note`i' `"`txt'"'
    di as txt "   " as res "note `i'" as txt ": " `"`txt'"'
end
