*! esrdiag 1.0.0  03oct2026  A. Araar
*! Corrected 02oct2026: after svy, subpop() the diagnostics work on the subpopulation.
*! Diagnostics of the selection equation after esreg (Section 4.4 of the paper):
*!   fit        -- log likelihood and pseudo-R2 of the selection equation on its
*!                 own (probit of D on the whole of Z), always reported;
*!   strength   -- LR test and incremental pseudo-R2 of the excluded instruments
*!                 in the probit, chi2(1) per instrument;
*!   variation  -- share of the variance of P(Z) not explained by X, and the VIF
*!                 of the Mills ratio lambda_j on X within each regime (how far
*!                 identification rests on the curvature of lambda alone);
*!   support    -- range of P(Z) by group, common support, weighted share of the
*!                 treated (ATT) and untreated (ATU) whose counterfactual is
*!                 extrapolated outside the other group's range, and the ATT/ATU
*!                 restricted to the common support next to the full ones.
*!
*!   esrdiag [, EST(name)]
*! Returns r(lr_excl) r(df_excl) r(p_excl) r(r2_full) r(r2_restr) r(r2_incr)
*!         r(ll_sel) r(ll0_sel) r(df_sel)
*!         r(varP_share) r(vif1) r(vif0) r(supp_lo) r(supp_hi) r(p_min1) r(p_max1)
*!         r(p_min0) r(p_max0) r(share_att) r(share_atu) r(att) r(atu) r(att_cs)
*!         r(atu_cs) r(instr) (chi2 per instrument) r(weak) r(form)
cap program drop esrdiag
program define esrdiag, rclass
    version 16
    syntax [, EST(name)]
    esreg_engine
    tempname hold
    _estimates hold `hold', restore copy
    _esreg_getest `est'

    * ---- what we need from e(), before any estimation command runs ----------------
    local y      "`e(depvar)'"
    local dv     "`e(treat)'"
    local xv     "`e(xvars)'"
    local zv     "`e(zvars)'"
    local method "`e(method)'"
    local wtype  "`e(wtype)'"
    local wexp0  "`e(wexp)'"
    local att    = e(att)
    local atu    = e(atu)
    local p1min  = e(p_min1)
    local p1max  = e(p_max1)
    local p0min  = e(p_min0)
    local p0max  = e(p_max0)
    local plo    = e(supp_lo)
    local phi    = e(supp_hi)
    tempvar smp w dd P l1 l0 zg
    * the fitted observations (the subpopulation after svy, subpop())
    _esreg_esample `smp'
    if ("`wtype'" != "") qui gen double `w' `wexp0' if `smp'
    else                 qui gen double `w' = 1 if `smp'
    local wexp  ""
    local wexpa ""
    if ("`wtype'" != "") {
        local wexp  "[`wtype' `wexp0']"
        local wexpa "[aweight `wexp0']"
    }
    qui predict double `dd' if `smp', effect
    qui predict double `zg' if `smp', xbsel
    qui gen double `P'  = normal(`zg') if `smp'
    qui gen double `l1' = exp(lnnormalden(`zg') - lnnormal(`zg'))  if `smp'
    qui gen double `l0' = exp(lnnormalden(`zg') - lnnormal(-`zg')) if `smp'
    local excl : list zv - xv
    local nexcl : word count `excl'
    * plain variables behind factor terms (2.region 3.region would be re-merged by
    * probit/regress with 2 as the base); __esr_fv* are dropped at the end
    _esreg_data if `smp'
    local xv_m "`r(x)'"
    local zv_m "`r(z)'"
    local excl_m : list zv_m - xv_m

    * ---- strength of the excluded instruments -----------------------------------------
    tempname I
    * the selection equation on its own: probit of D on the whole of Z.  Computed
    * unconditionally, so that its fit is reported even without an exclusion.
    qui probit `dv' `zv_m' `wexp' if `smp'
    local llf  = e(ll)
    local ll0f = e(ll_0)
    local r2f  = e(r2_p)
    local kzf  = e(df_m)
    if (`nexcl' > 0) {
        mat `I' = J(`nexcl', 3, .)
        local rn ""
        local i = 0
        foreach v of local excl_m {
            local ++i
            cap mat `I'[`i',1] = _b[`v']
            cap mat `I'[`i',2] = _se[`v']
            cap mat `I'[`i',3] = (_b[`v'] / _se[`v'])^2
            local vd : word `i' of `excl'
            local vn = subinstr("`vd'", ".", "_", .)
            local rn "`rn' `vn'"
        }
        mat colnames `I' = coef se chi2
        mat rownames `I' = `rn'
        qui probit `dv' `xv_m' `wexp' if `smp'
        local llr = e(ll)
        local r2r = e(r2_p)
        local lr  = 2 * (`llf' - `llr')
        local plr = chi2tail(`nexcl', `lr')
    }
    else {
        local lr = .
        local plr = .
        local r2r = .
    }

    * ---- variation: Var(P | X) / Var(P) and the VIF of the Mills ratios -----------------
    qui regress `P' `xv_m' `wexpa' if `smp'
    local varP = 1 - e(r2)
    qui regress `l1' `xv_m' `wexpa' if `smp' & `dv' == 1
    local vif1 = 1 / (1 - e(r2))
    qui regress `l0' `xv_m' `wexpa' if `smp' & `dv' == 0
    local vif0 = 1 / (1 - e(r2))
    cap drop __esr_fv*

    * ---- support and extrapolated shares -------------------------------------------------
    tempvar out
    qui gen byte `out' = cond(`dv' == 1, `P' < `p0min' | `P' > `p0max', `P' < `p1min' | `P' > `p1max') if `smp'
    qui summarize `out' [aw = `w'] if `smp' & `dv' == 1, meanonly
    local share_att = r(mean)
    qui summarize `out' [aw = `w'] if `smp' & `dv' == 0, meanonly
    local share_atu = r(mean)
    qui summarize `dd' [aw = `w'] if `smp' & `dv' == 1 & `out' == 0, meanonly
    local att_cs = r(mean)
    qui summarize `dd' [aw = `w'] if `smp' & `dv' == 0 & `out' == 0, meanonly
    local atu_cs = r(mean)

    * ---- display -------------------------------------------------------------------------
    di
    di as txt "Diagnostics of the selection equation (" as res "`method'" as txt " estimation of " ///
       as res "`y'" as txt " on treatment " as res "`dv'" as txt ")"
    di as txt "{hline 76}"
    di as txt "Fit of the selection equation (probit of " as res "`dv'" as txt " on the whole of Z)"
    if (`ll0f' < .) {
        local lr0 = 2 * (`llf' - `ll0f')
        di as txt "  log likelihood " as res %11.4f `llf' as txt "   intercept only " as res %11.4f `ll0f'
        di as txt "  pseudo-R2 " as res %6.4f `r2f' as txt "   LR chi2(" as res `kzf' as txt ") = " ///
           as res %8.2f `lr0' as txt "   Prob > chi2 = " as res %6.4f chi2tail(`kzf', `lr0')
    }
    else di as txt "  log likelihood " as res %11.4f `llf' as txt "   pseudo-R2 " as res %6.4f `r2f'
    di as txt "{hline 76}"
    di as txt "Strength of the excluded instruments" _col(40) as res "`excl'"
    if (`nexcl' > 0) {
        di as txt "  LR test, with vs without them:  chi2(" as res `nexcl' as txt ") = " ///
           as res %8.2f `lr' as txt "   Prob > chi2 = " as res %6.4f `plr'
        di as txt "  pseudo-R2 of the probit:  with " as res %6.4f `r2f' as txt "   without " as res %6.4f `r2r' ///
           as txt "   increment " as res %6.4f `r2f' - `r2r'
        di as txt "  per instrument (Wald chi2(1) in the full probit):"
        forvalues i = 1/`nexcl' {
            local v : word `i' of `excl'
            di as txt %26s "`v'" as res %10.4f `I'[`i',1] " (" %7.4f `I'[`i',2] ")   chi2 = " %8.2f `I'[`i',3]
        }
    }
    else di as txt "  none: the selection equation has no variable outside the outcome equations"
    di as txt "{hline 76}"
    di as txt "Variation of the score beyond X"
    di as txt "  Var(P | X) / Var(P) = " as res %6.4f `varP' as txt "   (share of the variance of P(Z) not explained by X)"
    di as txt "  VIF of lambda_1 on X among the treated   = " as res %8.2f `vif1'
    di as txt "  VIF of lambda_0 on X among the untreated = " as res %8.2f `vif0'
    di as txt "{hline 76}"
    di as txt "Support of P(Z)"
    di as txt "  treated   [" as res %6.4f `p1min' as txt ", " as res %6.4f `p1max' as txt "]" ///
       _col(40) as txt "untreated [" as res %6.4f `p0min' as txt ", " as res %6.4f `p0max' as txt "]"
    di as txt "  common support [" as res %6.4f `plo' as txt ", " as res %6.4f `phi' as txt "]"
    di as txt "  share of the treated outside the untreated range (ATT extrapolated):   " as res %6.4f `share_att'
    di as txt "  share of the untreated outside the treated range (ATU extrapolated):   " as res %6.4f `share_atu'
    di as txt "  ATT: full " as res %9.5g `att' as txt "   on the common support " as res %9.5g `att_cs'
    di as txt "  ATU: full " as res %9.5g `atu' as txt "   on the common support " as res %9.5g `atu_cs'
    di as txt "{hline 76}"
    * ---- warnings ------------------------------------------------------------------------
    local weak = 0
    local form = 0
    if (`nexcl' > 0 & `lr' < 10) {
        local weak = 1
        di as res "warning: weak instruments (LR chi2 < 10); the effects and kappa will be imprecise"
    }
    if (max(`vif1', `vif0') > 10) {
        local form = 1
        di as res "warning: the Mills ratio is nearly collinear with X (VIF > 10): identification rests"
        di as res "         on the curvature of lambda alone (identification by form)"
    }
    if (max(`share_att', `share_atu') > 0.10) {
        di as res "note: more than 10% of one group has no counterpart in the other: the corresponding"
        di as res "      effect extrapolates the parametric line outside the common support"
    }
    if (`weak' + `form' == 0 & max(`share_att', `share_atu') <= 0.10) di as txt "no warning"

    * ---- returns -------------------------------------------------------------------------
    if (`nexcl' > 0) return matrix instr = `I'
    return scalar lr_excl    = `lr'
    return scalar df_excl    = `nexcl'
    return scalar p_excl     = `plr'
    return scalar r2_full    = `r2f'
    return scalar r2_restr   = `r2r'
    return scalar r2_incr    = `r2f' - `r2r'
    return scalar ll_sel     = `llf'
    return scalar ll0_sel    = `ll0f'
    return scalar df_sel     = `kzf'
    return scalar varP_share = `varP'
    return scalar vif1       = `vif1'
    return scalar vif0       = `vif0'
    return scalar supp_lo    = `plo'
    return scalar supp_hi    = `phi'
    return scalar p_min1     = `p1min'
    return scalar p_max1     = `p1max'
    return scalar p_min0     = `p0min'
    return scalar p_max0     = `p0max'
    return scalar share_att  = `share_att'
    return scalar share_atu  = `share_atu'
    return scalar att        = `att'
    return scalar atu        = `atu'
    return scalar att_cs     = `att_cs'
    return scalar atu_cs     = `atu_cs'
    return scalar weak       = `weak'
    return scalar form       = `form'
    return local  excluded   "`excl'"
    return local  method     "`method'"
end
