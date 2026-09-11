*! _esreg_pwr.ado  v1.0  (esreg internal engine)  Araar A.  2026sep08
*! Percentile Weights Regression -- private copy for esrmte, semipar.
*! Derived from gepwreg v1.4 (Araar 2016, 2023; Araar 2026, Exploring
*! Heterogeneous Effects: Quantile Models and Percentile Weights Regression,
*! Zenodo, doi:10.5281/zenodo.20315684) with the target() bandwidth
*! (derivative-aware, n^(-1/7) pilot + empirical MSE search). All programs,
*! Mata functions and scalars are prefixed _esreg_pwr / _epwr_ so that this
*! file never collides with the user's own gepwreg installation.
*! Not meant to be called directly; use  esrmte, semipar.
/* ─────────────────────────────────────────────────────────────────────────────
   Syntax:
     gepwreg depvar [indepvars] [fw/aw/pw] [if] [in],
           [ PERcentile(#) CBAND(#) BAND(#) OPTbw BOOT(#) Level(#) noConstant
             TARget(varlist) TCONst(#) TGRid(numlist) ]
   target(v1 [v2]) : choose h by minimising the empirical MSE of the coefficient
     on v1 (or, if v2 = v1^2 is also given, of the derivative b1 + 2 q_tau b2 at
     the tau-quantile q_tau of the ranking variable). Pilot h0 = tconst *
     sigma_p * n^(-1/7) (tconst default 1.6), curvature from a pilot curve at
     band 2*h0, search over h0 x tgrid (default .5 .7 1 1.4 2 2.8 4).
     Use for slope targets such as the MTE from the propensity score.
   ───────────────────────────────────────────────────────────────────────── */
#delimit ;

capture program drop _esreg_pwr ;
program define _esreg_pwr, eclass properties(svyb) sortpreserve ;
version 16.0 ;

syntax varlist(min=1 numeric fv) [fw aw pw] [if] [in] [,
    PERcentile(real 0.5)
    CBAND(real 0.9)
    BAND(real 0)
    OPTbw
    SILverman
    BOOT(integer 0)
    Level(cilevel)
    noConstant
    VCE(string)
    RANKvar(varname)
    TARget(varlist numeric min=1 max=2)
    TCONst(real 1.6)
    TGRid(numlist >0)
    /* Internal options passed by svy: prefix */
    SVY
    noHEader
    noTABle
] ;

/* ── 0. Mark sample ──────────────────────────────────────────────────────── */
marksample touse ;
qui count if `touse' ;
if r(N) == 0 error 2000 ;

/* ── 1. Parse varlist ────────────────────────────────────────────────────── */
local depvar    : word 1 of `varlist' ;
local indepvars : list varlist - depvar ;

/* Step 1 : fvexpand then filter base (Nb.) and omitted (No.) categories
   0b.rururb → removed    1.rururb → kept
   This gives the display names AND the clean list for fvrevar             */
if "`indepvars'" != "" {;
    fvexpand `indepvars' if `touse' ;
    local indepvars_exp_raw `r(varlist)' ;
    local indepvars_exp ;
    foreach v of local indepvars_exp_raw {;
        if !regexm("`v'","^[0-9]+b\.") & !regexm("`v'","^[0-9]+o\.") {;
            local indepvars_exp `indepvars_exp' `v' ;
        } ;
    } ;
} ;
else local indepvars_exp "" ;

/* Step 2 : fvrevar on the FILTERED list → real temp vars for Mata
   fvrevar on  hhsize 1.rururb  →  hhsize __000002  (no base var)
   Dimensions of indepvars_mata always equal dimensions of indepvars_exp  */
if "`indepvars_exp'" != "" {;
    fvrevar `indepvars_exp' if `touse' ;
    local indepvars_mata `r(varlist)' ;
} ;
else local indepvars_mata "" ;

/* ── 1b. Ranking variable ───────────────────────────────────────────────── */
/* Default: rank on depvar. If rankvar() specified, rank on that variable.  */
local rank_var "`depvar'" ;
local rank_label "y (outcome)" ;
if "`rankvar'" != "" {;
    local rank_var "`rankvar'" ;
    local rank_label "`rankvar'" ;
    markout `touse' `rankvar' ;
    di as text "Ranking variable : `rankvar' (instead of `depvar')" ;
} ;

/* ── 2. Survey weights ───────────────────────────────────────────────────── */
/* Priority rule (highest to lowest) :                                       */
/*   1. Weights declared in the command  [pw=weight]  → used as-is          */
/*   2. Weights declared in svyset       [pw=weight]  → used automatically   */
/*   3. No weights anywhere              → unweighted (fw=1)                 */
/* If both command and svyset declare weights, command weights take priority  */
/* and a note is displayed so the user is aware.                             */
tempvar fw_var ;
qui gen double `fw_var' = 1 if `touse' ;
local wgt_name "(unweighted)" ;
if "`exp'" != "" {;
    /* ── Priority 1: user declared weights in the command ── */
    local wexpr = subinstr("`exp'","=","",1) ;
    /* Trim spaces from wexpr for clean comparison */
    local wexpr = trim("`wexpr'") ;
    qui replace `fw_var' = `wexpr' if `touse' ;
    local wgt_type = upper(substr("`weight'",1,2)) ;
    local wgt_name "`wgt_type'[`wexpr']" ;
    /* Check svyset — compare with command weights */
    qui svyset ;
    if "`r(wvar)'" != "" {;
        local svywvar = trim("`r(wvar)'") ;
        if "`svywvar'" == "`wexpr'" {;
            /* Same variable declared in both places — silently use it once */
            /* No message needed : consistent declaration */
        } ;
        else {;
            /* Different variables — command takes priority, warn user */
            di as text "Note: command weights [pw=`wexpr'] take priority over" ;
            di as text "      svyset weights [`svywvar']. Using [pw=`wexpr']." ;
        } ;
    } ;
} ;
else {;
    /* ── Priority 2: no command weights — check svyset ── */
    qui svyset ;
    if "`r(wvar)'" != "" {;
        qui replace `fw_var' = `r(wvar)' if `touse' ;
        local wgt_name "PW[`r(wvar)'] (svyset)" ;
    } ;
    /* else: Priority 3 — unweighted, wgt_name stays "(unweighted)" */
} ;

/* ── 3. Validate options ─────────────────────────────────────────────────── */
if `percentile' <= 0 | `percentile' >= 1 {;
    di as error "percentile() must be strictly between 0 and 1" ;
    error 198 ;
} ;
if `boot' < 0 {;
    di as error "boot() must be a non-negative integer" ;
    error 198 ;
} ;
/* ── Bandwidth method ───────────────────────────────────────────────────── */
/* MSE-optimal is the DEFAULT (paper Algorithm 1).  Override with:           */
/*   silverman   → Silverman rule  h = cband * min(sd,IQR/1.34) * n^(-1/5)    */
/*   band(#)     → manual fixed bandwidth                                     */
/* optbw is accepted as an explicit synonym for the (optimal) default.       */
local n_bwspec = (`band' != 0) + ("`silverman'" != "") + ("`optbw'" != "") + ("`target'" != "") ;
if `n_bwspec' > 1 {;
    di as error "specify at most one of band(), silverman, optbw, target()" ;
    error 198 ;
} ;
/* target() : derivative-aware bandwidth for one coefficient (v1.4) */
local jcoef = 0 ;
local jsq   = 0 ;
if "`tgrid'" == "" local tgrid ".5 .7 1 1.4 2 2.8 4" ;
if "`target'" != "" {;
    local t1 : word 1 of `target' ;
    local t2 : word 2 of `target' ;
    local jcoef : list posof "`t1'" in indepvars_exp ;
    if `jcoef' == 0 {;
        di as error "target(): `t1' is not among the independent variables" ;
        error 198 ;
    } ;
    if "`t2'" != "" {;
        local jsq : list posof "`t2'" in indepvars_exp ;
        if `jsq' == 0 {;
            di as error "target(): `t2' is not among the independent variables" ;
            error 198 ;
        } ;
    } ;
    if `tconst' <= 0 {;
        di as error "tconst() must be positive" ;
        error 198 ;
    } ;
} ;
if "`target'" != "" {;
    local do_optbw = 2 ;
    local bw_label "Target-MSE (n^-1/7 pilot)" ;
} ;
else if `band' != 0 {;
    local do_optbw = 0 ;
    local bw_label "Manual" ;
} ;
else if "`silverman'" != "" {;
    local do_optbw = 0 ;
    local bw_label "Silverman" ;
} ;
else {;
    local do_optbw = 1 ;          /* DEFAULT: MSE-optimal */
    local bw_label "MSE-optimal" ;
} ;
local addcons = ("`constant'" == "") ;

/* ── 3b. Parse vce() option for survey design ───────────────────────────── */
/* Detect svy: prefix (Stata passes vce(linearized) or sets svy option)    */
local do_svy = 0 ;
local svy_psu    "" ;
local svy_strata "" ;
if "`svy'" != "" {;
    local vce "linearized" ;
} ;
if "`vce'" != "" {;
    if lower("`vce'")=="svy" | lower("`vce'")=="survey" | lower("`vce'")=="linearized" {;
        /* Check svyset has been called */
        qui svyset ;
        if "`r(wvar)'" == "" & "`r(su1)'" == "" {;
            di as error "vce(svy) requires svyset to be called first" ;
            error 119 ;
        } ;
        local do_svy    = 1 ;
        local svy_psu    "`r(su1)'" ;     /* svyset stores PSU in r(su1)    */
        local svy_strata "`r(strata1)'" ;
        /* If no weights supplied by user, take them from svyset            */
        if "`exp'" == "" & "`r(wvar)'" != "" {;
            qui replace `fw_var' = `r(wvar)' if `touse' ;
            local wgt_name "PW[`r(wvar)'] (svyset)" ;
        } ;
        di as text "Survey design detected:" ;
        if "`svy_psu'" == "" {;
            di as text "  PSU    : (not declared — variance may be conservative)" ;
        } ;
        else {;
            di as text "  PSU    : `svy_psu'" ;
        } ;
        di as text "  Strata : `svy_strata'" ;
        /* Check minimum PSU per stratum */
        if "`svy_strata'" != "" & "`svy_psu'" != "" {;
            qui tab `svy_strata' if `touse' ;
            qui {;
                tempvar _npsu ;
                bysort `svy_strata' (`svy_psu') : gen `_npsu' = (_n==1) ;
                bysort `svy_strata' : replace `_npsu' = sum(`_npsu') ;
                bysort `svy_strata' : replace `_npsu' = `_npsu'[_N] ;
            } ;
            qui sum `_npsu' if `touse' ;
            local min_psu = r(min) ;
            if `min_psu' < 5 {;
                di as text "Warning: minimum PSU/stratum = `min_psu'." ;
                di as text "  Taylor SE may be unreliable. Recommend n_h >= 10." ;
            } ;
        } ;
    } ;
    else {;
        di as error "vce(`vce') not supported. Use vce(svy) for survey design." ;
        error 198 ;
    } ;
} ;

/* ── 4. Call Mata ────────────────────────────────────────────────────────── */
scalar _epwr_h_tmp = 0 ;
scalar _epwr_h0 = . ;
scalar _epwr_m2 = . ;
scalar _epwr_qt = . ;
mata: _epwr_main("`depvar'", "`indepvars_mata'", "`fw_var'", "`touse'",
                  `percentile', `cband', `band', `boot', `addcons',
                  `do_svy', "`svy_psu'", "`svy_strata'",
                  "`rank_var'", `do_optbw', `jcoef', `jsq', `tconst', "`tgrid'") ;

/* ── 5. Retrieve scalars ─────────────────────────────────────────────────── */
local h     = _epwr_h ;
local n_eff = _epwr_neff ;
local n_obs = _epwr_n ;
local tau   = `percentile' ;

/* ── 6. Name matrices ────────────────────────────────────────────────────── */
if `addcons' local vnames `indepvars_exp' _cons ;
else         local vnames `indepvars_exp' ;

matrix colnames _epwr_b  = `vnames' ;
matrix rownames _epwr_V  = `vnames' ;
matrix colnames _epwr_V  = `vnames' ;
matrix rownames _epwr_Vn = `vnames' ;
matrix colnames _epwr_Vn = `vnames' ;

/* ── 7. Post results ─────────────────────────────────────────────────────── */
ereturn post _epwr_b _epwr_V, esample(`touse') obs(`n_obs') ;
ereturn matrix V_naive = _epwr_Vn ;
if `do_svy' {;
    matrix rownames _epwr_Vs = `vnames' ;
    matrix colnames _epwr_Vs = `vnames' ;
    ereturn matrix V_svy = _epwr_Vs ;
} ;
if `boot' > 0 {;
    matrix rownames _epwr_Vb = `vnames' ;
    matrix colnames _epwr_Vb = `vnames' ;
    ereturn matrix V_boot = _epwr_Vb ;
} ;
ereturn scalar tau   = `tau' ;
ereturn scalar h     = `h' ;
ereturn local  bw_method "`bw_label'" ;
ereturn scalar N_eff = `n_eff' ;
if "`target'" != "" {;
    ereturn local  target   "`target'" ;
    ereturn scalar h0       = _epwr_h0 ;
    ereturn scalar m2_target = _epwr_m2 ;
    ereturn scalar q_tau    = _epwr_qt ;
} ;
ereturn scalar boot   = `boot' ;
ereturn scalar do_svy = `do_svy' ;
ereturn local  cmd      "_esreg_pwr" ;
ereturn local  rankvar  "`rank_var'" ;
ereturn local  ranklabel "`rank_label'" ;
ereturn local  wgt      "`wgt_name'" ;
ereturn local  cmdline  "_esreg_pwr `0'" ;
ereturn local  title    "Percentile Weights Regression" ;
if `do_svy' {;
    ereturn local  SE_type  "Taylor linearisation (survey design)" ;
    ereturn local  vce      "linearized" ;
} ;
else {;
    ereturn local  SE_type  "Analytical (Influence Function)" ;
    ereturn local  vce      "IF-corrected" ;
} ;

/* ── 8. Display ──────────────────────────────────────────────────────────── */
di "" ;
di as text "Percentile Weights Regression" ;
di as text "{hline 60}" ;
di as text %28s "Target percentile (tau)" " = " as result %6.4f `tau' ;
di as text %28s "Bandwidth (h)"           " = " as result %9.6f `h' as text " (`bw_label')" ;
di as text %28s "Weights"                 " = " as result "`wgt_name'" ;
di as text %28s "Ranking variable"        " = " as result "`rank_label'" ;
di as text %28s "Observations"            " = " as result %7.0f `n_obs' ;
di as text %28s "Kernel N_eff"            " = " as result %7.1f `n_eff' ;
if "`target'" != "" {;
    di as text %28s "Target coefficient"      " = " as result "`target'" ;
    di as text %28s "Pilot h0 (n^-1/7)"       " = " as result %9.6f _epwr_h0 ;
    di as text %28s "Curvature m''(tau)"      " = " as result %9.4f _epwr_m2 ;
} ;
if `boot' > 0 {;
    di as text %28s "Bootstrap replications" " = " as result %7.0f `boot' ;
} ;
di as text "{hline 60}" ;
if `do_svy' {;
    di as text "SE: Taylor linearisation under survey design (PSU + strata)" ;
} ;
else {;
    di as text "SE: IF-corrected analytical (Deville 1999)" ;
} ;
di as text "{hline 60}" ;
di "" ;
ereturn display, level(`level') ;
di "" ;
di as text "Naive WLS SE stored in e(V_naive)." ;
if `boot' > 0 di as text "Bootstrap SE stored in e(V_boot)." ;


end ;


/* ─────────────────────────────────────────────────────────────────────────────
   MATA ENGINE
   #delimit cr required: Mata uses { } and ; internally — must not be
   intercepted by Stata's #delimit ; mode
   ───────────────────────────────────────────────────────────────────────── */
#delimit cr

/* Drop previous definitions so the ado can be sourced more than once */
capture mata: mata drop _epwr_wp()
capture mata: mata drop _epwr_hopt()
capture mata: mata drop _epwr_wquant()
capture mata: mata drop _epwr_htarget()
capture mata: mata drop _wls()
capture mata: mata drop _revCumSum()
capture mata: mata drop _se_IF()
capture mata: mata drop _se_naive()
capture mata: mata drop _se_boot()
capture mata: mata drop _se_svy()
capture mata: mata drop _epwr_main()

mata:

/* ─────────────────────────────────────────────────────────────────────────────
   _epwr_wp()
   Compute gepwe kernel weights and percentile ranks.
   Returns n x 2 matrix : col 1 = w_i,  col 2 = pc_i
   Bandwidth h is written to Stata scalar _epwr_h_tmp

   Kernel (exact gepwe.ado formula):
     w_i = exp(-0.25 * ((pc_i - tau)/h)^2) / (h * sqrt(2*pi) * n)

   Silverman bandwidth on the percentile scale:
     h = cband * min(sd(pc), IQR(pc)/1.34) * n^(-1/5)
   Override with band_in > 0
   ───────────────────────────────────────────────────────────────────────── */
real matrix _epwr_wp(real colvector y,
                      real colvector fw,
                      real scalar    tau,
                      real scalar    cband,
                      real scalar    band_in,
                      real colvector zrank)
{
    /* zrank : alternative ranking variable z (rows=n) or empty (rows=0).
       If empty, rank on y (standard PWR).
       If provided, rank on z (generalised rankvar option).           */
    real scalar    n, h, tmp, q25, q75, sd_pc, i
    real colvector ord, fw_s, pc_s, pc_sorted, u, w_s, w, pc, z_use, z_s

    n         = rows(y)

    /* Use z if provided, else fall back to y */
    z_use     = (rows(zrank) == n) ? zrank : y

    ord       = order(z_use, 1)       /* sort on z (or y if default) */
    fw_s      = fw[ord]
    pc_s      = runningsum(fw_s) :/ sum(fw_s)

    /* Ties: assign each group of equal z the SAME rank = F_n(z) (paper eq 3),
       so p-hat depends only on the value, not on the sort order of ties.
       Removes run-to-run drift in beta and the MSE-optimal bandwidth.       */
    z_s       = z_use[ord]
    i         = n - 1
    while (i >= 1) {
        if (z_s[i] == z_s[i+1]) pc_s[i] = pc_s[i+1]
        i--
    }

    pc_sorted = sort(pc_s, 1)
    q25       = pc_sorted[ceil(0.25 * n)]
    q75       = pc_sorted[ceil(0.75 * n)]
    tmp       = (q75 - q25) / 1.34
    sd_pc     = sqrt(variance(pc_s))
    if (tmp > sd_pc) tmp = sd_pc
    h         = (band_in == 0) ? cband * tmp * n^(-0.2) : band_in

    u         = (pc_s :- tau) :/ h
    w_s       = exp(-0.25 :* u:^2) :/ (h * sqrt(2 * pi()) * n)

    w         = J(n, 1, 0)
    pc        = J(n, 1, 0)
    w[ord]    = w_s
    pc[ord]   = pc_s

    st_numscalar("_epwr_h_tmp", h)
    return((w, pc))
}


/* ─────────────────────────────────────────────────────────────────────────────
   _epwr_hopt()
   MSE-optimal bandwidth via two-step plug-in (see paper, eq. h*).

     h*(tau) = [ sigma2(tau) / (8 * n * sqrt(2*pi) * [mu''(tau)]^2) ]^(1/5)

   Algorithm (paper Algorithm 1):
     1. Pilot h_pilot = 2 * h_Silverman
     2. Estimate local-constant fit of y on percentile rank at
        tau-d, tau, tau+d using h_pilot  (d = 0.02)
     3. mu''(tau) ~ [mu(tau+d) - 2*mu(tau) + mu(tau-d)] / d^2
     4. sigma2(tau) = weighted residual variance at tau (from step 2)
     5. h* clipped to [0.3*h_Silverman, 5*h_Silverman] for stability
        (curvature near zero -> h* would otherwise diverge)

   NOTE: h* is computed using y (or zrank) curvature only -- a single
   scalar bandwidth is then applied to ALL regressors, consistent with
   how Silverman's rule is also depvar-based in this module.
   ───────────────────────────────────────────────────────────────────────── */
real scalar _epwr_hopt(real colvector y,
                           real colvector fw,
                           real scalar    tau,
                           real scalar    cband,
                           real scalar    n)
{
    real scalar    h_sil, h_pilot, d, mu_lo, mu_mid, mu_hi, curv
    real scalar    sigma2, h_star, lo_clip, hi_clip
    real matrix    wp_sil, wp_lo, wp_mid, wp_hi
    real colvector w_sil, w_lo, w_mid, w_hi

    /* Step 0: pilot Silverman bandwidth (band_in=0 triggers Silverman) */
    wp_sil = _epwr_wp(y, fw, tau, cband, 0, J(0,1,.))
    h_sil  = st_numscalar("_epwr_h_tmp")
    h_pilot = 2 * h_sil

    d = 0.02

    /* Step 1: local-constant fits (weighted mean of y) at tau-d, tau, tau+d */
    wp_lo  = _epwr_wp(y, fw, tau - d, cband, h_pilot, J(0,1,.))
    w_lo   = wp_lo[., 1]
    mu_lo  = sum(w_lo :* y) / sum(w_lo)

    wp_mid = _epwr_wp(y, fw, tau,     cband, h_pilot, J(0,1,.))
    w_mid  = wp_mid[., 1]
    mu_mid = sum(w_mid :* y) / sum(w_mid)

    wp_hi  = _epwr_wp(y, fw, tau + d, cband, h_pilot, J(0,1,.))
    w_hi   = wp_hi[., 1]
    mu_hi  = sum(w_hi :* y) / sum(w_hi)

    /* Step 2: curvature (second derivative) of local mean profile */
    curv = (mu_hi - 2*mu_mid + mu_lo) / d^2

    /* Step 3: local residual variance at tau (using pilot weights) */
    sigma2 = sum(w_mid :* (y :- mu_mid):^2) / sum(w_mid)

    /* Step 4: MSE-optimal h* (paper eq. h*); guard against curv ~ 0 */
    if (abs(curv) < 1e-8) curv = 1e-8
    h_star = (sigma2 / (8 * n * sqrt(2*pi()) * curv^2))^(0.2)

    /* Step 5: clip for stability (paper Algorithm 1, step 5) */
    lo_clip = 0.3 * h_sil
    hi_clip = 5.0 * h_sil
    if (h_star < lo_clip) h_star = lo_clip
    if (h_star > hi_clip) h_star = hi_clip

    return(h_star)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _epwr_wquant() : weighted tau-quantile of z
   ───────────────────────────────────────────────────────────────────────── */
real scalar _epwr_wquant(real colvector z, real colvector fw, real scalar tau)
{
    real colvector ord, cs
    real scalar    i, n
    n   = rows(z)
    ord = order(z, 1)
    cs  = runningsum(fw[ord]) :/ sum(fw)
    i   = 1
    while (i < n & cs[i] < tau) i++
    return(z[ord[i]])
}


/* ─────────────────────────────────────────────────────────────────────────────
   _epwr_htarget()  (v1.4)
   Derivative-aware bandwidth for ONE coefficient (or the derivative
   b_j + 2 q_tau b_jsq when the regressor and its square are both given).

   Why: the MSE-optimal rule targets a locally constant coefficient
   (variance ~ 1/(nh)). A slope in the ranking variable has variance
   ~ 1/(nh^3), so its optimal bandwidth is O(n^(-1/7)), wider by a
   factor 3-5 in practice. Used for the MTE from the propensity score.

   Algorithm:
     1. sigma_p (robust scale of ranks) from the Silverman routine;
        pilot  h0 = cpil * sigma_p * n^(-1/7)   (cpil default 1.6)
     2. pilot curve of the target at tau-d, tau, tau+d (band 2*h0, d = .05)
        -> curvature m''(tau) by second differences
     3. grid over h = c * h0, c in tgrid, clipped to [0.01, 0.30]:
        MSE(h) = h^4 * m''^2 + Var_IF(target; h)    (IF variance at each h)
     4. return argmin
   ───────────────────────────────────────────────────────────────────────── */
real scalar _epwr_htarget(real colvector y,
                             real matrix    X,
                             real colvector fw,
                             real scalar    tau,
                             real colvector zrank,
                             real scalar    n,
                             real scalar    cband,
                             real scalar    cpil,
                             string scalar  gridstr,
                             real scalar    jcoef,
                             real scalar    jsq)
{
    real scalar    h_sil, sig_p, h0, hp, d, tlo, thi, i, curv, k, q, qtau
    real scalar    h, best, hbest, mse, var_t
    real rowvector tt, grid
    real colvector m, w, pc, beta, a, z_use, z_ord
    real matrix    wp, V

    k     = cols(X)
    z_use = (rows(zrank) == n) ? zrank : y
    z_ord = order(z_use, 1)

    /* 1. pilot bandwidth in rank space */
    wp    = _epwr_wp(y, fw, tau, cband, 0, zrank)
    h_sil = st_numscalar("_epwr_h_tmp")
    sig_p = h_sil / (cband * n^(-0.2))
    h0    = cpil * sig_p * n^(-1/7)

    /* 2. curvature of the target profile from a wide pilot */
    hp  = 2 * h0
    d   = 0.05
    tlo = max((tau - d, 0.02))
    thi = min((tau + d, 0.98))
    tt  = (tlo, tau, thi)
    m   = J(3, 1, .)
    for (i = 1; i <= 3; i++) {
        wp   = _epwr_wp(y, fw, tt[i], cband, hp, zrank)
        w    = wp[., 1]
        beta = _wls(X, y, w)
        q    = _epwr_wquant(z_use, fw, tt[i])
        m[i] = beta[jcoef] + ((jsq > 0) ? 2 * q * beta[jsq] : 0)
    }
    curv = (m[3] - 2*m[2] + m[1]) / ((thi - tlo)/2)^2

    /* 3. grid search on the empirical MSE */
    qtau  = _epwr_wquant(z_use, fw, tau)
    a     = J(k, 1, 0)
    a[jcoef] = 1
    if (jsq > 0) a[jsq] = 2 * qtau
    grid  = strtoreal(tokens(gridstr))
    best  = .
    hbest = h0
    for (i = 1; i <= cols(grid); i++) {
        h = grid[i] * h0
        if (h < 0.01) h = 0.01
        if (h > 0.30) h = 0.30      /* cap: avoids near-global fits (n_eff ~ n) at no RMSE cost */
        wp    = _epwr_wp(y, fw, tau, cband, h, zrank)
        w     = wp[., 1]
        pc    = wp[., 2]
        beta  = _wls(X, y, w)
        V     = _se_IF(X, y, w, pc, h, tau, beta, z_ord)
        var_t = a' * V * a
        mse   = h^4 * curv^2 + var_t
        if (best == . | mse < best) {
            best  = mse
            hbest = h
        }
    }
    st_numscalar("_epwr_h0", h0)
    st_numscalar("_epwr_m2", curv)
    st_numscalar("_epwr_qt", qtau)
    return(hbest)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _wls() : beta = (X'WX)^{-1} X'Wy
   ───────────────────────────────────────────────────────────────────────── */
real colvector _wls(real matrix X, real colvector y, real colvector w)
{
    /* invsym() is more stable than lusolve() for near-singular X'WX:
       returns generalised inverse (zeros for collinear cols) instead
       of missing values, which would crash ereturn post.             */
    real matrix XtW, A
    XtW = (X :* w)'
    A   = XtW * X
    return(invsym(A) * (XtW * y))
}


/* ─────────────────────────────────────────────────────────────────────────────
   _revCumSum()
   Reverse cumulative row sum of matrix M :
     R[j,.] = M[j,.] + M[j+1,.] + ... + M[n,.]

   Algorithm :
     Forward cumsum column-by-column (quadrunningsum requires a vector).
     revCumSum[i,.] = totalSum - forwardCumSum[i-1,.]
     where forwardCumSum[0,.] = 0 (prepended row of zeros).
   Complexity : O(n*k), one while-loop over k columns
   ───────────────────────────────────────────────────────────────────────── */
real matrix _revCumSum(real matrix M)
{
    real scalar n, k, j
    real matrix fcs, total_row

    n         = rows(M)
    k         = cols(M)
    fcs       = J(n, k, 0)

    /* forward cumulative sum, one column at a time */
    j = 1
    while (j <= k) {
        fcs[., j] = quadrunningsum(M[., j])
        j++
    }

    /* total row sum (1 x k) */
    total_row = fcs[n, .]

    /* revCumSum[i,.] = total - fcs[i-1,.]  (fcs[0,.] = 0 by convention) */
    return(total_row :- (J(1, k, 0) \ fcs[1..(n-1), .]))
}


/* ─────────────────────────────────────────────────────────────────────────────
   _se_IF()
   Influence-function variance (Deville 1999 / Newey-McFadden 1994)

   Total influence of observation j :
     psi_j = w_j * x_j * e_j
             + (1/n) * sum_{i: y_i >= y_j}  dw_i * x_i * e_i

   Kernel derivative :
     dw_i = -0.5 * (pc_i - tau) / h^2 * w_i

   Variance :
     V = A^{-1} * (psi'psi / n) * A^{-1} / n,   A = X'WX / n
   ───────────────────────────────────────────────────────────────────────── */
real matrix _se_IF(real matrix  X,
                   real colvector y,
                   real colvector w,
                   real colvector pc,
                   real scalar    h,
                   real scalar    tau,
                   real colvector beta,
                   real colvector z_ord)
{
    /* z_ord: ordering vector based on ranking variable z.
       If rows(z_ord)=0, fall back to ordering by y (default). */
    real scalar    n, k
    real colvector e, dw, ord
    real matrix    A, Ainv, S_dir, dS, revCS, psi, V

    n          = rows(y)
    k          = cols(X)
    e          = y - X * beta
    A          = (X :* w)' * X / n
    Ainv       = invsym(A)   /* more stable than luinv for near-singular A */
    dw         = (-0.5 / h^2) :* (pc :- tau) :* w
    S_dir      = X :* (w :* e)
    dS         = X :* (dw :* e)
    /* Sort on z if provided, else sort on y */
    ord        = (rows(z_ord) == n) ? z_ord : order(y, 1)
    revCS      = J(n, k, 0)
    revCS[ord,] = _revCumSum(dS[ord,])
    psi        = S_dir + revCS :/ n
    V          = Ainv * (psi' * psi / n) * Ainv / n
    /* Warn if any diagonal of A is near zero (sparse category at tau) */
    if (min(diagonal(A)) < 1e-12) {
        printf("{txt}Warning: near-singular X'WX at tau=%g -- SE unreliable for sparse categories.\n", tau)
    }
    return(V)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _se_naive()
   Standard WLS variance assuming fixed weights.
   INCONSISTENT for PWR -- stored in e(V_naive) for reference only.
   ───────────────────────────────────────────────────────────────────────── */
real matrix _se_naive(real matrix  X,
                      real colvector y,
                      real colvector w,
                      real colvector beta)
{
    real scalar    n, k, s2
    real colvector e
    real matrix    Ainv

    n    = rows(y)
    k    = cols(X)
    e    = y - X * beta
    Ainv = invsym((X :* w)' * X / n)   /* stable for sparse categories */
    s2   = sum(w :* e:^2) / (n - k)
    return(s2 * Ainv / n)
}


/* ─────────────────────────────────────────────────────────────────────────────
   _se_boot()
   Bootstrap standard errors.

   If survey design is declared (do_svy=1):
     Cluster (PSU) bootstrap stratifié — equivalent to:
       bootstrap, strata(strata) cluster(psu) reps(B): gepwreg ...
     Algorithm:
       For each stratum h with n_h PSUs:
         Draw n_h PSUs WITH REPLACEMENT from stratum h
         Stack ALL observations from drawn PSUs (PSU drawn k times → k copies)
       Run gepwreg on the stacked sample with RECOMPUTED kernel weights
     This exactly replicates what Stata's bootstrap prefix does.

   If no survey design (do_svy=0):
     Pairs bootstrap — observations drawn with replacement individually.
     Kernel weights RECONSTRUCTED at every draw.
   ─────────────────────────────────────────────────────────────────────────── */
real matrix _se_boot(real matrix   X,
                     real colvector y,
                     real colvector fw,
                     real colvector w0,
                     real scalar    tau,
                     real scalar    cband,
                     real scalar    band_in,
                     real colvector beta0,
                     real scalar    B,
                     real scalar    do_svy,
                     string scalar  touse,
                     string scalar  svy_psu,
                     string scalar  svy_strata)
{
    real scalar    n, k, b, h_val, p_val, i_h, i_p, n_h, n_b
    real colvector idx, y_b, fw_b, w_b, drawn
    real colvector psu_vec, strata_vec, strata_ids, psu_in_h
    real colvector mask_h, mask_p, obs_idx, idx_b
    real matrix    BB, wp_b, X_b

    n      = rows(y)
    k      = rows(beta0)
    BB     = J(B, k, .)
    rseed(12345)

    if (do_svy & svy_psu != "" & svy_strata != "") {
        /* ── Cluster (PSU) bootstrap within strata ──────────────────────
           Mirrors: bootstrap, strata(strata) cluster(psu): gepwreg ...
           For each replication:
             1. For each stratum: draw n_h PSUs with replacement
             2. Stack observations from selected PSUs
             3. Recompute kernel weights on stacked sample
             4. Run WLS on stacked sample                                */
        st_view(psu_vec,    ., svy_psu,    touse)
        st_view(strata_vec, ., svy_strata, touse)
        strata_ids = uniqrows(strata_vec)

        b = 1
        while (b <= B) {
            /* Build index of selected observations for this replication */
            idx_b = J(0, 1, .)

            i_h = 1
            while (i_h <= rows(strata_ids)) {
                h_val    = strata_ids[i_h]
                mask_h   = (strata_vec :== h_val)
                psu_in_h = uniqrows(select(psu_vec, mask_h))
                n_h      = rows(psu_in_h)

                /* Draw n_h PSUs with replacement */
                drawn = ceil(n_h :* runiform(n_h, 1))

                /* Append observations from each drawn PSU */
                i_p = 1
                while (i_p <= n_h) {
                    p_val   = psu_in_h[drawn[i_p]]
                    mask_p  = mask_h :& (psu_vec :== p_val)
                    obs_idx = select((1::n), mask_p)
                    idx_b   = idx_b \ obs_idx
                    i_p++
                }
                i_h++
            }

            /* Bootstrap sample */
            n_b  = rows(idx_b)
            y_b  = y[idx_b]
            fw_b = fw[idx_b]
            X_b  = X[idx_b, .]

            /* Recompute kernel weights on bootstrap sample */
            wp_b    = _epwr_wp(y_b, fw_b, tau, cband, band_in, J(0,1,.))
            w_b     = wp_b[., 1]
            BB[b, ] = _wls(X_b, y_b, w_b)'
            b++
        }
    }
    else {
        /* ── Pairs bootstrap (no survey design declared) ─────────────── */
        b = 1
        while (b <= B) {
            idx     = ceil(n :* runiform(n, 1))
            y_b     = y[idx]
            fw_b    = fw[idx]
            wp_b    = _epwr_wp(y_b, fw_b, tau, cband, band_in, J(0,1,.))
            w_b     = wp_b[., 1]
            BB[b, ] = _wls(X[idx, ], y_b, w_b)'
            b++
        }
    }
    return(variance(BB))
}


/* ─────────────────────────────────────────────────────────────────────────────
   _se_svy()
   Taylor linearisation variance under stratified cluster design.
     V_svy = A^{-1} [ sum_h n_h/(n_h-1) sum_i (z_hi-z_bar_h)(z_hi-z_bar_h)' ] A^{-1} / n^2
   where z_hi = sum_{j in PSU(h,i)} psi_j  (k-vector of IF scores).
   ───────────────────────────────────────────────────────────────────────── */
real matrix _se_svy(real matrix  X,
                    real colvector y,
                    real colvector w,
                    real colvector pc,
                    real scalar    h,
                    real scalar    tau,
                    real colvector beta,
                    string scalar  touse,
                    string scalar  svy_psu,
                    string scalar  svy_strata,
                    real scalar    n,
                    real colvector z_ord)
{
    real scalar    k, n_h
    real matrix    psi, A, Ainv, meat, dev, Z_h, S_dir, dS, revCS
    real colvector psu_vec, strata_vec, strata_ids, psu_in_h
    real colvector z_bar, mask_h, mask_p, e, dw, ord
    real scalar    h_val, p_val, i_h, i_p

    k    = cols(X)
    A    = (X :* w)' * X / n
    Ainv = invsym(A)   /* stable for near-singular A (sparse PSU or category) */

    /* ── Raw IF scores psi (n x k) ───────────────────────────────────── */
    e     = y - X * beta
    dw    = (-0.5 / h^2) :* (pc :- tau) :* w
    S_dir = X :* (w :* e)
    dS    = X :* (dw :* e)
    ord   = (rows(z_ord) == n) ? z_ord : order(y, 1)
    revCS = J(n, k, 0)
    revCS[ord,] = _revCumSum(dS[ord,])
    psi   = S_dir + revCS :/ n

    /* ── Load PSU and strata identifiers ─────────────────────────────── */
    if (svy_psu != "") {
        st_view(psu_vec, ., svy_psu, touse)
    }
    else {
        /* No PSU declared: treat each observation as its own PSU */
        psu_vec = (1::n)
        printf("{txt}Note: no PSU declared in svyset; treating each obs as its own PSU.\n")
        printf("{txt}      Variance will be conservative. Recommend: svyset psu [pw=fw], strata(strata)\n")
    }

    if (svy_strata != "") {
        st_view(strata_vec, ., svy_strata, touse)
    }
    else {
        /* No strata: single stratum */
        strata_vec = J(n, 1, 1)
    }

    /* ── Taylor meat: sum over strata and PSUs ───────────────────────── */
    meat     = J(k, k, 0)
    strata_ids = uniqrows(strata_vec)

    i_h = 1
    while (i_h <= rows(strata_ids)) {
        h_val   = strata_ids[i_h]
        mask_h  = (strata_vec :== h_val)
        psu_in_h = uniqrows(select(psu_vec, mask_h))
        n_h     = rows(psu_in_h)

        if (n_h < 2) {
            /* Certainty stratum (single PSU): skip */
            i_h++
            continue
        }

        /* Collect z_hi = sum_{j in PSU i, stratum h} psi_j */
        Z_h = J(n_h, k, 0)
        i_p = 1
        while (i_p <= n_h) {
            p_val    = psu_in_h[i_p]
            mask_p   = mask_h :& (psu_vec :== p_val)
            Z_h[i_p,] = colsum(select(psi, mask_p))
            i_p++
        }

        /* Within-stratum deviation */
        z_bar = colsum(Z_h) :/ n_h       /* 1 x k */
        dev   = Z_h :- z_bar                 /* broadcast 1xk to n_h x k */
        meat  = meat + (n_h / (n_h - 1)) * (dev' * dev)
        i_h++
    }

    return(Ainv * meat * Ainv / n^2)
}

void _epwr_main(string scalar depvar,
                 string scalar indepvars,
                 string scalar fw_var,
                 string scalar touse,
                 real scalar   tau,
                 real scalar   cband,
                 real scalar   band_in,
                 real scalar   B,
                 real scalar   addcons,
                 real scalar   do_svy,
                 string scalar svy_psu,
                 string scalar svy_strata,
                 string scalar rank_var,
                 real scalar   do_optbw,
                 real scalar   jcoef,
                 real scalar   jsq,
                 real scalar   cpil,
                 string scalar tgrid)
{
    real colvector  y, fw, w, pc, beta, zrank, z_ord_v
    real matrix     X, wp, V_IF, V_naive, V_boot, V_svy
    real scalar     n, k, h, n_eff, band_use, z_for_h
    real colvector  zrank_for_h

    st_view(y,  ., depvar,  touse)
    st_view(fw, ., fw_var,  touse)
    n = rows(y)

    /* Load ranking variable (empty if same as depvar) */
    zrank = J(0, 1, .)
    if (rank_var != depvar & rank_var != "") {
        st_view(zrank, ., rank_var, touse)
    }

    if (indepvars != "") {
        st_view(X, ., tokens(indepvars), touse)
        if (addcons) X = (X, J(n, 1, 1))
    }
    else X = J(n, 1, 1)
    k = cols(X)

    /* MSE-optimal bandwidth: compute h* and override band_in           */
    /* Curvature/variance are based on the ranking variable (y or zrank) */
    band_use = band_in
    if (do_optbw == 1) {
        zrank_for_h = (rows(zrank) == n) ? zrank : y
        band_use = _epwr_hopt(zrank_for_h, fw, tau, cband, n)
    }
    else if (do_optbw == 2) {
        band_use = _epwr_htarget(y, X, fw, tau, zrank, n, cband, cpil,
                                    tgrid, jcoef, jsq)
    }

    wp    = _epwr_wp(y, fw, tau, cband, band_use, zrank)
    w     = wp[., 1]
    pc    = wp[., 2]
    h     = st_numscalar("_epwr_h_tmp")
    n_eff = sum(w)^2 / sum(w:^2)
    beta  = _wls(X, y, w)

    /* z-based ordering for IF indirect term */
    z_ord_v = (rows(zrank) == n) ? order(zrank, 1) : order(y, 1)

    /* Standard IF variance */
    V_IF    = _se_IF(X, y, w, pc, h, tau, beta, z_ord_v)
    V_naive = _se_naive(X, y, w, beta)

    /* Survey Taylor variance — overrides V_IF as main SE */
    if (do_svy) {
        V_svy = _se_svy(X, y, w, pc, h, tau, beta,
                        touse, svy_psu, svy_strata, n, z_ord_v)
        st_matrix("_epwr_Vs", V_svy)
        V_IF = V_svy
    }

    /* Bootstrap variance */
    if (B > 0) {
        V_boot = _se_boot(X, y, fw, w, tau, cband, band_in, beta, B,
                          do_svy, touse, svy_psu, svy_strata)
        st_matrix("_epwr_Vb", V_boot)
    }

    st_matrix("_epwr_b",  beta')
    st_matrix("_epwr_V",  V_IF)
    st_matrix("_epwr_Vn", V_naive)
    st_numscalar("_epwr_h",    h)
    st_numscalar("_epwr_neff", n_eff)
    st_numscalar("_epwr_n",    n)
}
end
#delimit ;

