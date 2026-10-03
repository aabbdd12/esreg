* mc_se.do <design> <first> <last> <n> -- Monte Carlo of the analytic standard
* errors of esreg on the perfect designs (the model holds exactly, every
* parameter strongly identified): replications first..last, sample size n, one
* seed per replication (base + replication number: the blocks can run in
* parallel processes and be combined, the results do not depend on the split).
*   N: paper A's baseline.  x, z, u ~ N(0,1); D = 1{0.2 + 0.3x + 0.9z + u > 0};
*      Y1 = 1 + 1.2x + w1, Y0 = 0.5 + 0.4x + w0, (w1, w0, u) jointly normal with
*      sigma1 = 1.3, sigma0 = 1, rho1 = .6, rho0 = -.3, cov(w1, w0) = .15:
*      kappa = 1.08.  Two-step and FIML: ATT, ATU, ATE, kappa; esrmte at
*      u = .1 .5 .9; esrtest, pdid (kappa_DD, two strata).
*   K: design K of the note.  As N with w1 = (0.78 + 0.4x)u + e1, w0 = -0.3u + e0
*      (e1, e0 independent normal, sd 1.04 and sqrt(.91)): kappa(x) = 1.08 + 0.4x.
*      Two-step with kappa(x): ATT, ATU, ATE, mean kappa, the slope of kappa(x);
*      esrmte at u = .1 .5 .9.
* Writes tests/mc/out/mc_<design>_<first>.dta, one row per replication (a failed
* estimation leaves its statistics missing); mc_summary.do combines the blocks.
args design first last n
version 16
clear all
set more off
* run from the folder that holds data/ and tests/: the root of the development
* project, or replication/ of the public repository (the package in ../src/)
local P = subinstr("`c(pwd)'", "\", "/", .)
local S "`P'/src"
capture confirm file "`S'/esreg.ado"
if (_rc) local S "`P'/../src"
local W "`P'/tests/mc/out"
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrmte esrtest {
    run "`S'/`f'.ado"
}
if ("`design'" == "N") {
    local base 20261002
    local vars "ts_att ts_atu ts_ate ts_kap ts_m10 ts_m50 ts_m90 kdd fi_att fi_atu fi_ate fi_kap fi_m10 fi_m50 fi_m90"
}
else if ("`design'" == "K") {
    local base 30261002
    local vars "kx_att kx_atu kx_ate kx_kap kx_slope kx_m10 kx_m50 kx_m90"
}
else {
    di as err "design: N or K"
    exit 198
}
local post ""
foreach v of local vars {
    local post "`post' `v' se_`v'"
}
tempname pf
postfile `pf' long rep double(`post') using "`W'/mc_`design'_`first'.dta", replace
timer clear 1
timer on 1
forvalues r = `first'/`last' {
    clear
    set seed `=`base' + `r''
    qui set obs `n'
    gen double x  = rnormal()
    gen double z  = rnormal()
    gen double u  = rnormal()
    gen double e1 = rnormal()
    gen double e0 = rnormal()
    gen byte d = (0.2 + 0.3*x + 0.9*z + u > 0)
    if ("`design'" == "N") {
        gen double w1 = 0.78*u + 1.04*e1
        gen double w0 = -0.3*u + 0.369230769*e1 + 0.879584356*e0
    }
    else {
        gen double w1 = (0.78 + 0.4*x)*u + 1.04*e1
        gen double w0 = -0.3*u + 0.953939201*e0
    }
    gen double y = cond(d, 1 + 1.2*x + w1, 0.5 + 0.4*x + w0)
    foreach v of local vars {
        local `v' = .
        local se_`v' = .
    }
    if ("`design'" == "N") {
        foreach route in ts fi {
            local mopt = cond("`route'" == "ts", "method(twostep)", "")
            cap qui esreg y x, select(d = x z) `mopt' nolog
            if (_rc == 0) {
                foreach s in att atu ate {
                    local `route'_`s' = e(`s')
                    local se_`route'_`s' = e(se_`s')
                }
                local `route'_kap = e(kappa)
                local se_`route'_kap = e(se_kappa)
                cap qui esrmte, at(.1 .5 .9)
                if (_rc == 0) {
                    tempname M
                    mat `M' = r(mte)
                    local j = 0
                    foreach p in 10 50 90 {
                        local ++j
                        local `route'_m`p' = `M'[`j', 2]
                        local se_`route'_m`p' = `M'[`j', 3]
                    }
                }
                if ("`route'" == "ts") {
                    cap qui esrtest, pdid
                    if (_rc == 0) {
                        local kdd = r(kappa_dd)
                        local se_kdd = r(se_kappa_dd)
                    }
                }
            }
        }
    }
    else {
        cap qui esreg y x, select(d = x z) method(twostep) kappa(x) nolog
        if (_rc == 0) {
            foreach s in att atu ate {
                local kx_`s' = e(`s')
                local se_kx_`s' = e(se_`s')
            }
            local kx_kap = e(kappa)
            local se_kx_kap = e(se_kappa)
            qui lincom [y_1]lambda_x - [y_0]lambda_x
            local kx_slope = r(estimate)
            local se_kx_slope = r(se)
            cap qui esrmte, at(.1 .5 .9)
            if (_rc == 0) {
                tempname M
                mat `M' = r(mte)
                local j = 0
                foreach p in 10 50 90 {
                    local ++j
                    local kx_m`p' = `M'[`j', 2]
                    local se_kx_m`p' = `M'[`j', 3]
                }
            }
        }
    }
    local pv "(`r')"
    foreach v of local vars {
        local pv "`pv' (``v'') (`se_`v'')"
    }
    post `pf' `pv'
    if (mod(`r' - `first' + 1, 50) == 0) di as txt "mc `design': replication `r' done"
}
postclose `pf'
timer off 1
timer list 1
