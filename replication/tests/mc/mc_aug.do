* mc_aug.do <design> <first> <last> <n> [routes] -- Monte Carlo of the augmented
* two-step (esreg, method(twostep) hermite(3)) against the line (method(twostep)):
* replications first..last, sample size n, one seed per replication.
* routes (default "ln h3"): ln the line, h3 hermite(3), r1 hermite(3 0) (the
* Hermite terms in the treated regime only), r0 hermite(0 3).  Another route on
* the same seeds regenerates the same samples: it is added to the stored ones
* (mc_summary.do <design> <route>) without rerunning them.
*   AN: paper A's baseline N (joint normal: E[w1 - w0 | u] = 1.08 u), where the
*       Hermite terms are zero: the cost of the augmented model.
*   AC: design C1 (normal errors, E[w0 | u] = -0.3 u, E[w1 | u] = 0.78 u +
*       0.54 (u^3 - 3u)), where hermite(3) is exactly specified and the line is
*       not: the benefit.
* Both: x, z, u ~ N(0,1); D = 1{0.2 + 0.3x + 0.9z + u > 0}; Y1 = 1 + 1.2x + w1,
* Y0 = 0.5 + 0.4x + w0.  For each route: ATT, ATU, ATE, kappa, the MTE at
* u = .1 .5 .9 (esrmte); the augmented route also dh2, dh3.
* Writes tests/mc/out/mc_<design>_<first>.dta, or mc_<design>-<routes>_<first>.dta
* for routes other than "ln h3" (mc_summary.do combines the blocks).
args design first last n routes
if ("`routes'" == "") local routes "ln h3"
local rtag = subinstr("`routes'", " ", "", .)
local fn = cond("`routes'" == "ln h3", "mc_`design'_`first'.dta", "mc_`design'-`rtag'_`first'.dta")
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
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p esrmte {
    run "`S'/`f'.ado"
}
if ("`design'" == "AN")      local base 40261002
else if ("`design'" == "AC") local base 50261002
else {
    di as err "design: AN or AC"
    exit 198
}
local vars ""
foreach route of local routes {
    if (!inlist("`route'", "ln", "h3", "r1", "r0")) {
        di as err "routes: ln h3 r1 r0"
        exit 198
    }
    local vars "`vars' `route'_att `route'_atu `route'_ate `route'_kap `route'_m10 `route'_m50 `route'_m90"
    if ("`route'" != "ln") local vars "`vars' `route'_dh2 `route'_dh3"
}
local post ""
foreach v of local vars {
    local post "`post' `v' se_`v'"
}
tempname pf
postfile `pf' long rep double(`post') using "`W'/`fn'", replace
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
    if ("`design'" == "AN") {
        gen double w1 = 0.78*u + 1.04*e1
        gen double w0 = -0.3*u + 0.369230769*e1 + 0.879584356*e0
    }
    else {
        gen double w0 = -0.3*u + 0.95*e0
        gen double w1 = w0 + 1.08*(u + 0.5*(u^3 - 3*u)) + 0.8*e1
    }
    gen double y = cond(d, 1 + 1.2*x + w1, 0.5 + 0.4*x + w0)
    foreach v of local vars {
        local `v' = .
        local se_`v' = .
    }
    foreach route of local routes {
        local hopt ""
        if ("`route'" == "h3") local hopt "hermite(3)"
        if ("`route'" == "r1") local hopt "hermite(3 0)"
        if ("`route'" == "r0") local hopt "hermite(0 3)"
        cap qui esreg y x, select(d = x z) method(twostep) `hopt' nolog
        if (_rc) continue
        foreach s in att atu ate {
            local `route'_`s' = e(`s')
            local se_`route'_`s' = e(se_`s')
        }
        local `route'_kap = e(kappa)
        local se_`route'_kap = e(se_kappa)
        cap qui esrmte, at(.1 .5 .9)
        if (_rc) continue
        tempname M CV VC
        mat `M' = r(mte)
        local j = 0
        foreach p in 10 50 90 {
            local ++j
            local `route'_m`p' = `M'[`j', 2]
            local se_`route'_m`p' = `M'[`j', 3]
        }
        if ("`route'" != "ln") {
            mat `CV' = r(curve)
            mat `VC' = r(V_curve)
            local `route'_dh2 = `CV'[1, 3]
            local se_`route'_dh2 = sqrt(`VC'[3, 3])
            local `route'_dh3 = `CV'[1, 4]
            local se_`route'_dh3 = sqrt(`VC'[4, 4])
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
