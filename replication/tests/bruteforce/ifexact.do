* ifexact.do -- the exact influence function of ATT, ATU, ATE and mean kappa for
* esreg (both routes), against the current standard errors (delta + sampling,
* added).  D_i = d(estimate)/d(w_i) = direct sampling term + J * d(theta)/d(w_i),
* with d(theta)/d(w_i) = V_oim * s_i (FIML) or -Jac^{-1} g_i (two-step).
* Writes ifexact_<method>.dta with the D_i of every observation.
args method vce
version 16
clear all
set more off
* run from the folder that holds data/ and tests/: the root of the development
* project, or replication/ of the public repository (the package in ../src/)
local P = subinstr("`c(pwd)'", "\", "/", .)
local S "`P'/src"
capture confirm file "`S'/esreg.ado"
if (_rc) local S "`P'/../src"
local W "`P'/tests/bruteforce/out"
foreach f in esreg esreg_lf1 esreg_engine esreg_mata _esreg_getest _esreg_data _esreg_effects_post _esreg_esample _esreg_ifcov esreg_p {
    run "`S'/`f'.ado"
}
use "`P'/data/esr_sim.dta", clear
keep in 1/1000
gen long id = _n
gen double one = 1
esreg y x, select(d = x z) method(`method') `vce' nolog
matrix E = e(effects)
matrix list E

mata:
void ifexact(string scalar method)
{
    real matrix X, Z, W1, W2, V, Jm, S, G, Jac, dth, D
    real colvector y, d, w, th, ev, dd, kap_i
    real rowvector est, se_cur, se_ex
    real scalar n, p, i, h, meth
    meth = (method == "fiml" ? 1 : 2)
    y = st_data(., "y"); d = st_data(., "d"); w = st_data(., "one")
    n = rows(y)
    X = st_data(., "x"), J(n, 1, 1)
    Z = st_data(., ("x", "z")), J(n, 1, 1)
    W1 = J(n, 1, 1); W2 = (meth == 1 ? J(n, 1, 1) : J(n, 0, .))
    th = st_matrix("e(b)")'
    // d theta/d w_i needs the inverse observed Hessian: e(V_modelbased) after vce(robust)
    V  = (rows(st_matrix("e(V_modelbased)")) ? st_matrix("e(V_modelbased)") : st_matrix("e(V)"))
    p  = rows(th)
    est = _esr_agg(meth, th, X, Z, d, W1, W2, w)
    kap_i = .
    dd = _esr_dd(meth, th, X, Z, d, W1, W2, kap_i)
    // Jacobian of the four estimates in theta (as in _esreg_effects)
    Jm = J(4, p, 0); h = 1e-6
    for (i = 1; i <= p; i++) {
        ev = J(p, 1, 0); ev[i] = h
        Jm[., i] = ((_esr_agg(meth, th + ev, X, Z, d, W1, W2, w) - _esr_agg(meth, th - ev, X, Z, d, W1, W2, w)) / (2*h))'
    }
    // d theta / d w_i
    if (meth == 1) {
        S = _esr_score_fiml(th, y, X, Z, d, W1, W1)
        dth = S * V'                      // V = inverse observed information (vce oim)
    }
    else {
        G = _esr_mom(th, y, X, Z, d, W1, w, 0)
        Jac = J(p, p, 0)
        for (i = 1; i <= p; i++) {
            ev = J(p, 1, 0); ev[i] = h
            Jac[., i] = ((colsum(_esr_mom(th + ev, y, X, Z, d, W1, w, 0)) - colsum(_esr_mom(th - ev, y, X, Z, d, W1, w, 0))) / (2*h))'
        }
        dth = -(G * luinv(Jac)')
    }
    // direct (sampling) terms of the four estimates
    D = J(n, 4, 0)
    D[., 1] = d:*(dd :- est[1]) / sum(w:*d)
    D[., 2] = (1:-d):*(dd :- est[2]) / sum(w:*(1:-d))
    D[., 3] = (dd :- est[3]) / sum(w)
    D[., 4] = (kap_i :- est[4]) / sum(w)
    D = D + dth * Jm'
    se_ex = sqrt(colsum((w:*D):^2))
    se_cur = st_matrix("E")[., 2]'
    printf("\n%s: ATT ATU ATE kappa\n", method)
    printf("  estimate     %9.5f %9.5f %9.5f %9.5f\n", est[1], est[2], est[3], est[4])
    printf("  se current   %9.5f %9.5f %9.5f %9.5f\n", se_cur[1], se_cur[2], se_cur[3], se_cur[4])
    printf("  se exact IF  %9.5f %9.5f %9.5f %9.5f\n", se_ex[1], se_ex[2], se_ex[3], se_ex[4])
    printf("  ratio cur/ex %9.4f %9.4f %9.4f %9.4f\n", se_cur[1]/se_ex[1], se_cur[2]/se_ex[2], se_cur[3]/se_ex[3], se_cur[4]/se_ex[4])
    st_store(., st_addvar("double", ("D_att", "D_atu", "D_ate", "D_kap")), D)
}
end
mata: ifexact("`method'")
keep id D_*
save "`W'/ifexact_`method'`vce'.dta", replace
