*! esreg_engine 0.4.0  11sep2026  A. Araar
*! Shared engine of the esreg family (esreg, predict, esrcurve, esrmte,
*! esrdiag, esrtest).
*!   Mata  : individual effects, weighted aggregation, delta-method variance,
*!           stacked-moment variance of the two-step route, cell means by group.
*!   Stata : esreg_engine     -- makes sure the Mata functions are defined
*!           _esreg_getest    -- selects the estimation to work from (est(), e(), _esreg)
*!           _esreg_data      -- plain variables (__esr_fv#) behind the five variable lists
*! The Mata block comes first so that -run- of this file (re)defines the
*! functions even when the program definitions below cannot be redefined.

cap mata: mata drop _esr_m1()
cap mata: mata drop _esr_m0()
cap mata: mata drop _esr_dd()
cap mata: mata drop _esr_agg()
cap mata: mata drop _esr_wvarmean()
cap mata: mata drop _esreg_effects()
cap mata: mata drop _esr_mom()
cap mata: mata drop _esreg_ts_var()
cap mata: mata drop _esr_cellmean()
cap mata: mata drop _esreg_cells()
cap mata: mata drop _esreg_ddvars()
cap mata: mata drop _esr_mk()
cap mata: mata drop _esreg_mte()
cap mata: mata drop _esr_score_fiml()
cap mata: mata drop _esreg_hausman()

* ===========================================================================
mata:
mata set matastrict off

real colvector _esr_m1(real colvector a) return(exp(lnnormalden(a) :- lnnormal(a)))
real colvector _esr_m0(real colvector a) return(exp(lnnormalden(a) :- lnnormal(-a)))

// ---- individual effects and kappa_i from theta -----------------------------
// method 1 (fiml): th = [b1 b0 g a1 a0 c1 c0], W1 = Ws, W2 = Wr
// method 2 (twostep): th = [g  b1 t1 b0 t0],  W1 = Wk   (t_j after b_j, as regress orders them:
//         columns x..., lambda_w..., lambda, _cons -> we reorder to X-part then lambda-part)
real colvector _esr_dd(real scalar method, real colvector th, real matrix X, real matrix Z,
                       real colvector d, real matrix W1, real matrix W2, real colvector kap_i)
{
    real scalar k, m, ps, pr, pk, i
    real colvector b1, b0, g, a1, a0, c1, c0, t1, t0, s1, s0, r1, r0, zg, l1, l0, mg
    k = cols(X); m = cols(Z)
    if (method == 1) {
        ps = cols(W1); pr = cols(W2); i = 0
        b1 = th[i+1..i+k]; i = i + k
        b0 = th[i+1..i+k]; i = i + k
        g  = th[i+1..i+m]; i = i + m
        a1 = th[i+1..i+ps]; i = i + ps
        a0 = th[i+1..i+ps]; i = i + ps
        c1 = th[i+1..i+pr]; i = i + pr
        c0 = th[i+1..i+pr]
        s1 = exp(W1*a1); s0 = exp(W1*a0); r1 = tanh(W2*c1); r0 = tanh(W2*c0)
        kap_i = r1:*s1 - r0:*s0
    }
    else {
        pk = cols(W1); i = 0
        g  = th[i+1..i+m]; i = i + m
        // regime 1 block as posted: x(1..k-1), lambda-terms(pk), _cons
        b1 = (k > 1 ? th[i+1..i+k-1] \ th[i+k+pk] : th[i+k+pk])
        t1 = th[i+k..i+k+pk-1]; i = i + k + pk
        b0 = (k > 1 ? th[i+1..i+k-1] \ th[i+k+pk] : th[i+k+pk])
        t0 = th[i+k..i+k+pk-1]
        kap_i = W1*(t1 - t0)
    }
    zg = Z*g; l1 = _esr_m1(zg); l0 = _esr_m0(zg); mg = X*(b1 - b0)
    return(d:*(mg + kap_i:*l1) + (1:-d):*(mg - kap_i:*l0))
}

real rowvector _esr_agg(real scalar method, real colvector th, real matrix X, real matrix Z,
                        real colvector d, real matrix W1, real matrix W2, real colvector w)
{
    real colvector dd, kap_i
    kap_i = .
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i)
    return((sum(w:*d:*dd)/sum(w:*d), sum(w:*(1:-d):*dd)/sum(w:*(1:-d)), sum(w:*dd)/sum(w), sum(w:*kap_i)/sum(w)))
}

real scalar _esr_wvarmean(real colvector v, real colvector w)
{
    real scalar mu
    mu = sum(w:*v)/sum(w)
    return(sum((w:*(v:-mu)):^2)/sum(w)^2)
}

void _esreg_effects(string scalar y, string scalar xl, string scalar zl, string scalar dv,
                     string scalar hsl, string scalar hrl, string scalar kl, string scalar wv,
                     string scalar touse, string scalar Ename)
{
    real matrix X, Z, W1, W2, V, Jm, Vp, E
    real colvector yv, d, w, th, dd, kap_i, P, P1, P0, zg, e, l1, l0, s1, s0, r1, r0
    real colvector rs1, rs0, e1, e0, c1, c0
    real matrix R1, R0
    real scalar v1, v0
    real rowvector est, samp
    real scalar method, p, i, h, k, m, ps, pr, pk, n
    string scalar mth
    mth = st_global("e(method)")
    method = (mth == "fiml" ? 1 : 2)
    yv = st_data(., y, touse); d = st_data(., dv, touse); w = st_data(., wv, touse)
    n = rows(yv)
    X = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse)), J(n, 1, 1)
    Z = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    if (method == 1) {
        W1 = (hsl == "" ? J(n, 0, .) : st_data(., tokens(hsl), touse)), J(n, 1, 1)
        W2 = (hrl == "" ? J(n, 0, .) : st_data(., tokens(hrl), touse)), J(n, 1, 1)
    }
    else {
        W1 = (kl == "" ? J(n, 0, .) : st_data(., tokens(kl), touse)), J(n, 1, 1)
        W2 = J(n, 0, .)
    }
    th = st_matrix("e(b)")'
    V  = st_matrix("e(V)")
    p  = rows(th)
    est = _esr_agg(method, th, X, Z, d, W1, W2, w)
    kap_i = .
    dd  = _esr_dd(method, th, X, Z, d, W1, W2, kap_i)
    Jm = J(4, p, 0)
    h = 1e-6
    for (i = 1; i <= p; i++) {
        e = J(p, 1, 0)
        e[i] = h
        Jm[., i] = ((_esr_agg(method, th + e, X, Z, d, W1, W2, w) - _esr_agg(method, th - e, X, Z, d, W1, W2, w)) / (2*h))'
    }
    Vp = Jm * V * Jm'
    samp = (_esr_wvarmean(select(dd, d), select(w, d)), _esr_wvarmean(select(dd, 1:-d), select(w, 1:-d)), _esr_wvarmean(dd, w), 0)
    E = est', sqrt(diagonal(Vp) + samp'), diagonal(Vp), samp'
    st_matrix(Ename, E)
    // sigma, rho, support, mean lambdas
    k = cols(X); m = cols(Z)
    if (method == 1) {
        ps = cols(W1); pr = cols(W2); i = 2*k + m
        s1 = exp(W1*th[i+1..i+ps]); s0 = exp(W1*th[i+ps+1..i+2*ps]); i = i + 2*ps
        r1 = tanh(W2*th[i+1..i+pr]); r0 = tanh(W2*th[i+pr+1..i+2*pr])
        zg = Z*th[2*k+1..2*k+m]
        st_numscalar("r(sigma1)", sum(w:*s1)/sum(w)); st_numscalar("r(sigma0)", sum(w:*s0)/sum(w))
        st_numscalar("r(rho1)",   sum(w:*r1)/sum(w)); st_numscalar("r(rho0)",   sum(w:*r0)/sum(w))
    }
    else {
        pk = cols(W1)
        zg = Z*th[1..m]
        l1 = _esr_m1(zg); l0 = _esr_m0(zg)
        // rho_j sigma_j (x) and the implied sigma_j (variance corrected for the truncation term)
        i = m
        c1 = th[i+1..i+k+pk]; i = i + k + pk
        c0 = th[i+1..i+k+pk]
        R1 = (k > 1 ? X[., 1..k-1] : J(n,0,.)), W1:*l1, J(n,1,1)
        R0 = (k > 1 ? X[., 1..k-1] : J(n,0,.)), -W1:*l0, J(n,1,1)
        e1 = yv - R1*c1; e0 = yv - R0*c0
        rs1 = W1*c1[k..k+pk-1]; rs0 = W1*c0[k..k+pk-1]
        v1 = sum(w:*d:*(e1:^2 + rs1:^2:*l1:*(l1 + zg)))/sum(w:*d)
        v0 = sum(w:*(1:-d):*(e0:^2 + rs0:^2:*l0:*(l0 - zg)))/sum(w:*(1:-d))
        st_numscalar("r(sigma1)", sqrt(v1)); st_numscalar("r(sigma0)", sqrt(v0))
        st_numscalar("r(rho1)", (sum(w:*d:*rs1)/sum(w:*d))/sqrt(v1))
        st_numscalar("r(rho0)", (sum(w:*(1:-d):*rs0)/sum(w:*(1:-d)))/sqrt(v0))
    }
    P = normal(zg); l1 = _esr_m1(zg); l0 = _esr_m0(zg)
    P1 = select(P, d)
    P0 = select(P, 1:-d)
    st_numscalar("r(p_min1)", min(P1))
    st_numscalar("r(p_max1)", max(P1))
    st_numscalar("r(p_min0)", min(P0))
    st_numscalar("r(p_max0)", max(P0))
    st_numscalar("r(supp_lo)", max((min(select(P, d)), min(select(P, 1:-d)))))
    st_numscalar("r(supp_hi)", min((max(select(P, d)), max(select(P, 1:-d)))))
    st_numscalar("r(ml1)", sum(w:*d:*l1)/sum(w:*d))
    st_numscalar("r(ml0)", sum(w:*(1:-d):*l0)/sum(w:*(1:-d)))
    st_numscalar("r(att)", E[1,1]); st_numscalar("r(se_att)", E[1,2])
    st_numscalar("r(atu)", E[2,1]); st_numscalar("r(se_atu)", E[2,2])
    st_numscalar("r(ate)", E[3,1]); st_numscalar("r(se_ate)", E[3,2])
    st_numscalar("r(kappa)", E[4,1]); st_numscalar("r(se_kappa)", E[4,2])
}

// ---- two-step: stacked moment conditions and sandwich variance ------------
real matrix _esr_mom(real colvector th, real colvector y, real matrix X, real matrix Z,
                     real colvector d, real matrix Wk, real colvector w, | real scalar ph)
{
    real scalar k, m, pk, i, n
    real colvector g, c1, c0, zg, l1, l0, e1, e0, gr
    real matrix R1, R0, H1, H0, HH1, HH0
    if (args() < 8) ph = 0
    n = rows(y)
    k = cols(X); m = cols(Z); pk = cols(Wk)
    g = th[1..m]; i = m
    c1 = th[i+1..i+k+pk+ph]; i = i + k + pk + ph
    c0 = th[i+1..i+k+pk+ph]
    zg = Z*g; l1 = _esr_m1(zg); l0 = _esr_m0(zg)
    H1 = J(n, 0, .)
    H0 = J(n, 0, .)
    if (ph > 0) {
        HH1 = (-zg:*l1, (zg:^2 :- 1):*l1)
        HH0 = ( zg:*l0, -(zg:^2 :- 1):*l0)
        H1 = HH1[., 1..ph]
        H0 = HH0[., 1..ph]
    }
    R1 = (k > 1 ? X[., 1..k-1] : J(n,0,.)), Wk:*l1, H1, J(n,1,1)
    R0 = (k > 1 ? X[., 1..k-1] : J(n,0,.)), -Wk:*l0, H0, J(n,1,1)
    e1 = y - R1*c1; e0 = y - R0*c0
    gr = d:*_esr_m1(zg) - (1:-d):*_esr_m0(zg)
    return((Z:*(w:*gr)), (R1:*(w:*d:*e1)), (R0:*(w:*(1:-d):*e0)))
}

void _esreg_ts_var(string scalar y, string scalar xl, string scalar zl, string scalar dv,
                    string scalar kl, string scalar wv, string scalar touse,
                    string scalar bname, string scalar Vname, | real scalar ph)
{
    if (args() < 10) ph = 0
    real matrix X, Z, Wk, G, S, Jac, Ginv
    real colvector yv, d, w, th, e
    real scalar n, p, i, h
    yv = st_data(., y, touse); d = st_data(., dv, touse); w = st_data(., wv, touse)
    n = rows(yv)
    X  = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse)), J(n, 1, 1)
    Z  = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    Wk = (kl == "" ? J(n, 0, .) : st_data(., tokens(kl), touse)), J(n, 1, 1)
    th = st_matrix(bname)'
    p = rows(th)
    G = _esr_mom(th, yv, X, Z, d, Wk, w, ph)
    S = G'G
    Jac = J(p, p, 0); h = 1e-6
    for (i = 1; i <= p; i++) {
        e = J(p, 1, 0); e[i] = h
        Jac[., i] = ((colsum(_esr_mom(th + e, yv, X, Z, d, Wk, w, ph)) - colsum(_esr_mom(th - e, yv, X, Z, d, Wk, w, ph))) / (2*h))'
    }
    Ginv = luinv(Jac)
    st_matrix(Vname, Ginv * S * Ginv')
}

// ---- individual effect and kappa_i written to two Stata variables ----------
void _esreg_ddvars(string scalar xl, string scalar zl, string scalar dv,
                    string scalar hsl, string scalar hrl, string scalar kl,
                    string scalar touse, string scalar ddvar, string scalar kvar)
{
    real matrix X, Z, W1, W2
    real colvector d, th, dd, kap_i
    real scalar method, n
    method = (st_global("e(method)") == "fiml" ? 1 : 2)
    d = st_data(., dv, touse)
    n = rows(d)
    X = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse)), J(n, 1, 1)
    Z = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    if (method == 1) {
        W1 = (hsl == "" ? J(n, 0, .) : st_data(., tokens(hsl), touse)), J(n, 1, 1)
        W2 = (hrl == "" ? J(n, 0, .) : st_data(., tokens(hrl), touse)), J(n, 1, 1)
    }
    else {
        W1 = (kl == "" ? J(n, 0, .) : st_data(., tokens(kl), touse)), J(n, 1, 1)
        W2 = J(n, 0, .)
    }
    th = st_matrix("e(b)")'
    kap_i = .
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i)
    st_store(., ddvar, touse, dd)
    st_store(., kvar, touse, kap_i)
}

// ---- cell means of the individual effect by group (treated / untreated / all)
// g: group index 1..G on the estimation sample.  Writes r(G x 7):
// columns  mean_1 se_1  mean_0 se_0  mean_all se_all  n_all  -- delta method on
// e(V) plus the sampling component within the cell (as in the ATT / ATU table)
real rowvector _esr_cellmean(real scalar method, real colvector th, real matrix X, real matrix Z,
                             real colvector d, real matrix W1, real matrix W2, real colvector w,
                             real colvector g, real scalar G)
{
    real colvector dd, kap_i, w1, w0
    real rowvector out
    real scalar q
    kap_i = .
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i)
    w1 = w:*d; w0 = w:*(1:-d)
    out = J(1, 3*G, .)
    for (q = 1; q <= G; q++) {
        out[q]       = sum((g:==q):*w1:*dd) / sum((g:==q):*w1)
        out[G+q]     = sum((g:==q):*w0:*dd) / sum((g:==q):*w0)
        out[2*G+q]   = sum((g:==q):*w:*dd)  / sum((g:==q):*w)
    }
    return(out)
}

void _esreg_cells(string scalar xl, string scalar zl, string scalar dv,
                   string scalar hsl, string scalar hrl, string scalar kl, string scalar wv,
                   string scalar gv, real scalar G, string scalar touse, string scalar Rname)
{
    real matrix X, Z, W1, W2, V, Jm, Vp, R
    real colvector d, w, g, th, dd, kap_i, e
    real rowvector est, samp
    real scalar method, n, p, i, h, q
    method = (st_global("e(method)") == "fiml" ? 1 : 2)
    d = st_data(., dv, touse); w = st_data(., wv, touse); g = st_data(., gv, touse)
    n = rows(d)
    X = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse)), J(n, 1, 1)
    Z = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    if (method == 1) {
        W1 = (hsl == "" ? J(n, 0, .) : st_data(., tokens(hsl), touse)), J(n, 1, 1)
        W2 = (hrl == "" ? J(n, 0, .) : st_data(., tokens(hrl), touse)), J(n, 1, 1)
    }
    else {
        W1 = (kl == "" ? J(n, 0, .) : st_data(., tokens(kl), touse)), J(n, 1, 1)
        W2 = J(n, 0, .)
    }
    th = st_matrix("e(b)")'
    V  = st_matrix("e(V)")
    p  = rows(th)
    est = _esr_cellmean(method, th, X, Z, d, W1, W2, w, g, G)
    Jm = J(3*G, p, 0)
    h = 1e-6
    for (i = 1; i <= p; i++) {
        e = J(p, 1, 0)
        e[i] = h
        Jm[., i] = ((_esr_cellmean(method, th + e, X, Z, d, W1, W2, w, g, G)
                   - _esr_cellmean(method, th - e, X, Z, d, W1, W2, w, g, G)) / (2*h))'
    }
    Vp = Jm * V * Jm'
    kap_i = .
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i)
    samp = J(1, 3*G, 0)
    for (q = 1; q <= G; q++) {
        samp[q]     = _esr_wvarmean(select(dd, (g:==q):*d),      select(w, (g:==q):*d))
        samp[G+q]   = _esr_wvarmean(select(dd, (g:==q):*(1:-d)), select(w, (g:==q):*(1:-d)))
        samp[2*G+q] = _esr_wvarmean(select(dd, (g:==q)),         select(w, (g:==q)))
    }
    R = J(G, 7, .)
    for (q = 1; q <= G; q++) {
        R[q, 1] = est[q];     R[q, 2] = sqrt(Vp[q, q] + samp[q])
        R[q, 3] = est[G+q];   R[q, 4] = sqrt(Vp[G+q, G+q] + samp[G+q])
        R[q, 5] = est[2*G+q]; R[q, 6] = sqrt(Vp[2*G+q, 2*G+q] + samp[2*G+q])
        R[q, 7] = sum((g:==q):*w)
    }
    st_matrix(Rname, R)
}
// ---- m = E[X(b1 - b0)] and kappa (means over the sample), joint delta-method
// covariance and the sampling variance of m; for the parametric MTE line.
// Writes Rname = (m, kappa, v_mm, v_kk, v_mk, samp_m)
real rowvector _esr_mk(real scalar method, real colvector th, real matrix X, real matrix Z,
                       real colvector d, real matrix W1, real matrix W2, real colvector w)
{
    real scalar k, m, ps, pr, pk, i
    real colvector b1, b0, t1, t0, s1, s0, r1, r0, kap_i, mg
    k = cols(X); m = cols(Z)
    if (method == 1) {
        ps = cols(W1); pr = cols(W2); i = 0
        b1 = th[i+1..i+k]; i = i + k
        b0 = th[i+1..i+k]; i = i + k + m
        s1 = exp(W1*th[i+1..i+ps]); s0 = exp(W1*th[i+ps+1..i+2*ps]); i = i + 2*ps
        r1 = tanh(W2*th[i+1..i+pr]); r0 = tanh(W2*th[i+pr+1..i+2*pr])
        kap_i = r1:*s1 - r0:*s0
    }
    else {
        pk = cols(W1); i = m
        b1 = (k > 1 ? th[i+1..i+k-1] \ th[i+k+pk] : th[i+k+pk])
        t1 = th[i+k..i+k+pk-1]; i = i + k + pk
        b0 = (k > 1 ? th[i+1..i+k-1] \ th[i+k+pk] : th[i+k+pk])
        t0 = th[i+k..i+k+pk-1]
        kap_i = W1*(t1 - t0)
    }
    mg = X*(b1 - b0)
    return((sum(w:*mg)/sum(w), sum(w:*kap_i)/sum(w)))
}

void _esreg_mte(string scalar xl, string scalar zl, string scalar dv,
                string scalar hsl, string scalar hrl, string scalar kl, string scalar wv,
                string scalar touse, string scalar Rname)
{
    real matrix X, Z, W1, W2, V, Jm, V2
    real colvector d, w, th, e, b1, b0, mg
    real rowvector est
    real scalar method, n, p, i, h, k, m, pk
    method = (st_global("e(method)") == "fiml" ? 1 : 2)
    d = st_data(., dv, touse); w = st_data(., wv, touse)
    n = rows(d)
    X = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse)), J(n, 1, 1)
    Z = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    if (method == 1) {
        W1 = (hsl == "" ? J(n, 0, .) : st_data(., tokens(hsl), touse)), J(n, 1, 1)
        W2 = (hrl == "" ? J(n, 0, .) : st_data(., tokens(hrl), touse)), J(n, 1, 1)
    }
    else {
        W1 = (kl == "" ? J(n, 0, .) : st_data(., tokens(kl), touse)), J(n, 1, 1)
        W2 = J(n, 0, .)
    }
    th = st_matrix("e(b)")'
    V  = st_matrix("e(V)")
    p  = rows(th)
    est = _esr_mk(method, th, X, Z, d, W1, W2, w)
    Jm = J(2, p, 0)
    h = 1e-6
    for (i = 1; i <= p; i++) {
        e = J(p, 1, 0)
        e[i] = h
        Jm[., i] = ((_esr_mk(method, th + e, X, Z, d, W1, W2, w) - _esr_mk(method, th - e, X, Z, d, W1, W2, w)) / (2*h))'
    }
    V2 = Jm * V * Jm'
    // sampling variance of the mean gain m
    k = cols(X); m = cols(Z)
    if (method == 1) {
        b1 = th[1..k]; b0 = th[k+1..2*k]
    }
    else {
        pk = cols(W1)
        b1 = (k > 1 ? th[m+1..m+k-1] \ th[m+k+pk] : th[m+k+pk])
        b0 = (k > 1 ? th[m+k+pk+1..m+2*k+pk-1] \ th[m+2*k+2*pk] : th[m+2*k+2*pk])
    }
    mg = X*(b1 - b0)
    st_matrix(Rname, (est[1], est[2], V2[1,1], V2[2,2], V2[1,2], _esr_wvarmean(mg, w)))
}
// ---- per-observation scores of the Gaussian ESR likelihood, layout [b1 b0 g a1 a0 c1 c0]
real matrix _esr_score_fiml(real colvector th, real colvector y, real matrix X, real matrix Z,
                            real colvector d, real matrix Ws, real matrix Wr)
{
    real scalar k, m, ps, pr, i, n, p
    real colvector b1, b0, g, a1, a0, c1, c0, zg, s1, s0, r1, r0, e1, e0, q1, q0, A1, A0, m1, m0, d1, d0
    real matrix S
    n = rows(y)
    k = cols(X)
    m = cols(Z)
    ps = cols(Ws)
    pr = cols(Wr)
    i = 0
    b1 = th[i+1..i+k]
    i = i + k
    b0 = th[i+1..i+k]
    i = i + k
    g  = th[i+1..i+m]
    i = i + m
    a1 = th[i+1..i+ps]
    i = i + ps
    a0 = th[i+1..i+ps]
    i = i + ps
    c1 = th[i+1..i+pr]
    i = i + pr
    c0 = th[i+1..i+pr]
    p = 2*k + m + 2*ps + 2*pr
    zg = Z*g
    s1 = exp(Ws*a1)
    s0 = exp(Ws*a0)
    r1 = tanh(Wr*c1)
    r0 = tanh(Wr*c0)
    e1 = (y - X*b1):/s1
    e0 = (y - X*b0):/s0
    q1 = sqrt(1 :- r1:^2)
    q0 = sqrt(1 :- r0:^2)
    A1 = (zg + r1:*e1):/q1
    A0 = (zg + r0:*e0):/q0
    m1 = _esr_m1(A1)
    m0 = _esr_m0(A0)
    d1 = d
    d0 = 1 :- d
    S = J(n, p, 0)
    i = 0
    S[., i+1..i+k]  = X:*(d1:*(e1 - m1:*r1:/q1):/s1)
    i = i + k
    S[., i+1..i+k]  = X:*(d0:*(e0 + m0:*r0:/q0):/s0)
    i = i + k
    S[., i+1..i+m]  = Z:*(d1:*m1:/q1 - d0:*m0:/q0)
    i = i + m
    S[., i+1..i+ps] = Ws:*(d1:*(e1:^2 :- 1 :- m1:*r1:*e1:/q1))
    i = i + ps
    S[., i+1..i+ps] = Ws:*(d0:*(e0:^2 :- 1 :+ m0:*r0:*e0:/q0))
    i = i + ps
    S[., i+1..i+pr] = Wr:*(d1:*m1:*(e1 + r1:*zg):/q1)
    i = i + pr
    S[., i+1..i+pr] = Wr:*(-d0:*m0:*(e0 + r0:*zg):/q0)
    return(S)
}

// ---- Hausman contrast FIML vs two-step with influence-function covariance ----------
// bF, VF: FIML e(b), e(V) (constant sigma, rho); bS: two-step e(b) (no kappa vars, no
// Hermite).  Writes Rname = (chi2, df, p, kappa_2s, kappa_fiml, kdiff, kse) and Qname =
// (q_2s, q_fiml, se_diff) with q = (gamma, b1, rho1sigma1, b0, rho0sigma0).
void _esreg_hausman(string scalar bF, string scalar VF, string scalar bS,
                    string scalar y, string scalar xl, string scalar zl, string scalar dv,
                    string scalar wv, string scalar touse, string scalar Rname, string scalar Qname,
                    | string scalar VDname)
{
    real matrix X, Z, W1, thF, thS, V_F, SF, psiF, G, Jac, Ginv, psiS, JF, JS, D, VD, Lk
    real colvector yv, d, w, qF, qS, diff, e
    real scalar n, k, m, pF, pS, i, h, a1, a0, c1, c0, rs1, rs0, qdim, r, j, W, dfree, kdiff, kse
    yv = st_data(., y, touse); d = st_data(., dv, touse); w = st_data(., wv, touse)
    n = rows(yv)
    X = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse)), J(n, 1, 1)
    Z = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    W1 = J(n, 1, 1)
    k = cols(X); m = cols(Z)
    thF = st_matrix(bF)'
    V_F = st_matrix(VF)
    thS = st_matrix(bS)'
    pF = rows(thF)
    pS = rows(thS)
    // influence functions
    SF = _esr_score_fiml(thF, yv, X, Z, d, W1, W1):*w
    psiF = SF * V_F
    G = _esr_mom(thS, yv, X, Z, d, W1, w, 0)
    Jac = J(pS, pS, 0)
    h = 1e-6
    for (i = 1; i <= pS; i++) {
        e = J(pS, 1, 0)
        e[i] = h
        Jac[., i] = ((colsum(_esr_mom(thS + e, yv, X, Z, d, W1, w, 0)) - colsum(_esr_mom(thS - e, yv, X, Z, d, W1, w, 0))) / (2*h))'
    }
    Ginv = luinv(Jac)
    psiS = -(G * Ginv')
    // q and its Jacobians;  FIML layout [b1(k) b0(k) g(m) a1 a0 c1 c0], two-step [g(m) b1(k-1) t1 cons | b0(k-1) t0 cons]
    a1 = thF[2*k+m+1]
    a0 = thF[2*k+m+2]
    c1 = thF[2*k+m+3]
    c0 = thF[2*k+m+4]
    rs1 = tanh(c1)*exp(a1)
    rs0 = tanh(c0)*exp(a0)
    qdim = m + 2*k + 2
    JF = J(qdim, pF, 0); JS = J(qdim, pS, 0)
    r = 0
    for (j = 1; j <= m; j++) {
        r++
        JF[r, 2*k+j] = 1
        JS[r, j] = 1
    }
    // b1: FIML columns 1..k (x..., cons); two-step columns m+1..m+k-1 (x...), cons at m+k+1
    for (j = 1; j <= k; j++) {
        r++
        JF[r, j] = 1
        JS[r, (j < k ? m + j : m + k + 1)] = 1
    }
    r++
    JF[r, 2*k+m+3] = (1 - tanh(c1)^2)*exp(a1)
    JF[r, 2*k+m+1] = rs1
    JS[r, m+k] = 1
    for (j = 1; j <= k; j++) {
        r++
        JF[r, k+j] = 1
        JS[r, (j < k ? m + k + 1 + j : m + 2*k + 2)] = 1
    }
    r++
    JF[r, 2*k+m+4] = (1 - tanh(c0)^2)*exp(a0)
    JF[r, 2*k+m+2] = rs0
    JS[r, m+2*k+1] = 1
    qF = JF * thF
    qS = JS * thS
    // the nonlinear entries of qF are rs1, rs0 themselves (JF*thF is only a linearization)
    qF[m+k+1] = rs1
    qF[m+2*k+2] = rs0
    D = JS * psiS' - JF * psiF'
    VD = D * D'
    diff = qS - qF
    W = diff' * pinv(VD) * diff
    dfree = rank(VD)
    Lk = J(1, qdim, 0); Lk[m+k+1] = 1; Lk[m+2*k+2] = -1
    kdiff = Lk * diff
    kse = sqrt(Lk * VD * Lk')
    st_matrix(Rname, (W, dfree, chi2tail(dfree, W), Lk*qS, Lk*qF, kdiff, kse))
    st_matrix(Qname, (qS, qF, sqrt(diagonal(VD))))
    if (args() >= 12) st_matrix(VDname, VD)
}

end


* ===========================================================================
* Stata-side helpers of the family
* ===========================================================================
cap program drop esreg_engine
program define esreg_engine
    version 16
    mata: st_numscalar("r(ok)", findexternal("_esreg_ddvars()") != NULL)
    if (r(ok) == 0) {
        findfile esreg_engine.ado
        cap noisily run "`r(fn)'"
        mata: st_numscalar("r(ok)", findexternal("_esreg_ddvars()") != NULL)
        if (r(ok) == 0) {
            di as err "the Mata functions of esreg_engine.ado could not be loaded"
            exit 601
        }
    }
end

* select the estimation to work from: est(name) > current e() if esreg > _esreg
* usage:  _esreg_getest [name]   (call after -_estimates hold-, see post-commands)
cap program drop _esreg_getest
program define _esreg_getest
    version 16
    args name
    if ("`name'" != "") {
        cap estimates restore `name'
        if (_rc) {
            di as err "stored estimation `name' not found"
            exit 301
        }
    }
    else if ("`e(cmd)'" != "esreg") {
        cap estimates restore _esreg
        if (_rc) {
            di as err "no esreg estimation in memory; run esreg first or give est(name)"
            exit 301
        }
        di as txt "(using the last esreg estimation, _esreg)"
    }
    if ("`e(cmd)'" != "esreg") {
        di as err "`name' is not an esreg estimation"
        exit 301
    }
    * the data in memory must still be the estimation sample
    qui count if e(sample)
    if (r(N) != e(N)) {
        di as err "the data in memory do not match the estimation sample (e(sample) has " r(N) " obs, e(N) = " e(N) ")"
        exit 459
    }
end

* the plain variable lists of the current esreg estimation, for Mata
* returns r(x) r(z) r(hs) r(hr) r(kap); variables __esr_fv* are created as needed
cap program drop _esreg_data
program define _esreg_data, rclass
    version 16
    syntax [if]
    cap drop __esr_fv*
    local j = 0
    local made_terms ""
    local made_vars ""
    foreach o in xvars zvars hetsigma hetrho kappavars {
        local out ""
        foreach v in `e(`o')' {
            cap confirm variable `v'
            if (_rc == 0) {
                local out "`out' `v'"
                continue
            }
            * the same factor term gets the same plain variable in every list
            local pos : list posof "`v'" in made_terms
            if (`pos' > 0) {
                local out "`out' `: word `pos' of `made_vars''"
                continue
            }
            local ++j
            local made_terms "`made_terms' `v'"
            local made_vars  "`made_vars' __esr_fv`j'"
            fvrevar `v' `if'
            local r "`r(varlist)'"
            qui summarize `r' `if', meanonly
            if (r(min) == r(max)) {
                local lev = substr("`v'", 1, strpos("`v'", ".") - 1)
                local var = substr("`v'", strpos("`v'", ".") + 1, .)
                qui gen double __esr_fv`j' = (`var' == `lev') `if'
            }
            else qui gen double __esr_fv`j' = `r' `if'
            local out "`out' __esr_fv`j'"
        }
        local `o' "`out'"
    }
    return local x   "`xvars'"
    return local z   "`zvars'"
    return local hs  "`hetsigma'"
    return local hr  "`hetrho'"
    return local kap "`kappavars'"
end

* ---------------------------------------------------------------------------
* recompute the effects, kappa, support and ancillary results from the current
* e(b) and e(V): after the svy prefix has replaced e(V) by its linearized
* variance, or on request (esreg, effects); reposts them in e()
cap program drop _esreg_effects_post
program define _esreg_effects_post, eclass
    version 16
    tempvar smp w
    qui gen byte `smp' = e(sample)
    if ("`e(wtype)'" != "") {
        cap qui gen double `w' `e(wexp)' if `smp'
        if (_rc) {
            * under svy the weight passed to esreg was a temporary variable: use the design weight
            cap qui svyset
            if ("`r(wvar)'" != "") qui gen double `w' = `r(wvar)' if `smp'
            else                   qui gen double `w' = 1 if `smp'
        }
    }
    else qui gen double `w' = 1 if `smp'
    _esreg_data if `smp'
    local xl  "`r(x)'"
    local zl  "`r(z)'"
    local hsl "`r(hs)'"
    local hrl "`r(hr)'"
    local kl  "`r(kap)'"
    tempname E A S L
    mata: _esreg_effects("`e(depvar)'", "`xl'", "`zl'", "`e(treat)'", "`hsl'", "`hrl'", "`kl'", "`w'", "`smp'", "`E'")
    cap drop __esr_fv*
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
    local ev = cond("`e(prefix)'" == "svy", "svy", "model")
    ereturn local eff_vce "`ev'"
end
