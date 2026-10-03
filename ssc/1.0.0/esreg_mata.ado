*! esreg_mata.ado 1.0.0  03oct2026  A. Araar
*! Mata code of the esreg family, in its own file: an ado-file that is autoloaded
*! does not execute a mata: block placed before the program definitions, so the
*! loader (esreg_engine) compiles this file explicitly with -run-, under the
*! session's Stata version.  Functions: _esr_m1 _esr_m0 _esr_dd _esr_agg
*! _esr_wvarmean _esreg_effects _esr_mom _esreg_ts_var _esr_cellmean _esreg_cells
*! _esreg_ddvars _esr_mk _esreg_mte _esr_score_fiml _esreg_hausman _esr_psi
*! _esr_clsum _esr_ifcov_m _esr_wls _esr_pdid_mom _esreg_pdid.
*! Corrected 02oct2026: the standard errors of the effects (ATT, ATU, ATE, mean kappa), of the
*! MTE line and of the cells of esrcurve add the covariance between the sampling
*! part and the parameter part, from the influence function of the whole
*! procedure (_esr_psi; checked by brute force); the two-step variance by cluster
*! and by survey design; the Hausman contrast with the observed Hessian; the
*! pseudo-DiD of esrtest with the X-adjusted Mills contrasts and the influence
*! function of the whole procedure (_esreg_pdid).
*! Corrected 03oct2026: the Hermite controls of the augmented two-step by regime,
*! ph = (ph1, ph0) (a scalar is the same number in both regimes; _esr_phv, _esr_dh).

cap mata: mata drop _esr_psi()
cap mata: mata drop _esr_clsum()
cap mata: mata drop _esr_ifcov_m()
cap mata: mata drop _esr_wls()
cap mata: mata drop _esr_pdid_mom()
cap mata: mata drop _esreg_pdid()
cap mata: mata drop _esr_m1()
cap mata: mata drop _esr_m0()
cap mata: mata drop _esr_ph()
cap mata: mata drop _esr_phv()
cap mata: mata drop _esr_dh()
cap mata: mata drop _esr_ts_blocks()
cap mata: mata drop _esr_hterm()
cap mata: mata drop _esr_tmom()
cap mata: mata drop _esr_tvar()
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

// ---- the numbers of Hermite controls of the current two-step estimation by
// regime, (e(k_h1), e(k_h0)), or e(k_h) in both when the regime counts are absent;
// (0, 0) without hermite, after method(fiml), or in an estimation stored without it
real rowvector _esr_ph()
{
    real matrix p, p1, p0
    if (st_global("e(method)") != "twostep") return((0, 0))
    p1 = st_numscalar("e(k_h1)")
    p0 = st_numscalar("e(k_h0)")
    if (rows(p1) == 1 & rows(p0) == 1) {
        if (p1 < . & p0 < .) return((p1, p0))
    }
    p = st_numscalar("e(k_h)")
    if (rows(p) == 0) return((0, 0))
    if (p >= .) return((0, 0))
    return((p, p))
}

// ---- (ph1, ph0) from a scalar (the same number in both regimes) or a pair
real rowvector _esr_phv(real rowvector ph)
{
    if (cols(ph) == 0) return((0, 0))
    if (cols(ph) == 1) return((ph[1], ph[1]))
    return((ph[1], ph[2]))
}

// ---- the differences h_1k - h_0k of the Hermite coefficients of the two regimes,
// a term absent from a regime counting as zero (max(ph1, ph0) entries)
real colvector _esr_dh(real colvector hc1, real colvector hc0, real rowvector phv)
{
    real colvector a1, a0
    real scalar pm
    pm = max(phv)
    a1 = J(pm, 1, 0)
    a0 = J(pm, 1, 0)
    if (phv[1] > 0) a1[1..phv[1]] = hc1
    if (phv[2] > 0) a0[1..phv[2]] = hc0
    return(a1 - a0)
}

// ---- the blocks of the two-step theta: [g (m) | regime 1: x (k-1), lambda terms
// (pk), Hermite controls (ph1), _cons | regime 0: the same with ph0]; b_j = (x
// coefficients \ _cons), t_j the lambda terms (rho_j sigma_j, by the kappa()
// design), hc_j the coefficients of the Hermite controls of regime j
void _esr_ts_blocks(real colvector th, real scalar k, real scalar m, real scalar pk,
                    real rowvector ph, b1, t1, hc1, b0, t0, hc0)
{
    real scalar i, nb1, nb0
    real rowvector phv
    phv = _esr_phv(ph)
    nb1 = k + pk + phv[1]
    nb0 = k + pk + phv[2]
    i = m
    b1  = (k > 1 ? th[i+1..i+k-1] \ th[i+nb1] : th[i+nb1])
    t1  = th[i+k..i+k+pk-1]
    hc1 = (phv[1] > 0 ? th[i+k+pk..i+k+pk+phv[1]-1] : J(0, 1, .))
    i = i + nb1
    b0  = (k > 1 ? th[i+1..i+k-1] \ th[i+nb0] : th[i+nb0])
    t0  = th[i+k..i+k+pk-1]
    hc0 = (phv[2] > 0 ? th[i+k+pk..i+k+pk+phv[2]-1] : J(0, 1, .))
}

// ---- E[H(u) | D, Z] for the Hermite controls H = (u^2 - 1, u^3 - 3u) (the first
// ph): (-zg l1, (zg^2 - 1) l1) when treated (u > -zg), (zg l0, -(zg^2 - 1) l0) when
// not (u < -zg): the regressors of the two regimes, and the terms of the effects
real matrix _esr_hterm(real colvector d, real colvector zg, real colvector l1,
                       real colvector l0, real scalar ph)
{
    real matrix H1, H0
    H1 = (-zg:*l1, (zg:^2 :- 1):*l1)
    H0 = ( zg:*l0, -(zg:^2 :- 1):*l0)
    return(d:*H1[., 1..ph] + (1:-d):*H0[., 1..ph])
}

// ---- the moments E[u^r | u > a] (upper = 1) or E[u^r | u < a] (upper = 0) of the
// standard normal, r = 0..6 (n x 7): M_r = a^(r-1) lambda + (r-1) M_(r-2), lambda =
// phi(a)/(1 - Phi(a)) above, -phi(a)/Phi(a) below
real matrix _esr_tmom(real colvector a, real scalar upper)
{
    real matrix M
    real colvector lam
    real scalar r
    M = J(rows(a), 7, .)
    lam = (upper ? exp(lnnormalden(a) :- lnnormal(-a)) : -exp(lnnormalden(a) :- lnnormal(a)))
    M[., 1] = J(rows(a), 1, 1)
    M[., 2] = lam
    for (r = 2; r <= 6; r++) M[., r+1] = a:^(r-1):*lam + (r-1)*M[., r-1]
    return(M)
}

// ---- the variance of g(u) = c u + h_2 (u^2 - 1) + h_3 (u^3 - 3u) under the
// truncation whose moments are M (_esr_tmom), by observation (c a column, hc the
// ph Hermite coefficients): E[g^2] - E[g]^2 from the coefficients of g in the
// powers of u, (-h_2, c - 3 h_3, h_2, h_3)
real colvector _esr_tvar(real matrix M, real colvector c, real colvector hc, real scalar ph)
{
    real scalar h2, h3, r, s
    real matrix A
    real colvector Eg, Eg2
    h2 = (ph >= 1 ? hc[1] : 0)
    h3 = (ph >= 2 ? hc[2] : 0)
    A = (J(rows(c), 1, -h2), c :- 3*h3, J(rows(c), 1, h2), J(rows(c), 1, h3))
    Eg = rowsum(A :* M[., 1..4])
    Eg2 = J(rows(c), 1, 0)
    for (r = 1; r <= 4; r++) {
        for (s = 1; s <= 4; s++) Eg2 = Eg2 + A[., r]:*A[., s]:*M[., r+s-1]
    }
    return(Eg2 - Eg:^2)
}

// ---- individual effects and kappa_i from theta -----------------------------
// method 1 (fiml): th = [b1 b0 g a1 a0 c1 c0], W1 = Ws, W2 = Wr
// method 2 (twostep): th = [g  b1 t1 b0 t0],  W1 = Wk   (t_j after b_j, as regress orders them:
//         columns x..., lambda_w..., lambda, [Hermite controls,] _cons; see _esr_ts_blocks)
// With ph = (ph1, ph0) Hermite controls (two-step, hermite()), the conditional
// effect adds dh' E[H(u) | D, Z], dh = h_1 - h_0 (a term absent from a regime is
// zero there): E[w1 - w0 | u] = kappa(x) u + sum_k (h1k - h0k) H_k(u), kappa_i
// staying the linear coefficient (Cov(w1 - w0, u), the Hermite terms being
// orthogonal to u)
real colvector _esr_dd(real scalar method, real colvector th, real matrix X, real matrix Z,
                       real colvector d, real matrix W1, real matrix W2, real colvector kap_i,
                       | real rowvector ph)
{
    real scalar k, m, ps, pr, pk, i
    real colvector b1, b0, g, a1, a0, c1, c0, t1, t0, s1, s0, r1, r0, zg, l1, l0, mg, hc1, hc0, out
    real rowvector phv
    if (args() < 9) ph = 0
    phv = _esr_phv(ph)
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
        pk = cols(W1)
        g  = th[1..m]
        // regime blocks as posted: x(1..k-1), lambda-terms(pk), Hermite(ph_j), _cons
        _esr_ts_blocks(th, k, m, pk, phv, b1, t1, hc1, b0, t0, hc0)
        kap_i = W1*(t1 - t0)
    }
    zg = Z*g; l1 = _esr_m1(zg); l0 = _esr_m0(zg); mg = X*(b1 - b0)
    out = d:*(mg + kap_i:*l1) + (1:-d):*(mg - kap_i:*l0)
    if (method == 2 & max(phv) > 0) out = out + _esr_hterm(d, zg, l1, l0, max(phv))*_esr_dh(hc1, hc0, phv)
    return(out)
}

real rowvector _esr_agg(real scalar method, real colvector th, real matrix X, real matrix Z,
                        real colvector d, real matrix W1, real matrix W2, real colvector w,
                        | real rowvector ph)
{
    real colvector dd, kap_i
    if (args() < 9) ph = 0
    kap_i = .
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i, ph)
    return((sum(w:*d:*dd)/sum(w:*d), sum(w:*(1:-d):*dd)/sum(w:*(1:-d)), sum(w:*dd)/sum(w), sum(w:*kap_i)/sum(w)))
}

real scalar _esr_wvarmean(real colvector v, real colvector w)
{
    real scalar mu
    mu = sum(w:*v)/sum(w)
    return(sum((w:*(v:-mu)):^2)/sum(w)^2)
}

// ---- the influence of each observation on theta, weighted: Psi (n x p) with
// theta_hat - theta ~ sum_i Psi_i'.  FIML: -H^{-1} w_i s_i, H the observed Hessian
// of the weighted log-likelihood (central differences of the analytic score);
// two-step: -Jac^{-1} w_i g_i, Jac the Jacobian of the stacked moments.  It does
// not depend on the vce of the estimation (oim, robust, cluster, svy).
real matrix _esr_psi(real scalar method, real colvector th, real colvector y, real matrix X,
                     real matrix Z, real colvector d, real matrix W1, real matrix W2,
                     real colvector w, | real rowvector ph)
{
    real matrix G, H
    real colvector ev
    real scalar p, i, h
    if (args() < 10) ph = 0
    p = rows(th)
    h = 1e-6
    H = J(p, p, 0)
    if (method == 1) {
        G = _esr_score_fiml(th, y, X, Z, d, W1, W2) :* w
        for (i = 1; i <= p; i++) {
            ev = J(p, 1, 0)
            ev[i] = h
            H[., i] = ((colsum(_esr_score_fiml(th + ev, y, X, Z, d, W1, W2) :* w)
                      - colsum(_esr_score_fiml(th - ev, y, X, Z, d, W1, W2) :* w)) / (2*h))'
        }
    }
    else {
        G = _esr_mom(th, y, X, Z, d, W1, w, ph)
        for (i = 1; i <= p; i++) {
            ev = J(p, 1, 0)
            ev[i] = h
            H[., i] = ((colsum(_esr_mom(th + ev, y, X, Z, d, W1, w, ph))
                      - colsum(_esr_mom(th - ev, y, X, Z, d, W1, w, ph))) / (2*h))'
        }
    }
    return(-G * luinv(H)')
}

// ---- sums of the rows of U within the groups of c (numeric codes)
real matrix _esr_clsum(real matrix U, real colvector c)
{
    real colvector o
    real matrix Uo, info, Us
    real scalar g
    o = order(c, 1)
    Uo = U[o, .]
    info = panelsetup(c[o], 1)
    Us = J(rows(info), cols(U), 0)
    for (g = 1; g <= rows(info); g++) Us[g, .] = colsum(panelsubmatrix(Uo, g, info))
    return(Us)
}

// ---- covariance of the sums of influence functions held in Stata variables:
// independent observations (cross products), or by cluster (sums within the
// clusters, G/(G-1) as ml does); sets the local ncl to the number of clusters
void _esr_ifcov_m(string scalar vl, string scalar touse, string scalar cl, string scalar Cname)
{
    real matrix U, Us
    real scalar G
    U = st_data(., tokens(vl), touse)
    if (cl == "") {
        st_matrix(Cname, cross(U, U))
        st_local("ncl", ".")
    }
    else {
        Us = _esr_clsum(U, st_data(., cl, touse))
        G = rows(Us)
        st_matrix(Cname, cross(Us, Us) * G / (G - 1))
        st_local("ncl", strofreal(G))
    }
}

// The effects ATT, ATU, ATE and the mean kappa.  Ename (4 x 5): est, se, var_param
// (delta method on e(V)), var_samp (sampling part), cov_ps (twice the covariance of
// the two parts), with the sampling and covariance terms for independent
// observations.  With ifvars (8 new variable names), the influence functions are
// stored, U (sampling part, 4) then P (parameter part, 4), so that the caller can
// recompute var_samp and cov_ps by cluster or by survey design (_esreg_ifcov).
void _esreg_effects(string scalar y, string scalar xl, string scalar zl, string scalar dv,
                     string scalar hsl, string scalar hrl, string scalar kl, string scalar wv,
                     string scalar touse, string scalar Ename, | string scalar ifvars)
{
    real matrix X, Z, W1, W2, V, Jm, Vp, E, U, PP
    real colvector yv, d, w, th, dd, kap_i, P, P1, P0, zg, ev, l1, l0, s1, s0, r1, r0
    real colvector rs1, rs0, e1, e0, c1, c0, b1, t1, hc1, b0, t0, hc0
    real matrix R1, R0, H1, H0
    real scalar v1, v0, nb1, nb0
    real rowvector est, samp, cov, ph
    real scalar method, p, i, h, k, m, ps, pr, pk, n
    string scalar mth
    mth = st_global("e(method)")
    method = (mth == "fiml" ? 1 : 2)
    ph = _esr_ph()
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
    est = _esr_agg(method, th, X, Z, d, W1, W2, w, ph)
    kap_i = .
    dd  = _esr_dd(method, th, X, Z, d, W1, W2, kap_i, ph)
    Jm = J(4, p, 0)
    h = 1e-6
    for (i = 1; i <= p; i++) {
        ev = J(p, 1, 0)
        ev[i] = h
        Jm[., i] = ((_esr_agg(method, th + ev, X, Z, d, W1, W2, w, ph) - _esr_agg(method, th - ev, X, Z, d, W1, W2, w, ph)) / (2*h))'
    }
    Vp = Jm * V * Jm'
    // influence functions of the four estimates: the sampling part U (the average
    // over the units) and the parameter part P = Psi J' (the estimation of theta);
    // the variance of the sum is var_param + var_samp + 2 cov, the covariance
    // being zero only by accident (the sampling part depends on (d, z) as the
    // probit score does)
    U = (w:*d:*(dd :- est[1])/sum(w:*d), w:*(1:-d):*(dd :- est[2])/sum(w:*(1:-d)),
         w:*(dd :- est[3])/sum(w), w:*(kap_i :- est[4])/sum(w))
    PP = _esr_psi(method, th, yv, X, Z, d, W1, W2, w, ph) * Jm'
    samp = colsum(U:^2)
    cov = 2*colsum(U:*PP)
    if (args() > 10) {
        if (ifvars != "") {
            (void) st_addvar("double", tokens(ifvars))
            st_store(., tokens(ifvars), touse, (U, PP))
        }
    }
    E = est', sqrt(diagonal(Vp) + samp' + cov'), diagonal(Vp), samp', cov'
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
        nb1 = k + pk + ph[1]
        nb0 = k + pk + ph[2]
        i = m
        c1 = th[i+1..i+nb1]; i = i + nb1
        c0 = th[i+1..i+nb0]
        H1 = (ph[1] > 0 ? (-zg:*l1, (zg:^2 :- 1):*l1)[., 1..ph[1]] : J(n, 0, .))
        H0 = (ph[2] > 0 ? ( zg:*l0, -(zg:^2 :- 1):*l0)[., 1..ph[2]] : J(n, 0, .))
        R1 = (k > 1 ? X[., 1..k-1] : J(n,0,.)), W1:*l1, H1, J(n,1,1)
        R0 = (k > 1 ? X[., 1..k-1] : J(n,0,.)), -W1:*l0, H0, J(n,1,1)
        e1 = yv - R1*c1; e0 = yv - R0*c0
        rs1 = W1*c1[k..k+pk-1]; rs0 = W1*c0[k..k+pk-1]
        if (max(ph) == 0) {
            v1 = sum(w:*d:*(e1:^2 + rs1:^2:*l1:*(l1 + zg)))/sum(w:*d)
            v0 = sum(w:*(1:-d):*(e0:^2 + rs0:^2:*l0:*(l0 - zg)))/sum(w:*(1:-d))
        }
        else {
            // with the Hermite controls, E[w_j | u] = g_j(u) = rs_j u + sum_k h_jk H_k(u):
            // sigma_j^2 = E[resid^2 | D = j] - Var(g_j(u) | truncation) + Var(g_j(u)),
            // Var(g_j(u)) = rs_j^2 + 2 h_j2^2 + 6 h_j3^2 (orthogonal Hermite polynomials);
            // a regime without Hermite controls has h_j = 0 (the linear formula)
            _esr_ts_blocks(th, k, m, pk, ph, b1, t1, hc1, b0, t0, hc0)
            v1 = sum(w:*d:*(e1:^2 - _esr_tvar(_esr_tmom(-zg, 1), rs1, hc1, ph[1])
                     + rs1:^2 :+ (ph[1] >= 1 ? 2*hc1[1]^2 : 0) + (ph[1] >= 2 ? 6*hc1[2]^2 : 0)))/sum(w:*d)
            v0 = sum(w:*(1:-d):*(e0:^2 - _esr_tvar(_esr_tmom(-zg, 0), rs0, hc0, ph[2])
                     + rs0:^2 :+ (ph[2] >= 1 ? 2*hc0[1]^2 : 0) + (ph[2] >= 2 ? 6*hc0[2]^2 : 0)))/sum(w:*(1:-d))
        }
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
                     real colvector d, real matrix Wk, real colvector w, | real rowvector ph)
{
    real scalar k, m, pk, i, n
    real colvector g, c1, c0, zg, l1, l0, e1, e0, gr
    real matrix R1, R0, H1, H0, HH1, HH0
    real rowvector phv
    if (args() < 8) ph = 0
    phv = _esr_phv(ph)
    n = rows(y)
    k = cols(X); m = cols(Z); pk = cols(Wk)
    g = th[1..m]; i = m
    c1 = th[i+1..i+k+pk+phv[1]]; i = i + k + pk + phv[1]
    c0 = th[i+1..i+k+pk+phv[2]]
    zg = Z*g; l1 = _esr_m1(zg); l0 = _esr_m0(zg)
    H1 = J(n, 0, .)
    H0 = J(n, 0, .)
    if (max(phv) > 0) {
        HH1 = (-zg:*l1, (zg:^2 :- 1):*l1)
        HH0 = ( zg:*l0, -(zg:^2 :- 1):*l0)
        if (phv[1] > 0) H1 = HH1[., 1..phv[1]]
        if (phv[2] > 0) H0 = HH0[., 1..phv[2]]
    }
    R1 = (k > 1 ? X[., 1..k-1] : J(n,0,.)), Wk:*l1, H1, J(n,1,1)
    R0 = (k > 1 ? X[., 1..k-1] : J(n,0,.)), -Wk:*l0, H0, J(n,1,1)
    e1 = y - R1*c1; e0 = y - R0*c0
    gr = d:*_esr_m1(zg) - (1:-d):*_esr_m0(zg)
    return((Z:*(w:*gr)), (R1:*(w:*d:*e1)), (R0:*(w:*(1:-d):*e0)))
}

// The variance of the two-step estimator: V = sum_i Psi_i Psi_i' with Psi_i =
// -Jac^{-1} w_i g_i (= Jac^{-1} S Jac^{-1}', the sandwich of the stacked moments).
// With cl (numeric cluster codes), the sums of Psi within the clusters, times
// G/(G-1) as ml does, and the local esr_ncl = G.  With psivars (p new names), Psi
// is stored, for the survey-design variance computed by the caller (_esreg_ifcov).
void _esreg_ts_var(string scalar y, string scalar xl, string scalar zl, string scalar dv,
                    string scalar kl, string scalar wv, string scalar touse,
                    string scalar bname, string scalar Vname, | real rowvector ph,
                    string scalar cl, string scalar psivars)
{
    real matrix X, Z, Wk, Psi, Us
    real colvector yv, d, w, th
    real scalar n, G
    if (args() < 10) ph = 0
    if (args() < 11) cl = ""
    if (args() < 12) psivars = ""
    yv = st_data(., y, touse); d = st_data(., dv, touse); w = st_data(., wv, touse)
    n = rows(yv)
    X  = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse)), J(n, 1, 1)
    Z  = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    Wk = (kl == "" ? J(n, 0, .) : st_data(., tokens(kl), touse)), J(n, 1, 1)
    th = st_matrix(bname)'
    Psi = _esr_psi(2, th, yv, X, Z, d, Wk, J(n, 0, .), w, ph)
    if (cl == "") {
        st_matrix(Vname, cross(Psi, Psi))
    }
    else {
        Us = _esr_clsum(Psi, st_data(., cl, touse))
        G = rows(Us)
        st_matrix(Vname, cross(Us, Us) * G / (G - 1))
        st_local("esr_ncl", strofreal(G))
    }
    if (psivars != "") {
        (void) st_addvar("double", tokens(psivars))
        st_store(., tokens(psivars), touse, Psi)
    }
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
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i, _esr_ph())
    st_store(., ddvar, touse, dd)
    st_store(., kvar, touse, kap_i)
}

// ---- cell means of the individual effect by group (treated / untreated / all)
// g: group index 1..G on the estimation sample.  Writes r(G x 7):
// columns  mean_1 se_1  mean_0 se_0  mean_all se_all  n_all  -- delta method on
// e(V) plus the sampling component within the cell (as in the ATT / ATU table)
real rowvector _esr_cellmean(real scalar method, real colvector th, real matrix X, real matrix Z,
                             real colvector d, real matrix W1, real matrix W2, real colvector w,
                             real colvector g, real scalar G, | real rowvector ph)
{
    real colvector dd, kap_i, w1, w0
    real rowvector out
    real scalar q
    if (args() < 11) ph = 0
    kap_i = .
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i, ph)
    w1 = w:*d; w0 = w:*(1:-d)
    out = J(1, 3*G, .)
    for (q = 1; q <= G; q++) {
        out[q]       = sum((g:==q):*w1:*dd) / sum((g:==q):*w1)
        out[G+q]     = sum((g:==q):*w0:*dd) / sum((g:==q):*w0)
        out[2*G+q]   = sum((g:==q):*w:*dd)  / sum((g:==q):*w)
    }
    return(out)
}

// touse is the estimation sample (the influence function of theta is computed on
// it); g is missing outside the analysis sample.  Rname (G x 7): mean_1 se_1
// mean_0 se_0 mean_all se_all n_all, the standard errors with the sampling part,
// the parameter part and their covariance (independent observations); with ifvars
// (6G new names) the influence functions U (3G) then P (3G) are stored for a
// recomputation by cluster or survey design, and Vname receives the 3G variances
// of the parameter part
void _esreg_cells(string scalar y, string scalar xl, string scalar zl, string scalar dv,
                   string scalar hsl, string scalar hrl, string scalar kl, string scalar wv,
                   string scalar gv, real scalar G, string scalar touse, string scalar Rname,
                   | string scalar ifvars, string scalar Vname)
{
    real matrix X, Z, W1, W2, V, Jm, Vp, R, U, PP
    real colvector yv, d, w, g, th, dd, kap_i, ev, iq, w1, w0
    real rowvector est, samp, cov, ph
    real scalar method, n, p, i, h, q
    method = (st_global("e(method)") == "fiml" ? 1 : 2)
    ph = _esr_ph()
    yv = st_data(., y, touse)
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
    est = _esr_cellmean(method, th, X, Z, d, W1, W2, w, g, G, ph)
    Jm = J(3*G, p, 0)
    h = 1e-6
    for (i = 1; i <= p; i++) {
        ev = J(p, 1, 0)
        ev[i] = h
        Jm[., i] = ((_esr_cellmean(method, th + ev, X, Z, d, W1, W2, w, g, G, ph)
                   - _esr_cellmean(method, th - ev, X, Z, d, W1, W2, w, g, G, ph)) / (2*h))'
    }
    Vp = Jm * V * Jm'
    kap_i = .
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i, ph)
    // influence functions of the 3G cell means: sampling part U, parameter part P
    w1 = w:*d
    w0 = w:*(1:-d)
    U = J(n, 3*G, 0)
    for (q = 1; q <= G; q++) {
        iq = (g:==q)
        U[., q]     = iq:*w1:*(dd :- est[q])     / sum(iq:*w1)
        U[., G+q]   = iq:*w0:*(dd :- est[G+q])   / sum(iq:*w0)
        U[., 2*G+q] = iq:*w:*(dd :- est[2*G+q])  / sum(iq:*w)
    }
    PP = _esr_psi(method, th, yv, X, Z, d, W1, W2, w, ph) * Jm'
    samp = colsum(U:^2)
    cov = 2*colsum(U:*PP)
    if (args() > 12) {
        if (ifvars != "") {
            (void) st_addvar("double", tokens(ifvars))
            st_store(., tokens(ifvars), touse, (U, PP))
            st_matrix(Vname, diagonal(Vp)')
        }
    }
    R = J(G, 7, .)
    for (q = 1; q <= G; q++) {
        R[q, 1] = est[q];     R[q, 2] = sqrt(Vp[q, q] + samp[q] + cov[q])
        R[q, 3] = est[G+q];   R[q, 4] = sqrt(Vp[G+q, G+q] + samp[G+q] + cov[G+q])
        R[q, 5] = est[2*G+q]; R[q, 6] = sqrt(Vp[2*G+q, 2*G+q] + samp[2*G+q] + cov[2*G+q])
        R[q, 7] = sum((g:==q):*w)
    }
    st_matrix(Rname, R)
}
// ---- m = E[X(b1 - b0)] and kappa (means over the sample), for the parametric
// MTE line (_esreg_mte); with ph = (ph1, ph0) Hermite controls (two-step,
// hermite()), also the differences h_1k - h_0k of their coefficients (zero for a
// term absent from a regime): the MTE curve m + kappa v + sum_k (h_1k - h_0k) H_k(v)
real rowvector _esr_mk(real scalar method, real colvector th, real matrix X, real matrix Z,
                       real colvector d, real matrix W1, real matrix W2, real colvector w,
                       | real rowvector ph)
{
    real scalar k, m, ps, pr, pk, i
    real colvector b1, b0, t1, t0, s1, s0, r1, r0, kap_i, mg, hc1, hc0
    real rowvector phv
    if (args() < 9) ph = 0
    phv = _esr_phv(ph)
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
        pk = cols(W1)
        _esr_ts_blocks(th, k, m, pk, phv, b1, t1, hc1, b0, t0, hc0)
        kap_i = W1*(t1 - t0)
    }
    mg = X*(b1 - b0)
    if (method == 2 & max(phv) > 0) return((sum(w:*mg)/sum(w), sum(w:*kap_i)/sum(w), _esr_dh(hc1, hc0, phv)'))
    return((sum(w:*mg)/sum(w), sum(w:*kap_i)/sum(w)))
}

// m and kappa are means over the analysis sample (av = 1, within touse, the
// estimation sample on which the influence function of theta is computed); with
// ph Hermite controls the curve has q = 2 + ph coefficients b = (m, kappa,
// h_12 - h_02[, h_13 - h_03]), MTE(v) = m + kappa v + sum_k b_(k+1) H_k(v).
// Rname = (m, kappa, v_mm, v_kk, v_mk): the covariance of (m, kappa) with the
// parameter part (delta method on e(V)), the sampling part and their covariance
// (independent observations); with ifvars (2q new names) the influence functions
// U_1..U_q P_1..P_q are stored (U = 0 for the Hermite differences, which have no
// sampling part) and Vname receives the q x q parameter part, for a
// recomputation by cluster or survey design; Cname receives (b', V), q x (q + 1),
// V the covariance of b for independent observations
void _esreg_mte(string scalar y, string scalar xl, string scalar zl, string scalar dv,
                string scalar hsl, string scalar hrl, string scalar kl, string scalar wv,
                string scalar touse, string scalar av, string scalar Rname,
                | string scalar ifvars, string scalar Vname, string scalar Cname)
{
    real matrix X, Z, W1, W2, V, Jm, V2, U, PP, C, VT
    real colvector yv, d, w, wf, th, ev, b1, b0, t1, t0, hc1, hc0, mg, dd, kap_i
    real rowvector est, ph
    real scalar method, n, p, i, h, k, m, pk, q
    method = (st_global("e(method)") == "fiml" ? 1 : 2)
    ph = _esr_ph()
    q = 2 + max(ph)
    yv = st_data(., y, touse)
    d = st_data(., dv, touse); wf = st_data(., wv, touse)
    // the weight of the means: zero outside the analysis sample
    w = wf :* st_data(., av, touse)
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
    est = _esr_mk(method, th, X, Z, d, W1, W2, w, ph)
    Jm = J(q, p, 0)
    h = 1e-6
    for (i = 1; i <= p; i++) {
        ev = J(p, 1, 0)
        ev[i] = h
        Jm[., i] = ((_esr_mk(method, th + ev, X, Z, d, W1, W2, w, ph) - _esr_mk(method, th - ev, X, Z, d, W1, W2, w, ph)) / (2*h))'
    }
    V2 = Jm * V * Jm'
    // sampling variance of the mean gain m and of the mean kappa
    k = cols(X); m = cols(Z)
    if (method == 1) {
        b1 = th[1..k]; b0 = th[k+1..2*k]
    }
    else {
        pk = cols(W1)
        _esr_ts_blocks(th, k, m, pk, ph, b1, t1, hc1, b0, t0, hc0)
    }
    mg = X*(b1 - b0)
    kap_i = .
    dd = _esr_dd(method, th, X, Z, d, W1, W2, kap_i, ph)
    U = (w:*(mg :- est[1])/sum(w), w:*(kap_i :- est[2])/sum(w), J(n, max(ph), 0))
    PP = _esr_psi(method, th, yv, X, Z, d, W1, W2, wf, ph) * Jm'
    C = cross((U, PP), (U, PP))
    VT = V2 + C[1..q, 1..q] + C[1..q, q+1..2*q] + C[q+1..2*q, 1..q]
    if (args() > 11) {
        if (ifvars != "") {
            (void) st_addvar("double", tokens(ifvars))
            st_store(., tokens(ifvars), touse, (U, PP))
            st_matrix(Vname, V2)
        }
    }
    if (args() > 13) {
        if (Cname != "") st_matrix(Cname, (est', VT))
    }
    st_matrix(Rname, (est[1], est[2], VT[1,1], VT[2,2], VT[1,2]))
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
// bF: FIML e(b) (constant sigma, rho); bS: two-step e(b) (no kappa vars, no Hermite).
// The influence function of the FIML uses the observed Hessian (_esr_psi), whatever
// the vce of the FIML estimation (with pweights or vce(robust), e(V) is a sandwich
// and e(V)*s_i would not be the influence function).  Writes Rname = (chi2, df, p,
// kappa_2s, kappa_fiml, kdiff, kse) and Qname = (q_2s, q_fiml, se_diff) with
// q = (gamma, b1, rho1sigma1, b0, rho0sigma0).
void _esreg_hausman(string scalar bF, string scalar bS,
                    string scalar y, string scalar xl, string scalar zl, string scalar dv,
                    string scalar wv, string scalar touse, string scalar Rname, string scalar Qname,
                    | string scalar VDname)
{
    real matrix X, Z, W1, thF, thS, psiF, psiS, JF, JS, D, VD, Lk
    real colvector yv, d, w, qF, qS, diff
    real scalar n, k, m, pF, pS, a1, a0, c1, c0, rs1, rs0, qdim, r, j, W, dfree, kdiff, kse
    yv = st_data(., y, touse); d = st_data(., dv, touse); w = st_data(., wv, touse)
    n = rows(yv)
    X = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse)), J(n, 1, 1)
    Z = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    W1 = J(n, 1, 1)
    k = cols(X); m = cols(Z)
    thF = st_matrix(bF)'
    thS = st_matrix(bS)'
    pF = rows(thF)
    pS = rows(thS)
    // influence functions
    psiF = _esr_psi(1, thF, yv, X, Z, d, W1, W1, w)
    psiS = _esr_psi(2, thS, yv, X, Z, d, W1, J(n, 0, .), w, 0)
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
    if (args() >= 11) st_matrix(VDname, VD)
}


// ---- pseudo-DiD on the selection index (esrtest, pdid) ---------------------------
// Weighted least squares
real colvector _esr_wls(real matrix X, real colvector y, real colvector w)
{
    real matrix XtW
    XtW = (X :* w)'
    return(invsym(XtW * X) * (XtW * y))
}

// The estimating equations of the pseudo-DiD, one row per observation (weighted).
// theta = (gamma [m], c_1..c_{G-1}, g_1..g_G [kx+2 each: x, d, 1],
//          b1, b0, a1, a0 [kx+2 each: x, T, 1]):
//   the probit score; the weighted shares of the index below the cut points
//   (levels tau, the shares realised by the strata); in each stratum q the
//   regression of y on (x, d, 1), whose coefficient of d is the gap; in the two
//   extreme strata, for each group j, the regressions of y and of lambda_j(v) on
//   (x, T, 1), T the top-stratum indicator: C_j = b_j[T], pi_j = a_j[T].
// h = 0: the strata by indicators (the estimator); h > 0: the indicators replaced
// by normal((c - v)/h), for the Jacobian (the derivative of the expected moments
// in gamma and in the cut points, which the step functions do not have).
real matrix _esr_pdid_mom(real colvector th, real colvector y, real matrix X, real matrix Z,
                          real colvector d, real colvector w, real colvector tau,
                          real scalar G, real scalar h)
{
    real scalar m, kx, kq, kr, n, i, q
    real colvector g, c, v, T, B, S, l1, l0, gr, b1, b0, a1, a0, gq
    real matrix F, M, out, Q, R
    n = rows(y); m = cols(Z); kx = cols(X)
    kq = kx + 2; kr = kx + 2
    i = 0
    g = th[i+1..i+m]; i = i + m
    c = th[i+1..i+G-1]; i = i + G - 1
    v = Z*g
    F = J(n, G-1, .)
    for (q = 1; q <= G-1; q++) F[., q] = (h > 0 ? normal((c[q] :- v):/h) : (v :<= c[q]))
    M = J(n, G, .)
    M[., 1] = F[., 1]
    for (q = 2; q <= G-1; q++) M[., q] = F[., q] - F[., q-1]
    M[., G] = 1 :- F[., G-1]
    T = M[., G]; B = M[., 1]; S = T + B
    l1 = exp(lnnormalden(v) - lnnormal(v))
    l0 = exp(lnnormalden(v) - lnnormal(-v))
    gr = d:*l1 - (1:-d):*l0
    out = Z :* (w:*gr)
    for (q = 1; q <= G-1; q++) out = out, w:*(F[., q] :- tau[q])
    Q = X, d, J(n, 1, 1)
    for (q = 1; q <= G; q++) {
        gq = th[i+1..i+kq]; i = i + kq
        out = out, Q :* (w:*M[., q]:*(y - Q*gq))
    }
    R = X, T, J(n, 1, 1)
    b1 = th[i+1..i+kr]; i = i + kr
    b0 = th[i+1..i+kr]; i = i + kr
    a1 = th[i+1..i+kr]; i = i + kr
    a0 = th[i+1..i+kr]
    out = out, R :* (w:*d:*S:*(y - R*b1))
    out = out, R :* (w:*(1:-d):*S:*(y - R*b0))
    out = out, R :* (w:*d:*S:*(l1 - R*a1))
    out = out, R :* (w:*(1:-d):*S:*(l0 - R*a0))
    return(out)
}

// The pseudo-DiD and its influence functions.  s: the strata 1..G (xtile of the
// probit index, weighted), gname: the probit e(b).  Rname (2 x (8 + G)): row 1
// the estimates, row 2 their standard errors for independent observations, in the
// order C1 C0 pi1 pi0 rho1sigma1 rho0sigma0 kappa_DD dd gap_1..gap_G, with
// rho1sigma1 = C1/pi1, rho0sigma0 = -C0/pi0, kappa_DD = C1/pi1 + C0/pi0 and dd the
// double difference of the gaps (top - bottom).  With ifvars (8 + G new names) the
// influence functions are stored, for a variance by cluster or survey design.
// The influence function of theta is -Jac^{-1} m_i (just-identified estimating
// equations); in Jac, the derivatives in gamma and in the cut points are those of
// the kernel-smoothed moment sums (bandwidth 1.06 sd(v) n^(-1/5)), the others are
// exact.
void _esreg_pdid(string scalar y, string scalar xl, string scalar zl, string scalar dv,
                 string scalar wv, string scalar sv, string scalar touse, real scalar G,
                 string scalar gname, string scalar Rname, | string scalar ifvars,
                 string scalar Tname)
{
    real scalar n, m, kx, kq, kr, q, p, k, h, hk, eps, base, iC1, iC0, ip1, ip0, ns, mv
    real scalar W, a, b, ma, mb, tq, t, sdq
    real colvector yv, d, w, s, gam, v, c, tau, th, T, B, S, l1, l0, ev, est
    real matrix X, Z, Q, R, m0, Jac, Psi, A, IF, D
    yv = st_data(., y, touse); d = st_data(., dv, touse); w = st_data(., wv, touse)
    s = st_data(., sv, touse)
    n = rows(yv)
    X = (xl == "" ? J(n, 0, .) : st_data(., tokens(xl), touse))
    Z = (zl == "" ? J(n, 0, .) : st_data(., tokens(zl), touse)), J(n, 1, 1)
    m = cols(Z); kx = cols(X); kq = kx + 2; kr = kx + 2
    gam = st_matrix(gname)'
    v = Z*gam
    // cut points between the strata (midpoints) and the shares they realise
    c = J(G-1, 1, .); tau = J(G-1, 1, .)
    for (q = 1; q <= G-1; q++) {
        c[q]   = (max(select(v, s:<=q)) + min(select(v, s:>q))) / 2
        tau[q] = sum(w:*(s:<=q)) / sum(w)
    }
    // the estimates (weighted least squares on the realised strata)
    Q = X, d, J(n, 1, 1)
    th = gam \ c
    for (q = 1; q <= G; q++) th = th \ _esr_wls(Q, yv, w:*(s:==q))
    T = (s:==G); B = (s:==1); S = T + B
    R = X, T, J(n, 1, 1)
    l1 = exp(lnnormalden(v) - lnnormal(v))
    l0 = exp(lnnormalden(v) - lnnormal(-v))
    th = th \ _esr_wls(R, yv, w:*d:*S) \ _esr_wls(R, yv, w:*(1:-d):*S) \
              _esr_wls(R, l1, w:*d:*S) \ _esr_wls(R, l0, w:*(1:-d):*S)
    // the influence function of theta
    m0 = _esr_pdid_mom(th, yv, X, Z, d, w, tau, G, 0)
    mv = sum(w:*v) / sum(w)
    h  = 1.06 * sqrt(sum(w:*(v :- mv):^2) / sum(w)) * n^(-0.2)
    p  = rows(th)
    Jac = J(p, p, 0)
    eps = 1e-6
    // the columns of gamma and of the cut points from the kernel-smoothed moments
    // (the step functions have no derivative there); the columns of the
    // regression coefficients from the moments themselves, which are linear in
    // them (smoothing the strata there would turn the top indicator into a
    // fractional regressor and bias the Jacobian)
    for (k = 1; k <= p; k++) {
        ev = J(p, 1, 0)
        ev[k] = eps
        hk = (k <= m + G - 1 ? h : 0)
        Jac[., k] = (colsum(_esr_pdid_mom(th + ev, yv, X, Z, d, w, tau, G, hk))
                   - colsum(_esr_pdid_mom(th - ev, yv, X, Z, d, w, tau, G, hk)))' / (2*eps)
    }
    Psi = -m0 * luinv(Jac)'
    Psi = Psi :- mean(Psi)
    // the statistics and their gradients in theta
    base = m + (G - 1) + G*kq
    iC1 = base + kx + 1
    iC0 = base + kr + kx + 1
    ip1 = base + 2*kr + kx + 1
    ip0 = base + 3*kr + kx + 1
    ns = 8 + G
    est = J(ns, 1, .)
    A = J(ns, p, 0)
    est[1] = th[iC1];  A[1, iC1] = 1
    est[2] = th[iC0];  A[2, iC0] = 1
    est[3] = th[ip1];  A[3, ip1] = 1
    est[4] = th[ip0];  A[4, ip0] = 1
    est[5] = th[iC1]/th[ip1]
    A[5, iC1] = 1/th[ip1];  A[5, ip1] = -th[iC1]/th[ip1]^2
    est[6] = -th[iC0]/th[ip0]
    A[6, iC0] = -1/th[ip0]; A[6, ip0] = th[iC0]/th[ip0]^2
    est[7] = est[5] - est[6]
    A[7, .] = A[5, .] - A[6, .]
    for (q = 1; q <= G; q++) {
        k = m + (G - 1) + (q - 1)*kq + kx + 1
        est[8 + q] = th[k]
        A[8 + q, k] = 1
    }
    est[8] = est[8 + G] - est[9]
    A[8, .] = A[8 + G, .] - A[9, .]
    IF = Psi * A'
    if (args() > 10) {
        if (ifvars != "") {
            (void) st_addvar("double", tokens(ifvars))
            st_store(., tokens(ifvars), touse, IF)
        }
    }
    st_matrix(Rname, (est' \ sqrt(colsum(IF:^2))))
    st_numscalar("r(pdid_h)", h)
    // the boundaries of the strata: can a resample move a block of tied values
    // of the index across them?  Tname ((G-1) x 5): realised share below the
    // boundary, its target q/G, the shares of the last value below and of the
    // first value above, and a flag (1 when a block above 1% of the sample can
    // move: the target lies within two standard deviations of a break of the
    // empirical distribution next to the boundary)
    if (args() > 11) {
        D = J(G-1, 5, .)
        W = sum(w)
        for (q = 1; q <= G-1; q++) {
            a   = max(select(v, s:<=q))
            b   = min(select(v, s:>q))
            ma  = sum(w:*(v:==a)) / W
            mb  = sum(w:*(v:==b)) / W
            tq  = sum(w:*(s:<=q)) / W
            t   = q / G
            sdq = sqrt(t*(1 - t)/n)
            D[q, .] = (tq, t, ma, mb, ((t - (tq - ma) < 2*sdq & ma > .01) | (tq - t < 2*sdq & mb > .01)))
        }
        st_matrix(Tname, D)
    }
}

end


