"""
Specification tests after the ESR (Section 4 of the esreg paper).

wald_kappa_x(r)  -- constant selection on gains:  H0: kappa(x) = kappa.
    two-step route (kappa(x) = W(t1 - t0)):  H0: (t1 - t0)_{-1} = 0, exact null;
    FIML with hetsigma/hetrho:                H0: all non-constant coefficients of
        ln sigma_j and atanh rho_j are zero (homogeneous regime laws), a
        sufficient condition for kappa(x) constant.
Returns dict(chi2, df, p, slopes, se_slopes, names).
"""
import numpy as np
from scipy.stats import chi2


def wald_kappa_x(r):
    fit = r.fit
    th, V = fit["theta"], fit["V"]
    lay = fit["layout"]
    if fit["method"] == "twostep":
        pk = lay.pk
        if pk < 2:
            raise ValueError("no kappa() variables: nothing to test")
        i1 = np.arange(lay.t1.start, lay.t1.stop - 1)      # non-constant entries of t1
        i0 = np.arange(lay.t0.start, lay.t0.stop - 1)
        R = np.zeros((pk - 1, len(th)))
        R[np.arange(pk - 1), i1] = 1.0
        R[np.arange(pk - 1), i0] = -1.0
        names = [n.split("lambda1*")[1] for n in r.names if "lambda1*" in n][:-1]
    else:
        ps, pr = lay.ps, lay.pr
        idx = []
        for sl, p in ((lay.a1, ps), (lay.a0, ps), (lay.c1, pr), (lay.c0, pr)):
            idx += list(range(sl.start, sl.stop - 1))
        if not idx:
            raise ValueError("no hetsigma()/hetrho() variables: nothing to test")
        R = np.zeros((len(idx), len(th)))
        R[np.arange(len(idx)), idx] = 1.0
        names = [r.names[i] for i in idx]
    d = R @ th
    Vd = R @ V @ R.T
    W = float(d @ np.linalg.solve(Vd, d))
    df = R.shape[0]
    return dict(chi2=W, df=df, p=chi2.sf(W, df), slopes=d, se_slopes=np.sqrt(np.diag(Vd)), names=names)


def wald_sufficiency(df, y, x, d, z, keep=None, weights=None):
    """
    Index sufficiency (Section 4.2):  H0: E[omega_j | X, Z, u] = E[omega_j | X, u].
    The excluded instruments other than `keep` (default: the strongest in the
    probit) are added to the regime regressions of the two-step route; under H0
    their coefficients are zero in both regimes.  Wald on the stacked variance,
    joint over the two regimes; also each regime separately.
    Returns dict(chi2, df, p, chi2_1, df_1, p_1, chi2_0, df_0, p_0, tested, keep, coefs, se).
    """
    from .api import esreg
    from .probit import fit_probit
    import numpy as np
    from scipy.stats import chi2 as _chi2
    excl = [v for v in z if v not in x]
    if len(excl) == 0:
        raise ValueError("no excluded instrument in z")
    if keep is None:
        Z = np.c_[df[list(z)].to_numpy(float), np.ones(len(df))]
        dv = df[d].to_numpy() > 0
        w = df[weights].to_numpy(float) if weights else None
        pr = fit_probit(Z, dv, w)
        g, V = pr["gamma"], pr["V"]
        stats = {v: (g[list(z).index(v)] ** 2) / V[list(z).index(v), list(z).index(v)] for v in excl}
        keep = max(stats, key=stats.get)
    tested = [v for v in excl if v != keep]
    if len(tested) == 0:
        # single instrument: identification by the curvature of lambda alone
        tested = [keep]
        keep = None
    r = esreg(df, y, list(x) + tested, d, list(z), method="twostep", weights=weights)
    lay = r.fit["layout"]
    th, V = r.fit["theta"], r.fit["V"]
    k = len(x) + len(tested) + 1
    i1 = [lay.b1.start + len(x) + j for j in range(len(tested))]
    i0 = [lay.b0.start + len(x) + j for j in range(len(tested))]
    def wald(idx):
        R = np.zeros((len(idx), len(th))); R[np.arange(len(idx)), idx] = 1.0
        dlt = R @ th; Vd = R @ V @ R.T
        W = float(dlt @ np.linalg.solve(Vd, dlt)); return W, len(idx), _chi2.sf(W, len(idx))
    W, dfree, p = wald(i1 + i0)
    W1, d1, p1 = wald(i1)
    W0, d0, p0 = wald(i0)
    return dict(chi2=W, df=dfree, p=p, chi2_1=W1, df_1=d1, p_1=p1, chi2_0=W0, df_0=d0, p_0=p0,
                tested=tested, keep=keep, coefs=th[i1 + i0], se=np.sqrt(np.diag(V))[i1 + i0], fit=r)


def wald_hermite(df, y, x, d, z, weights=None, order=2):
    """
    Linearity of E[omega_j | u] in u (Assumption A2, Section 4.5).  The Hermite
    controls E[u^2 - 1 | D, Z] and E[u^3 - 3u | D, Z] are added beside the Mills
    ratio in each regime; under A2 their coefficients are zero.  Wald on the
    stacked variance, joint and by regime.  order = 1 (quadratic only) or 2.
    """
    from .api import esreg
    import numpy as np
    from scipy.stats import chi2 as _chi2
    r = esreg(df, y, list(x), d, list(z), method="twostep", weights=weights, hermite=order)
    lay = r.fit["layout"]; th, V = r.fit["theta"], r.fit["V"]
    i1 = list(range(lay.h1.start, lay.h1.stop)); i0 = list(range(lay.h0.start, lay.h0.stop))
    def wald(idx):
        R = np.zeros((len(idx), len(th))); R[np.arange(len(idx)), idx] = 1.0
        dlt = R @ th; Vd = R @ V @ R.T
        W = float(dlt @ np.linalg.solve(Vd, dlt)); return W, len(idx), _chi2.sf(W, len(idx))
    W, dfree, p = wald(i1 + i0); W1, d1, p1 = wald(i1); W0, d0, p0 = wald(i0)
    return dict(chi2=W, df=dfree, p=p, chi2_1=W1, df_1=d1, p_1=p1, chi2_0=W0, df_0=d0, p_0=p0,
                coefs=th[i1 + i0], se=np.sqrt(np.diag(V))[i1 + i0], fit=r)


def hausman_normal(df, y, x, d, z, weights=None):
    """
    FIML against two-step (Assumption A3, Section 4.5): Hausman contrast on the
    common parameters q = (gamma, beta_1, rho1 sigma1, beta_0, rho0 sigma0) with the
    covariance of the difference built from the influence functions of the two
    estimators (always positive semi-definite, valid under H0 and H1).
    Returns chi2, df, p on q, and the kappa contrast (kappa_2S - kappa_FIML) with
    its standard error.
    """
    from .api import esreg
    import numpy as np
    from scipy.stats import chi2 as _chi2, norm as _norm
    rF = esreg(df, y, list(x), d, list(z), method="fiml", weights=weights)
    rS = esreg(df, y, list(x), d, list(z), method="twostep", weights=weights)
    fF, fS = rF.fit, rS.fit
    layF, layS = fF["layout"], fS["layout"]
    k, m = layF.k, layF.m
    # influence functions: psi_i (n x p) with theta_hat - theta ~ mean(psi_i)
    psiF = fF["scores"] @ fF["V_oim"]                 # theta_F - theta ~ sum_i s_i (-H)^{-1}
    psiS = -(fS["moments"] @ fS["Ginv"].T)            # theta_S - theta ~ -Jac^{-1} sum_i g_i
    # q from FIML: (g, b1, tanh(c1)exp(a1), b0, tanh(c0)exp(a0)); Jacobian rows
    thF = fF["theta"]
    a1, a0, c1, c0 = thF[layF.a1][0], thF[layF.a0][0], thF[layF.c1][0], thF[layF.c0][0]
    rs1, rs0 = np.tanh(c1) * np.exp(a1), np.tanh(c0) * np.exp(a0)
    qdim = m + 2 * k + 2
    JF = np.zeros((qdim, layF.p)); JS = np.zeros((qdim, layS.p))
    r = 0
    for j in range(m): JF[r, layF.g.start + j] = 1; JS[r, layS.g.start + j] = 1; r += 1
    for j in range(k): JF[r, layF.b1.start + j] = 1; JS[r, layS.b1.start + j] = 1; r += 1
    JF[r, layF.c1.start] = (1 - np.tanh(c1) ** 2) * np.exp(a1); JF[r, layF.a1.start] = rs1; JS[r, layS.t1.start] = 1; r += 1
    for j in range(k): JF[r, layF.b0.start + j] = 1; JS[r, layS.b0.start + j] = 1; r += 1
    JF[r, layF.c0.start] = (1 - np.tanh(c0) ** 2) * np.exp(a0); JF[r, layF.a0.start] = rs0; JS[r, layS.t0.start] = 1; r += 1
    thS = fS["theta"]
    qF = np.r_[thF[layF.g], thF[layF.b1], rs1, thF[layF.b0], rs0]
    qS = np.r_[thS[layS.g], thS[layS.b1], thS[layS.t1], thS[layS.b0], thS[layS.t0]]
    D = JS @ psiS.T - JF @ psiF.T                     # qdim x n, influence of q_S - q_F
    VD = D @ D.T
    diff = qS - qF
    W = float(diff @ np.linalg.pinv(VD) @ diff)
    dfree = int(np.linalg.matrix_rank(VD))
    # kappa contrast
    Lk = np.zeros(qdim); Lk[m + k] = 1; Lk[m + 2 * k + 1] = -1
    kdiff = Lk @ diff; kse = np.sqrt(Lk @ VD @ Lk)
    return dict(chi2=W, df=dfree, p=_chi2.sf(W, dfree), q_2s=qS, q_fiml=qF, V_diff=VD,
                kappa_2s=Lk @ qS, kappa_fiml=Lk @ qF, kappa_diff=kdiff, kappa_se=kse,
                kappa_z=kdiff / kse, kappa_p=2 * _norm.sf(abs(kdiff / kse)), rF=rF, rS=rS)


def _xtile(v, w, nq):
    """Strata 0..nq-1 as Stata's xtile [aw = w]: the cut points are the weighted
    percentiles of _pctile (the average of two order statistics when the
    cumulated weight hits the level exactly), values equal to a cut point go to
    the lower stratum."""
    import numpy as np
    o = np.argsort(v, kind="stable")
    vs, cw = v[o], np.cumsum(w[o])
    W = cw[-1]
    cuts = []
    for q in range(1, nq):
        P = W * q / nq
        i = int(np.searchsorted(cw, P, side="left"))     # first cw >= P
        cuts.append((vs[i] + vs[i + 1]) / 2 if cw[i] == P else vs[i])
    cuts = np.array(cuts)
    return np.sum(v[:, None] > cuts[None, :], axis=1)


def pseudo_did(df, y, x, d, z, nq=2, weights=None):
    """
    Pseudo-DiD on the selection index (Section 4.6).  Strata of the probit
    index v = Z gamma (nq weighted quantile groups); within each stratum the
    X-adjusted gap between treated and untreated (y on X, D).  In the two extreme
    strata and within each group j, the regression of y on (X, T), T the top
    stratum, gives the contrast C_j; the same regression with the Mills ratio
    lambda_j(v) for y gives pi_j, the X-adjusted contrast of lambda_j.  Under A1
    (E[eps_j | u] linear) C_j = rho_j sigma_j pi_j for the treated and
    C_0 = -rho_0 sigma_0 pi_0 for the untreated, whatever the distribution of X
    within the strata: rho1 sigma1 = C1/pi1, rho0 sigma0 = -C0/pi0,
    kappa_dd = C1/pi1 + C0/pi0.  The signs of C1 and C0 need only A2.

    Standard errors: the influence function of the whole procedure, stacked
    just-identified estimating equations (probit score; weighted shares of the
    index below the cut points; the stratum regressions; the four contrast
    regressions); IF = -Jac^{-1} m_i, with the derivatives in gamma and in the
    cut points taken on kernel-smoothed moments (normal kernel, bandwidth
    1.06 sd(v) n^(-1/5)) and the others exact.  se_*_ols: the conventional OLS
    standard errors with the strata taken as given (for comparison).
    """
    import numpy as np
    from scipy.stats import norm
    from .probit import fit_probit
    from .twostep import mills
    n = len(df)
    dv = df[d].to_numpy() > 0
    D = dv.astype(float)
    w = df[weights].to_numpy(float) if weights else np.ones(n)
    Xr = df[list(x)].to_numpy(float) if x else np.zeros((n, 0))
    Z = np.c_[df[list(z)].to_numpy(float), np.ones(n)]
    yv = df[y].to_numpy(float)
    g = fit_probit(Z, dv, w)["gamma"]
    v = Z @ g
    l1, l0 = mills(v)
    s = _xtile(v, w, nq)
    G = int(s.max()) + 1
    m, kx = Z.shape[1], Xr.shape[1]
    kq = kr = kx + 2
    one = np.ones(n)

    def wls(yy, XX, ww):
        A = XX * ww[:, None]
        return np.linalg.solve(A.T @ XX, A.T @ yy)

    def wls_ols_se(yy, XX, ww):
        sel = ww > 0
        b = wls(yy[sel], XX[sel], ww[sel])
        e = yy[sel] - XX[sel] @ b
        A = XX[sel] * np.sqrt(ww[sel])[:, None]
        Vb = np.linalg.inv(A.T @ A) * np.sum(ww[sel] * e ** 2) / (sel.sum() - XX.shape[1])
        return np.sqrt(np.diag(Vb))

    # cut points (midpoints between the strata) and the shares they realise
    c = np.array([(v[s <= q].max() + v[s > q].min()) / 2 for q in range(G - 1)])
    tau = np.array([w[s <= q].sum() / w.sum() for q in range(G - 1)])
    Q = np.c_[Xr, D, one]
    T, B = (s == G - 1).astype(float), (s == 0).astype(float)
    S = T + B
    R = np.c_[Xr, T, one]
    th = [g, c]
    th += [wls(yv, Q, w * (s == q)) for q in range(G)]
    th += [wls(yv, R, w * D * S), wls(yv, R, w * (1 - D) * S),
           wls(l1, R, w * D * S), wls(l0, R, w * (1 - D) * S)]
    th = np.concatenate(th)
    p = th.size

    def mom(t, h):
        i = 0
        gg = t[i:i + m]; i += m
        cc = t[i:i + G - 1]; i += G - 1
        vv = Z @ gg
        F = np.column_stack([norm.cdf((cc[q] - vv) / h) if h > 0 else (vv <= cc[q]).astype(float)
                             for q in range(G - 1)])
        M = np.column_stack([F[:, 0]] + [F[:, q] - F[:, q - 1] for q in range(1, G - 1)] + [1 - F[:, G - 2]])
        TT, SS = M[:, G - 1], M[:, G - 1] + M[:, 0]
        a1_, a0_ = mills(vv)
        cols = [Z * (w * (D * a1_ - (1 - D) * a0_))[:, None]]
        cols += [(w * (F[:, q] - tau[q]))[:, None] for q in range(G - 1)]
        for q in range(G):
            gq = t[i:i + kq]; i += kq
            cols.append(Q * (w * M[:, q] * (yv - Q @ gq))[:, None])
        RR = np.c_[Xr, TT, one]
        b1_, b0_, a1c, a0c = (t[i + j * kr:i + (j + 1) * kr] for j in range(4))
        cols.append(RR * (w * D * SS * (yv - RR @ b1_))[:, None])
        cols.append(RR * (w * (1 - D) * SS * (yv - RR @ b0_))[:, None])
        cols.append(RR * (w * D * SS * (a1_ - RR @ a1c))[:, None])
        cols.append(RR * (w * (1 - D) * SS * (a0_ - RR @ a0c))[:, None])
        return np.hstack(cols)

    m0 = mom(th, 0.0)
    mv = np.sum(w * v) / w.sum()
    h = 1.06 * np.sqrt(np.sum(w * (v - mv) ** 2) / w.sum()) * n ** (-0.2)
    eps = 1e-6
    Jac = np.empty((p, p))
    for k in range(p):
        e = np.zeros(p); e[k] = eps
        hk = h if k < m + G - 1 else 0.0
        Jac[:, k] = (mom(th + e, hk).sum(0) - mom(th - e, hk).sum(0)) / (2 * eps)
    Psi = -m0 @ np.linalg.inv(Jac).T
    Psi -= Psi.mean(0)

    # the statistics and their gradients in theta
    base = m + (G - 1) + G * kq
    iC1, iC0 = base + kx, base + kr + kx
    ip1, ip0 = base + 2 * kr + kx, base + 3 * kr + kx
    igap = [m + (G - 1) + q * kq + kx for q in range(G)]
    C1, C0, pi1, pi0 = th[iC1], th[iC0], th[ip1], th[ip0]
    A = np.zeros((8 + G, p))
    A[0, iC1] = A[1, iC0] = A[2, ip1] = A[3, ip0] = 1
    A[4, iC1], A[4, ip1] = 1 / pi1, -C1 / pi1 ** 2
    A[5, iC0], A[5, ip0] = -1 / pi0, C0 / pi0 ** 2
    A[6] = A[4] - A[5]
    for q in range(G):
        A[8 + q, igap[q]] = 1
    A[7] = A[8 + G - 1] - A[8]
    se = np.sqrt(((Psi @ A.T) ** 2).sum(0))
    rs1, rs0 = C1 / pi1, -C0 / pi0
    gaps = th[igap]

    rows = []
    for q in range(G):
        mq = s == q
        rows.append(dict(stratum=q + 1, n1=int(dv[mq].sum()), n0=int((~dv[mq]).sum()),
                         v_mean=np.average(v[mq], weights=w[mq]), p_mean=norm.cdf(np.average(v[mq], weights=w[mq])),
                         gap=gaps[q], se_gap=se[8 + q], se_gap_ols=wls_ols_se(yv, Q, w * mq)[kx],
                         lam1=np.average(l1[mq & dv], weights=w[mq & dv]),
                         lam0=np.average(l0[mq & ~dv], weights=w[mq & ~dv])))
    return dict(table=rows, C1=C1, se_C1=se[0], C0=C0, se_C0=se[1], pi1=pi1, se_pi1=se[2],
                pi0=pi0, se_pi0=se[3], rhosig1=rs1, se_rhosig1=se[4], rhosig0=rs0, se_rhosig0=se[5],
                kappa_dd=rs1 - rs0, se_kappa_dd=se[6], dd=gaps[-1] - gaps[0], se_dd=se[7],
                dlam1=rows[-1]["lam1"] - rows[0]["lam1"], dlam0=rows[-1]["lam0"] - rows[0]["lam0"],
                se_C1_ols=wls_ols_se(yv, R, w * D * S)[kx], se_C0_ols=wls_ols_se(yv, R, w * (1 - D) * S)[kx],
                bandwidth=h, nq=G)


def link_test(df, d, z, order=2, add=None, weights=None):
    """
    Specification of the selection index (Section 4.1), on the probit alone.
    Default (add=None): Pregibon-type link test.  v = Z gamma_hat from the probit
    of D on Z; a second probit of D on (v, v^2, ..., v^order); H0: the coefficients
    of the powers 2..order are zero.  Wald chi2(order - 1) on the observed
    information.  Rejection: the index is not linear in Z (missing polynomials or
    interactions) or the link is not the normal cdf.
    add = list of columns: RESET-type variant.  The probit of D on Z and the added
    terms; H0: their coefficients are zero.  Wald chi2(#add) and the LR test of
    the nested probits.
    Returns dict(chi2, df, p, coefs, se, names[, lr, p_lr, ll_full, ll_restr]).
    """
    from .api import _design
    from .probit import fit_probit
    from scipy.stats import chi2 as _chi2
    n = len(df)
    w = np.ones(n) if weights is None else df[weights].to_numpy(float)
    dd = df[d].to_numpy(bool)
    Z = _design(df, list(z))
    pr = fit_probit(Z, dd, w)
    if add is None:
        v = Z @ pr["gamma"]
        Zl = np.c_[np.column_stack([v ** j for j in range(1, order + 1)]), np.ones(n)]
        fl = fit_probit(Zl, dd, w)
        idx = list(range(1, order))
        names = [f"v^{j + 1}" for j in idx]
        coefs, V = fl["gamma"], fl["V"]
        out = dict(v_coef=fl["gamma"][0], v_se=np.sqrt(fl["V"][0, 0]))
    else:
        Zf = np.c_[Z[:, :-1], df[list(add)].to_numpy(float), np.ones(n)]
        fl = fit_probit(Zf, dd, w)
        idx = list(range(len(z), len(z) + len(add)))
        names = list(add)
        coefs, V = fl["gamma"], fl["V"]
        lr = 2 * (fl["ll"] - pr["ll"])
        out = dict(lr=lr, p_lr=_chi2.sf(lr, len(add)), ll_full=fl["ll"], ll_restr=pr["ll"])
    b = coefs[idx]
    Vb = V[np.ix_(idx, idx)]
    W = float(b @ np.linalg.solve(Vb, b))
    out.update(chi2=W, df=len(idx), p=_chi2.sf(W, len(idx)), coefs=b, se=np.sqrt(np.diag(Vb)),
               names=names, gamma=pr["gamma"])
    return out


def gamma_contrast(df, y, x, d, z, weights=None):
    """
    What trivariate normality imposes on the selection equation (Section 4.1):
    the gamma block of the Hausman contrast of hausman_normal.  The probit
    estimate of gamma is consistent under A1 alone; the FIML estimate uses the
    outcomes through the correlations and is consistent under A1-A3.  H0: A3;
    Wald chi2 on the gamma block of the influence-function covariance of the
    difference.  Returns dict(chi2, df, p, g_probit, g_fiml, diff, se_diff, hausman).
    """
    from scipy.stats import chi2 as _chi2
    h = hausman_normal(df, y, x, d, z, weights=weights)
    m = h["rF"].fit["layout"].m
    diff = (h["q_2s"] - h["q_fiml"])[:m]
    VD = h["V_diff"][:m, :m]
    W = float(diff @ np.linalg.pinv(VD) @ diff)
    dfree = int(np.linalg.matrix_rank(VD))
    return dict(chi2=W, df=dfree, p=_chi2.sf(W, dfree), g_probit=h["q_2s"][:m], g_fiml=h["q_fiml"][:m],
                diff=diff, se_diff=np.sqrt(np.diag(VD)), hausman=h)
