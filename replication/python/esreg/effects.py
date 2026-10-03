"""
Treatment effects after the ESR (either route): the four conditional
expectations, ATT / ATU / ATE, kappa, and the common support of P(Z).
Standard errors from the influence function of the whole procedure (esreg
0.5.0): the parameter part (delta method on the parameter covariance), the
sampling part of the averaged individual effects, and twice their covariance,
se^2 = J V J' + sum U^2 + 2 sum U P, with P = Psi J' and Psi the influence of each
observation on theta (FIML: V_oim w s_i; two-step: -Jac^{-1} w g_i).
"""
import numpy as np
from scipy.stats import norm
from .twostep import mills, hermite_controls


def _components(th, fit, X, Z, d):
    lay = fit["layout"]
    if fit["method"] == "fiml":
        b1, b0, g, a1, a0, c1, c0 = lay.unpack(th)
        s1, s0 = np.exp(fit["Ws"] @ a1), np.exp(fit["Ws"] @ a0)
        r1, r0 = np.tanh(fit["Wr"] @ c1), np.tanh(fit["Wr"] @ c0)
        kap_i = r1 * s1 - r0 * s0
    else:
        g, b1, t1, b0, t0 = th[lay.g], th[lay.b1], th[lay.t1], th[lay.b0], th[lay.t0]
        kap_i = fit["Wk"] @ (t1 - t0)
    zg = Z @ g
    l1, l0 = mills(zg)
    mg = X @ (b1 - b0)
    return b1, b0, g, kap_i, zg, l1, l0, mg


def individual_effects(th, fit, X, Z, d):
    b1, b0, g, kap_i, zg, l1, l0, mg = _components(th, fit, X, Z, d)
    dd = np.where(d, mg + kap_i * l1, mg - kap_i * l0)
    if fit["method"] == "twostep" and fit["layout"].ph:
        # E[omega_1 - omega_0 | u] = kappa u + sum_k (h_1k - h_0k) H_k(u): the
        # Hermite terms of the conditional effect, E[H_k(u) | D, Z]; a term absent
        # from a regime (Hermite controls by regime) counts as zero there
        lay = fit["layout"]
        ph = lay.ph
        a1 = np.zeros(ph); a1[:lay.ph1] = th[lay.h1]
        a0 = np.zeros(ph); a0[:lay.ph0] = th[lay.h0]
        dh = a1 - a0
        H1, H0 = hermite_controls(zg, l1, l0)
        dd = dd + np.where(d, H1[:, :ph] @ dh, H0[:, :ph] @ dh)
    return dd, kap_i


def aggregate(th, fit, X, Z, d, w):
    dd, kap_i = individual_effects(th, fit, X, Z, d)
    att = np.average(dd[d], weights=w[d])
    atu = np.average(dd[~d], weights=w[~d])
    ate = np.average(dd, weights=w)
    kap = np.average(kap_i, weights=w)
    return np.array([att, atu, ate, kap]), dd


def treatment_effects(fit, X, Z, d, w=None, h=1e-6):
    n = len(d)
    w = np.ones(n) if w is None else np.asarray(w, float)
    th, V = fit["theta"], fit["V"]
    est, dd = aggregate(th, fit, X, Z, d, w)
    p = len(th)
    J = np.zeros((4, p))
    for i in range(p):
        e = np.zeros(p); e[i] = h
        J[:, i] = (aggregate(th + e, fit, X, Z, d, w)[0] - aggregate(th - e, fit, X, Z, d, w)[0]) / (2 * h)
    Vp = J @ V @ J.T
    b1, b0, g, kap_i, zg, l1, l0, mg = _components(th, fit, X, Z, d)
    # influence functions: the sampling part U (the averages over the units) and
    # the parameter part P = Psi J' (the estimation of theta); their covariance is
    # zero only by accident (the sampling part depends on (d, z) as the probit
    # score does)
    if fit["method"] == "fiml":
        Psi = fit["scores"] @ fit["V_oim"]
    else:
        Psi = -fit["moments"] @ fit["Ginv"].T
    U = np.zeros((n, 4))
    U[:, 0] = np.where(d, w * (dd - est[0]), 0.0) / w[d].sum()
    U[:, 1] = np.where(~d, w * (dd - est[1]), 0.0) / w[~d].sum()
    U[:, 2] = w * (dd - est[2]) / w.sum()
    U[:, 3] = w * (kap_i - est[3]) / w.sum()
    Pm = Psi @ J.T
    samp = (U ** 2).sum(0)
    cov = 2 * (U * Pm).sum(0)
    se = np.sqrt(np.diag(Vp) + samp + cov)
    P = norm.cdf(zg)
    return dict(names=["ATT", "ATU", "ATE", "kappa"], est=est, se=se, z=est / se,
                var_param=np.diag(Vp), var_samp=samp, cov_ps=cov, dd=dd, kappa_i=kap_i, P=P,
                support=(P[d].min(), P[d].max(), P[~d].min(), P[~d].max()),
                common_support=(max(P[d].min(), P[~d].min()), min(P[d].max(), P[~d].max())),
                mean_l1_treated=np.average(l1[d], weights=w[d]),
                mean_l0_untreated=np.average(l0[~d], weights=w[~d]))
