"""
Two-step estimation of the ESR (probit, then OLS by regime on X and the Mills
ratio), with the exact variance of the whole procedure from the stacked moment
conditions (probit score; regime-1 normal equations; regime-0 normal equations).

Regime regressions
    D = 1 :  y = X b1 + (Wk t1) * lambda1 + e,   lambda1 = phi(Zg)/Phi(Zg)
    D = 0 :  y = X b0 - (Wk t0) * lambda0 + e,   lambda0 = phi(Zg)/(1 - Phi(Zg))
so that rho_j sigma_j (x) = Wk t_j and kappa(x) = Wk (t1 - t0).
Wk is a column of ones by default (constant kappa).

Parameter vector
    theta = [g (m), b1 (k), t1 (pk), b0 (k), t0 (pk)]
"""
import numpy as np
from scipy.stats import norm
from scipy.special import log_ndtr
from .probit import fit_probit, probit_score_obs


def _ph_pair(ph):
    """(ph1, ph0), the numbers of Hermite controls of the treated and the untreated
    regime, from a number (the same in both) or a pair."""
    if np.isscalar(ph):
        return int(ph), int(ph)
    return int(ph[0]), int(ph[1])


class Layout2:
    def __init__(self, k, m, pk=1, ph=0):
        self.k, self.m, self.pk = k, m, pk
        self.ph1, self.ph0 = _ph_pair(ph)
        self.ph = max(self.ph1, self.ph0)
        i = 0
        self.g = slice(i, i + m); i += m
        self.b1 = slice(i, i + k); i += k
        self.t1 = slice(i, i + pk); i += pk
        self.h1 = slice(i, i + self.ph1); i += self.ph1
        self.b0 = slice(i, i + k); i += k
        self.t0 = slice(i, i + pk); i += pk
        self.h0 = slice(i, i + self.ph0); i += self.ph0
        self.p = i


def hermite_controls(zg, l1, l0):
    """E[u^2 - 1 | D, Z] and E[u^3 - 3u | D, Z] under u ~ N(0,1): the controls that a
    quadratic / cubic E[omega_j | u] adds to the Mills ratio.  Returns (H1, H0), n x 2."""
    H1 = np.c_[-zg * l1, (zg ** 2 - 1.0) * l1]
    H0 = np.c_[zg * l0, -(zg ** 2 - 1.0) * l0]
    return H1, H0


def _sigma2_hermite(e, rs, hc, lo, hi, w):
    """sigma_j^2 with the Hermite controls: mean of e^2 - Var(g(u) | lo < u < hi)
    + Var(g(u)), g(u) = rs u + h2 (u^2 - 1) + h3 (u^3 - 3u), u ~ N(0,1); the moments
    of the truncated normal from scipy (truncnorm), Var(g(u)) = rs^2 + 2 h2^2 + 6 h3^2."""
    from scipy.stats import truncnorm
    h2 = hc[0] if len(hc) >= 1 else 0.0
    h3 = hc[1] if len(hc) >= 2 else 0.0
    a = np.broadcast_to(lo, e.shape).astype(float)
    b = np.broadcast_to(hi, e.shape).astype(float)
    M = np.column_stack([np.ones_like(e)] + [truncnorm.moment(r, a, b) for r in range(1, 7)])
    A = np.column_stack([np.full_like(e, -h2), rs - 3 * h3, np.full_like(e, h2), np.full_like(e, h3)])
    Eg = (A * M[:, :4]).sum(1)
    Eg2 = sum(A[:, r] * A[:, s] * M[:, r + s] for r in range(4) for s in range(4))
    vt = Eg2 - Eg ** 2
    return np.average(e ** 2 - vt + rs ** 2 + 2 * h2 ** 2 + 6 * h3 ** 2, weights=w)


def mills(zg):
    l1 = np.exp(norm.logpdf(zg) - log_ndtr(zg))
    l0 = np.exp(norm.logpdf(zg) - log_ndtr(-zg))
    return l1, l0


def _regime_designs(g, X, Z, Wk, ph):
    ph1, ph0 = _ph_pair(ph)
    zg = Z @ g
    l1, l0 = mills(zg)
    R1 = np.c_[X, Wk * l1[:, None]]
    R0 = np.c_[X, -Wk * l0[:, None]]
    if ph1 or ph0:
        H1, H0 = hermite_controls(zg, l1, l0)
        if ph1:
            R1 = np.c_[R1, H1[:, :ph1]]
        if ph0:
            R0 = np.c_[R0, H0[:, :ph0]]
    return R1, R0, l1, l0, zg


def moments_obs(th, lay, y, X, Z, d, Wk, w):
    """n x p matrix of per-observation moment contributions g_i(theta)."""
    g = th[lay.g]
    R1, R0, l1, l0, zg = _regime_designs(g, X, Z, Wk, (lay.ph1, lay.ph0))
    e1 = y - R1 @ th[lay.b1.start:lay.h1.stop]
    e0 = y - R0 @ th[lay.b0.start:lay.h0.stop]
    d1 = d.astype(float); d0 = 1.0 - d1
    G = np.zeros((len(y), lay.p))
    G[:, lay.g] = probit_score_obs(g, Z, d, w)
    G[:, lay.b1.start:lay.h1.stop] = R1 * (w * d1 * e1)[:, None]
    G[:, lay.b0.start:lay.h0.stop] = R0 * (w * d0 * e0)[:, None]
    return G


def fit_twostep(y, X, Z, d, w=None, Wk=None, h=1e-6, hermite=0):
    """hermite = 0, 1 or 2: number of Hermite controls (u^2 - 1, u^3 - 3u) added to
    each regime regression beside the Mills ratio (Section 4.5 test of A2); or a
    pair (treated, untreated), a number for each regime (Stata's hermite(#1 #0)
    with # = the number + 1, 0 for none)."""
    n, k = X.shape
    m = Z.shape[1]
    w = np.ones(n) if w is None else np.asarray(w, float)
    Wk = np.ones((n, 1)) if Wk is None else Wk
    lay = Layout2(k, m, Wk.shape[1], hermite)

    pr = fit_probit(Z, d, w)
    g = pr["gamma"]
    R1, R0, l1, l0, zg = _regime_designs(g, X, Z, Wk, hermite)
    sw = np.sqrt(w)
    c1 = np.linalg.lstsq(R1[d] * sw[d, None], y[d] * sw[d], rcond=None)[0]
    c0 = np.linalg.lstsq(R0[~d] * sw[~d, None], y[~d] * sw[~d], rcond=None)[0]
    th = np.zeros(lay.p)
    th[lay.g] = g
    pk = Wk.shape[1]
    th[lay.b1], th[lay.t1], th[lay.h1] = c1[:k], c1[k:k + pk], c1[k + pk:]
    th[lay.b0], th[lay.t0], th[lay.h0] = c0[:k], c0[k:k + pk], c0[k + pk:]

    # stacked sandwich: V = Ginv S Ginv',  G = d(sum g_i)/d theta',  S = sum g_i g_i'
    Gm = moments_obs(th, lay, y, X, Z, d, Wk, w)
    S = Gm.T @ Gm
    Jac = np.zeros((lay.p, lay.p))
    for j in range(lay.p):
        e = np.zeros(lay.p); e[j] = h
        Jac[:, j] = (moments_obs(th + e, lay, y, X, Z, d, Wk, w).sum(0)
                     - moments_obs(th - e, lay, y, X, Z, d, Wk, w).sum(0)) / (2 * h)
    Ginv = np.linalg.inv(Jac)
    V = Ginv @ S @ Ginv.T

    # implied sigma_j (residual variance corrected for the truncation term) and rho_j
    e1 = y[d] - R1[d] @ c1
    e0 = y[~d] - R0[~d] @ c0
    rs1 = Wk[d] @ c1[k:k + pk]
    rs0 = Wk[~d] @ c0[k:k + pk]
    if lay.ph == 0:
        var1 = np.average(e1 ** 2, weights=w[d]) + np.average(rs1 ** 2 * l1[d] * (l1[d] + zg[d]), weights=w[d])
        var0 = np.average(e0 ** 2, weights=w[~d]) + np.average(rs0 ** 2 * l0[~d] * (l0[~d] - zg[~d]), weights=w[~d])
    else:
        # E[omega_j | u] = g_j(u) = rs_j u + h_j2 (u^2 - 1) + h_j3 (u^3 - 3u):
        # sigma_j^2 = E[resid^2 | D = j] - Var(g_j(u) | truncation) + Var(g_j(u))
        var1 = _sigma2_hermite(e1, rs1, c1[k + pk:], -zg[d], np.inf, w[d])
        var0 = _sigma2_hermite(e0, rs0, c0[k + pk:], -np.inf, -zg[~d], w[~d])
    s1, s0 = np.sqrt(var1), np.sqrt(var0)
    out = dict(method="twostep", theta=th, V=V, layout=lay, probit=pr, n=n, k=k, m=m,
               sigma1=s1, sigma0=s0,
               rho1=np.average(rs1, weights=w[d]) / s1, rho0=np.average(rs0, weights=w[~d]) / s0,
               rhosig1=Wk @ c1[k:k + pk], rhosig0=Wk @ c0[k:k + pk], Wk=Wk,
               moments=Gm, Ginv=Ginv, hermite=hermite)
    return out
