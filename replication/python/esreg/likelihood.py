"""
Gaussian endogenous switching regression: per-observation log-likelihood and
analytic score.

Model
    Y_1 = X b1 + w1,   Y_0 = X b0 + w0,   D = 1{Z g + u > 0}
    (w1, w0, u) trivariate normal, Var(u) = 1,
    ln sigma_j = Ws a_j,   atanh rho_j = Wr c_j      (Ws, Wr = column of ones by default)

Parameter vector (this order is the one esreg.ado will use)
    theta = [b1 (k), b0 (k), g (m), a1 (ps), a0 (ps), c1 (pr), c0 (pr)]

Observation i contributes, with e_j = (y - X b_j)/s_j and
a_j = (Z g + r_j e_j)/sqrt(1 - r_j^2):
    D = 1 : log phi(e1) - ln s1 + log Phi( a1)
    D = 0 : log phi(e0) - ln s0 + log Phi(-a0)
"""
import numpy as np
from scipy.stats import norm
from scipy.special import log_ndtr


class Layout:
    """Index bookkeeping for theta."""

    def __init__(self, k, m, ps=1, pr=1):
        self.k, self.m, self.ps, self.pr = k, m, ps, pr
        i = 0
        self.b1 = slice(i, i + k); i += k
        self.b0 = slice(i, i + k); i += k
        self.g = slice(i, i + m); i += m
        self.a1 = slice(i, i + ps); i += ps
        self.a0 = slice(i, i + ps); i += ps
        self.c1 = slice(i, i + pr); i += pr
        self.c0 = slice(i, i + pr); i += pr
        self.p = i

    def unpack(self, th):
        return (th[self.b1], th[self.b0], th[self.g],
                th[self.a1], th[self.a0], th[self.c1], th[self.c0])


def _mills_pos(a):
    """phi(a)/Phi(a), stable."""
    return np.exp(norm.logpdf(a) - log_ndtr(a))


def loglik_obs(th, lay, y, X, Z, d, Ws, Wr):
    """Vector of per-observation log-likelihood contributions."""
    b1, b0, g, a1, a0, c1, c0 = lay.unpack(th)
    zg = Z @ g
    s1, s0 = np.exp(Ws @ a1), np.exp(Ws @ a0)
    r1, r0 = np.tanh(Wr @ c1), np.tanh(Wr @ c0)
    e1 = (y - X @ b1) / s1
    e0 = (y - X @ b0) / s0
    q1 = np.sqrt(1.0 - r1 ** 2)
    q0 = np.sqrt(1.0 - r0 ** 2)
    A1 = (zg + r1 * e1) / q1
    A0 = (zg + r0 * e0) / q0
    ll1 = norm.logpdf(e1) - np.log(s1) + log_ndtr(A1)
    ll0 = norm.logpdf(e0) - np.log(s0) + log_ndtr(-A0)
    return np.where(d, ll1, ll0)


def score_obs(th, lay, y, X, Z, d, Ws, Wr):
    """n x p matrix of per-observation scores (analytic)."""
    b1, b0, g, a1, a0, c1, c0 = lay.unpack(th)
    n = len(y)
    zg = Z @ g
    s1, s0 = np.exp(Ws @ a1), np.exp(Ws @ a0)
    r1, r0 = np.tanh(Wr @ c1), np.tanh(Wr @ c0)
    e1 = (y - X @ b1) / s1
    e0 = (y - X @ b0) / s0
    q1 = np.sqrt(1.0 - r1 ** 2)
    q0 = np.sqrt(1.0 - r0 ** 2)
    A1 = (zg + r1 * e1) / q1
    A0 = (zg + r0 * e0) / q0
    m1 = _mills_pos(A1)       # phi(A1)/Phi(A1)
    m0 = _mills_pos(-A0)      # phi(A0)/Phi(-A0)
    d1 = d.astype(float)
    d0 = 1.0 - d1

    S = np.zeros((n, lay.p))
    # regime coefficients
    S[:, lay.b1] = X * (d1 * (e1 - m1 * r1 / q1) / s1)[:, None]
    S[:, lay.b0] = X * (d0 * (e0 + m0 * r0 / q0) / s0)[:, None]
    # selection coefficients
    S[:, lay.g] = Z * (d1 * m1 / q1 - d0 * m0 / q0)[:, None]
    # ln sigma
    S[:, lay.a1] = Ws * (d1 * (e1 ** 2 - 1.0 - m1 * r1 * e1 / q1))[:, None]
    S[:, lay.a0] = Ws * (d0 * (e0 ** 2 - 1.0 + m0 * r0 * e0 / q0))[:, None]
    # atanh rho
    S[:, lay.c1] = Wr * (d1 * m1 * (e1 + r1 * zg) / q1)[:, None]
    S[:, lay.c0] = Wr * (-d0 * m0 * (e0 + r0 * zg) / q0)[:, None]
    return S


def hessian_from_score(th, lay, y, X, Z, d, Ws, Wr, w, h=1e-5):
    """Hessian of the weighted log-likelihood by central differences of the analytic score."""
    p = lay.p
    H = np.zeros((p, p))
    for i in range(p):
        e = np.zeros(p); e[i] = h
        gp = w @ score_obs(th + e, lay, y, X, Z, d, Ws, Wr)
        gm = w @ score_obs(th - e, lay, y, X, Z, d, Ws, Wr)
        H[:, i] = (gp - gm) / (2 * h)
    return 0.5 * (H + H.T)
