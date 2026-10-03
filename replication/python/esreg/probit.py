"""Weighted probit with analytic score and Hessian (first step of the two-step route)."""
import numpy as np
from scipy.stats import norm
from scipy.special import log_ndtr
from scipy.optimize import minimize


def _mills(a):
    return np.exp(norm.logpdf(a) - log_ndtr(a))


def probit_loglik(g, Z, d, w):
    zg = Z @ g
    return float(w @ np.where(d, log_ndtr(zg), log_ndtr(-zg)))


def probit_score_obs(g, Z, d, w):
    zg = Z @ g
    gr = np.where(d, _mills(zg), -_mills(-zg))          # generalized residual
    return Z * (w * gr)[:, None]


def probit_hessian(g, Z, d, w):
    zg = Z @ g
    m1, m0 = _mills(zg), _mills(-zg)
    lam = np.where(d, m1 * (m1 + zg), m0 * (m0 - zg))    # minus second derivative
    return -(Z * (w * lam)[:, None]).T @ Z


def fit_probit(Z, d, w=None, tol=1e-10, maxit=100):
    n, m = Z.shape
    w = np.ones(n) if w is None else w
    g = minimize(lambda g_: -probit_loglik(g_, Z, d, w),
                 np.zeros(m), jac=lambda g_: -probit_score_obs(g_, Z, d, w).sum(0),
                 method="BFGS").x
    for _ in range(maxit):                                # Newton polish
        s = probit_score_obs(g, Z, d, w).sum(0)
        H = probit_hessian(g, Z, d, w)
        step = np.linalg.solve(H, s)
        g = g - step
        if np.max(np.abs(step)) < tol:
            break
    H = probit_hessian(g, Z, d, w)
    return dict(gamma=g, V=np.linalg.inv(-H), ll=probit_loglik(g, Z, d, w),
                score_obs=probit_score_obs(g, Z, d, w), H=H)
