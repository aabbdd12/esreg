"""Full-information maximum likelihood for the Gaussian ESR."""
import numpy as np
from scipy.optimize import minimize
from .likelihood import Layout, loglik_obs, score_obs, hessian_from_score
from .probit import fit_probit


def starting_values(y, X, Z, d, w, lay):
    pr = fit_probit(Z, d, w)
    sw = np.sqrt(w)
    b1 = np.linalg.lstsq(X[d] * sw[d, None], y[d] * sw[d], rcond=None)[0]
    b0 = np.linalg.lstsq(X[~d] * sw[~d, None], y[~d] * sw[~d], rcond=None)[0]
    th = np.zeros(lay.p)
    th[lay.b1], th[lay.b0], th[lay.g] = b1, b0, pr["gamma"]
    th[lay.a1.start] = np.log(np.sqrt(np.average((y[d] - X[d] @ b1) ** 2, weights=w[d])))
    th[lay.a0.start] = np.log(np.sqrt(np.average((y[~d] - X[~d] @ b0) ** 2, weights=w[~d])))
    return th, pr


def fit_fiml(y, X, Z, d, w=None, Ws=None, Wr=None, vce="oim", tol=1e-9, maxit=200, verbose=False):
    """
    Returns dict with theta, V, layout, loglik, iterations, and the derived
    sigma_j, rho_j (per observation when heterogeneous).
    vce: "oim" (inverse observed information), "opg", or "robust" (sandwich).
    """
    n, k = X.shape
    m = Z.shape[1]
    w = np.ones(n) if w is None else np.asarray(w, float)
    Ws = np.ones((n, 1)) if Ws is None else Ws
    Wr = np.ones((n, 1)) if Wr is None else Wr
    lay = Layout(k, m, Ws.shape[1], Wr.shape[1])
    args = (lay, y, X, Z, d, Ws, Wr)

    f = lambda th: -float(w @ loglik_obs(th, *args))
    jac = lambda th: -(w @ score_obs(th, *args))

    th0, pr = starting_values(y, X, Z, d, w, lay)
    res = minimize(f, th0, jac=jac, method="BFGS", options={"gtol": 1e-6, "maxiter": 5000})
    th = res.x
    # Newton polish on the analytic score with a finite-difference Hessian
    it = 0
    for it in range(1, maxit + 1):
        g = w @ score_obs(th, *args)
        H = hessian_from_score(th, *args, w)
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
        th_new = th - step
        if f(th_new) > f(th) + 1e-8:          # safeguard: halve the step
            for _ in range(20):
                step *= 0.5
                th_new = th - step
                if f(th_new) <= f(th):
                    break
        th = th_new
        if verbose:
            print(f"iter {it}: ll = {-f(th):.6f}  max|step| = {np.max(np.abs(step)):.2e}")
        if np.max(np.abs(step)) < tol:
            break

    H = hessian_from_score(th, *args, w)
    S = score_obs(th, *args) * w[:, None]
    Vinv = -H
    V_oim = np.linalg.inv(Vinv)
    if vce == "oim":
        V = V_oim
    elif vce == "opg":
        V = np.linalg.inv(S.T @ S)
    else:  # robust sandwich
        V = V_oim @ (S.T @ S) @ V_oim

    b1, b0, g, a1, a0, c1, c0 = lay.unpack(th)
    out = dict(method="fiml", theta=th, V=V, layout=lay, loglik=-f(th), iterations=it,
               scores=S, V_oim=V_oim,
               converged=np.max(np.abs(w @ score_obs(th, *args))) < 1e-3,
               probit=pr, n=n, k=k, m=m, vce=vce,
               sigma1=np.exp(Ws @ a1), sigma0=np.exp(Ws @ a0),
               rho1=np.tanh(Wr @ c1), rho0=np.tanh(Wr @ c0),
               Ws=Ws, Wr=Wr)
    # log-likelihood of the independent-equations model (rho1 = rho0 = 0) for the LR test
    th_ind = th.copy(); th_ind[lay.c1] = 0; th_ind[lay.c0] = 0
    lay_c = lay
    fr = lambda tt: -float(w @ loglik_obs(_restore(tt, th_ind, lay_c), lay_c, y, X, Z, d, Ws, Wr))
    free = np.r_[np.arange(lay.c1.start)]
    r = minimize(fr, th_ind[free], method="BFGS", options={"gtol": 1e-6})
    out["loglik_indep"] = -r.fun
    out["lr_indep"] = 2 * (out["loglik"] - out["loglik_indep"])
    out["lr_df"] = lay.pr * 2
    return out


def _restore(free_vals, template, lay):
    th = template.copy()
    th[:lay.c1.start] = free_vals
    return th
