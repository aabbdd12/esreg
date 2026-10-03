"""
Monte Carlo for Paper A, Section 5: ONE design with switches.

Baseline (N): x, z, u ~ N(0,1); D = 1{0.2 + 0.3x + 0.9z + u > 0};
  Y1 = 1.0 + 1.2x + w1, Y0 = 0.5 + 0.4x + w0, (w1, w0, u) jointly normal with
  s1 = 1.3, s0 = 1, r1 = .6, r0 = -.3, cov10 = .15  ->  kappa = 1.08.
Two periods with permanent unobservables only: Y_pre = 0.5 + 0.4x + w0 (nobody treated),
Y_post = D*Y1 + (1-D)*Y0 + delta_t, delta_t = 0.5 (SIG_V = 0: no transitory shock, so that the
post cross-section is exactly esr_sim / esr_nonlin / esr_skew).

Switches (one assumption removed at a time):
  H   : kappa = 0 and homogeneous effect (w1 = w0, r = .5, effect 2)
  C1  : E[w1 - w0 | u] nonlinear (Hermite cubic), normal errors
  C2  : linear gain, skewed regime errors (standardized lognormal)
  Z   : exclusion restriction fails: z enters both outcomes with coefficient -.3
  T   : differential trend: delta_t = 0.5 + 0.5 u  (post period only)
  M   : non-additive time change: Y0_post = 2 (0.5 + 0.4x + w0) + 0.5
Estimators on the post cross-section: naive, IPW on x (ATT), ESR FIML, ESR two-step,
single-equation ETR (common error, ML); on the two periods: DiD, CiC (Athey-Imbens).
MTE curves (ESR line, two-step line, PWR on the score) for N, C1, C2.
"""
import sys, json, time, numpy as np
from multiprocessing import Pool
from scipy.stats import norm
from scipy.optimize import minimize
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from esr_fiml_reference import _nll, _hessian

B1 = np.array([1.0, 1.2]); B0 = np.array([0.5, 0.4]); G = np.array([0.2, 0.3, 0.9])
KAPPA = 1.08; C3 = 0.5; SKEW = 0.5; DELTA_T = 0.5; SIG_V = 0.0
THETA_T = 0.5; A_M = 2.0; GZ_Y = -0.3
SWITCHES = ["N", "H", "C1", "C2", "Z", "T", "M"]

def h3(u): return u**3 - 3*u

def skew(rng, n, s):
    e = np.exp(SKEW*rng.standard_normal(n)); e = (e - e.mean())/e.std(); return s*e

def dgp(rng, n, sw):
    x = rng.standard_normal(n); z = rng.standard_normal(n); u = rng.standard_normal(n)
    d = (G[0] + G[1]*x + G[2]*z + u > 0)
    X = np.c_[np.ones(n), x]; Z = np.c_[np.ones(n), x, z]
    b1, b0 = B1.copy(), B0.copy()
    if sw == "H":
        b1 = np.array([2.5, 0.4])
        w0 = 0.5*u + np.sqrt(1-0.25)*rng.standard_normal(n); w1 = w0.copy()
        gain_u = lambda uu: 0*uu
    elif sw == "C1":
        w0 = -0.3*u + 0.95*rng.standard_normal(n)
        w1 = w0 + KAPPA*(u + C3*h3(u)) + 0.8*rng.standard_normal(n)
        gain_u = lambda uu: KAPPA*(uu + C3*h3(uu))
    elif sw == "C2":
        w0 = -0.3*u + skew(rng, n, 0.95)
        w1 = w0 + KAPPA*u + skew(rng, n, 0.8)
        gain_u = lambda uu: KAPPA*uu
    else:   # N, Z, T, M : joint normal as esr_sim
        S = np.array([[1.69, 0.15, 0.6*1.3], [0.15, 1.0, -0.3], [0.6*1.3, -0.3, 1.0]])
        mu = np.array([0.78, -0.3])[:, None]*u[None, :]
        Sc = S[:2, :2] - np.outer(S[:2, 2], S[2, :2])
        e = mu + np.linalg.cholesky(Sc) @ rng.standard_normal((2, n))
        w1, w0 = e[0], e[1]
        gain_u = lambda uu: KAPPA*uu
    y1 = X @ b1 + w1; y0 = X @ b0 + w0
    if sw == "Z":
        y1 = y1 + GZ_Y*z; y0 = y0 + GZ_Y*z
    ypre = X @ B0 + w0 + SIG_V*rng.standard_normal(n)        # untreated, period 1
    if sw == "Z": ypre = ypre + GZ_Y*z
    trend = DELTA_T + (THETA_T*u if sw == "T" else 0.0)
    if sw == "M":
        y0 = A_M*y0                                          # non-additive time change
    vpost = SIG_V*rng.standard_normal(n)                     # transitory error, period 2
    y1p = y1 + trend + vpost; y0p = y0 + trend + vpost
    y = np.where(d, y1p, y0p)
    return dict(x=x, z=z, u=u, d=d, X=X, Z=Z, y=y, ypre=ypre, delta=y1p-y0p, gain_u=gain_u)

# ---------------------------------------------------------------- estimators
def probit(Z, d):
    return minimize(lambda g_: -np.where(d, norm.logcdf(Z@g_), norm.logcdf(-Z@g_)).sum(),
                    np.zeros(Z.shape[1]), method="BFGS").x

def unpack(th, k, m):
    b1, b0, g = th[:k], th[k:2*k], th[2*k:2*k+m]
    s1, s0 = np.exp(th[2*k+m]), np.exp(th[2*k+m+1]); r1, r0 = np.tanh(th[2*k+m+2]), np.tanh(th[2*k+m+3])
    return b1, b0, g, s1, s0, r1, r0

def effects_from(th, X, Z, d):
    k, m = X.shape[1], Z.shape[1]
    b1, b0, g, s1, s0, r1, r0 = unpack(th, k, m)
    zg = Z @ g; l1 = norm.pdf(zg)/norm.cdf(zg); l0 = norm.pdf(zg)/(1-norm.cdf(zg))
    kap = r1*s1 - r0*s0; mg = X @ (b1-b0)
    dd = np.where(d, mg + kap*l1, mg - kap*l0)
    return np.array([dd[d].mean(), dd[~d].mean(), dd.mean(), kap]), dd

def fit_esr(y, X, Z, d, se=False):
    k, m = X.shape[1], Z.shape[1]
    g0 = probit(Z, d)
    b1 = np.linalg.lstsq(X[d], y[d], rcond=None)[0]; b0 = np.linalg.lstsq(X[~d], y[~d], rcond=None)[0]
    st = np.r_[b1, b0, g0, np.log(np.std(y[d]-X[d]@b1)), np.log(np.std(y[~d]-X[~d]@b0)), 0.0, 0.0]
    th = minimize(_nll, st, args=(y, X, Z, d), method="L-BFGS-B", options={"maxiter": 3000}).x
    eff, dd = effects_from(th, X, Z, d)
    b1, b0, g, s1, s0, r1, r0 = unpack(th, k, m)
    out = dict(att=eff[0], atu=eff[1], ate=eff[2], kappa=eff[3], m0=(b1-b0)[0],
               r1=r1, r0=r0, s1=s1, s0=s0, gz=g[2], p=norm.cdf(Z@g))
    if se:
        H = _hessian(_nll, th, (y, X, Z, d))
        V = np.linalg.inv(H)
        # gradient of the effects w.r.t. theta (central differences)
        p = len(th); Jm = np.zeros((4, p)); hh = 1e-5
        for i in range(p):
            e = np.zeros(p); e[i] = hh
            Jm[:, i] = (effects_from(th+e, X, Z, d)[0] - effects_from(th-e, X, Z, d)[0])/(2*hh)
        Vp = Jm @ V @ Jm.T                                  # parameter uncertainty
        n1, n0, n = d.sum(), (~d).sum(), len(d)
        samp = np.array([dd[d].var()/n1, dd[~d].var()/n0, dd.var()/n, 0.0])   # sampling component
        out["se_delta"] = np.sqrt(np.diag(Vp) + samp)
        out["se_disp"] = np.sqrt(samp)                       # dispersion of predicted effects only
    return out

def twostep(y, X, Z, d):
    g = probit(Z, d); zg = Z@g; l1 = norm.pdf(zg)/norm.cdf(zg); l0 = norm.pdf(zg)/(1-norm.cdf(zg))
    b1 = np.linalg.lstsq(np.c_[X[d], l1[d]], y[d], rcond=None)[0]
    b0 = np.linalg.lstsq(np.c_[X[~d], -l0[~d]], y[~d], rcond=None)[0]
    kap = b1[-1] - b0[-1]; mg = X @ (b1[:-1]-b0[:-1])
    dd = np.where(d, mg + kap*l1, mg - kap*l0)
    return dict(att=dd[d].mean(), atu=dd[~d].mean(), ate=dd.mean(), kappa=kap, m0=(b1[:-1]-b0[:-1])[0], p=norm.cdf(zg))

def etr(y, X, Z, d):
    """single-equation endogenous treatment: y = X b + delta d + e, (e,u) bivariate normal, ML"""
    k, m = X.shape[1], Z.shape[1]
    def expand(t):     # t = [b(k), delta, g(m), ln s, atanh r]
        b, dl, g, ls, ar = t[:k], t[k], t[k+1:k+1+m], t[k+1+m], t[k+2+m]
        b1 = b.copy(); b1[0] += dl
        return np.r_[b1, b, g, ls, ls, ar, ar]
    g0 = probit(Z, d)
    bo = np.linalg.lstsq(np.c_[X, d], y, rcond=None)[0]
    st = np.r_[bo[:k], bo[k], g0, np.log(np.std(y - np.c_[X, d]@bo)), 0.0]
    t = minimize(lambda t_: _nll(expand(t_), y, X, Z, d), st, method="L-BFGS-B", options={"maxiter": 3000}).x
    return dict(att=t[k], rho=np.tanh(t[k+2+m]))

def ipw_att(y, x, d):
    from scipy.special import expit
    Xl = np.c_[np.ones(len(x)), x]
    f = lambda b: -(d*np.log(expit(Xl@b)+1e-12) + (~d)*np.log(1-expit(Xl@b)+1e-12)).sum()
    b = minimize(f, np.zeros(2), method="BFGS").x; ps = expit(Xl@b)
    w0 = ps/(1-ps); w0 = w0[~d]/w0[~d].sum()
    return y[d].mean() - (w0*y[~d]).sum()

def did(y, ypre, d):
    return (y[d].mean() - ypre[d].mean()) - (y[~d].mean() - ypre[~d].mean())

def cic(y, ypre, d):
    """Athey-Imbens changes-in-changes ATT: E[Y11] - E[F01^{-1}(F00(Y10))]"""
    y00 = np.sort(ypre[~d]); y01 = np.sort(y[~d]); y10 = ypre[d]
    F = np.searchsorted(y00, y10, side="right")/len(y00)
    F = np.clip(F, 0.5/len(y00), 1-0.5/len(y00))
    cf = np.quantile(y01, F)
    return y[d].mean() - cf.mean()

def pwr_mte(y, x, p, taus, c0=2.5):
    n = len(y); r = (np.argsort(np.argsort(p)) + 0.5)/n
    h = c0*np.std(r)*n**(-0.2)
    xc = x - x.mean(); W = np.c_[np.ones(n), xc, xc*p, p, p**2]
    out = []
    for t in taus:
        w = np.exp(-(r-t)**2/(4*h*h)); w /= w.sum()
        Ww = W*w[:, None]; beta = np.linalg.solve(W.T @ Ww, Ww.T @ y)
        pt = np.quantile(p, t); neff = 1/np.sum(w**2)
        lo, hi = np.quantile(p, max(t-h*np.sqrt(2), 0)), np.quantile(p, min(t+h*np.sqrt(2), 1))
        out.append((beta[3] + 2*beta[4]*pt, neff, hi-lo))
    return np.array(out)

# ---------------------------------------------------------------- one switch
def run_switch(args):
    sw, reps, n, seed = args
    rng = np.random.default_rng(seed)
    taus = np.arange(0.10, 0.901, 0.05)
    rows = []; mte = []
    t0 = time.time()
    for r in range(reps):
        D = dgp(rng, n, sw); y, x, d, X, Z = D["y"], D["x"], D["d"], D["X"], D["Z"]
        es = fit_esr(y, X, Z, d); ts = twostep(y, X, Z, d); et = etr(y, X, Z, d)
        row = dict(att_true=D["delta"][d].mean(), atu_true=D["delta"][~d].mean(), ate_true=D["delta"].mean(),
                   naive=y[d].mean()-y[~d].mean(), ipw=ipw_att(y, x, d),
                   esr_att=es["att"], esr_atu=es["atu"], esr_ate=es["ate"], esr_kappa=es["kappa"],
                   esr_r1=es["r1"], esr_r0=es["r0"], esr_gz=es["gz"],
                   ts_att=ts["att"], ts_atu=ts["atu"], ts_ate=ts["ate"], ts_kappa=ts["kappa"],
                   etr_att=et["att"], etr_rho=et["rho"],
                   did=did(y, D["ypre"], d), cic=cic(y, D["ypre"], d))
        rows.append(row)
        if sw in ("N", "C1", "C2"):
            pt = np.quantile(ts["p"], taus)
            truth = (B1-B0)[0] + D["gain_u"](norm.ppf(1-pt))
            line = es["m0"] + es["kappa"]*norm.ppf(1-pt); line2 = ts["m0"] + ts["kappa"]*norm.ppf(1-pt)
            mte.append(np.c_[pt, truth, line, line2, pwr_mte(y, x, ts["p"], taus)])
        if r % 10 == 0: print(sw, r, f"{time.time()-t0:.0f}s", flush=True)
    big = dgp(np.random.default_rng(seed+999), 2_000_000, sw)
    truth = dict(att=float(big["delta"][big["d"]].mean()), atu=float(big["delta"][~big["d"]].mean()),
                 ate=float(big["delta"].mean()))
    return sw, dict(rows=rows, mte=np.array(mte).tolist() if mte else None, truth=truth)

def run_coverage(args):
    """switch H at n, reps: coverage of 95% intervals, delta-method vs dispersion-only"""
    n, reps, seed = args
    rng = np.random.default_rng(seed); out = []
    for r in range(reps):
        D = dgp(rng, n, "H"); es = fit_esr(D["y"], D["X"], D["Z"], D["d"], se=True)
        out.append(dict(att=es["att"], atu=es["atu"], ate=es["ate"], kappa=es["kappa"],
                        se_delta=es["se_delta"].tolist(), se_disp=es["se_disp"].tolist()))
        if r % 20 == 0: print("H", n, r, flush=True)
    return n, out

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    n = 10000
    jobs = [(sw, reps, n, 20260908 + i) for i, sw in enumerate(SWITCHES)]
    with Pool(2) as pool:
        res = dict(pool.map(run_switch, jobs))
        cov = dict(pool.map(run_coverage, [(1000, 2*reps, 777), (5000, 2*reps, 778)]))
    res["coverage"] = {str(k): v for k, v in cov.items()}
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "mc_common.json"), "w"))
    print("saved")
