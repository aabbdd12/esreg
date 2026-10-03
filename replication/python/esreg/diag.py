"""
Diagnostics of the selection equation (Section 4.4): strength of the excluded
instruments, collinearity of the Mills ratios with X, common support and the
share of the effects extrapolated outside it.
diagnostics(r, df, x, z, weights=None) with r an esreg Results.
"""
import numpy as np
from scipy.stats import norm, chi2
from .probit import fit_probit, probit_loglik
from .twostep import mills


def _wls_r2(y, X, w):
    sw = np.sqrt(w)
    b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    e = y - X @ b
    mu = np.average(y, weights=w)
    return 1.0 - np.sum(w * e ** 2) / np.sum(w * (y - mu) ** 2)


def diagnostics(r, df, x, z, weights=None):
    meta, fit = r.meta, r.fit
    n = len(df)
    y, d = meta["y"], meta["d"]
    dv = df[d].to_numpy() > 0
    w = df[weights].to_numpy(float) if weights else np.ones(n)
    X = np.c_[df[list(x)].to_numpy(float), np.ones(n)]
    Z = np.c_[df[list(z)].to_numpy(float), np.ones(n)]
    excl = [v for v in z if v not in x]
    # --- strength: probit with and without the excluded instruments
    pf = fit_probit(Z, dv, w)
    Xz = np.c_[df[list(x)].to_numpy(float), np.ones(n)]
    pr = fit_probit(Xz, dv, w)
    ll0 = float(w @ np.where(dv, np.log(np.average(dv, weights=w)), np.log(1 - np.average(dv, weights=w))))
    r2_full = 1 - pf["ll"] / ll0
    r2_restr = 1 - pr["ll"] / ll0
    lr = 2 * (pf["ll"] - pr["ll"])
    g, V = pf["gamma"], pf["V"]
    per = {v: (g[list(z).index(v)] ** 2) / V[list(z).index(v), list(z).index(v)] for v in excl}
    # --- score, variation of P not explained by X, VIF of the Mills ratios
    if fit["method"] == "fiml":
        lay = fit["layout"]; gam = fit["theta"][lay.g]
    else:
        lay = fit["layout"]; gam = fit["theta"][lay.g]
    zg = Z @ gam
    P = norm.cdf(zg)
    l1, l0 = mills(zg)
    varP_share = 1 - _wls_r2(P, X, w)
    vif1 = 1 / (1 - _wls_r2(l1[dv], X[dv], w[dv]))
    vif0 = 1 / (1 - _wls_r2(l0[~dv], X[~dv], w[~dv]))
    # --- support and extrapolated shares
    p1min, p1max = P[dv].min(), P[dv].max()
    p0min, p0max = P[~dv].min(), P[~dv].max()
    lo, hi = max(p1min, p0min), min(p1max, p0max)
    out_t = (P[dv] < p0min) | (P[dv] > p0max)        # treated without untreated counterpart
    out_u = (P[~dv] < p1min) | (P[~dv] > p1max)      # untreated without treated counterpart
    share_att = np.sum(w[dv] * out_t) / np.sum(w[dv])
    share_atu = np.sum(w[~dv] * out_u) / np.sum(w[~dv])
    dd = r.eff["dd"]
    atu_cs = np.average(dd[~dv][~out_u], weights=w[~dv][~out_u]) if (~out_u).any() else np.nan
    att_cs = np.average(dd[dv][~out_t], weights=w[dv][~out_t]) if (~out_t).any() else np.nan
    return dict(lr_excl=lr, df_excl=len(excl), p_excl=chi2.sf(lr, len(excl)),
                r2_full=r2_full, r2_restr=r2_restr, r2_incr=r2_full - r2_restr,
                per_instrument=per, varP_share=varP_share, vif1=vif1, vif0=vif0,
                support=(lo, hi), p_range1=(p1min, p1max), p_range0=(p0min, p0max),
                share_att_extrap=share_att, share_atu_extrap=share_atu,
                att=r.eff["est"][0], atu=r.eff["est"][1], att_cs=att_cs, atu_cs=atu_cs)
