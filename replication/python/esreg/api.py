"""
User-level entry point.

    from esreg import esreg
    r = esreg(df, y="y", x=["x"], d="d", z=["x", "z"], method="fiml")   # or "twostep"
    r.summary()

Design matrices get a constant appended last (Stata convention).
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
from .fiml import fit_fiml
from .twostep import fit_twostep
from .effects import treatment_effects


def _design(df, cols):
    M = df[cols].to_numpy(float) if cols else np.empty((len(df), 0))
    return np.c_[M, np.ones(len(df))]


class Results:
    def __init__(self, fit, eff, names, meta):
        self.fit, self.eff, self.names, self.meta = fit, eff, names, meta

    @property
    def theta(self):
        return self.fit["theta"]

    @property
    def se(self):
        return np.sqrt(np.diag(self.fit["V"]))

    def table(self):
        return pd.DataFrame({"coef": self.theta, "se": self.se, "z": self.theta / self.se},
                            index=self.names)

    def effects(self):
        e = self.eff
        return pd.DataFrame({"est": e["est"], "se": e["se"], "z": e["z"],
                             "p": 2 * norm.sf(np.abs(e["z"]))}, index=e["names"])

    def summary(self, digits=6):
        f, m = self.fit, self.meta
        hdr = "Endogenous switching regression" + (" -- full information ML" if f["method"] == "fiml"
                                                   else " -- two-step (probit + OLS with Mills ratios)")
        print(hdr)
        print(f"Outcome: {m['y']}   Treatment: {m['d']}   n = {f['n']}   "
              f"treated = {int(m['d_sum'])}")
        if f["method"] == "fiml":
            print(f"Log likelihood = {f['loglik']:.4f}   vce = {f['vce']}   "
                  f"LR test rho1 = rho0 = 0: chi2({f['lr_df']}) = {f['lr_indep']:.2f}")
        print("-" * 72)
        print(f"{'':>22s}{'Coef.':>12s}{'Std. err.':>12s}{'z':>10s}{'P>|z|':>8s}")
        print("-" * 72)
        se = self.se
        cur = None
        for nm, b, s in zip(self.names, self.theta, se):
            eq = nm.split(":")[0]
            if eq != cur:
                print(eq); cur = eq
            zz = b / s
            print(f"{nm.split(':')[1]:>22s}{b:12.{digits}g}{s:12.{digits}g}{zz:10.2f}{2*norm.sf(abs(zz)):8.3f}")
        print("-" * 72)
        s1, s0, r1, r0 = f["sigma1"], f["sigma0"], f["rho1"], f["rho0"]
        s1 = float(np.mean(s1)); s0 = float(np.mean(s0)); r1 = float(np.mean(r1)); r0 = float(np.mean(r0))
        print(f"sigma_1 = {s1:.6g}   sigma_0 = {s0:.6g}   rho_1 = {r1:.6g}   rho_0 = {r0:.6g}"
              + ("   (means; heterogeneous)" if f["method"] == "fiml" and (f["Ws"].shape[1] > 1 or f["Wr"].shape[1] > 1) else ""))
        print()
        print("Treatment effects (std. err.: influence function of the whole procedure)")
        print("-" * 72)
        print(f"{'':>8s}{'Estimate':>13s}{'Std. err.':>12s}{'z':>10s}{'P>|z|':>8s}{'[95% conf. int.]':>22s}")
        e = self.eff
        for nm, b, s, zz in zip(e["names"], e["est"], e["se"], e["z"]):
            print(f"{nm:>8s}{b:13.6g}{s:12.6g}{zz:10.2f}{2*norm.sf(abs(zz)):8.3f}"
                  f"{b-1.96*s:11.5g}{b+1.96*s:11.5g}")
        print("-" * 72)
        cs = e["common_support"]
        print(f"kappa = rho1*sigma1 - rho0*sigma0;  ATT - ATU = kappa*(mean lambda1 + mean lambda0)")
        print(f"Variance of ATT: parameters {e['var_param'][0]:.3g} + sampling {e['var_samp'][0]:.3g}"
              f" + covariance {e['cov_ps'][0]:.3g}")
        print(f"Common support of P(Z): [{cs[0]:.4f}, {cs[1]:.4f}]")


def esreg(df, y, x, d, z, method="fiml", weights=None, hetsigma=None, hetrho=None,
           kappa=None, vce="oim", verbose=False, hermite=0):
    """
    df       : pandas DataFrame
    y, d     : outcome and binary treatment column names
    x, z     : lists of regressors of the outcome equations and of the selection equation
    method   : "fiml" or "twostep"
    weights  : column name of frequency/sampling weights (optional)
    hetsigma, hetrho : lists of columns entering ln sigma_j and atanh rho_j (fiml only)
    kappa    : list of columns W entering kappa(x) = W(t1 - t0) in the two-step route
    hermite  : two-step, the number of Hermite controls (u^2 - 1, u^3 - 3u) in each
               regime, 0, 1 or 2, or a pair (treated, untreated)
    """
    df = df.dropna(subset=[y, d] + list(x) + list(z)).reset_index(drop=True)
    yv = df[y].to_numpy(float)
    dv = df[d].to_numpy() > 0
    X, Z = _design(df, list(x)), _design(df, list(z))
    w = df[weights].to_numpy(float) if weights else np.ones(len(df))
    xn = list(x) + ["_cons"]; zn = list(z) + ["_cons"]
    if method == "fiml":
        Ws = _design(df, list(hetsigma)) if hetsigma else None
        Wr = _design(df, list(hetrho)) if hetrho else None
        fit = fit_fiml(yv, X, Z, dv, w, Ws, Wr, vce=vce, verbose=verbose)
        sn = (list(hetsigma) if hetsigma else []) + ["_cons"]
        rn = (list(hetrho) if hetrho else []) + ["_cons"]
        names = ([f"{y}_1:{v}" for v in xn] + [f"{y}_0:{v}" for v in xn] + [f"{d}:{v}" for v in zn]
                 + [f"lnsigma_1:{v}" for v in sn] + [f"lnsigma_0:{v}" for v in sn]
                 + [f"atanhrho_1:{v}" for v in rn] + [f"atanhrho_0:{v}" for v in rn])
    elif method == "twostep":
        Wk = _design(df, list(kappa)) if kappa else None
        fit = fit_twostep(yv, X, Z, dv, w, Wk, hermite=hermite)
        kn = (list(kappa) if kappa else []) + ["_cons"]
        lay = fit["layout"]
        hn1 = ["h2", "h3"][:lay.ph1]
        hn0 = ["h2", "h3"][:lay.ph0]
        names = ([f"{d}:{v}" for v in zn] + [f"{y}_1:{v}" for v in xn] + [f"{y}_1:lambda1*{v}" for v in kn]
                 + [f"{y}_1:{v}" for v in hn1]
                 + [f"{y}_0:{v}" for v in xn] + [f"{y}_0:-lambda0*{v}" for v in kn] + [f"{y}_0:{v}" for v in hn0])
    else:
        raise ValueError("method must be 'fiml' or 'twostep'")
    eff = treatment_effects(fit, X, Z, dv, w)
    meta = dict(y=y, d=d, x=list(x), z=list(z), d_sum=dv.sum(), method=method)
    return Results(fit, eff, names, meta)
