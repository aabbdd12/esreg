"""
Monte Carlo for the test of constant selection on gains (Section 4.3):
size under H0 (kappa constant) and power under design K, kappa(x) = kappa + s*x,
two-step Wald on (t1 - t0)_x with the stacked variance; FIML Wald on the
hetsigma/hetrho coefficients.  Baseline design of paper A otherwise.
usage: python3 mc_kappa_x.py [reps]
"""
import sys, time, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg import esreg
from esreg.tests import wald_kappa_x

def dgp(rng, n, slope):
    x = rng.normal(size=n); z = rng.normal(size=n); u = rng.normal(size=n)
    d = 0.2 + 0.3*x + 0.9*z + u > 0
    # regime errors: omega_j = rs_j(x) u + e_j with Var(e_j) chosen so that sigma_j ~ (1.3, 1)
    rs1 = 0.78 + slope*x            # rho1*sigma1 at x = 0: 0.6*1.3 = 0.78
    rs0 = -0.30                     # rho0*sigma0 = -0.3
    w1 = rs1*u + rng.normal(scale=np.sqrt(1.3**2 - 0.78**2), size=n)
    w0 = rs0*u + rng.normal(scale=np.sqrt(1.0 - 0.09), size=n)
    y = np.where(d, 1.0 + 1.2*x + w1, 0.5 + 0.4*x + w0)
    return pd.DataFrame(dict(y=y, d=d.astype(int), x=x, z=z))

def run(reps, n, slope, seed, fiml=False):
    rng = np.random.default_rng(seed)
    rej = 0; est = []; t0 = time.time()
    for r in range(reps):
        df = dgp(rng, n, slope)
        if fiml:
            res = esreg(df, "y", ["x"], "d", ["x", "z"], method="fiml", hetsigma=["x"], hetrho=["x"])
        else:
            res = esreg(df, "y", ["x"], "d", ["x", "z"], method="twostep", kappa=["x"])
        t = wald_kappa_x(res)
        rej += t["p"] < 0.05
        est.append(t["slopes"][0] if not fiml else np.nan)
    return rej/reps, np.nanmean(est), np.nanstd(est), time.time() - t0

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    print(f"{'route':8s} {'n':>6s} {'slope':>6s} {'rej@5%':>7s} {'mean slope':>11s} {'sd':>7s} {'sec':>6s}")
    for n in (2000, 5000):
        for slope in (0.0, 0.1, 0.2, 0.4):
            rj, m, s, tt = run(reps, n, slope, 20260911 + n)
            print(f"{'twostep':8s} {n:6d} {slope:6.2f} {rj:7.3f} {m:11.3f} {s:7.3f} {tt:6.0f}")
    for slope in (0.0, 0.2, 0.4):
        rj, m, s, tt = run(reps // 2, 2000, slope, 777, fiml=True)
        print(f"{'fiml':8s} {2000:6d} {slope:6.2f} {rj:7.3f} {'':>11s} {'':>7s} {tt:6.0f}")
