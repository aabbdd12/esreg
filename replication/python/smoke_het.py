"""Smoke test of the heterogeneous options (not yet validated against Stata):
hetsigma()/hetrho() in the likelihood and kappa() in the two-step route, on a
design where kappa rises with x."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg import esreg
rng = np.random.default_rng(7); n = 20000
x = rng.normal(size=n); z = rng.normal(size=n); u = rng.normal(size=n)
d = 0.2 + 0.3*x + 0.9*z + u > 0
# rho1*sigma1 = 0.6 + 0.4x (kappa(x) = 0.9 + 0.4x), rho0*sigma0 = -0.3
w1 = (0.6 + 0.4*x)*u + rng.normal(size=n); w0 = -0.3*u + rng.normal(size=n)
y = np.where(d, 1 + 1.2*x + w1, 0.5 + 0.4*x + w0)
df = pd.DataFrame(dict(y=y, d=d.astype(int), x=x, z=z))
print("true kappa(x) = 0.9 + 0.4x ; mean over sample =", (0.9 + 0.4*x).mean())
r = esreg(df, "y", ["x"], "d", ["x", "z"], method="twostep", kappa=["x"]); r.summary()
print("\nkappa(x) slope (t1 - t0 on x) =", r.theta[r.fit['layout'].t1][0] - r.theta[r.fit['layout'].t0][0])
r = esreg(df, "y", ["x"], "d", ["x", "z"], method="fiml", hetsigma=["x"], hetrho=["x"]); r.summary()
