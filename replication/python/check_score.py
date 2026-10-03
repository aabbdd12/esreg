"""Check the analytic score against finite differences of the log-likelihood."""
import numpy as np, pandas as pd, sys
sys.path.insert(0, ".")
from esreg.likelihood import Layout, loglik_obs, score_obs
rng = np.random.default_rng(1)
n = 400
x = rng.normal(size=n); z = rng.normal(size=n); u = rng.normal(size=n)
d = (0.2 + 0.3*x + 0.9*z + u > 0)
y = np.where(d, 1 + 1.2*x + 0.6*u + rng.normal(size=n), 0.5 + 0.4*x - 0.3*u + rng.normal(size=n))
X = np.c_[x, np.ones(n)]; Z = np.c_[x, z, np.ones(n)]
Ws = np.c_[x, np.ones(n)]; Wr = np.c_[z, np.ones(n)]      # heterogeneous to exercise every block
lay = Layout(2, 3, 2, 2)
th = rng.normal(scale=0.3, size=lay.p)
w = rng.uniform(0.5, 1.5, size=n)
f = lambda t: w @ loglik_obs(t, lay, y, X, Z, d, Ws, Wr)
g_an = w @ score_obs(th, lay, y, X, Z, d, Ws, Wr)
g_num = np.array([(f(th + h*e) - f(th - h*e))/(2*h) for h in [1e-6] for e in np.eye(lay.p)])
print("max |analytic - numerical| =", np.max(np.abs(g_an - g_num)))
print("max relative              =", np.max(np.abs(g_an - g_num)/np.maximum(1, np.abs(g_num))))
