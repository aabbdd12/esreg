"""Design K sample for the Stata validation of esrtest, kappa: kappa(x) = 1.08 + 0.4 x,
n = 5,000, seed 2026; reference values from the Python library."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg import esreg
from esreg.tests import wald_kappa_x
from mc_kappa_x import dgp
from paths import DATA
rng = np.random.default_rng(2026)
df = dgp(rng, 5000, 0.4)
df["kappa_true"] = 1.08 + 0.4*df["x"]
df.to_stata(DATA + "esr_kx.dta", write_index=False, version=118,
            data_label="Design K: kappa(x) = 1.08 + 0.4x, n = 5000, seed 2026")
r = esreg(df, "y", ["x"], "d", ["x", "z"], method="twostep", kappa=["x"]); r.summary()
t = wald_kappa_x(r); print("two-step: slope", t["slopes"], "se", t["se_slopes"], "chi2", t["chi2"], "df", t["df"], "p", t["p"])
r2 = esreg(df, "y", ["x"], "d", ["x", "z"], method="fiml", hetsigma=["x"], hetrho=["x"]); r2.summary()
t2 = wald_kappa_x(r2); print("fiml: chi2", t2["chi2"], "df", t2["df"], "p", t2["p"]); print(t2["names"], t2["slopes"], t2["se_slopes"])
# constant-kappa fits for comparison
r3 = esreg(df, "y", ["x"], "d", ["x", "z"], method="twostep"); print("constant kappa two-step:", r3.effects().loc["kappa"].to_dict())
