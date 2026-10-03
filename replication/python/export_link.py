"""Sample for the Stata validation of esrtest, spec (Section 4.1): the selection
index has a quadratic term in z that the working probit on (x, z) omits.
D = 1{0.2 + 0.3x + 0.9z + 0.5(z^2 - 1) + u > 0}, outcomes and (w1, w0, u) as esr_sim
(kappa = 1.08).  esr_link.dta carries z2 = z^2.  Reference values from Python."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg.tests import link_test, gamma_contrast
from esreg import esreg
from paths import DATA as OUT
rng = np.random.default_rng(2026)
n = 5000
x = rng.standard_normal(n); z = rng.standard_normal(n); u = rng.standard_normal(n)
d = (0.2 + 0.3 * x + 0.9 * z + 0.5 * (z ** 2 - 1) + u > 0)
S = np.array([[1.69, 0.15, 0.6 * 1.3], [0.15, 1.0, -0.3], [0.6 * 1.3, -0.3, 1.0]])
mu = np.array([0.78, -0.3])[:, None] * u[None, :]
Sc = S[:2, :2] - np.outer(S[:2, 2], S[2, :2])
e = mu + np.linalg.cholesky(Sc) @ rng.standard_normal((2, n))
y1 = 1.0 + 1.2 * x + e[0]; y0 = 0.5 + 0.4 * x + e[1]
y = np.where(d, y1, y0)
df = pd.DataFrame(dict(y=y, d=d.astype(int), x=x, z=z, z2=z ** 2))
df.to_stata(OUT + "esr_link.dta", write_index=False, version=118,
            data_label="Selection index with an omitted quadratic in z: D = 1{.2+.3x+.9z+.5(z^2-1)+u>0}; n = 5000, seed 2026")

def show(t, lab):
    print(f"{lab}: chi2({t['df']}) = {t['chi2']:.4f}  p = {t['p']:.2e}  coefs {np.round(t['coefs'], 5)} se {np.round(t['se'], 5)}"
          + (f"  | LR = {t['lr']:.4f} p = {t['p_lr']:.2e}" if 'lr' in t else ""))

for name in ("esr_link", "esr_sim"):
    dd = df if name == "esr_link" else pd.read_stata(OUT + name + ".dta")
    print(f"---- {name}")
    show(link_test(dd, "d", ["x", "z"]), "link test order 2, probit on (x z)")
    show(link_test(dd, "d", ["x", "z"], order=3), "link test order 3, probit on (x z)")
    if name == "esr_link":
        show(link_test(dd, "d", ["x", "z"], add=["z2"]), "RESET-type, add z2")
        show(link_test(dd, "d", ["x", "z", "z2"]), "link test order 2, probit on (x z z2)")
    for zz in (["x", "z"], ["x", "z", "z2"]) if name == "esr_link" else (["x", "z"],):
        g = gamma_contrast(dd, "y", ["x"], "d", zz)
        print(f"gamma contrast, select({' '.join(zz)}): chi2({g['df']}) = {g['chi2']:.4f} p = {g['p']:.3f}")
        print("   probit ", np.round(g["g_probit"], 5)); print("   fiml   ", np.round(g["g_fiml"], 5))
        print("   se diff", np.round(g["se_diff"], 5))
        h = g["hausman"]
        print(f"   full Hausman chi2({h['df']}) = {h['chi2']:.2f}; kappa 2S {h['kappa_2s']:.4f} FIML {h['kappa_fiml']:.4f}")
    if name == "esr_link":
        for zz in (["x", "z"], ["x", "z", "z2"]):
            for meth in ("fiml", "twostep"):
                r = esreg(dd, "y", ["x"], "d", zz, method=meth)
                ef = r.effects()
                print(f"   esreg {meth:7s} select({' '.join(zz):6s}): ATT {ef.loc['ATT','est']:.4f}  kappa {ef.loc['kappa','est']:.4f} ({ef.loc['kappa','se']:.4f})")
