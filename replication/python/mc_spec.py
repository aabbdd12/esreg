"""Monte Carlo for Section 4.1: the link test of the probit index (order 2) and the
gamma contrast (probit vs FIML, the gamma block of the Hausman statistic) on the
designs N (index linear in Z, joint normal: size), Q (omitted quadratic .5(z^2-1) in
the index: power of the link test, A1 holds so the contrast is a size check on a
misspecified index), C2 (skewed regime errors, A3 fails: power of the contrast,
size of the link test).  n = 5000."""
import sys, time, numpy as np, pandas as pd
sys.path.insert(0, ".")
import mc_common as M
from esreg.tests import link_test, gamma_contrast

def dgp(rng, n, des):
    if des != "Q":
        D = M.dgp(rng, n, des)
        return pd.DataFrame(dict(y=D["y"], d=D["d"].astype(int), x=D["x"], z=D["z"]))
    x = rng.standard_normal(n); z = rng.standard_normal(n); u = rng.standard_normal(n)
    d = (0.2 + 0.3 * x + 0.9 * z + 0.5 * (z ** 2 - 1) + u > 0)
    S = np.array([[1.69, 0.15, 0.6 * 1.3], [0.15, 1.0, -0.3], [0.6 * 1.3, -0.3, 1.0]])
    mu = np.array([0.78, -0.3])[:, None] * u[None, :]
    Sc = S[:2, :2] - np.outer(S[:2, 2], S[2, :2])
    e = mu + np.linalg.cholesky(Sc) @ rng.standard_normal((2, n))
    y = np.where(d, 1.0 + 1.2 * x + e[0], 0.5 + 0.4 * x + e[1])
    return pd.DataFrame(dict(y=y, d=d.astype(int), x=x, z=z))

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    n = 5000
    print(f"{'design':7s} {'link':>6s} {'link3':>6s} {'gamma':>6s} {'Hausman':>8s} {'mean k2S':>8s} {'mean kF':>8s} {'sec':>5s}")
    for des in ("N", "Q", "C2"):
        rng = np.random.default_rng(20260911); t0 = time.time()
        rl = rl3 = rg = rh = 0; k2 = []; kf = []
        for r in range(reps):
            df = dgp(rng, n, des)
            rl += link_test(df, "d", ["x", "z"])["p"] < 0.05
            rl3 += link_test(df, "d", ["x", "z"], order=3)["p"] < 0.05
            try:
                g = gamma_contrast(df, "y", ["x"], "d", ["x", "z"])
                rg += g["p"] < 0.05; rh += g["hausman"]["p"] < 0.05
                k2.append(g["hausman"]["kappa_2s"]); kf.append(g["hausman"]["kappa_fiml"])
            except Exception:
                pass
        print(f"{des:7s} {rl/reps:6.3f} {rl3/reps:6.3f} {rg/reps:6.3f} {rh/reps:8.3f} {np.mean(k2):8.3f} {np.mean(kf):8.3f} {time.time()-t0:5.0f}", flush=True)
