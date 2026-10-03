"""Monte Carlo for Section 4.5: the Hausman FIML-vs-two-step contrast (A3) and the
Hermite-control test (A2) on the designs N (joint normal), C1 (gain nonlinear in u,
A2 fails) and C2 (skewed regime errors, A3 fails, A2 holds) of paper A, n = 5000."""
import sys, time, numpy as np, pandas as pd
sys.path.insert(0, ".")
sys.path.insert(0, ".")
import mc_common as M
from esreg.tests import hausman_normal, wald_hermite

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    n = 5000
    print(f"{'design':7s} {'Hausman q':>9s} {'kappa z':>8s} {'Hermite':>8s} {'H reg1':>7s} {'H reg0':>7s} {'mean kdiff':>10s} {'sec':>5s}")
    for des in ("N", "C2", "C1"):
        rng = np.random.default_rng(20260911); t0 = time.time()
        rh = rk = rw = r1 = r0 = 0; kd = []
        for r in range(reps):
            D = M.dgp(rng, n, des)
            df = pd.DataFrame(dict(y=D["y"], d=D["d"].astype(int), x=D["x"], z=D["z"]))
            try:
                h = hausman_normal(df, "y", ["x"], "d", ["x", "z"])
                rh += h["p"] < 0.05; rk += h["kappa_p"] < 0.05; kd.append(h["kappa_diff"])
            except Exception as e:
                pass
            w = wald_hermite(df, "y", ["x"], "d", ["x", "z"])
            rw += w["p"] < 0.05; r1 += w["p_1"] < 0.05; r0 += w["p_0"] < 0.05
        print(f"{des:7s} {rh/reps:9.3f} {rk/reps:8.3f} {rw/reps:8.3f} {r1/reps:7.3f} {r0/reps:7.3f} {np.mean(kd):10.3f} {time.time()-t0:5.0f}")
