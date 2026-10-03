"""Monte Carlo for Section 4.6: the pseudo-DiD on the selection index.  Designs of
paper A: N (kappa = 1.08), H (kappa = 0 with rho1 sigma1 = rho0 sigma0 = 0.5), C2
(skewed), C1 (nonlinear gain); n = 5000.  kappa_dd with nq = 2 and 4 strata
against the two-step kappa; rejection of kappa = 0 at 5%."""
import sys, time, numpy as np, pandas as pd
sys.path.insert(0, ".")
import mc_common as M
from esreg.tests import pseudo_did
from esreg import esreg
if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    n = 5000
    print(f"{'design':6s} {'nq':>3s} {'kdd mean':>8s} {'sd':>6s} {'rej k=0':>7s} {'C1<0':>5s} {'C0>0':>5s} | {'2S kappa':>8s} {'sd':>6s}")
    for des in ("N", "H", "C2", "C1"):
        for nq in (2, 4):
            rng = np.random.default_rng(20260911); kd = []; k2 = []; rej = 0; c1n = 0; c0p = 0
            for r in range(reps):
                D = M.dgp(rng, n, des)
                df = pd.DataFrame(dict(y=D["y"], d=D["d"].astype(int), x=D["x"], z=D["z"]))
                t = pseudo_did(df, "y", ["x"], "d", ["x", "z"], nq=nq)
                kd.append(t["kappa_dd"]); rej += abs(t["kappa_dd"] / t["se_kappa_dd"]) > 1.96
                c1n += t["C1"] / t["se_C1"] < -1.96; c0p += t["C0"] / t["se_C0"] > 1.96
                if nq == 2:
                    k2.append(esreg(df, "y", ["x"], "d", ["x", "z"], method="twostep").eff["est"][3])
            kd = np.array(kd)
            print(f"{des:6s} {nq:3d} {kd.mean():8.3f} {kd.std():6.3f} {rej/reps:7.3f} {c1n/reps:5.2f} {c0p/reps:5.2f} | "
                  + (f"{np.mean(k2):8.3f} {np.std(k2):6.3f}" if nq == 2 else ""))
