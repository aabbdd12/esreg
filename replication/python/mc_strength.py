"""Monte Carlo for Section 4.4: how the strength diagnostics track the precision of
kappa as the excluded instrument weakens (gamma_z in .9 .45 .2 0), n = 5000.
Reports the mean diagnostics and the sd / RMSE of kappa-hat by route."""
import sys, time, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg import esreg
from esreg.diag import diagnostics

def dgp(rng, n, gz):
    x = rng.normal(size=n); z = rng.normal(size=n); u = rng.normal(size=n)
    d = 0.2 + 0.3*x + gz*z + u > 0
    w1 = 0.78*u + rng.normal(scale=np.sqrt(1.3**2 - 0.78**2), size=n)
    w0 = -0.30*u + rng.normal(scale=np.sqrt(1.0 - 0.09), size=n)
    y = np.where(d, 1.0 + 1.2*x + w1, 0.5 + 0.4*x + w0)
    return pd.DataFrame(dict(y=y, d=d.astype(int), x=x, z=z))

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    n = 5000
    print(f"{'gz':>4s} {'LR_z':>8s} {'dR2':>6s} {'VarP|X':>7s} {'VIF1':>7s} {'VIF0':>7s} {'ATU extr':>8s} | {'2S: sd k':>8s} {'rmse':>6s} | {'FIML: sd k':>10s} {'rmse':>6s} {'sec':>5s}")
    for gz in (0.9, 0.45, 0.2, 0.0):
        rng = np.random.default_rng(20260911); t0 = time.time()
        D = []; k2 = []; kf = []
        for r in range(reps):
            df = dgp(rng, n, gz)
            r2 = esreg(df, "y", ["x"], "d", ["x", "z"], method="twostep")
            dg = diagnostics(r2, df, ["x"], ["x", "z"])
            D.append([dg["lr_excl"], dg["r2_incr"], dg["varP_share"], dg["vif1"], dg["vif0"], dg["share_atu_extrap"]])
            k2.append(r2.eff["est"][3])
            try:
                rf = esreg(df, "y", ["x"], "d", ["x", "z"], method="fiml")
                kf.append(rf.eff["est"][3])
            except Exception:
                kf.append(np.nan)
        D = np.array(D); k2 = np.array(k2); kf = np.array(kf)
        rm = lambda k: np.sqrt(np.nanmean((k - 1.08)**2))
        print(f"{gz:4.2f} {D[:,0].mean():8.1f} {D[:,1].mean():6.3f} {D[:,2].mean():7.3f} {D[:,3].mean():7.1f} {D[:,4].mean():7.1f} {D[:,5].mean():8.3f} | "
              f"{k2.std():8.3f} {rm(k2):6.3f} | {np.nanstd(kf):10.3f} {rm(kf):6.3f} {time.time()-t0:5.0f}")
