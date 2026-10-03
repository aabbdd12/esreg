"""
Validation of the Python esreg against movestay (SJ build 2.0.0) and msat 1.3.2
on the four datasets of the replication package.  Reference values come from the
Stata logs (msat_validate.log, msat13_test.log).
"""
import sys, time, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg import esreg

from paths import DATA

def show(title, r, ref):
    print("=" * 78); print(title); print("=" * 78)
    r.summary()
    e = r.effects()
    print("\nReference (Stata):")
    for k, v in ref.items(): print(f"   {k:>8s}: {v}")
    print()

# ---- esr_sim, FIML -------------------------------------------------------
df = pd.read_stata(DATA + "esr_sim.dta")
t = time.time(); r = esreg(df, "y", ["x"], "d", ["x", "z"], method="fiml"); t = time.time() - t
show(f"esr_sim  FIML  ({t:.1f}s)", r,
     {"movestay": "y_1: x 1.205903 (.0121553), _cons .9831363 (.0189707); y_0: x .3865108 (.0109584), _cons .4803906 (.0196117); "
                  "d: x .292636 (.0103871), z .9081515 (.0128427), _cons .2042513 (.0102002); lns1 .2549284 (.0083776), lns2 -.0094671 (.0082558); "
                  "sigma_1 1.290369, sigma_0 .9905775, rho_1 .6038142, rho_0 -.3129654; ll -39940.291; LR 787.01",
      "msat": "ATT 1.17919 (.0342352)  ATU -.343526 (.0380632)  ATE .508128 (.0278977)  kappa 1.08916 (.0373076); "
              "var comp ATT: param .00111 samp .0000584; support [0.0213, 0.9910]"})

# ---- esr_sim, two-step ----------------------------------------------------
r2 = esreg(df, "y", ["x"], "d", ["x", "z"], method="twostep")
show("esr_sim  two-step", r2, {"msat twostep": "ATT 1.17991  ATU -.287578  kappa 1.04355"})

# ---- esr_annex2, FIML -----------------------------------------------------
da = pd.read_stata(DATA + "esr_annex2.dta")
da["region"] = da["region"].astype(int)
for j in (2, 3): da[f"reg{j}"] = (da["region"] == j).astype(float)
r3 = esreg(da, "income", ["educ"], "treatment", ["educ", "reg2", "reg3", "inst"], method="fiml")
show("esr_annex2  FIML", r3,
     {"movestay": "income_1: educ .5039501 (.0052076) _cons 62.8872 (.0187406); income_0: educ .503302 (.0073309) _cons 60.9076 (.0297132); "
                  "sel: educ -.0034487 reg2 .2872393 reg3 .3312872 inst .7008038 _cons .3195647; sigma_1 .3566222 sigma_0 .366249 rho_1 .5304626 rho_0 .5189426; ll -4229.3298; LR 172.38",
      "msat": "ATT 1.98115 (.0291947)  ATU 1.98233 (.0267532)  ATE 1.98155 (.0230826)  kappa -.000887422 (.0263384); support [0.0421, 0.9754]"})

# ---- esr_skew: FIML vs two-step, with the bootstrap SEs of the two-step ----
ds = pd.read_stata(DATA + "esr_skew.dta")
r4 = esreg(ds, "y", ["x"], "d", ["x", "z"], method="fiml")
show("esr_skew  FIML", r4, {"msat": "ATT 2.33027  ATU -1.63347  kappa 2.82489"})
r5 = esreg(ds, "y", ["x"], "d", ["x", "z"], method="twostep")
show("esr_skew  two-step (stacked-moment SEs vs bootstrap-100 in Stata)", r5,
     {"msat twostep": "ATT 1.25885  ATU -.344916  kappa 1.16692;  bootstrap SEs: ATT .0539  ATU .0747  ATE .0446  kappa .0742"})

# ---- esr_nonlin two-step --------------------------------------------------
dn = pd.read_stata(DATA + "esr_nonlin.dta")
r6 = esreg(dn, "y", ["x"], "d", ["x", "z"], method="twostep")
show("esr_nonlin  two-step", r6, {"msat twostep": "ATT 1.08803  ATU -.210066  kappa .87816"})
