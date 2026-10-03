"""Weighted validation sample for esreg (pweights, hsize(), the svy prefix): esr_annex2
with a sampling weight wt = exp(.4 N(0,1)) and a household size hsize in 1..6 (seed 2026),
saved as esr_wt.dta.  Reference values from Python: FIML with robust (sandwich) variance
under pweights, two-step with the stacked variance, both with wt and with wt*hsize.
Under a pweight-only design, svy's linearized variance is n/(n-1) times the sandwich."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg import esreg
from paths import DATA as OUT
rng = np.random.default_rng(2026)
df = pd.read_stata(OUT + "esr_annex2.dta", convert_categoricals=False)
n = len(df)
df["wt"] = np.exp(0.4 * rng.standard_normal(n))
df["hsize"] = rng.integers(1, 7, n).astype(float)
df["wth"] = df.wt * df.hsize
for r in (2, 3):
    df[f"reg{r}"] = (df.region == r).astype(float)
df.to_stata(OUT + "esr_wt.dta", write_index=False, version=118,
            data_label="esr_annex2 with sampling weight wt and household size hsize (seed 2026)")
x = ["educ"]; z = ["educ", "reg2", "reg3", "inst"]
np.set_printoptions(precision=5, suppress=True)
for w in ("wt", "wth"):
    print("=" * 70, "weight", w)
    rF = esreg(df, "income", x, "treatment", z, method="fiml", weights=w, vce="robust")
    t = rF.table()
    print("FIML robust: ll", round(rF.fit["loglik"], 4))
    print(t.round(6).to_string())
    print(rF.effects().round(5).to_string())
    print("svy (pweight only) se factor sqrt(n/(n-1)) =", round(np.sqrt(n / (n - 1)), 6))
    rS = esreg(df, "income", x, "treatment", z, method="twostep", weights=w)
    print("two-step:"); print(rS.table().round(6).to_string()); print(rS.effects().round(5).to_string())
