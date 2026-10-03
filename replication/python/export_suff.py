"""Samples for the Stata validation of esrtest, suff (Section 4.2): two excluded
instruments z1 z2, n = 5000, seed 2026.  esr_suff.dta: sufficiency holds;
esr_dh.dta: double hurdle (a second index).  Reference values from Python."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg.tests import wald_sufficiency
from mc_sufficiency import dgp
from paths import DATA as OUT
for des, name, lab in (("S", "esr_suff", "Two instruments, index sufficiency holds; n = 5000, seed 2026"),
                       ("DH", "esr_dh", "Double hurdle: D = 1{.3+.3x+.9z1+e1>0 & .6+.3x+.9z2+e2>0}, omega ~ e1; n = 5000, seed 2026")):
    rng = np.random.default_rng(2026)
    df = dgp(rng, 5000, des)
    df.to_stata(OUT + name + ".dta", write_index=False, version=118, data_label=lab)
    for keep in ("z1", None):
        t = wald_sufficiency(df, "y", ["x"], "d", ["x", "z1", "z2"], keep=keep)
        print(f"{name} keep={t['keep']} tested={t['tested']}: coefs {t['coefs']} se {t['se']} | joint chi2({t['df']}) = {t['chi2']:.4f} p={t['p']:.2e} | "
              f"reg1 chi2({t['df_1']}) = {t['chi2_1']:.4f} | reg0 chi2({t['df_0']}) = {t['chi2_0']:.4f}")
    t = wald_sufficiency(df, "y", ["x"], "d", ["x", "z1", "z2"])   # default keep = strongest
    print("   default keep:", t["keep"])
