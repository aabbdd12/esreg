"""
Monte Carlo for the index-sufficiency test (Section 4.2).  Two excluded
instruments z1 (kept) and z2 (tested).  Designs:
  S   sufficiency holds (baseline of paper A with two instruments)
  Z   z2 enters the outcomes directly (-0.3 z2): exclusion fails for z2
  DH  double hurdle, D = 1{a1 + z1 + e1 > 0 and a2 + z2 + e2 > 0}, omega_j
      correlated with e1 only: a second index, P(Z) is not a single index
  U3  sufficiency holds but E[omega_j | u] is cubic (A2 fails): what the test
      picks up when the failure is in the shape, not in the instrument
usage: python3 mc_sufficiency.py [reps]
"""
import sys, time, numpy as np, pandas as pd
sys.path.insert(0, ".")
from esreg.tests import wald_sufficiency

def dgp(rng, n, design):
    x = rng.normal(size=n); z1 = rng.normal(size=n); z2 = rng.normal(size=n)
    if design == "DH":
        e1 = rng.normal(size=n); e2 = rng.normal(size=n)
        d = (0.3 + 0.3*x + 0.9*z1 + e1 > 0) & (0.6 + 0.3*x + 0.9*z2 + e2 > 0)
        u = e1
    else:
        u = rng.normal(size=n)
        d = 0.2 + 0.3*x + 0.7*z1 + 0.6*z2 + u > 0
    if design == "U3":
        h = (u**3 - 3*u)/np.sqrt(6)          # standardized Hermite, orthogonal to u
        w1 = 0.78*u + 0.6*h + rng.normal(scale=np.sqrt(1.3**2 - 0.78**2 - 0.36), size=n)
        w0 = -0.30*u + rng.normal(scale=np.sqrt(1.0 - 0.09), size=n)
    else:
        w1 = 0.78*u + rng.normal(scale=np.sqrt(1.3**2 - 0.78**2), size=n)
        w0 = -0.30*u + rng.normal(scale=np.sqrt(1.0 - 0.09), size=n)
    y1 = 1.0 + 1.2*x + w1; y0 = 0.5 + 0.4*x + w0
    if design == "Z":
        y1 = y1 - 0.3*z2; y0 = y0 - 0.3*z2
    y = np.where(d, y1, y0)
    return pd.DataFrame(dict(y=y, d=d.astype(int), x=x, z1=z1, z2=z2))

def run(reps, n, design, seed):
    rng = np.random.default_rng(seed); rej = rej1 = rej0 = 0; c = []; t0 = time.time()
    for r in range(reps):
        df = dgp(rng, n, design)
        t = wald_sufficiency(df, "y", ["x"], "d", ["x", "z1", "z2"], keep="z1")
        rej += t["p"] < 0.05; rej1 += t["p_1"] < 0.05; rej0 += t["p_0"] < 0.05
        c.append(t["coefs"])
    c = np.array(c)
    return rej/reps, rej1/reps, rej0/reps, c.mean(0), time.time() - t0

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    print(f"{'design':7s} {'n':>6s} {'joint':>6s} {'reg.1':>6s} {'reg.0':>6s} {'coef z2 (1)':>12s} {'coef z2 (0)':>12s} {'sec':>5s}")
    for n in (2000, 5000):
        for des in ("S", "Z", "DH", "U3"):
            rj, r1, r0, cm, tt = run(reps, n, des, 20260911 + n)
            print(f"{des:7s} {n:6d} {rj:6.3f} {r1:6.3f} {r0:6.3f} {cm[0]:12.3f} {cm[1]:12.3f} {tt:5.0f}")
