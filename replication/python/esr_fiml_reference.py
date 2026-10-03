#!/usr/bin/env python3
# =====================================================================
#  ESR (Endogenous Switching Regression) -- implementation de reference
#
#  1. simule un echantillon a parametres CONNUS
#  2. sauvegarde les donnees en .dta (contre-epreuve movestay sous Stata)
#  3. estime le modele par FIML et compare aux vraies valeurs
#  4. ecrit un .do compagnon pour rejouer l'estimation sous Stata
#
#  A. Araar -- projet "Selection sur niveaux, tendances et gains"
# =====================================================================

import os
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize

# ---------------------------------------------------------------- config
CONFIG = dict(
    n      = 20_000,     # taille d'echantillon
    seed   = 2024,
    b1     = [1.0, 1.2],      # regime traite      : Y1 = b1_0 + b1_x * x + e1
    b0     = [0.5, 0.4],      # regime non traite  : Y0 = b0_0 + b0_x * x + e0
    g      = [0.2, 0.3, 0.9], # selection          : D* = g_0 + g_x * x + g_z * z + u
    sigma1 = 1.30,
    sigma0 = 1.00,
    rho1   = 0.60,       # corr(e1, u)
    rho0   = -0.30,      # corr(e0, u)   -> kappa = r1*s1 - r0*s0 > 0 : avantage comparatif
    cov10  = 0.15,       # cov(e1, e0) : NON IDENTIFIE, n'affecte ni ATT ni ATU ni ATE
    outdir = ".",
    stem   = "esr_sim",
)

# --------------------------------------------------------------- 1. DGP
def simulate(cfg):
    rng = np.random.default_rng(cfg["seed"])
    n = cfg["n"]
    s1, s0, r1, r0 = cfg["sigma1"], cfg["sigma0"], cfg["rho1"], cfg["rho0"]

    S = np.array([[s1**2,      cfg["cov10"], r1 * s1],
                  [cfg["cov10"], s0**2,      r0 * s0],
                  [r1 * s1,     r0 * s0,     1.0    ]])   # var(u) = 1 : NORMALISATION
    ev = np.linalg.eigvalsh(S)
    if ev.min() <= 0:
        raise ValueError(f"matrice de covariance non definie positive (vp min = {ev.min():.4f}) : "
                         "reduire |rho1|, |rho0| ou cov10")

    e1, e0, u = rng.multivariate_normal(np.zeros(3), S, size=n).T
    x = rng.normal(0, 1, n)
    z = rng.normal(0, 1, n)                      # restriction d'exclusion : z absent de X

    X = np.c_[np.ones(n), x]
    Z = np.c_[np.ones(n), x, z]
    b1, b0, g = map(np.asarray, (cfg["b1"], cfg["b0"], cfg["g"]))

    y1 = X @ b1 + e1
    y0 = X @ b0 + e0
    D  = (Z @ g + u > 0)
    Y  = np.where(D, y1, y0)

    df = pd.DataFrame({
        "y": Y, "d": D.astype(np.int8), "x": x, "z": z,
        "delta_true": y1 - y0,        # gain individuel : inconnaissable en pratique
        "y1_true": y1, "y0_true": y0, # contrefactuels  : pour validation seulement
    })
    truth = dict(
        b1_0=b1[0], b1_x=b1[1], b0_0=b0[0], b0_x=b0[1],
        g_0=g[0], g_x=g[1], g_z=g[2],
        sigma1=s1, sigma0=s0, rho1=r1, rho0=r0,
        kappa=r1 * s1 - r0 * s0,
        ATT=(y1 - y0)[D].mean(), ATU=(y1 - y0)[~D].mean(), ATE=(y1 - y0).mean(),
    )
    return df, X, Z, Y, D, truth

# ------------------------------------------------------- 2. FIML de l'ESR
def _nll(th, Y, X, Z, D):
    """Log-vraisemblance negative. th = [b1, b0, g, ln s1, ln s0, atanh r1, atanh r0]."""
    k, m = X.shape[1], Z.shape[1]
    b1, b0, g = th[:k], th[k:2*k], th[2*k:2*k+m]
    s1, s0 = np.exp(th[2*k+m]),   np.exp(th[2*k+m+1])
    r1, r0 = np.tanh(th[2*k+m+2]), np.tanh(th[2*k+m+3])

    zg = Z @ g
    e1 = (Y - X @ b1) / s1
    e0 = (Y - X @ b0) / s0
    a1 = (zg + r1 * e1) / np.sqrt(1 - r1**2)
    a0 = (zg + r0 * e0) / np.sqrt(1 - r0**2)

    ll = np.where(D,
                  norm.logpdf(e1) - np.log(s1) + norm.logcdf( a1),
                  norm.logpdf(e0) - np.log(s0) + norm.logcdf(-a0))
    return -ll.sum()


def _hessian(f, th, args, h=1e-4):
    """Hessienne par differences finies centrees."""
    p = len(th); H = np.zeros((p, p))
    for i in range(p):
        for j in range(i, p):
            ei = np.zeros(p); ei[i] = h
            ej = np.zeros(p); ej[j] = h
            H[i, j] = H[j, i] = (f(th+ei+ej, *args) - f(th+ei-ej, *args)
                                 - f(th-ei+ej, *args) + f(th-ei-ej, *args)) / (4*h*h)
    return H


def fit_esr(Y, X, Z, D):
    k, m = X.shape[1], Z.shape[1]

    # valeurs initiales : probit pour gamma, MCO par regime pour beta et sigma
    g0 = minimize(lambda g_: -np.where(D, norm.logcdf(Z@g_), norm.logcdf(-Z@g_)).sum(),
                  np.zeros(m), method="BFGS").x
    b1_0 = np.linalg.lstsq(X[D],  Y[D],  rcond=None)[0]
    b0_0 = np.linalg.lstsq(X[~D], Y[~D], rcond=None)[0]
    st = np.r_[b1_0, b0_0, g0,
               np.log(np.std(Y[D]  - X[D]  @ b1_0)),
               np.log(np.std(Y[~D] - X[~D] @ b0_0)), 0.0, 0.0]

    res = minimize(_nll, st, args=(Y, X, Z, D), method="BFGS", options={"maxiter": 5000})
    th = res.x

    # ecarts-types : Hessienne numerique + methode delta pour sigma = exp(.) et rho = tanh(.)
    try:
        se_th = np.sqrt(np.diag(np.linalg.inv(_hessian(_nll, th, (Y, X, Z, D)))))
    except np.linalg.LinAlgError:
        se_th = np.full(len(th), np.nan)

    s1, s0 = np.exp(th[2*k+m]),   np.exp(th[2*k+m+1])
    r1, r0 = np.tanh(th[2*k+m+2]), np.tanh(th[2*k+m+3])
    names = ["b1_0", "b1_x", "b0_0", "b0_x", "g_0", "g_x", "g_z"]
    est = dict(zip(names, th[:2*k+m]))
    se  = dict(zip(names, se_th[:2*k+m]))
    est.update(sigma1=s1, sigma0=s0, rho1=r1, rho0=r0)
    se.update(sigma1=s1 * se_th[2*k+m],        sigma0=s0 * se_th[2*k+m+1],
              rho1=(1-r1**2) * se_th[2*k+m+2], rho0=(1-r0**2) * se_th[2*k+m+3])

    est["kappa"] = r1 * s1 - r0 * s0
    se["kappa"]  = np.nan
    return est, se, th, res


def effects(est, X, Z, D):
    """ATT et ATU a partir des parametres estimes."""
    b1 = np.array([est["b1_0"], est["b1_x"]])
    b0 = np.array([est["b0_0"], est["b0_x"]])
    zg = Z @ np.array([est["g_0"], est["g_x"], est["g_z"]])
    l1 = norm.pdf(zg) / norm.cdf(zg)
    l0 = norm.pdf(zg) / (1 - norm.cdf(zg))
    k  = est["kappa"]
    att = (X[D]  @ (b1 - b0) + k * l1[D]).mean()
    atu = (X[~D] @ (b1 - b0) - k * l0[~D]).mean()
    ate = D.mean() * att + (1 - D.mean()) * atu
    return dict(ATT=att, ATU=atu, ATE=ate)

# ---------------------------------------------------- 3. export Stata
def export_stata(df, truth, cfg):
    os.makedirs(cfg["outdir"], exist_ok=True)
    dta = os.path.join(cfg["outdir"], cfg["stem"] + ".dta")
    do  = os.path.join(cfg["outdir"], cfg["stem"] + ".do")

    df.to_stata(
        dta, write_index=False, version=117,
        data_label="ESR simule -- parametres connus (voir .do compagnon)",
        variable_labels={
            "y": "Resultat observe", "d": "Traitement (1=adopte)",
            "x": "Covariable (dans X et dans Z)", "z": "Restriction d'exclusion (Z seul)",
            "delta_true": "Gain individuel vrai (non observable)",
            "y1_true": "Y1 contrefactuel vrai", "y0_true": "Y0 contrefactuel vrai",
        },
    )

    lines = [
        "* " + "=" * 66,
        "* Contre-epreuve Stata de l'ESR -- donnees a parametres CONNUS",
        "* Genere par esr_fiml_reference.py",
        "* " + "=" * 66,
        "*",
        "* Vraies valeurs :",
    ]
    for k_ in ("b1_0","b1_x","b0_0","b0_x","g_0","g_x","g_z",
               "sigma1","sigma0","rho1","rho0","kappa","ATT","ATU","ATE"):
        lines.append(f"*   {k_:<8} = {truth[k_]:+.6f}")
    lines += [
        "*",
        "* Normalisation : var(u) = 1. cov(e1,e0) n'est PAS identifie et n'affecte",
        "* ni ATT, ni ATU, ni ATE -- seulement la variance des gains individuels.",
        "*",
        f'use "{cfg["stem"]}.dta", clear',
        "",
        "* ssc install movestay",
        "* Verifier la syntaxe exacte avec -help movestay- (versions variables).",
        "movestay y x, select(d = x z)",
        "",
        "* Comparer aux vraies valeurs ci-dessus :",
        "*   rho_1 -> rho1 , rho_2 -> rho0 , sigma_1 -> sigma1 , sigma_2 -> sigma0",
        "*   ATTENTION a la convention de signe de u selon la version : si les rho",
        "*   sortent avec le signe oppose, c'est la convention, pas une erreur.",
        "",
        "* Contrefactuels par regime :",
        "* mspredict yhat1, yc1",
        "* mspredict yhat0, yc2",
        "",
        "* Test emboite ETR c. ESR (kappa = 0) :",
        "* etregress y x, treat(d = x z)",
        "",
        "* Verification directe, possible ici seulement car les vrais gains sont sauvegardes :",
        "summarize delta_true if d==1   // ATT vrai",
        "summarize delta_true if d==0   // ATU vrai",
        "summarize delta_true           // ATE vrai",
    ]
    with open(do, "w") as f:
        f.write("\n".join(lines) + "\n")
    return dta, do

# ------------------------------------------------------------- 4. main
if __name__ == "__main__":
    cfg = CONFIG
    df, X, Z, Y, D, truth = simulate(cfg)
    dta, do = export_stata(df, truth, cfg)

    est, se, th, res = fit_esr(Y, X, Z, D)
    eff = effects(est, X, Z, D)

    print(f"\nn = {cfg['n']:,}   part traitee = {D.mean():.3f}   "
          f"log-vraisemblance = {-res.fun:,.1f}   |gradient| = {np.linalg.norm(res.jac):.5f}")
    print("(le drapeau de convergence de BFGS porte sur son propre seuil de gradient "
          "numerique ;\n la norme du gradient ci-dessus est le critere pertinent)\n")

    print(f"{'parametre':>10} {'vrai':>10} {'estime':>10} {'ec.-type':>10} {'ecart':>9}")
    print("-" * 53)
    for k_ in ("b1_0","b1_x","b0_0","b0_x","g_0","g_x","g_z",
               "sigma1","sigma0","rho1","rho0","kappa"):
        s = f"{se[k_]:>10.4f}" if np.isfinite(se.get(k_, np.nan)) else f"{'--':>10}"
        print(f"{k_:>10} {truth[k_]:>10.4f} {est[k_]:>10.4f} {s} {est[k_]-truth[k_]:>+9.4f}")

    print("-" * 53)
    for k_ in ("ATT", "ATU", "ATE"):
        print(f"{k_:>10} {truth[k_]:>10.4f} {eff[k_]:>10.4f} {'':>10} {eff[k_]-truth[k_]:>+9.4f}")

    print(f"\nEcrit : {dta}\n        {do}")
