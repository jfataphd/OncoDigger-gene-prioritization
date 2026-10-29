"""
factor_diagnostics.py -- which factor of R = W * ln(1+E) * ln(1+S) actually orders the genes? (exploratory, reads features_raw.csv)
For each cancer, over the 100 highest-scoring candidate genes: Spearman correlation of each factor with R, the correlation between W and S
(the 'support counted twice' question), and the share of the variance of ln R carried by each log-factor. Output: factor_diagnostics.csv
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
HERE = Path(__file__).resolve().parent
f = pd.read_csv(HERE / "features_raw.csv")
f["Be"] = np.maximum(f.B, f.S); f["E"] = (f.S / f.N) / (f.Be / f.C)
f["lE"] = np.log1p(f.E); f["lS"] = np.log1p(f.S); f["R"] = f.W * f.lE * f.lS
rows = []
for c, g in f.groupby("cancer", sort=False):
    g = g.sort_values("R", ascending=False)
    for label, h in (("top100", g.head(100)), ("top10", g.head(10))):
        lr = np.log(h.R)
        parts = {"W": np.log(h.W), "lnE": np.log(h.lE), "lnS": np.log(h.lS)}
        cov = {k: np.cov(v, lr)[0, 1] / np.var(lr, ddof=1) for k, v in parts.items()}
        rows.append({"cancer": c, "set": label, "rho_R_W": spearmanr(h.R, h.W)[0], "rho_R_S": spearmanr(h.R, h.S)[0], "rho_R_lnE": spearmanr(h.R, h.lE)[0],
                     "rho_W_S": spearmanr(h.W, h.S)[0], "rho_S_E": spearmanr(h.S, h.E)[0],
                     "var_share_W": cov["W"], "var_share_lnE": cov["lnE"], "var_share_lnS": cov["lnS"]})
d = pd.DataFrame(rows); d.to_csv(HERE / "factor_diagnostics.csv", index=False)
print(d.groupby("set").mean(numeric_only=True).round(3).T.to_string())
