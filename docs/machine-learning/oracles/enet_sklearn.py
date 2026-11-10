"""Reference values for ElasticNetRegression from scikit-learn ElasticNet and enet_path.

Run from the repository root after enet_glmnet.R:  python docs/machine-learning/oracles/enet_sklearn.py
Writes docs/machine-learning/oracles/enet-sklearn.json.

Parameter mapping. glmnet (and ElasticNetRegression) standardize the response before fitting, so a
glmnet solution at penalty lambda and mixing alpha equals scikit-learn fitted to
    Xs = (X - mean) / sd_pop(X),   ys = (y - mean) / sd_pop(y)
with alpha_sklearn = lambda / sd_pop(y) and l1_ratio = alpha, back-transformed by
    beta_j = b_j * sd_pop(y) / sd_pop(X_j),   intercept = mean(y) - sum_j beta_j * mean(X_j).
For alpha = 1 the response scaling cancels: alpha_sklearn = lambda on Xs with the raw centered y.
"""
import json
import pathlib

import numpy as np
import sklearn
from sklearn.linear_model import ElasticNet, enet_path

root = pathlib.Path(__file__).resolve().parents[3]
here = pathlib.Path(__file__).resolve().parent
d = np.loadtxt(root / "docs/example-data/prostate-cancer.csv", delimiter=",", comments="#", skiprows=4)
d = d[d[:, 9] == 1]
X, y = d[:, :8], d[:, 8]
glm = json.loads((here / "enet-glmnet.json").read_text())

mx, sx = X.mean(axis=0), X.std(axis=0)
my, sy = y.mean(), y.std()
Xs, ys = (X - mx) / sx, (y - my) / sy
INDEXES = [5, 20, 40, 60]


def back(b):
    beta = b * sy / sx
    return beta, my - beta @ mx


out = {"sklearn": sklearn.__version__, "numpy": np.__version__, "indexes": INDEXES}
for name, alpha in (("lasso", 1.0), ("enet", 0.5)):
    lam = np.array(glm[name]["lambda"])
    rows = []
    for k in INDEXES:
        model = ElasticNet(alpha=lam[k] / sy, l1_ratio=alpha, tol=1e-30, max_iter=1_000_000).fit(Xs, ys)
        beta, a0 = back(model.coef_)
        rows.append({"lambda": lam[k], "intercept": a0, "beta": beta.tolist()})
    out[name] = rows

# The same points from one call to enet_path on the standardized data.
lam = np.array(glm["enet"]["lambda"])
_, coefs, _ = enet_path(Xs, ys, l1_ratio=0.5, alphas=lam / sy, tol=1e-30, max_iter=1_000_000)
out["enet_path"] = [back(coefs[:, k])[0].tolist() for k in INDEXES]

# Lasso without response scaling: alpha_sklearn = lambda on standardized X and centered y.
lam = np.array(glm["lasso"]["lambda"])
model = ElasticNet(alpha=lam[20], l1_ratio=1.0, tol=1e-30, max_iter=1_000_000).fit(Xs, y - my)
out["lasso_unscaled_y_index20"] = (model.coef_ / sx).tolist()

(here / "enet-sklearn.json").write_text(json.dumps(out, indent=1))
print(here / "enet-sklearn.json")
