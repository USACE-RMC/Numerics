"""Reference values for PrincipalComponentAnalysis from scikit-learn PCA.

Run from the repository root:  python docs/machine-learning/oracles/pca_oracle.py
Writes pca-sklearn.json and the random test matrix pca-random-matrix.csv next to this script.

Correlation basis: scikit-learn has no correlation option, so the centered data are divided by the
sample standard deviation (ddof = 1) before fitting, which gives the eigenvalues of the correlation
matrix (R prcomp with scale. = TRUE).
"""
import json
import pathlib

import numpy as np
import sklearn
from sklearn.decomposition import PCA

root = pathlib.Path(__file__).resolve().parents[3]
iris = np.loadtxt(root / "docs/example-data/iris-dataset.csv", delimiter=",", comments="#", skiprows=4)[:, :4]

# Random 12 x 5 matrix with correlated columns, printed in full so the test can embed it.
rng = np.random.default_rng(20261105)
A = rng.normal(size=(12, 5)) @ np.array([[1.0, 0.5, 0.2, 0.0, 0.1],
                                          [0.0, 1.0, 0.4, 0.3, 0.0],
                                          [0.0, 0.0, 1.0, 0.6, 0.2],
                                          [0.0, 0.0, 0.0, 1.0, 0.5],
                                          [0.0, 0.0, 0.0, 0.0, 1.0]]) + np.array([10.0, -3.0, 0.5, 100.0, 0.0])
A = np.round(A, 6)


def fit(X, basis):
    Z = X - X.mean(axis=0)
    if basis == "correlation":
        Z = Z / X.std(axis=0, ddof=1)
    pca = PCA(svd_solver="full").fit(Z)
    return {
        "explained_variance": pca.explained_variance_.tolist(),
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "singular_values": pca.singular_values_.tolist(),
        "components": pca.components_.tolist(),
        "scores_first3": pca.transform(Z)[:3].tolist(),
    }


out = {
    "sklearn": sklearn.__version__,
    "numpy": np.__version__,
    "random_matrix": A.tolist(),
    "iris_covariance": fit(iris, "covariance"),
    "iris_correlation": fit(iris, "correlation"),
    "random_covariance": fit(A, "covariance"),
    "random_correlation": fit(A, "correlation"),
    # Wide case (p > n): the first four rows; the last component has zero variance.
    "wide_covariance": fit(A[:4], "covariance"),
}
np.savetxt(pathlib.Path(__file__).with_name("pca-random-matrix.csv"), A, delimiter=",", fmt="%.6f")
path = pathlib.Path(__file__).with_name("pca-sklearn.json")
path.write_text(json.dumps(out, indent=1))
print(path)
