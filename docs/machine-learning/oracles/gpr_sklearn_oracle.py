"""scikit-learn oracle values for GaussianProcessRegression.

Writes gpr-sklearn.json with:
  * the fixture: 20 two-dimensional training points, 5 test points (rounded to 6 decimals,
    so the rounded values are the data);
  * for each kernel (RBF, Matern 3/2, Matern 5/2, Matern 1/2 = exponential), with fixed
    hyperparameters sigma_f^2 = 1.3, length scales (0.7, 1.4), noise 0.05 and a zero mean
    (normalize_y=False): the predictive mean and latent variance at the test points, the log
    marginal likelihood at theta and its gradient;
  * a fitted Matern 5/2 ARD model with an estimated noise level (WhiteKernel), using 50
    optimizer restarts: the optimal theta and log marginal likelihood.

Usage:
    python gpr_sklearn_oracle.py <out.json>
"""
import json
import sys

import numpy as np
import sklearn
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern, WhiteKernel

SF2, ELL, NOISE = 1.3, [0.7, 1.4], 0.05


def fixture():
    rng = np.random.RandomState(20261105)
    x = np.round(rng.uniform(0.0, 5.0, size=(20, 2)), 6)
    y = np.round(np.sin(x[:, 0]) + 0.5 * np.cos(2.0 * x[:, 1]) + 0.1 * rng.standard_normal(20), 6)
    xt = np.round(rng.uniform(0.0, 5.0, size=(5, 2)), 6)
    return x, y, xt


def base_kernel(name, ell, bounds="fixed"):
    if name == "SquaredExponential":
        return RBF(length_scale=ell, length_scale_bounds=bounds)
    nu = {"Matern32": 1.5, "Matern52": 2.5, "Exponential": 0.5}[name]
    return Matern(length_scale=ell, length_scale_bounds=bounds, nu=nu)


def main(out):
    x, y, xt = fixture()
    result = {
        "sklearn": sklearn.__version__,
        "numpy": np.__version__,
        "x": x.tolist(), "y": y.tolist(), "xTest": xt.tolist(),
        "signalVariance": SF2, "lengthScales": ELL, "noiseVariance": NOISE,
        "fixed": {},
    }
    for name in ["SquaredExponential", "Matern32", "Matern52", "Exponential"]:
        kernel = ConstantKernel(SF2, "fixed") * base_kernel(name, ELL)
        gp = GaussianProcessRegressor(kernel=kernel, alpha=NOISE, optimizer=None, normalize_y=False).fit(x, y)
        mean, std = gp.predict(xt, return_std=True)
        # The white-noise form has the same likelihood and exposes the noise derivative.
        k2 = ConstantKernel(SF2, (1e-5, 1e5)) * base_kernel(name, ELL, (1e-5, 1e5)) + WhiteKernel(NOISE, (1e-8, 1e5))
        gp2 = GaussianProcessRegressor(kernel=k2, alpha=0.0, optimizer=None).fit(x, y)
        lml, grad = gp2.log_marginal_likelihood(gp2.kernel_.theta, eval_gradient=True)
        result["fixed"][name] = {
            "mean": mean.tolist(), "variance": (std ** 2).tolist(),
            "theta": gp2.kernel_.theta.tolist(), "logMarginalLikelihood": float(lml), "gradient": grad.tolist(),
        }
    kernel = ConstantKernel(1.0, (1e-5, 1e5)) * Matern([1.0, 1.0], (1e-3, 1e3), nu=2.5) + WhiteKernel(0.1, (1e-8, 1e2))
    gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0, n_restarts_optimizer=50, random_state=0).fit(x, y)
    result["fitted"] = {
        "kernel": "Matern52", "theta": gp.kernel_.theta.tolist(),
        "logMarginalLikelihood": float(gp.log_marginal_likelihood_value_),
    }
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
