# Machine learning

`pynumerics.machine_learning` wraps the models in `Numerics.MachineLearning`. Each class holds one
Numerics model in its `net` attribute and mirrors that model's members. The constructor takes the
training data, `train` (or `fit` for the elastic net) runs the Numerics estimator, and `predict` and the
other queries call the matching Numerics members. Every algorithm, default, tolerance and random seed
comes from Numerics.

Features are NumPy arrays with one row per observation and one column per feature. Each call copies
them to .NET in one bulk transfer and returns NumPy arrays. A 1-D feature array counts as a single
feature column, as in the Numerics `double[]` overloads. Three methods instead read a 1-D array as one
observation, and their docstrings say so: `ElasticNetRegression.predict`,
`PrincipalComponentAnalysis.transform` and `Kriging.predict_with_variance`.

| Area | Classes |
| --- | --- |
| Linear models | `GeneralizedLinearModel`, `ElasticNetRegression` |
| Trees and ensembles | `DecisionTree`, `DecisionNode`, `RandomForest` |
| Instance and probabilistic classifiers | `KNearestNeighbors`, `NaiveBayes` |
| Gaussian processes and kriging | `GaussianProcessRegression`, `Kriging` |
| Clustering | `KMeans`, `GaussianMixtureModel`, `JenksNaturalBreaks`, `JenksCluster` |
| Dimension reduction | `PrincipalComponentAnalysis` |
| Enums | `GaussianProcessKernel`, `GaussianProcessMean`, `VariogramModel`, `KrigingType`, `PrincipalComponentBasis` |

Models that draw random numbers take a seed. A positive seed makes the results reproducible, and zero
or a negative seed seeds the generator from the clock.

## Generalized linear models

```python
import numpy as np
from pynumerics.machine_learning import GeneralizedLinearModel

rng = np.random.default_rng(1)
x = rng.normal(size=(200, 2))
p = 1.0 / (1.0 + np.exp(-(0.5 + x @ [1.0, -0.8])))
y = (rng.random(200) < p).astype(float)

model = GeneralizedLinearModel(x, y, link="logit")
model.train()
print(model.parameters)                 # intercept first
print(model.parameter_standard_errors)
print(model.aic, model.bic)          # the normal family uses the MLE variance SSE / n, as statsmodels
print(model.predict(x[:3]))             # mean response
print(model.predict(x[:3], alpha=0.1))  # columns lower, mean, upper
print("\n".join(model.summary()))
```

The available links are `"identity"` (normal), `"log"` (Poisson), and `"logit"`, `"probit"` and
`"complementary_log_log"` (Bernoulli). If you set `use_robust_se = True` before training, the
standard errors are HC0 sandwich estimates.

## Elastic net

```python
import numpy as np
from pynumerics.machine_learning import ElasticNetRegression

rng = np.random.default_rng(2)
x = rng.normal(size=(100, 6))
y = x @ [1.5, 0.0, -0.7, 0.0, 2.0, 0.0] + rng.normal(size=100)

model = ElasticNetRegression(x, y)
model.alpha = 0.5                 # 1 is the lasso, 0 is ridge
model.cross_validate(folds=5, seed=12345)
best = model.lambda_min_index
print(model.lambda_min, model.get_coefficients(best))
print(model.selected_variables(best))
print(model.predict(x[:3], best))
```

The `coefficients` attribute holds the whole path, with shape (penalties, p).

## Trees, forests, k-nearest neighbors and naive Bayes

```python
import numpy as np
from pynumerics.machine_learning import DecisionTree, KNearestNeighbors, NaiveBayes, RandomForest

rng = np.random.default_rng(3)
x = rng.normal(size=(150, 3))
y = x[:, 0] - 2.0 * x[:, 1] + 0.1 * rng.normal(size=150)
labels = (x[:, 0] > 0).astype(float)

tree = DecisionTree(x, y, seed=12345)
tree.max_depth = 6
tree.train()
print(tree.predict(x[:3]), tree.root.feature_index)

forest = RandomForest(x, y, seed=12345)
forest.number_of_trees = 100
forest.train()
print(forest.predict(x[:3]))      # columns lower, median, upper, mean

knn = KNearestNeighbors(x, y, 5)
print(knn.predict(x[:3]), knn.get_neighbors(x[:3]))

bayes = NaiveBayes(x, labels)
bayes.train()
print(bayes.predict(x[:5]), bayes.priors)
```

For classification, set `is_regression = False` on the tree, the forest or the nearest-neighbor model,
and pass the class labels as numbers. If a model is untrained, or its input has the wrong number of
columns, `predict` returns `None`, because the Numerics member returns null.

## Gaussian process regression and kriging

```python
import numpy as np
from pynumerics.machine_learning import GaussianProcessRegression, Kriging

x = np.linspace(0.0, 6.0, 30)
y = np.sin(x)
gp = GaussianProcessRegression(x, y, kernel="squared_exponential", mean="constant")
gp.train()
mean, variance = gp.predict_with_variance(np.array([1.5, 7.0]))
print(gp.log_marginal_likelihood, gp.length_scales, mean, variance)

rng = np.random.default_rng(4)
sites = rng.uniform(0.0, 10.0, size=(20, 2))
values = np.sin(sites[:, 0]) + 0.1 * sites[:, 1]
krig = Kriging(sites, values, "exponential", nugget=0.0, partial_sill=1.0, range=3.0)
print(krig.predict_with_variance(np.array([[5.0, 5.0], [1.0, 9.0]])))
print(krig.semivariance([0.0, 1.0, 5.0]))
```

`Kriging` also accepts a custom `distance(a, b)`, given as a Python callable, as C# source compiled
with `pynumerics.callbacks`, or as a .NET `Func<double[], double[], double>`.

## Clustering and principal components

```python
import numpy as np
from pynumerics.machine_learning import (GaussianMixtureModel, JenksNaturalBreaks, KMeans,
                                         PrincipalComponentAnalysis)

rng = np.random.default_rng(5)
data = np.vstack([rng.normal(c, 0.5, size=(50, 2)) for c in ([0, 0], [4, 4], [0, 5])])

kmeans = KMeans(data, 3)
kmeans.train(seed=12345)
print(kmeans.means, np.bincount(kmeans.labels))

mixture = GaussianMixtureModel(data, 3)
mixture.train(seed=12345)
print(mixture.weights, mixture.log_likelihood, mixture.sigmas.shape)

breaks = JenksNaturalBreaks(rng.normal(size=40), 3)
print(breaks.breaks, breaks.goodness_of_variance_fit)

pca = PrincipalComponentAnalysis(data, basis="correlation")
pca.train()
print(pca.explained_variance_ratio, pca.components)
print(pca.transform(data[:2]))
```

`KMeans.initialize(x, k, seed)` returns the starting centroids that `train` uses with the same seed.
The PCA axes follow the scikit-learn sign convention: the element of largest magnitude in each axis is
positive.
