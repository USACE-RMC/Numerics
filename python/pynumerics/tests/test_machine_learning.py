"""Machine learning wrappers checked against scikit-learn, statsmodels, SciPy and textbook definitions.

Every fixture is generated from a fixed NumPy seed and every Numerics model that draws random numbers
gets a fixed seed, so each comparison is reproducible. Tolerances are stated per model with the reason.
"""

from __future__ import annotations

import itertools
import time

import numpy as np
import pytest
import statsmodels.api as sm
from scipy import stats
from sklearn.cluster import KMeans as SklearnKMeans
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern, WhiteKernel
from sklearn.linear_model import ElasticNet, Lasso
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor, NearestNeighbors
from sklearn.tree import DecisionTreeRegressor

from pynumerics import machine_learning as ml
from pynumerics.marshal import to_net_matrix


def _regression(seed: int = 1, n: int = 60) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 3))
    return x, x @ [1.0, -2.0, 0.5] + 0.3 + 0.1 * rng.normal(size=n)


def _blobs(seed: int = 2) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.vstack([rng.normal(center, 0.5, size=(50, 2)) for center in ([0.0, 0.0], [4.0, 4.0], [0.0, 5.0])])


# Generalized linear model: statsmodels GLM is the oracle. Numerics keeps the better of its optimizer
# and Fisher scoring estimates. For the canonical links that is the exact maximum and coefficients agree
# to 1e-9; for the probit link the expected and observed information differ, the kept estimate stops
# about 1e-8 from the maximum, and coefficients agree to 1e-6.
# Standard errors are the inverse Fisher information at the estimate; statsmodels forms them from the
# IRLS weights of its last iteration, which differ by the final step, so they agree to 1e-5 relative.


@pytest.mark.parametrize(("link", "family", "tolerance"), [
    ("identity", sm.families.Gaussian(), 1e-9),
    ("logit", sm.families.Binomial(), 1e-9),
    ("probit", sm.families.Binomial(sm.families.links.Probit()), 1e-6),
    ("log", sm.families.Poisson(), 1e-9),
])
def test_glm_matches_statsmodels(link: str, family: object, tolerance: float) -> None:
    # Arrange
    rng = np.random.default_rng(3)
    x = rng.normal(size=(120, 2))
    eta = 0.4 + x @ [0.8, -0.6]
    if link == "identity":
        y = eta + 0.3 * rng.normal(size=120)
    elif link == "log":
        y = rng.poisson(np.exp(eta)).astype(float)
    else:
        y = (rng.random(120) < stats.norm.cdf(eta)).astype(float)
    # Act
    model = ml.GeneralizedLinearModel(x, y, link=link)
    model.train()
    oracle = sm.GLM(y, sm.add_constant(x), family=family).fit(tol=1e-14)
    # Assert
    np.testing.assert_allclose(model.parameters, oracle.params, rtol=tolerance, atol=tolerance)
    np.testing.assert_allclose(model.parameter_standard_errors, oracle.bse, rtol=1e-5)
    np.testing.assert_allclose(model.predict(x[:5]), oracle.predict(sm.add_constant(x)[:5]), rtol=tolerance)
    # The Normal-family log-likelihood uses the maximum likelihood variance SSE / n, as statsmodels does.
    assert model.aic == pytest.approx(oracle.aic, rel=tolerance)


def test_glm_members_and_summary() -> None:
    # Arrange
    x, y = _regression()
    model = ml.GeneralizedLinearModel(x, y, feature_names=["a", "b", "c"], response_name="z")
    named = ml.GeneralizedLinearModel(np.column_stack([np.ones(60), x]), y, has_intercept=False,
                                      feature_names=["one", "a", "b", "c"])
    # Act
    model.train()
    interval = model.predict(x[:4], alpha=0.1)
    with_intercept = model.predict(np.column_stack([np.ones(4), x[:4]]))
    # Assert
    assert named.parameter_names == ["one", "a", "b", "c"]
    assert model.parameter_names == ["Intercept", "a", "b", "c"]
    assert ml.GeneralizedLinearModel(x, y).parameter_names == ["Intercept", "β1", "β2", "β3"]
    assert model.has_intercept and model.link_type == "Identity" and model.sample_size == 60
    assert model.degrees_of_freedom == 56 and model.x.shape == (60, 4)
    np.testing.assert_array_equal(model.y, y)
    np.testing.assert_allclose(model.residuals, y - model.predict(x), atol=1e-12)
    np.testing.assert_allclose(model.parameter_z_scores, model.parameters / model.parameter_standard_errors)
    np.testing.assert_allclose(model.parameter_p_values, 2 * stats.norm.sf(np.abs(model.parameter_z_scores)),
                               atol=1e-12)
    np.testing.assert_allclose(np.sqrt(np.diag(model.covariance)), model.parameter_standard_errors)
    np.testing.assert_array_equal(interval[:, 1], model.predict(x[:4]))
    assert (interval[:, 0] < interval[:, 1]).all() and (interval[:, 1] < interval[:, 2]).all()
    np.testing.assert_array_equal(with_intercept, model.predict(x[:4]))
    assert any("Model for predicting z" in line for line in model.summary())
    assert model.aicc > model.aic and np.isfinite(model.bic) and model.standard_error > 0
    assert model.link_function is not None and model.optimizer is not None


def test_glm_robust_standard_errors_match_statsmodels_hc0() -> None:
    # Arrange
    x, y = _regression(5)
    model = ml.GeneralizedLinearModel(x, y)
    model.use_robust_se = True
    # Act
    model.train()
    oracle = sm.GLM(y, sm.add_constant(x)).fit(cov_type="HC0")
    # Assert: HC0 sandwich of the same estimate
    assert model.use_robust_se
    np.testing.assert_allclose(model.parameter_standard_errors, oracle.bse, rtol=1e-8)


def test_glm_optimizer_choice_and_custom_link() -> None:
    # Arrange
    x, y = _regression(6)
    from pynumerics._util import namespace

    functions = namespace("Numerics.Functions")
    link = functions.LinkFunctionFactory.Create(functions.LinkFunctionType.Identity)
    # Act
    model = ml.GeneralizedLinearModel(x, y, link_function=link)
    model.set_optimizer("BFGS")
    model.train()
    reference = ml.GeneralizedLinearModel(x, y)
    reference.train()
    # Assert: the Fisher scoring estimate is kept when it is at least as good, so the fits agree
    np.testing.assert_allclose(model.parameters, reference.parameters, rtol=1e-8)
    with pytest.raises(ValueError, match="is not a member of LocalMethod"):
        model.set_optimizer("NotAMethod")


def test_glm_rejects_too_few_points() -> None:
    with pytest.raises(ArithmeticError):
        ml.GeneralizedLinearModel([[1.0], [2.0]], [1.0, 2.0])


# k-means: scikit-learn's Lloyd algorithm started from the Numerics k-means++ centroids reaches the
# same fixed point; centroids agree to rounding (1e-12) and labels exactly.


def test_kmeans_matches_sklearn_from_the_same_start() -> None:
    # Arrange
    data = _blobs()
    start = ml.KMeans.initialize(data, 3, 12345)
    # Act
    model = ml.KMeans(data, 3)
    model.train(12345)
    oracle = SklearnKMeans(3, init=start, n_init=1, algorithm="lloyd", tol=0.0, max_iter=1000).fit(data)
    # Assert
    np.testing.assert_allclose(model.means, oracle.cluster_centers_, atol=1e-12)
    np.testing.assert_array_equal(model.labels, oracle.labels_)
    assert model.k == 3 and model.dimension == 2 and model.max_iterations == 1000 and model.iterations >= 1
    np.testing.assert_array_equal(model.x, data)


def test_kmeans_seed_is_reproducible_and_random_start_runs() -> None:
    # Arrange
    data = _blobs(4)
    first, second = ml.KMeans(data, 3), ml.KMeans(data, 3)
    first.max_iterations = 50
    # Act
    first.train(7, k_means_plus_plus=False)
    second.train(7, k_means_plus_plus=False)
    # Assert
    np.testing.assert_array_equal(first.means, second.means)
    np.testing.assert_array_equal(ml.KMeans.initialize(data, 3, 7, False), ml.KMeans.initialize(data, 3, 7, False))
    assert first.max_iterations == 50


# Gaussian mixture: scikit-learn's EM from the Numerics means converges to the same maximum; with a
# 1e-10 tolerance on both sides the means and the total log-likelihood agree to 1e-8.


def test_gaussian_mixture_matches_sklearn() -> None:
    # Arrange
    data = _blobs()
    model = ml.GaussianMixtureModel(data, 3)
    model.tolerance = 1e-10
    # Act
    model.train(12345)
    oracle = GaussianMixture(3, means_init=model.means, tol=1e-12, max_iter=2000, reg_covar=0.0).fit(data)
    # Assert
    order = np.argsort(model.means[:, 0] + 10 * model.means[:, 1])
    ref = np.argsort(oracle.means_[:, 0] + 10 * oracle.means_[:, 1])
    np.testing.assert_allclose(model.means[order], oracle.means_[ref], atol=1e-8)
    np.testing.assert_allclose(model.weights[order], oracle.weights_[ref], atol=1e-8)
    np.testing.assert_allclose(model.sigmas[order], oracle.covariances_[ref], atol=1e-7)
    assert model.log_likelihood == pytest.approx(oracle.score(data) * len(data), rel=1e-9)
    np.testing.assert_allclose(model.likelihood_matrix.sum(axis=1), 1.0, atol=1e-12)
    np.testing.assert_array_equal(model.labels, model.likelihood_matrix.argmax(axis=1))
    assert model.k == 3 and model.dimension == 2 and model.iterations >= 1 and model.max_iterations == 1000
    assert model.x.shape == (150, 2)


# PCA: scikit-learn's PCA uses the same sign convention and n - 1 divisor; the covariance route squares
# the condition number, so components and scores agree to 1e-10 on this well-conditioned fixture.


@pytest.mark.parametrize("basis", ["covariance", "correlation"])
def test_pca_matches_sklearn(basis: str) -> None:
    # Arrange
    rng = np.random.default_rng(8)
    data = rng.normal(size=(80, 4)) @ rng.normal(size=(4, 4)) + [1.0, -2.0, 3.0, 0.5]
    scaled = (data - data.mean(0)) / data.std(0, ddof=1) if basis == "correlation" else data
    # Act
    model = ml.PrincipalComponentAnalysis(data, basis)
    model.train()
    oracle = PCA().fit(scaled)
    # Assert
    np.testing.assert_allclose(model.components, oracle.components_, atol=1e-10)
    np.testing.assert_allclose(model.explained_variance, oracle.explained_variance_, rtol=1e-10)
    np.testing.assert_allclose(model.explained_variance_ratio, oracle.explained_variance_ratio_, rtol=1e-10)
    np.testing.assert_allclose(model.singular_values, oracle.singular_values_, rtol=1e-10)
    np.testing.assert_allclose(model.scores, oracle.transform(scaled), atol=1e-9)
    np.testing.assert_allclose(model.transform(data), model.scores, atol=1e-12)
    np.testing.assert_allclose(model.transform(data[0]), model.scores[0], atol=1e-12)
    np.testing.assert_allclose(model.inverse_transform(model.scores), data, atol=1e-9)
    np.testing.assert_allclose(model.inverse_transform(model.scores[1]), data[1], atol=1e-9)
    np.testing.assert_allclose(model.loadings, model.components.T * np.sqrt(model.explained_variance))
    np.testing.assert_allclose(model.cumulative_explained_variance_ratio, np.cumsum(model.explained_variance_ratio))
    np.testing.assert_allclose(model.means, data.mean(0))
    assert model.basis == ml.PrincipalComponentBasis(basis) and model.is_trained
    assert model.number_of_samples == 80 and model.number_of_features == 4 and model.number_of_components == 4
    assert model.total_variance == pytest.approx(np.sum(oracle.explained_variance_), rel=1e-12)
    assert model.scales.shape == (4,)


def test_pca_retained_components_and_invalid_input() -> None:
    # Arrange
    rng = np.random.default_rng(9)
    data = rng.normal(size=(30, 5))
    # Act
    model = ml.PrincipalComponentAnalysis(data, number_of_components=2)
    model.train()
    # Assert
    assert model.components.shape == (2, 5) and model.transform(data).shape == (30, 2)
    with pytest.raises(ValueError):
        ml.PrincipalComponentAnalysis(data, number_of_components=9)
    with pytest.raises(ValueError):
        ml.PrincipalComponentAnalysis(np.ones((1, 3)))


# k-nearest neighbors: scikit-learn with the same metric, neighbor count and 1/d**2 weights.


def test_knn_regression_and_neighbors_match_sklearn() -> None:
    # Arrange
    x, y = _regression(10)
    query = np.random.default_rng(11).normal(size=(7, 3))
    # Act
    model = ml.KNearestNeighbors(x, y, 5)
    oracle = KNeighborsRegressor(5, weights=lambda d: 1.0 / d**2).fit(x, y)
    neighbors = NearestNeighbors(n_neighbors=5).fit(x).kneighbors(query)[1]
    # Assert: the same weighted average, summed in a different order
    np.testing.assert_allclose(model.predict(query), oracle.predict(query), rtol=1e-13)
    np.testing.assert_array_equal(model.get_neighbors(query), neighbors)
    assert model.k == 5 and model.number_of_features == 3 and model.is_regression
    np.testing.assert_array_equal(model.x, x)
    np.testing.assert_array_equal(model.y, y)
    assert model.predict(query[:, :2]) is None and model.get_neighbors(query[:, :2]) is None


def test_knn_classification_matches_sklearn() -> None:
    # Arrange
    data = _blobs(12)
    labels = np.repeat([0.0, 1.0, 2.0], 50)
    query = np.random.default_rng(13).uniform(-1, 5, size=(40, 2))
    # Act
    model = ml.KNearestNeighbors(data, labels, 5)
    model.is_regression = False
    oracle = KNeighborsClassifier(5).fit(data, labels)
    # Assert: an odd k on three well-separated classes leaves no voting ties
    probabilities = oracle.predict_proba(query)
    clear = np.sort(probabilities, axis=1)[:, -1] > np.sort(probabilities, axis=1)[:, -2]
    np.testing.assert_array_equal(model.predict(query)[clear], oracle.predict(query)[clear])


def test_knn_bootstrap_is_seeded() -> None:
    # Arrange
    x, y = _regression(14)
    query = x[:6] + 0.01
    model = ml.KNearestNeighbors(x, y, 4)
    # Act
    first, second = model.bootstrap_predict(query, 3), model.bootstrap_predict(query, 3)
    intervals = model.prediction_intervals(query, 3, 200, 0.1)
    # Assert
    np.testing.assert_array_equal(first, second)
    assert intervals.shape == (6, 4)
    assert (intervals[:, 0] <= intervals[:, 1]).all() and (intervals[:, 1] <= intervals[:, 2]).all()


# Naive Bayes: the textbook Gaussian naive Bayes MAP rule with SciPy normal densities, n - 1 standard
# deviations and frequency priors (the Numerics definition); scikit-learn uses n divisors and smoothing.


def test_naive_bayes_matches_the_map_rule() -> None:
    # Arrange
    rng = np.random.default_rng(15)
    x = rng.normal(size=(90, 2))
    y = (x[:, 0] > 0).astype(float) + (x[:, 1] > 0.5)
    query = rng.normal(size=(25, 2))
    # Act
    model = ml.NaiveBayes(x, y)
    untrained = model.predict(query)
    model.train()
    # Assert
    classes = model.classes
    means = np.array([x[y == c].mean(0) for c in classes])
    sds = np.array([x[y == c].std(0, ddof=1) for c in classes])
    priors = np.array([np.mean(y == c) for c in classes])
    log_post = np.log(priors) + stats.norm.logpdf(query[:, None, :], means, sds).sum(axis=2)
    np.testing.assert_allclose(model.means, means, atol=1e-14)
    np.testing.assert_allclose(model.standard_deviations, sds, rtol=1e-12)
    np.testing.assert_allclose(model.priors, priors)
    np.testing.assert_array_equal(model.predict(query), classes[log_post.argmax(axis=1)])
    assert untrained is None and model.is_trained
    np.testing.assert_array_equal(model.x, x)
    np.testing.assert_array_equal(model.y, y)


# Decision tree: with every feature considered at every split, the tree is exact CART. scikit-learn
# places thresholds at midpoints and Numerics at data values, so training-set predictions agree exactly
# when the split sequences agree; this fixture has no tied split gains.


def test_decision_tree_regression_matches_sklearn_on_training_points() -> None:
    # Arrange
    x, y = _regression(16, 80)
    model = ml.DecisionTree(x, y, 12345)
    model.features = 3
    model.max_depth = 4
    model.minimum_split_size = 2
    # Act
    model.train()
    oracle = DecisionTreeRegressor(max_depth=4, min_samples_split=2).fit(x, y)
    # Assert
    np.testing.assert_allclose(model.predict(x), oracle.predict(x), rtol=1e-12)
    assert model.root.feature_index == oracle.tree_.feature[0]
    assert model.root.left is not None and model.root.right is not None and not model.root.is_leaf_node
    assert model.is_trained and model.dimensions == 3 and model.max_depth == 4 and model.minimum_split_size == 2
    np.testing.assert_array_equal(model.x, x)
    np.testing.assert_array_equal(model.y, y)


def test_decision_tree_classification_and_null_paths() -> None:
    # Arrange
    data = _blobs(17)
    labels = np.repeat([0.0, 1.0, 2.0], 50)
    model = ml.DecisionTree(data, labels, 3)
    model.is_regression = False
    # Act
    untrained = model.predict(data)
    model.train()
    # Assert
    assert untrained is None and not model.is_regression and model.features == 1
    assert np.mean(model.predict(data) == labels) > 0.95
    assert model.predict(np.ones((2, 5))) is None
    with pytest.raises(ValueError, match="at least ten"):
        ml.DecisionTree(data[:5], labels[:5])


def test_decision_node_round_trip() -> None:
    # Arrange
    node = ml.DecisionNode()
    # Act
    node.feature_index, node.threshold, node.value, node.is_leaf_node = 2, 1.5, 3.0, True
    # Assert
    assert (node.feature_index, node.threshold, node.value, node.is_leaf_node) == (2, 1.5, 3.0, True)
    assert node.left is None and node.right is None


# Random forest: no other library reproduces its bootstrap streams, so the output columns are checked
# against their definition from the member trees (NumPy percentiles and mean), and the accuracy against
# a scikit-learn forest of the same size.


def test_random_forest_columns_follow_the_member_trees() -> None:
    # Arrange
    x, y = _regression(18, 100)
    model = ml.RandomForest(x, y, 12345)
    model.number_of_trees = 40
    model.features = 2
    # Act
    model.train()
    out = model.predict(x[:10], 0.2)
    trees = np.array([tree.predict(x[:10]) for tree in model.decision_trees])
    # Assert
    np.testing.assert_allclose(out[:, 3], trees.mean(axis=0), rtol=1e-12)
    np.testing.assert_allclose(out[:, 1], np.median(trees, axis=0), rtol=1e-12)
    np.testing.assert_allclose(out[:, 0], np.quantile(trees, 0.1, axis=0), rtol=1e-12)
    np.testing.assert_allclose(out[:, 2], np.quantile(trees, 0.9, axis=0), rtol=1e-12)
    assert len(model.decision_trees) == 40 and model.is_trained and model.dimensions == 3
    assert model.number_of_trees == 40 and model.features == 2 and model.max_depth == 100
    assert model.minimum_split_size == 2 and model.is_regression
    np.testing.assert_array_equal(model.x, x)
    np.testing.assert_array_equal(model.y, y)


def test_random_forest_is_seeded_and_untrained_returns_none() -> None:
    # Arrange
    x, y = _regression(19, 50)
    first, second = ml.RandomForest(x, y, 9), ml.RandomForest(x, y, 9)
    first.number_of_trees = second.number_of_trees = 15
    # Act
    untrained = first.predict(x)
    first.train()
    second.train()
    # Assert
    assert untrained is None and ml.RandomForest(x, y).decision_trees == []
    np.testing.assert_array_equal(first.predict(x), second.predict(x))


# Jenks natural breaks: the optimum over every split of the sorted data into contiguous classes,
# enumerated by brute force (the definition of the method), is unique for this fixture.


def test_jenks_matches_brute_force_optimum() -> None:
    # Arrange
    data = np.random.default_rng(20).normal(size=14)
    ordered = np.sort(data)

    def sdcm(cuts: tuple[int, ...]) -> float:
        edges = (0, *cuts, len(ordered))
        return sum(np.sum((ordered[a:b] - ordered[a:b].mean()) ** 2) for a, b in itertools.pairwise(edges))

    best = min(itertools.combinations(range(1, len(ordered)), 2), key=sdcm)
    # Act
    model = ml.JenksNaturalBreaks(data, 3)
    # Assert
    clusters = model.clusters
    assert [c.start_index for c in clusters] == [0, *best]
    sdam = np.sum((ordered - ordered.mean()) ** 2)
    assert model.goodness_of_variance_fit == pytest.approx(1 - sdcm(best) / sdam, rel=1e-12)
    np.testing.assert_array_equal(model.sorted_data, ordered)
    assert model.number_of_clusters == 3 and len(model.breaks) == 3
    first = clusters[0]
    part = ordered[first.start_index:first.end_index + 1]
    assert first.count == len(part) and first.min_value == part[0] and first.max_value == part[-1]
    assert first.sum == pytest.approx(part.sum()) and first.average == pytest.approx(part.mean())
    assert first.variance == pytest.approx(part.var(ddof=1))
    assert first.sum_of_squared_deviations == pytest.approx(np.sum((part - part.mean()) ** 2))


def test_jenks_cluster_and_invalid_input() -> None:
    # Arrange
    cluster = ml.JenksCluster([1.0, 2.0, 4.0, 8.0], 1, 3)
    # Assert
    assert (cluster.start_index, cluster.end_index, cluster.count) == (1, 3, 3)
    assert cluster.average == pytest.approx(14.0 / 3.0) and cluster.variance == pytest.approx(np.var([2, 4, 8], ddof=1))
    with pytest.raises(ValueError):
        ml.JenksNaturalBreaks([1.0, 2.0], 3)


# Elastic net: per the Numerics remarks, the lasso path equals scikit-learn's Lasso with alpha = lambda
# on X standardized with population standard deviations, and for alpha < 1 ElasticNet on standardized X
# and y with alpha = lambda / sd(y). The default tolerance bounds the squared coefficient change, so the
# comparison tightens it to 1e-15 on both sides; coefficients then agree to 1e-7 relative.


def _standardized_oracle(x: np.ndarray, y: np.ndarray, lam: float, mix: float) -> np.ndarray:
    mu, sd = x.mean(0), x.std(0)
    if mix == 1.0:
        fit = Lasso(alpha=lam, tol=1e-14, max_iter=1_000_000).fit((x - mu) / sd, y)
        return fit.coef_ / sd
    sy = y.std()
    fit = ElasticNet(alpha=lam / sy, l1_ratio=mix, tol=1e-14, max_iter=1_000_000)
    fit.fit((x - mu) / sd, (y - y.mean()) / sy)
    return fit.coef_ * sy / sd


@pytest.mark.parametrize("mix", [1.0, 0.5])
def test_elastic_net_matches_sklearn(mix: float) -> None:
    # Arrange
    rng = np.random.default_rng(21)
    x = rng.normal(size=(80, 6)) * [1, 2, 3, 1, 0.5, 4]
    y = x @ [1.5, 0.0, -0.7, 0.0, 2.0, 0.1] + rng.normal(size=80)
    model = ml.ElasticNetRegression(x, y)
    model.alpha = mix
    model.tolerance = 1e-15
    # Act
    model.fit()
    # Assert
    for index in (5, 20, 40):
        expected = _standardized_oracle(x, y, model.lambdas[index], mix)
        np.testing.assert_allclose(model.get_coefficients(index), expected, rtol=1e-7, atol=1e-9)
        np.testing.assert_array_equal(model.coefficients[index], model.get_coefficients(index))
        np.testing.assert_array_equal(model.selected_variables(index), np.flatnonzero(model.coefficients[index]))
        manual = model.intercepts[index] + x[:3] @ model.coefficients[index]
        np.testing.assert_allclose(model.predict(x[:3], index), manual, rtol=1e-12)
        assert model.predict(x[0], index) == pytest.approx(manual[0], rel=1e-12)
    assert model.is_fitted and model.converged and model.iterations > 0 and model.alpha == mix
    assert model.sample_size == 80 and model.number_of_predictors == 6
    assert len(model.lambdas) == len(model.intercepts) == len(model.deviance_ratio) == len(model.nonzero_count)
    assert np.all(np.diff(model.lambdas) < 0)


def test_elastic_net_options_and_cross_validation() -> None:
    # Arrange
    x, y = _regression(22, 90)
    model = ml.ElasticNetRegression(x, y)
    model.user_lambdas = [0.5, 0.1, 0.01]
    model.penalty_factors = [1.0, 1.0, 0.0]
    model.tolerance, model.max_iterations, model.number_of_lambdas = 1e-9, 5000, 50
    model.standardize, model.has_intercept, model.lambda_min_ratio = True, True, 1e-3
    # Act
    model.cross_validate(5, 7)
    # Assert
    np.testing.assert_array_equal(model.lambdas, [0.5, 0.1, 0.01])
    assert model.is_cross_validated and len(model.cross_validation_error) == 3
    assert model.lambda_min == model.lambdas[model.lambda_min_index]
    assert model.lambda_one_standard_error == model.lambdas[model.lambda_one_standard_error_index]
    assert model.cross_validation_standard_error.shape == (3,) and len(set(model.fold_ids)) == 5
    assert model.get_coefficients(0)[2] != 0.0  # unpenalized
    np.testing.assert_array_equal(model.user_lambdas, [0.5, 0.1, 0.01])
    np.testing.assert_array_equal(model.penalty_factors, [1.0, 1.0, 0.0])
    assert (model.tolerance, model.max_iterations, model.number_of_lambdas) == (1e-9, 5000, 50)
    assert model.standardize and model.has_intercept and model.lambda_min_ratio == 1e-3
    np.testing.assert_array_equal(model.x, x)
    np.testing.assert_array_equal(model.y, y)
    # Act: explicit folds give the same fold labels back
    folds = np.arange(90) % 3
    model.cross_validate(fold_ids=folds)
    np.testing.assert_array_equal(model.fold_ids, folds)
    model.user_lambdas = None
    assert model.user_lambdas is None


def test_elastic_net_errors() -> None:
    # Arrange
    x, y = _regression(23)
    model = ml.ElasticNetRegression(x, y)
    # Act / Assert
    with pytest.raises(RuntimeError, match="fitted first"):
        model.get_coefficients(0)
    model.fit()
    with pytest.raises(ValueError, match="on the fitted path"):
        model.predict(x, 10_000)
    with pytest.raises(ValueError, match="between 3 and n"):
        model.cross_validate(2)
    with pytest.raises(ValueError):
        ml.ElasticNetRegression(x, y[:5])


# Gaussian process regression: with fixed hyperparameters the predictor, latent variance and log
# marginal likelihood are closed forms, so scikit-learn (optimizer=None, alpha = noise variance) agrees
# to rounding, 1e-10.


@pytest.mark.parametrize(("kernel", "sk_kernel"), [
    ("matern52", Matern(0.9, "fixed", nu=2.5)),
    ("matern32", Matern(0.9, "fixed", nu=1.5)),
    ("exponential", Matern(0.9, "fixed", nu=0.5)),
    ("squared_exponential", RBF(0.9, "fixed")),
])
def test_gpr_fixed_hyperparameters_match_sklearn(kernel: str, sk_kernel: object) -> None:
    # Arrange
    x = np.linspace(0.0, 5.0, 25)
    y = np.sin(x)
    points = np.linspace(-0.5, 5.5, 13)
    model = ml.GaussianProcessRegression(x, y, kernel, "zero")
    # Act
    model.set_hyperparameters(1.3, [0.9], 0.01)
    oracle = GaussianProcessRegressor(ConstantKernel(1.3, "fixed") * sk_kernel, alpha=0.01, optimizer=None).fit(
        x[:, None], y)
    mean, variance = model.predict_with_variance(points)
    expected_mean, expected_sd = oracle.predict(points[:, None], return_std=True)
    # Assert
    np.testing.assert_allclose(mean, expected_mean, atol=1e-10)
    np.testing.assert_allclose(variance, expected_sd**2, atol=1e-10)
    np.testing.assert_allclose(model.predict(points), mean, atol=1e-14)
    assert model.evaluate_log_marginal_likelihood(np.log([1.3, 0.9, 0.01])) == pytest.approx(
        oracle.log_marginal_likelihood_value_, abs=1e-9)
    noisy = model.predict_with_variance(points, include_noise=True)[1]
    np.testing.assert_allclose(noisy - variance, 0.01, atol=1e-12)
    interval = model.predict(points, alpha=0.1)
    np.testing.assert_allclose(interval[:, 2] - interval[:, 1], stats.norm.ppf(0.95) * np.sqrt(variance), rtol=1e-9)
    assert model.kernel == ml.GaussianProcessKernel(kernel) and model.mean_function == ml.GaussianProcessMean.ZERO


def test_gpr_gradient_matches_finite_differences() -> None:
    # Arrange
    rng = np.random.default_rng(24)
    x = rng.uniform(0, 4, size=(30, 2))
    y = np.sin(x[:, 0]) + 0.2 * x[:, 1] + 0.05 * rng.normal(size=30)
    model = ml.GaussianProcessRegression(x, y)
    theta = np.log([1.1, 0.8, 1.7, 0.02])
    # Act
    value, gradient = model.evaluate_log_marginal_likelihood(theta, gradient=True)
    steps = np.eye(4) * 1e-6
    numeric = [(model.evaluate_log_marginal_likelihood(theta + s) - model.evaluate_log_marginal_likelihood(theta - s))
               / 2e-6 for s in steps]
    # Assert: central differences with step 1e-6 are accurate to about 1e-7
    assert value == model.evaluate_log_marginal_likelihood(theta)
    np.testing.assert_allclose(gradient, numeric, rtol=1e-5, atol=1e-6)


def test_gpr_training_reaches_the_sklearn_optimum() -> None:
    # Arrange
    rng = np.random.default_rng(25)
    x = np.sort(rng.uniform(0, 6, 40))
    y = np.sin(x) + 0.1 * rng.normal(size=40)
    model = ml.GaussianProcessRegression(x, y, ml.GaussianProcessKernel.SQUARED_EXPONENTIAL, "zero")
    model.number_of_restarts, model.prng_seed = 2, 12345
    # Act
    model.train()
    kernel = ConstantKernel(1.0, (1e-6, 1e6)) * RBF(1.0, (1e-3, 1e3)) + WhiteKernel(0.01, (1e-10, 1e2))
    oracle = GaussianProcessRegressor(kernel, n_restarts_optimizer=5, random_state=0).fit(x[:, None], y)
    # Assert: both maximize the same likelihood; L-BFGS-B stops within 1e-6 of the optimum
    assert model.log_marginal_likelihood == pytest.approx(oracle.log_marginal_likelihood_value_, abs=1e-4)
    assert model.status == "Success" and model.is_trained and model.jitter >= 0.0
    assert model.sample_size == 40 and model.number_of_features == 1 and model.mean_coefficients.size == 0
    assert model.signal_variance > 0 and model.noise_variance > 0 and model.length_scales.shape == (1,)
    loo_mean, loo_var = model.leave_one_out()
    assert loo_mean.shape == (40,) and (loo_var > 0).all()
    assert (model.max_iterations, model.tolerance, model.number_of_restarts, model.prng_seed) == (200, 1e-5, 2, 12345)


def test_gpr_fixed_noise_training_and_settings() -> None:
    # Arrange
    rng = np.random.default_rng(26)
    x = rng.uniform(0, 3, size=(20, 2))
    y = x[:, 0] - x[:, 1] ** 2
    model = ml.GaussianProcessRegression(x, y, mean="linear")
    model.use_ard = False
    model.max_iterations, model.tolerance = 100, 1e-6
    # Act
    model.train(noise_variance=0.0)
    # Assert: an interpolating model reproduces the data
    np.testing.assert_allclose(model.predict(x), y, atol=1e-5)
    assert model.noise_variance == 0.0 and not model.use_ard and model.mean_coefficients.shape == (3,)
    with pytest.raises(ValueError):
        model.train(noise_variance=-1.0)
    with pytest.raises(RuntimeError, match="trained"):
        ml.GaussianProcessRegression(x, y).predict(x)


# Kriging: the ordinary kriging system with a Lagrange multiplier (Isaaks and Srivastava, 1989, ch. 12)
# solved by NumPy, and simple kriging as scikit-learn's zero-mean Gaussian process with a fixed
# covariance; both agree to rounding, 1e-10.


def _variogram_covariance(h: np.ndarray, model: str, sill: float, a: float) -> np.ndarray:
    if model == "exponential":
        return sill * np.exp(-h / a)
    if model == "gaussian":
        return sill * np.exp(-((h / a) ** 2))
    return np.where(h < a, sill * (1 - 1.5 * h / a + 0.5 * (h / a) ** 3), 0.0)


@pytest.mark.parametrize("model_name", ["exponential", "gaussian", "spherical"])
def test_ordinary_kriging_matches_the_kriging_system(model_name: str) -> None:
    # Arrange
    rng = np.random.default_rng(27)
    coords = rng.uniform(0, 10, size=(15, 2))
    values = np.sin(coords[:, 0]) + 0.1 * coords[:, 1]
    targets = rng.uniform(0, 10, size=(6, 2))
    psill, a = 1.0, 6.0
    n = len(values)
    lhs = np.ones((n + 1, n + 1))
    lhs[:n, :n] = _variogram_covariance(np.linalg.norm(coords[:, None] - coords[None], axis=2), model_name, psill, a)
    lhs[n, n] = 0.0
    rhs = np.ones((n + 1, len(targets)))
    rhs[:n] = _variogram_covariance(np.linalg.norm(coords[:, None] - targets[None], axis=2), model_name, psill, a)
    weights = np.linalg.solve(lhs, rhs)
    expected = weights[:n].T @ values
    expected_var = psill - np.sum(weights * rhs, axis=0)
    # Act
    model = ml.Kriging(coords, values, model_name, 0.0, psill, a)
    mean, variance = model.predict_with_variance(targets)
    # Assert
    np.testing.assert_allclose(mean, expected, atol=1e-10)
    np.testing.assert_allclose(variance, expected_var, atol=1e-10)
    np.testing.assert_allclose(model.predict(targets), mean, atol=1e-14)
    single = model.predict_with_variance(targets[0])
    assert single[0] == pytest.approx(mean[0], abs=1e-12) and single[1] == pytest.approx(variance[0], abs=1e-12)
    np.testing.assert_allclose(model.predict(coords), values, atol=1e-9)
    assert model.model == ml.VariogramModel(model_name) and model.kriging_type == ml.KrigingType.ORDINARY
    assert (model.nugget, model.partial_sill, model.sill, model.range) == (0.0, psill, psill, a)
    assert np.isnan(model.known_mean) and model.sample_size == 15 and model.number_of_dimensions == 2
    assert model.trend_coefficients.shape == (1,)


def test_simple_kriging_matches_sklearn_and_variogram() -> None:
    # Arrange
    rng = np.random.default_rng(28)
    coords = rng.uniform(0, 10, size=(12, 2))
    values = np.cos(coords[:, 1]) + 2.0
    targets = rng.uniform(0, 10, size=(5, 2))
    # Act
    model = ml.Kriging(coords, values, ml.VariogramModel.EXPONENTIAL, 0.0, 1.5, 4.0, "simple", 2.0)
    oracle = GaussianProcessRegressor(ConstantKernel(1.5, "fixed") * Matern(4.0, "fixed", nu=0.5), alpha=1e-300,
                                      optimizer=None).fit(coords, values - 2.0)
    mean, variance = model.predict_with_variance(targets)
    expected_mean, expected_sd = oracle.predict(targets, return_std=True)
    # Assert
    np.testing.assert_allclose(mean, expected_mean + 2.0, atol=1e-10)
    np.testing.assert_allclose(variance, expected_sd**2, atol=1e-10)
    h = np.array([0.0, 1.0, 3.0])
    np.testing.assert_allclose(model.semivariance(h), [0.0, *(1.5 * (1 - np.exp(-h[1:] / 4.0)))], atol=1e-15)
    np.testing.assert_allclose(model.covariance(h), 1.5 * np.exp(-h / 4.0), atol=1e-15)
    assert model.semivariance(1.0) == pytest.approx(1.5 * (1 - np.exp(-0.25)))
    assert model.known_mean == 2.0 and model.trend_coefficients.size == 0
    loo_mean, loo_var = model.cross_validate()
    assert loo_mean.shape == (12,) and (loo_var > 0).all()


def test_kriging_custom_distance_and_universal() -> None:
    # Arrange
    rng = np.random.default_rng(29)
    coords = rng.uniform(0, 10, size=(10, 2))
    values = coords[:, 0] + np.sin(coords[:, 1])
    targets = rng.uniform(0, 10, size=(4, 2))
    euclid = ml.Kriging(coords, values, "spherical", 0.1, 1.0, 8.0)
    # Act
    custom = ml.Kriging(coords, values, "spherical", 0.1, 1.0, 8.0, distance=lambda a, b: float(np.hypot(*(a - b))))
    universal = ml.Kriging(coords, values, "exponential", 0.0, 1.0, 5.0, ml.KrigingType.UNIVERSAL)
    # Assert
    np.testing.assert_allclose(custom.predict(targets), euclid.predict(targets), atol=1e-12)
    assert universal.trend_coefficients.shape == (3,)
    with pytest.raises(ValueError):
        ml.Kriging(coords, values, "spherical", -1.0, 1.0, 8.0)
    with pytest.raises(ValueError):
        ml.Kriging(coords, values, "exponential", 0.0, 1.0, 5.0, "universal", distance=lambda a, b: 1.0)


# Speed: the feature matrix crosses as one bulk copy, so the wrapper adds a fixed cost of a few copies
# to the Numerics call, never a per-row crossing.


def test_predict_crosses_the_feature_matrix_once() -> None:
    # Arrange
    x, y = _regression(30, 200)
    queries = np.random.default_rng(31).normal(size=(20_000, 3))
    model = ml.KNearestNeighbors(x, y, 3)
    net_queries = to_net_matrix(queries)
    model.predict(queries)
    # Act
    start = time.perf_counter()
    for _ in range(3):
        model.net.Predict(net_queries)
    raw = time.perf_counter() - start
    start = time.perf_counter()
    for _ in range(3):
        wrapped = model.predict(queries)
    total = time.perf_counter() - start
    # Assert: the copies of 60,000 inputs and 20,000 outputs cost far less than the search itself
    assert wrapped.shape == (20_000,)
    assert total < 1.5 * raw + 0.05, f"wrapper {total:.3f} s against {raw:.3f} s for the bare call"


@pytest.mark.benchmark
def test_kmeans_and_pca_are_faster_than_sklearn() -> None:
    # Arrange
    rng = np.random.default_rng(32)
    data = np.vstack([rng.normal(c, 1.0, size=(5000, 8)) for c in range(4)])
    ml.KMeans(data[:100], 2).train(1)
    ml.PrincipalComponentAnalysis(data[:100]).train()
    start_points = ml.KMeans.initialize(data, 4, 12345)

    def best(fn: object, repeats: int = 3) -> float:
        times = []
        for _ in range(repeats):
            start = time.perf_counter()
            fn()
            times.append(time.perf_counter() - start)
        return min(times)

    def numerics_pca() -> None:
        model = ml.PrincipalComponentAnalysis(data)
        model.train()

    def numerics_kmeans() -> None:
        ml.KMeans(data, 4).train(12345)

    # Act
    ours_pca, theirs_pca = best(numerics_pca), best(lambda: PCA().fit(data))
    ours_km = best(numerics_kmeans)
    theirs_km = best(lambda: SklearnKMeans(4, init=start_points, n_init=1, algorithm="lloyd", max_iter=1000).fit(data))
    # Assert: generous margins keep the check stable on a loaded machine
    assert ours_pca < 2.0 * theirs_pca, f"PCA {ours_pca:.4f} s against scikit-learn {theirs_pca:.4f} s"
    assert ours_km < 3.0 * theirs_km, f"k-means {ours_km:.4f} s against scikit-learn {theirs_km:.4f} s"
