"""Supervised and unsupervised models from ``Numerics.MachineLearning``.

Each class wraps one Numerics model and mirrors its members: the constructor takes the training
data, ``train`` (or ``fit``) runs the Numerics estimator, and ``predict`` and the other queries call
the Numerics members. Feature matrices cross to .NET once per call as a ``double[,]`` (rows are
observations, columns are features) and results come back as NumPy arrays. A one-dimensional
feature array is a single feature column, as in the Numerics ``double[]`` overloads, except where a
member documents it as one observation. Every estimator, default, seed and tolerance is the Numerics
one.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import cached_overload, loop, namespace, net_enum
from pynumerics.enums import NetEnum
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _numerics_matrix,
    _numerics_vector,
    _optional_doubles,
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_ints,
    to_numpy_matrix,
)


class GaussianProcessKernel(NetEnum):
    """Covariance functions of Gaussian process regression (``Numerics.MachineLearning.GaussianProcessKernel``)."""

    SQUARED_EXPONENTIAL = "SquaredExponential"
    MATERN32 = "Matern32"
    MATERN52 = "Matern52"
    EXPONENTIAL = "Exponential"


class GaussianProcessMean(NetEnum):
    """Mean functions of Gaussian process regression (``Numerics.MachineLearning.GaussianProcessMean``)."""

    ZERO = "Zero"
    CONSTANT = "Constant"
    LINEAR = "Linear"


class VariogramModel(NetEnum):
    """Isotropic variogram models for kriging (``Numerics.MachineLearning.VariogramModel``)."""

    SPHERICAL = "Spherical"
    EXPONENTIAL = "Exponential"
    GAUSSIAN = "Gaussian"


class KrigingType(NetEnum):
    """Kriging variants (``Numerics.MachineLearning.KrigingType``)."""

    SIMPLE = "Simple"
    ORDINARY = "Ordinary"
    UNIVERSAL = "Universal"


class PrincipalComponentBasis(NetEnum):
    """Bases of principal component analysis (``Numerics.MachineLearning.PrincipalComponentBasis``)."""

    COVARIANCE = "Covariance"
    CORRELATION = "Correlation"


def _net_type(name: str) -> Any:
    """Return a ``Numerics.MachineLearning`` type by name."""
    return getattr(namespace("Numerics.MachineLearning"), name)


def _columns(x: Any) -> np.ndarray:
    """Return features as a float64 matrix; a 1-D array becomes one feature column."""
    array = np.asarray(x, dtype=np.float64)
    if array.ndim == 1:
        array = array.reshape(-1, 1)
    if array.ndim != 2:
        raise ValueError(f"Expected a 1-D or 2-D feature array; got shape {array.shape}.")
    return array


def _net_columns(x: Any) -> Any:
    """Copy features into a ``System.Double[,]`` (see :func:`_columns`)."""
    return to_net_matrix(_columns(x))


def _doubles(value: Any) -> np.ndarray | None:
    """Copy a .NET ``double[]`` or ``Vector`` to NumPy; a null reference gives None."""
    if value is None:
        return None
    return to_numpy_doubles(value.Array if hasattr(value, "Header") else value)


def _matrix(value: Any) -> np.ndarray | None:
    """Copy a .NET ``double[,]`` or ``Matrix`` to NumPy; a null reference gives None."""
    if value is None:
        return None
    return to_numpy_matrix(value.Array if hasattr(value, "NumberOfRows") else value)


def _ints(value: Any) -> np.ndarray | None:
    """Copy a .NET ``int[]`` to NumPy; a null reference gives None."""
    return None if value is None else to_numpy_ints(value)


def _pair(result: Any) -> tuple[np.ndarray, np.ndarray]:
    """Copy a .NET ``(double[] Mean, double[] Variance)`` tuple to two NumPy arrays."""
    return to_numpy_doubles(result.Item1), to_numpy_doubles(result.Item2)


def _member(name: str, convert: Callable[[Any], Any] | None, doc: str,
            to_net: Callable[[Any], Any] | None = None) -> property:
    """Return a property reading (and, with ``to_net``, writing) the .NET member ``name`` of ``self.net``."""
    def fget(self: Any) -> Any:
        value = getattr(self.net, name)
        return value if convert is None else convert(value)

    if to_net is None:
        return property(fget, doc=doc)

    def fset(self: Any, value: Any) -> None:
        setattr(self.net, name, to_net(value))

    return property(fget, fset, doc=doc)


class _NetModel:
    """Base of the wrappers: holds the .NET object in ``net`` and wraps an existing one."""

    net: Any

    @classmethod
    def _wrap(cls, net: Any) -> Any:
        """Return a wrapper around an existing .NET object."""
        out = cls.__new__(cls)
        out.net = net
        return out


class DecisionNode(_NetModel):
    """A node of a decision tree (``Numerics.MachineLearning.DecisionNode``).

    Attributes
    ----------
    net : Numerics.MachineLearning.DecisionNode
        The wrapped .NET object.
    """

    def __init__(self) -> None:
        self.net = _net_type("DecisionNode")()

    feature_index = _member("FeatureIndex", int, "Index of the feature the node splits on; -1 for a leaf.", int)
    threshold = _member("Threshold", float, "Split threshold; NaN for a leaf.", float)
    value = _member("Value", float, "Leaf prediction; NaN for an internal node.", float)
    is_leaf_node = _member("IsLeafNode", bool, "Whether the node is a leaf.", bool)

    @property
    def left(self) -> DecisionNode | None:
        """Left child, or None."""
        node = self.net.Left
        return None if node is None else DecisionNode._wrap(node)

    @property
    def right(self) -> DecisionNode | None:
        """Right child, or None."""
        node = self.net.Right
        return None if node is None else DecisionNode._wrap(node)


class DecisionTree(_NetModel):
    """Classification and regression tree (``Numerics.MachineLearning.DecisionTree``).

    Parameters
    ----------
    x : array_like
        Training features, shape (n, d); a 1-D array is one feature.
    y : array_like
        Training responses (class labels as numbers for classification), length n, at least 10.
    seed : int, optional
        Mersenne Twister seed for the feature subsets; zero or negative seeds from the clock.

    Attributes
    ----------
    net : Numerics.MachineLearning.DecisionTree
        The wrapped .NET object.

    Notes
    -----
    Set :attr:`is_regression`, :attr:`max_depth`, :attr:`minimum_split_size` and :attr:`features`
    before :meth:`train`. Each split considers :attr:`features` randomly chosen features (default
    ``max(1, d - 1)``).
    """

    def __init__(self, x: Any, y: Any, seed: int = -1) -> None:
        self.net = _net_type("DecisionTree")(_net_columns(x), to_net_doubles(y), int(seed))

    def train(self) -> None:
        """Grow the tree (``Train``)."""
        self.net.Train()

    def predict(self, x: Any) -> np.ndarray | None:
        """Return the tree prediction for each row (``Predict``).

        Parameters
        ----------
        x : array_like
            Features, shape (m, d); a 1-D array is one feature.

        Returns
        -------
        numpy.ndarray or None
            The m predictions, or None when the tree is untrained or ``x`` has a different number of
            columns (the Numerics member returns null).
        """
        return _doubles(self.net.Predict(_net_columns(x)))

    minimum_split_size = _member("MinimumSplitSize", int, "Minimum number of rows to split a node (default 2).", int)
    max_depth = _member("MaxDepth", int, "Maximum depth of the tree (default 100).", int)
    features = _member("Features", int, "Number of features considered at each split.", int)
    is_regression = _member("IsRegression", bool, "Regression (True, default) or classification.", bool)
    dimensions = _member("Dimensions", int, "Number of feature columns d.")
    is_trained = _member("IsTrained", bool, "Whether the tree has been trained.")
    x = _member("X", _matrix, "Training features, shape (n, d).")
    y = _member("Y", _doubles, "Training responses.")

    @property
    def root(self) -> DecisionNode:
        """Root node of the tree."""
        return DecisionNode._wrap(self.net.Root)


class RandomForest(_NetModel):
    """Random forest of decision trees on bootstrap samples (``Numerics.MachineLearning.RandomForest``).

    Parameters
    ----------
    x : array_like
        Training features, shape (n, d); a 1-D array is one feature.
    y : array_like
        Training responses, length n, at least 10.
    seed : int, optional
        Mersenne Twister seed of the forest; zero or negative seeds from the clock.

    Attributes
    ----------
    net : Numerics.MachineLearning.RandomForest
        The wrapped .NET object.

    Notes
    -----
    Trees are trained in parallel, each from a seed drawn from the forest generator, so a seeded
    forest is deterministic.
    """

    def __init__(self, x: Any, y: Any, seed: int = -1) -> None:
        self.net = _net_type("RandomForest")(_net_columns(x), to_net_doubles(y), int(seed))

    def train(self) -> None:
        """Train every tree (``Train``)."""
        self.net.Train()

    def predict(self, x: Any, alpha: float = 0.1) -> np.ndarray | None:
        """Return the forest prediction and its interval across trees (``Predict``).

        Parameters
        ----------
        x : array_like
            Features, shape (m, d); a 1-D array is one feature.
        alpha : float, optional
            One minus the interval coverage; default 0.1 (a 90% interval).

        Returns
        -------
        numpy.ndarray or None
            Shape (m, 4) with columns lower, median, upper and mean of the tree predictions (floored
            for classification), or None when the forest is untrained or ``x`` has a different number
            of columns.
        """
        return _matrix(self.net.Predict(_net_columns(x), float(alpha)))

    number_of_trees = _member("NumberOfTrees", int, "Number of trees (default 1000).", int)
    minimum_split_size = _member("MinimumSplitSize", int, "Minimum number of rows to split a node (default 2).", int)
    max_depth = _member("MaxDepth", int, "Maximum depth of each tree (default 100).", int)
    features = _member("Features", int, "Number of features considered at each split.", int)
    is_regression = _member("IsRegression", bool, "Regression (True, default) or classification.", bool)
    dimensions = _member("Dimensions", int, "Number of feature columns d.")
    is_trained = _member("IsTrained", bool, "Whether the forest has been trained.")
    x = _member("X", _matrix, "Training features, shape (n, d).")
    y = _member("Y", _doubles, "Training responses.")

    @property
    def decision_trees(self) -> list[DecisionTree]:
        """The trained trees; empty before training."""
        trees = self.net.DecisionTrees
        return [] if trees is None else [DecisionTree._wrap(tree) for tree in trees]


class KNearestNeighbors(_NetModel):
    """The k-nearest neighbors regression and classification (``Numerics.MachineLearning.KNearestNeighbors``).

    Parameters
    ----------
    x : array_like
        Training features, shape (n, d); a 1-D array is one feature.
    y : array_like
        Training responses, length n.
    k : int
        Number of neighbors.

    Attributes
    ----------
    net : Numerics.MachineLearning.KNearestNeighbors
        The wrapped .NET object.

    Notes
    -----
    Distances are Euclidean, ties broken by the lower training index. Regression averages the
    neighbors with inverse squared distance weights (exact matches take all the weight);
    classification returns the most common neighbor class.
    """

    def __init__(self, x: Any, y: Any, k: int) -> None:
        self.net = _net_type("KNearestNeighbors")(_net_columns(x), to_net_doubles(y), int(k))

    def get_neighbors(self, x: Any) -> np.ndarray | None:
        """Return the training indexes of the k nearest neighbors of each row (``GetNeighbors``).

        Parameters
        ----------
        x : array_like
            Query features, shape (m, d); a 1-D array is one feature.

        Returns
        -------
        numpy.ndarray or None
            Shape (m, k) int32, nearest first, or None when ``x`` has a different number of columns.
        """
        rows = _columns(x)
        found = _ints(self.net.GetNeighbors(to_net_matrix(rows)))
        return None if found is None else found.reshape(rows.shape[0], -1)

    def predict(self, x: Any) -> np.ndarray | None:
        """Return the k-nearest-neighbor prediction of each row (``Predict``).

        Parameters
        ----------
        x : array_like
            Query features, shape (m, d); a 1-D array is one feature.

        Returns
        -------
        numpy.ndarray or None
            The m predictions, or None when ``x`` has a different number of columns.
        """
        return _doubles(self.net.Predict(_net_columns(x)))

    def bootstrap_predict(self, x: Any, seed: int = -1) -> np.ndarray | None:
        """Return the prediction from one bootstrap resample of the training data (``BootstrapPredict``).

        Parameters
        ----------
        x : array_like
            Query features, shape (m, d); a 1-D array is one feature.
        seed : int, optional
            Mersenne Twister seed of the resample; zero or negative seeds from the clock.

        Returns
        -------
        numpy.ndarray or None
            The m predictions, or None when ``x`` has a different number of columns.
        """
        return _doubles(self.net.BootstrapPredict(_net_columns(x), int(seed)))

    def prediction_intervals(self, x: Any, seed: int = -1, realizations: int = 1000,
                             alpha: float = 0.1) -> np.ndarray | None:
        """Return bootstrap prediction intervals (``PredictionIntervals``).

        Parameters
        ----------
        x : array_like
            Query features, shape (m, d); a 1-D array is one feature.
        seed : int, optional
            Mersenne Twister seed; zero or negative seeds from the clock.
        realizations : int, optional
            Number of bootstrap resamples; default 1000.
        alpha : float, optional
            One minus the interval coverage; default 0.1.

        Returns
        -------
        numpy.ndarray or None
            Shape (m, 4) with columns lower, median, upper and mean.
        """
        return _matrix(self.net.PredictionIntervals(_net_columns(x), int(seed), int(realizations), float(alpha)))

    k = _member("K", int, "Number of neighbors.")
    is_regression = _member("IsRegression", bool, "Regression (True, default) or classification.", bool)
    number_of_features = _member("NumberOfFeatures", int, "Number of feature columns d.")
    x = _member("X", _matrix, "Training features, shape (n, d).")
    y = _member("Y", _doubles, "Training responses.")


class NaiveBayes(_NetModel):
    """Gaussian naive Bayes classifier (``Numerics.MachineLearning.NaiveBayes``).

    Parameters
    ----------
    x : array_like
        Training features, shape (n, d); a 1-D array is one feature.
    y : array_like
        Class labels as numbers, length n.

    Attributes
    ----------
    net : Numerics.MachineLearning.NaiveBayes
        The wrapped .NET object.

    Notes
    -----
    Class-conditional standard deviations use the n - 1 divisor and are floored at
    ``sqrt(1e-9 * max variance)``; priors are the class frequencies.
    """

    def __init__(self, x: Any, y: Any) -> None:
        self.net = _net_type("NaiveBayes")(_net_columns(x), to_net_doubles(y))

    def train(self) -> None:
        """Estimate the priors and class-conditional normals (``Train``)."""
        self.net.Train()

    def predict(self, x: Any) -> np.ndarray | None:
        """Return the maximum a posteriori class of each row (``Predict``).

        Parameters
        ----------
        x : array_like
            Features, shape (m, d); a 1-D array is one feature.

        Returns
        -------
        numpy.ndarray or None
            The m predicted classes, or None when untrained or ``x`` has a different number of columns.
        """
        return _doubles(self.net.Predict(_net_columns(x)))

    classes = _member("Classes", _doubles, "The distinct class labels.")
    means = _member("Means", _matrix, "Class-conditional feature means, shape (classes, d).")
    standard_deviations = _member("StandardDeviations", _matrix,
                                  "Class-conditional feature standard deviations, shape (classes, d).")
    priors = _member("Priors", _doubles, "Class prior probabilities.")
    is_trained = _member("IsTrained", bool, "Whether the classifier has been trained.")
    x = _member("X", _matrix, "Training features, shape (n, d).")
    y = _member("Y", _doubles, "Training class labels.")


def _glm_predict(arity: int) -> Any:
    """Return ``GeneralizedLinearModel.Predict(Matrix)`` or ``Predict(Matrix, double)``, bound once."""
    def factory() -> Any:
        import System  # noqa: PLC0415 - the runtime is running

        matrix = namespace("Numerics.Mathematics.LinearAlgebra").Matrix
        signature = (matrix,) if arity == 1 else (matrix, System.Double)
        return _net_type("GeneralizedLinearModel").Predict.Overloads[signature]

    return cached_overload(f"GeneralizedLinearModel.Predict/{arity}", factory)


class GeneralizedLinearModel(_NetModel):
    """Generalized linear model fitted by maximum likelihood (``Numerics.MachineLearning.GeneralizedLinearModel``).

    Parameters
    ----------
    x : array_like
        Predictors, shape (n, p); a 1-D array is one predictor.
    y : array_like
        Responses, length n (0/1 for the binary links, counts for the log link).
    has_intercept : bool, optional
        Whether an intercept column is added; default True.
    link : str, optional
        ``LinkFunctionType`` member name: ``"identity"`` (normal, default), ``"log"`` (Poisson),
        ``"logit"``, ``"probit"`` or ``"complementary_log_log"`` (Bernoulli). With ``link_function``
        it selects the family and the optimizer start.
    link_function : Numerics.Functions.ILinkFunction, optional
        A custom .NET link function instance.
    feature_names : sequence of str, optional
        Predictor names for :meth:`summary` (``Matrix.Header``). With an added intercept column the
        coefficient names are ``"Intercept"`` followed by these names; without names they are
        ``β1 ... βp``.
    response_name : str, optional
        Response name for :meth:`summary` (``Vector.Header``).

    Attributes
    ----------
    net : Numerics.MachineLearning.GeneralizedLinearModel
        The wrapped .NET object.
    """

    def __init__(self, x: Any, y: Any, has_intercept: bool = True, link: str = "identity",
                 link_function: Any = None, feature_names: Any = None, response_name: str | None = None) -> None:
        functions = namespace("Numerics.Functions")
        net_type = _net_type("GeneralizedLinearModel")
        matrix, vector = _numerics_matrix(_columns(x)), _numerics_vector(y)
        if feature_names is not None:
            import System  # noqa: PLC0415 - the runtime is running

            matrix.Header = System.Array[System.String]([str(name) for name in feature_names])
        if response_name is not None:
            vector.Header = str(response_name)
        link_type = net_enum(functions.LinkFunctionType, link)
        if link_function is None:
            self.net = net_type(matrix, vector, bool(has_intercept), link_type)
        else:
            self.net = net_type(matrix, vector, link_function, bool(has_intercept), link_type)

    def set_optimizer(self, method: str = "NelderMead") -> None:
        """Choose the local optimizer used by :meth:`train` (``SetOptimizer``).

        Parameters
        ----------
        method : str, optional
            ``LocalMethod`` member name; default ``"NelderMead"``.
        """
        optimization = namespace("Numerics.Mathematics.Optimization")
        self.net.SetOptimizer(net_enum(optimization.LocalMethod, method))

    def train(self) -> None:
        """Fit the parameters and compute the diagnostics (``Train``).

        The optimizer result is kept unless the Fisher scoring estimate has a higher log-likelihood.
        """
        self.net.Train()

    def predict(self, x: Any, alpha: float | None = None) -> np.ndarray:
        """Return the mean response, optionally with a confidence interval (``Predict``).

        Parameters
        ----------
        x : array_like
            Predictors, shape (m, p), or (m, p + 1) with the intercept column already included.
        alpha : float, optional
            When given, one minus the confidence level of the interval for the mean response.

        Returns
        -------
        numpy.ndarray
            The m mean responses, or with ``alpha`` shape (m, 3) with columns lower, mean and upper
            (formed on the link scale and mapped through the inverse link).
        """
        matrix = _numerics_matrix(_columns(x))
        if alpha is None:
            return to_numpy_doubles(_glm_predict(1)(self.net, matrix))
        return to_numpy_matrix(_glm_predict(2)(self.net, matrix, float(alpha)))

    def summary(self) -> list[str]:
        """Return the text summary table (``Summary``).

        Returns
        -------
        list of str
            The lines of the report.
        """
        return [str(line) for line in self.net.Summary()]

    has_intercept = _member("HasIntercept", bool, "Whether the model has an intercept.")
    use_robust_se = _member("UseRobustSE", bool, "Use HC0 sandwich standard errors (set before training).", bool)
    parameters = _member("Parameters", _doubles, "Fitted coefficients, intercept first; None before training.")
    parameter_standard_errors = _member("ParameterStandardErrors", _doubles, "Coefficient standard errors.")
    parameter_z_scores = _member("ParameterZScores", _doubles, "Coefficient z statistics.")
    parameter_p_values = _member("ParameterPValues", _doubles, "Two-sided normal p-values of the coefficients.")
    covariance = _member("Covariance", _matrix, "Coefficient covariance matrix.")
    residuals = _member("Residuals", _doubles, "Response residuals y - mu.")
    standard_error = _member("StandardError", float, "Residual standard error sqrt(SSE / (n - p)).")
    sample_size = _member("SampleSize", int, "Number of observations n.")
    degrees_of_freedom = _member("DegreesOfFreedom", int, "Residual degrees of freedom n - p.")
    aic = _member("AIC", float, "Akaike information criterion; the normal family evaluates the likelihood at "
                                "the maximum likelihood variance SSE / n, as statsmodels GLM does.")
    aicc = _member("AICc", float, "Small-sample corrected AIC.")
    bic = _member("BIC", float, "Bayesian information criterion.")
    x = _member("X", _matrix, "Design matrix, including the intercept column when present.")
    y = _member("Y", _doubles, "Responses.")
    link_type = _member("LinkType", str, "``LinkFunctionType`` member name of the model.")
    link_function = _member("LinkFunction", None, "The .NET ``ILinkFunction`` instance.")
    optimizer = _member("Optimizer", None, "The .NET ``Optimizer`` used by :meth:`train`.")

    @property
    def parameter_names(self) -> list[str]:
        """Coefficient names, ``"Intercept"`` first when present."""
        names = self.net.ParameterNames
        return [] if names is None else [str(name) for name in names]


class ElasticNetRegression(_NetModel):
    """Elastic net regression path by coordinate descent (``Numerics.MachineLearning.ElasticNetRegression``).

    Parameters
    ----------
    x : array_like
        Predictors, shape (n, p); a 1-D array is one predictor.
    y : array_like
        Responses, length n.

    Attributes
    ----------
    net : Numerics.MachineLearning.ElasticNetRegression
        The wrapped .NET object.

    Notes
    -----
    Follows glmnet: for each penalty lambda the method minimizes
    ``1/(2n) sum (y - b0 - x b)**2 + lambda sum v_j ((1 - alpha)/2 b_j**2 + alpha |b_j|)`` with the
    predictors and the response standardized. Set the options before :meth:`fit`.
    """

    def __init__(self, x: Any, y: Any) -> None:
        self.net = _net_type("ElasticNetRegression")(_net_columns(x), to_net_doubles(y))

    def fit(self) -> None:
        """Fit the coefficient path (``Fit``)."""
        self.net.Fit()

    def cross_validate(self, folds: int = 10, seed: int = 12345, fold_ids: Any = None) -> None:
        """Run k-fold cross-validation over the path, fitting it first if needed (``CrossValidate``).

        Parameters
        ----------
        folds : int, optional
            Number of folds, at least 3; default 10. Ignored with ``fold_ids``.
        seed : int, optional
            Seed of the random fold assignment; default 12345. Ignored with ``fold_ids``.
        fold_ids : array_like of int, optional
            The fold label of each observation, at least three distinct labels.
        """
        if fold_ids is None:
            self.net.CrossValidate(int(folds), int(seed))
        else:
            self.net.CrossValidate(to_net_ints(fold_ids))

    def get_coefficients(self, index: int) -> np.ndarray:
        """Return the coefficients, without the intercept, at a penalty on the path (``GetCoefficients``).

        Parameters
        ----------
        index : int
            Penalty index.

        Returns
        -------
        numpy.ndarray
            The p coefficients on the original scale.
        """
        return to_numpy_doubles(self.net.GetCoefficients(int(index)))

    def selected_variables(self, index: int) -> np.ndarray:
        """Return the indexes of the nonzero coefficients at a penalty (``SelectedVariables``).

        Parameters
        ----------
        index : int
            Penalty index.

        Returns
        -------
        numpy.ndarray
            Predictor indexes in increasing order, int32.
        """
        return to_numpy_ints(self.net.SelectedVariables(int(index)))

    def predict(self, x: Any, index: int) -> np.ndarray | float:
        """Return the prediction at a penalty on the path (``Predict``).

        Parameters
        ----------
        x : array_like
            Predictors, shape (m, p), or a 1-D array of the p values of one observation.
        index : int
            Penalty index.

        Returns
        -------
        numpy.ndarray or float
            The m predictions, or one prediction for a 1-D ``x``.
        """
        array = np.asarray(x, dtype=np.float64)
        if array.ndim == 1:
            return float(self.net.Predict(to_net_doubles(array), int(index)))
        return to_numpy_doubles(self.net.Predict(to_net_matrix(array), int(index)))

    alpha = _member("Alpha", float, "Mixing parameter in [0, 1]: 1 is the lasso (default), 0 ridge.", float)
    number_of_lambdas = _member("NumberOfLambdas", int, "Number of penalties on the default path (100).", int)
    lambda_min_ratio = _member("LambdaMinRatio", float,
                               "Smallest penalty as a fraction of lambda_max; NaN selects the glmnet default.", float)
    user_lambdas = _member("UserLambdas", _doubles, "User penalty sequence, or None for the default path.",
                           _optional_doubles)
    standardize = _member("Standardize", bool, "Whether the predictors are standardized (default True).", bool)
    has_intercept = _member("HasIntercept", bool, "Whether an intercept is fitted (default True).", bool)
    tolerance = _member("Tolerance", float, "Coordinate descent convergence tolerance (default 1e-7).", float)
    max_iterations = _member("MaxIterations", int, "Maximum number of coordinate passes (default 100000).", int)
    penalty_factors = _member("PenaltyFactors", _doubles, "Per-predictor penalty factors, or None for ones.",
                              _optional_doubles)
    sample_size = _member("SampleSize", int, "Number of observations n.")
    number_of_predictors = _member("NumberOfPredictors", int, "Number of predictors p.")
    is_fitted = _member("IsFitted", bool, "Whether the path has been fitted.")
    converged = _member("Converged", bool, "Whether every penalty converged.")
    iterations = _member("Iterations", int, "Total number of coordinate passes.")
    lambdas = _member("Lambdas", _doubles, "Penalties on the fitted path.")
    intercepts = _member("Intercepts", _doubles, "Intercept at each penalty.")
    coefficients = _member("Coefficients", _matrix, "Coefficients on the original scale, shape (penalties, p).")
    deviance_ratio = _member("DevianceRatio", _doubles, "Fraction of deviance explained at each penalty.")
    nonzero_count = _member("NonzeroCount", _ints, "Number of nonzero coefficients at each penalty.")
    is_cross_validated = _member("IsCrossValidated", bool, "Whether cross-validation has run.")
    fold_ids = _member("FoldIds", _ints, "Fold label of each observation used by cross-validation.")
    cross_validation_error = _member("CrossValidationError", _doubles, "Mean squared prediction error per penalty.")
    cross_validation_standard_error = _member("CrossValidationStandardError", _doubles,
                                              "Standard error of the cross-validation error per penalty.")
    lambda_min_index = _member("LambdaMinIndex", int, "Index of the penalty with the smallest error; -1 before.")
    lambda_one_standard_error_index = _member("LambdaOneStandardErrorIndex", int,
                                              "Index of the largest penalty within one standard error; -1 before.")
    lambda_min = _member("LambdaMin", float, "Penalty with the smallest cross-validation error.")
    lambda_one_standard_error = _member("LambdaOneStandardError", float,
                                        "Largest penalty within one standard error of the minimum.")
    x = _member("X", _matrix, "Predictors, shape (n, p).")
    y = _member("Y", _doubles, "Responses.")


def _gpr_predict_interval() -> Any:
    """Return ``GaussianProcessRegression.Predict(Matrix, double, bool)``, bound once."""
    def factory() -> Any:
        import System  # noqa: PLC0415 - the runtime is running

        matrix = namespace("Numerics.Mathematics.LinearAlgebra").Matrix
        signature = (matrix, System.Double, System.Boolean)
        return _net_type("GaussianProcessRegression").Predict.Overloads[signature]

    return cached_overload("GaussianProcessRegression.Predict(Matrix, double, bool)", factory)


class GaussianProcessRegression(_NetModel):
    """Exact Gaussian process regression (``Numerics.MachineLearning.GaussianProcessRegression``).

    Parameters
    ----------
    x : array_like
        Training features, shape (n, d); a 1-D array is one feature.
    y : array_like
        Training responses, length n.
    kernel : GaussianProcessKernel or str, optional
        Covariance function; default Matern 5/2.
    mean : GaussianProcessMean or str, optional
        Mean function, estimated by generalized least squares; default constant.

    Attributes
    ----------
    net : Numerics.MachineLearning.GaussianProcessRegression
        The wrapped .NET object.

    Notes
    -----
    The covariance is ``signal_variance * rho(r) + noise_variance * I`` with the scaled distance
    ``r`` (one length scale per feature when :attr:`use_ard`). :meth:`train` maximizes the log
    marginal likelihood with L-BFGS-B; :meth:`set_hyperparameters` conditions on fixed values.
    """

    def __init__(self, x: Any, y: Any, kernel: GaussianProcessKernel | str = GaussianProcessKernel.MATERN52,
                 mean: GaussianProcessMean | str = GaussianProcessMean.CONSTANT) -> None:
        net_type = _net_type("GaussianProcessRegression")
        ml = namespace("Numerics.MachineLearning")
        self.net = net_type(_net_columns(x), to_net_doubles(y), net_enum(ml.GaussianProcessKernel, kernel),
                            net_enum(ml.GaussianProcessMean, mean))

    def train(self, noise_variance: float | None = None) -> None:
        """Estimate the hyperparameters by maximum marginal likelihood (``Train``).

        Parameters
        ----------
        noise_variance : float, optional
            When given, the noise variance is held at this value (zero interpolates); otherwise it is
            estimated with the other hyperparameters.
        """
        if noise_variance is None:
            self.net.Train()
        else:
            self.net.Train(float(noise_variance))

    def set_hyperparameters(self, signal_variance: float, length_scales: Any, noise_variance: float) -> None:
        """Set the hyperparameters without optimization and condition on the data (``SetHyperparameters``).

        Parameters
        ----------
        signal_variance : float
            Signal variance, positive.
        length_scales : float or array_like
            One shared length scale, or one per feature (which turns :attr:`use_ard` on).
        noise_variance : float
            Noise variance, non-negative.
        """
        self.net.SetHyperparameters(float(signal_variance), to_net_doubles(length_scales), float(noise_variance))

    def evaluate_log_marginal_likelihood(self, theta: Any, gradient: bool = False) -> float | tuple[float, np.ndarray]:
        """Return the log marginal likelihood at log hyperparameters (``EvaluateLogMarginalLikelihood``).

        Parameters
        ----------
        theta : array_like
            ``(ln signal_variance, ln l_1, ..., ln l_q, ln noise_variance)`` with q = d under ARD, else 1.
        gradient : bool, optional
            Whether to also return the gradient.

        Returns
        -------
        float or tuple of (float, numpy.ndarray)
            The log marginal likelihood (NaN when the covariance cannot be factored), and with
            ``gradient`` its gradient with respect to ``theta``.
        """
        net_theta = to_net_doubles(theta)
        if not gradient:
            return float(self.net.EvaluateLogMarginalLikelihood(net_theta, None))
        import System  # noqa: PLC0415 - the runtime is running

        out = System.Array[System.Double](int(net_theta.Length))
        value = float(self.net.EvaluateLogMarginalLikelihood(net_theta, out))
        return value, to_numpy_doubles(out)

    def predict(self, x: Any, alpha: float | None = None, include_noise: bool = False) -> np.ndarray:
        """Return the predictive mean, optionally with a central interval (``Predict``).

        Parameters
        ----------
        x : array_like
            Prediction points, shape (m, d); a 1-D array is one feature.
        alpha : float, optional
            When given, one minus the interval coverage.
        include_noise : bool, optional
            With ``alpha``, whether the interval covers a new noisy observation; default False.

        Returns
        -------
        numpy.ndarray
            The m means, or with ``alpha`` shape (m, 3) with columns lower, mean and upper.
        """
        if alpha is None:
            return to_numpy_doubles(self.net.Predict(_net_columns(x)))
        matrix = _numerics_matrix(_columns(x))
        return to_numpy_matrix(_gpr_predict_interval()(self.net, matrix, float(alpha), bool(include_noise)))

    def predict_with_variance(self, x: Any, include_noise: bool = False) -> tuple[np.ndarray, np.ndarray]:
        """Return the predictive mean and variance (``PredictWithVariance``).

        Parameters
        ----------
        x : array_like
            Prediction points, shape (m, d); a 1-D array is one feature.
        include_noise : bool, optional
            Whether the variance is that of a new observation (True) or of the latent function.

        Returns
        -------
        tuple of numpy.ndarray
            The m means and the m variances.
        """
        return _pair(self.net.PredictWithVariance(_net_columns(x), bool(include_noise)))

    def leave_one_out(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the leave-one-out predictive mean and variance of each response (``LeaveOneOut``).

        Returns
        -------
        tuple of numpy.ndarray
            The n means and the n variances (the variances include the noise).
        """
        return _pair(self.net.LeaveOneOut())

    use_ard = _member("UseARD", bool, "One length scale per feature (default True).", bool)
    number_of_restarts = _member("NumberOfRestarts", int, "Random optimizer restarts beyond the two data starts.", int)
    prng_seed = _member("PRNGSeed", int, "Seed of the random restarts (default 12345).", int)
    max_iterations = _member("MaxIterations", int, "Maximum optimizer iterations per start (default 200).", int)
    tolerance = _member("Tolerance", float, "Projected gradient tolerance (default 1e-5).", float)
    sample_size = _member("SampleSize", int, "Number of observations n.")
    number_of_features = _member("NumberOfFeatures", int, "Number of features d.")
    signal_variance = _member("SignalVariance", float, "Signal variance.")
    length_scales = _member("LengthScales", _doubles, "Length scales.")
    noise_variance = _member("NoiseVariance", float, "Noise variance.")
    jitter = _member("Jitter", float, "Diagonal jitter added to factor the covariance.")
    mean_coefficients = _member("MeanCoefficients", _doubles, "Generalized least squares mean coefficients.")
    log_marginal_likelihood = _member("LogMarginalLikelihood", float, "Log marginal likelihood after training.")
    status = _member("Status", str, "``OptimizationStatus`` member name of the last training.")
    is_trained = _member("IsTrained", bool, "Whether the model is conditioned on the data.")

    @property
    def kernel(self) -> GaussianProcessKernel:
        """Covariance function."""
        return GaussianProcessKernel(str(self.net.Kernel))

    @property
    def mean_function(self) -> GaussianProcessMean:
        """Mean function."""
        return GaussianProcessMean(str(self.net.MeanFunction))


class Kriging(_NetModel):
    """Simple, ordinary and universal kriging with a variogram model (``Numerics.MachineLearning.Kriging``).

    Parameters
    ----------
    coordinates : array_like
        Observation locations, shape (n, d); a 1-D array is one coordinate.
    values : array_like
        Observed values, length n.
    model : VariogramModel or str
        Variogram model.
    nugget : float
        Nugget c0, non-negative.
    partial_sill : float
        Partial sill c, positive.
    range : float
        Range parameter a (gstat ``vgm`` convention), positive.
    kriging_type : KrigingType or str, optional
        Kriging variant; default ordinary.
    mean : float, optional
        Known mean for simple kriging; default 0.
    distance : callable, str or System.Delegate, optional
        Distance between two locations, ``distance(a, b) -> float``, as a Python callable, C# source
        or a ``Func<double[], double[], double>``; default Euclidean.

    Attributes
    ----------
    net : Numerics.MachineLearning.Kriging
        The wrapped .NET object.
    """

    def __init__(self, coordinates: Any, values: Any, model: VariogramModel | str, nugget: float,
                 partial_sill: float, range: float, kriging_type: KrigingType | str = KrigingType.ORDINARY,
                 mean: float = 0.0, distance: Any = None) -> None:
        net_type = _net_type("Kriging")
        ml = namespace("Numerics.MachineLearning")
        head = (_net_columns(coordinates), to_net_doubles(values), net_enum(ml.VariogramModel, model),
                float(nugget), float(partial_sill), float(range))
        tail = (net_enum(ml.KrigingType, kriging_type), float(mean))
        if distance is None:
            self.net = net_type(*head, *tail)
        else:
            function = callbacks.to_delegate(distance, "Func<double[],double[],double>")
            self.net = net_type(*head, function, *tail)
            callbacks.apply_parallel_policy(self.net, function)

    def semivariance(self, h: Any) -> Any:
        """Return the variogram at lag distance ``h`` (``Semivariance``).

        Parameters
        ----------
        h : float or array_like
            Lag distances, non-negative.

        Returns
        -------
        float or numpy.ndarray
            ``gamma(h)``, zero at h = 0.
        """
        return loop(h, self.net.Semivariance)

    def covariance(self, h: Any) -> Any:
        """Return the covariance at lag distance ``h`` (``Covariance``).

        Parameters
        ----------
        h : float or array_like
            Lag distances, non-negative.

        Returns
        -------
        float or numpy.ndarray
            ``C(h)``, the sill at h = 0.
        """
        return loop(h, self.net.Covariance)

    def predict(self, x: Any) -> np.ndarray:
        """Return the kriging predictor at each location (``Predict``).

        Parameters
        ----------
        x : array_like
            Locations, shape (m, d); a 1-D array is one coordinate.

        Returns
        -------
        numpy.ndarray
            The m predictions.
        """
        return to_numpy_doubles(self.net.Predict(_net_columns(x)))

    def predict_with_variance(self, x: Any) -> tuple[np.ndarray, np.ndarray] | tuple[float, float]:
        """Return the kriging predictor and kriging variance (``PredictWithVariance``).

        Parameters
        ----------
        x : array_like
            Locations, shape (m, d), or a 1-D array holding the d coordinates of one location.

        Returns
        -------
        tuple
            The m predictions and m variances, or a (mean, variance) pair of floats for one location.
        """
        array = np.asarray(x, dtype=np.float64)
        if array.ndim == 1:
            result = self.net.PredictWithVariance(to_net_doubles(array))
            return float(result.Item1), float(result.Item2)
        return _pair(self.net.PredictWithVariance(to_net_matrix(array)))

    def cross_validate(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the leave-one-out prediction and variance at each observation (``CrossValidate``).

        Returns
        -------
        tuple of numpy.ndarray
            The n predictions and the n variances.
        """
        return _pair(self.net.CrossValidate())

    nugget = _member("Nugget", float, "Nugget c0.")
    partial_sill = _member("PartialSill", float, "Partial sill c.")
    sill = _member("Sill", float, "Sill c0 + c.")
    range = _member("Range", float, "Range parameter a.")
    known_mean = _member("KnownMean", float, "Known mean of simple kriging; NaN otherwise.")
    sample_size = _member("SampleSize", int, "Number of observations n.")
    number_of_dimensions = _member("NumberOfDimensions", int, "Number of coordinates d.")
    trend_coefficients = _member("TrendCoefficients", _doubles, "Estimated trend coefficients (empty for simple).")

    @property
    def model(self) -> VariogramModel:
        """Variogram model."""
        return VariogramModel(str(self.net.Model))

    @property
    def kriging_type(self) -> KrigingType:
        """Kriging variant."""
        return KrigingType(str(self.net.Type))


class KMeans(_NetModel):
    """The k-means clustering by Lloyd iterations (``Numerics.MachineLearning.KMeans``).

    Parameters
    ----------
    x : array_like
        Data, shape (n, d); a 1-D array is one feature.
    k : int
        Number of clusters.

    Attributes
    ----------
    net : Numerics.MachineLearning.KMeans
        The wrapped .NET object.
    """

    def __init__(self, x: Any, k: int) -> None:
        self.net = _net_type("KMeans")(_net_columns(x), int(k))

    def train(self, seed: int = -1, k_means_plus_plus: bool = True) -> None:
        """Cluster the data (``Train``); iterations stop when no label changes.

        Parameters
        ----------
        seed : int, optional
            Mersenne Twister seed of the initialization; zero or negative seeds from the clock.
        k_means_plus_plus : bool, optional
            The k-means++ seeding (default) or k distinct random rows.
        """
        self.net.Train(int(seed), bool(k_means_plus_plus))

    @staticmethod
    def initialize(x: Any, k: int, seed: int = -1, k_means_plus_plus: bool = True) -> np.ndarray:
        """Return the initial centroids :meth:`train` starts from (``Initialize``).

        Parameters
        ----------
        x : array_like
            Data, shape (n, d); a 1-D array is one feature.
        k : int
            Number of clusters.
        seed : int, optional
            Mersenne Twister seed; zero or negative seeds from the clock.
        k_means_plus_plus : bool, optional
            The k-means++ seeding (default) or k distinct random rows.

        Returns
        -------
        numpy.ndarray
            Shape (k, d).
        """
        net_type = _net_type("KMeans")
        matrix = _numerics_matrix(_columns(x))
        return to_numpy_matrix(net_type.Initialize(matrix, int(k), int(seed), bool(k_means_plus_plus)))

    k = _member("K", int, "Number of clusters.")
    dimension = _member("Dimension", int, "Number of features d.")
    max_iterations = _member("MaxIterations", int, "Maximum number of iterations (default 1000).", int)
    iterations = _member("Iterations", int, "Number of iterations of the last training.")
    means = _member("Means", _matrix, "Cluster centroids, shape (k, d).")
    labels = _member("Labels", _ints, "Cluster label of each row.")
    x = _member("X", _matrix, "Data, shape (n, d).")


class GaussianMixtureModel(_NetModel):
    """Gaussian mixture model fitted by expectation maximization (``Numerics.MachineLearning.GaussianMixtureModel``).

    Parameters
    ----------
    x : array_like
        Data, shape (n, d); a 1-D array is one feature.
    k : int
        Number of components.

    Attributes
    ----------
    net : Numerics.MachineLearning.GaussianMixtureModel
        The wrapped .NET object.

    Notes
    -----
    The components start from a k-means fit with the same seed, with full covariance matrices.
    """

    def __init__(self, x: Any, k: int) -> None:
        self.net = _net_type("GaussianMixtureModel")(_net_columns(x), int(k))

    def train(self, seed: int = -1, k_means_plus_plus: bool = True) -> None:
        """Fit the mixture (``Train``).

        Parameters
        ----------
        seed : int, optional
            Seed of the k-means initialization; zero or negative seeds from the clock.
        k_means_plus_plus : bool, optional
            The k-means++ seeding (default) or random rows.
        """
        self.net.Train(int(seed), bool(k_means_plus_plus))

    @property
    def sigmas(self) -> np.ndarray | None:
        """Component covariance matrices, shape (k, d, d)."""
        sigmas = self.net.Sigmas
        return None if sigmas is None else np.stack([to_numpy_matrix(s.Array) for s in sigmas])

    k = _member("K", int, "Number of components.")
    dimension = _member("Dimension", int, "Number of features d.")
    max_iterations = _member("MaxIterations", int, "Maximum number of EM iterations (default 1000).", int)
    tolerance = _member("Tolerance", float, "Log-likelihood convergence tolerance (default 1e-8).", float)
    iterations = _member("Iterations", int, "Number of EM iterations of the last training.")
    means = _member("Means", _matrix, "Component means, shape (k, d).")
    weights = _member("Weights", _doubles, "Component weights.")
    labels = _member("Labels", _ints, "Most likely component of each row.")
    likelihood_matrix = _member("LikelihoodMatrix", _matrix, "Component membership probabilities, shape (n, k).")
    log_likelihood = _member("LogLikelihood", float, "Log-likelihood of the fitted mixture.")
    x = _member("X", _matrix, "Data, shape (n, d).")


class JenksCluster(_NetModel):
    """One class of a Jenks natural breaks classification (``Numerics.MachineLearning.JenksCluster``).

    Parameters
    ----------
    data : array_like
        Sorted data.
    start_index : int
        Index of the first value of the cluster.
    end_index : int
        Index of the last value of the cluster (inclusive).

    Attributes
    ----------
    net : Numerics.MachineLearning.JenksCluster
        The wrapped .NET object.
    """

    def __init__(self, data: Any, start_index: int, end_index: int) -> None:
        self.net = _net_type("JenksCluster")(to_net_doubles(data), int(start_index), int(end_index))

    start_index = _member("StartIndex", int, "Index of the first value.")
    end_index = _member("EndIndex", int, "Index of the last value (inclusive).")
    count = _member("Count", int, "Number of values.")
    min_value = _member("MinValue", float, "Smallest value.")
    max_value = _member("MaxValue", float, "Largest value.")
    sum = _member("Sum", float, "Sum of the values.")
    average = _member("Average", float, "Mean of the values.")
    variance = _member("Variance", float, "Sample variance (n - 1 divisor; zero for one value).")
    sum_of_squared_deviations = _member("SumOfSquaredDeviations", float, "Sum of squared deviations from the mean.")


class JenksNaturalBreaks(_NetModel):
    """Jenks natural breaks optimization (``Numerics.MachineLearning.JenksNaturalBreaks``).

    Parameters
    ----------
    data : array_like
        Values to classify.
    number_of_clusters : int
        Number of classes.
    is_data_sorted : bool, optional
        Whether ``data`` is already sorted ascending; default False.

    Attributes
    ----------
    net : Numerics.MachineLearning.JenksNaturalBreaks
        The wrapped .NET object.

    Notes
    -----
    The classification is computed by the constructor: it minimizes the total within-class sum of
    squared deviations over contiguous classes of the sorted data.
    """

    def __init__(self, data: Any, number_of_clusters: int, is_data_sorted: bool = False) -> None:
        self.net = _net_type("JenksNaturalBreaks")(to_net_doubles(data), int(number_of_clusters), bool(is_data_sorted))

    sorted_data = _member("SortedData", _doubles, "The data sorted ascending.")
    number_of_clusters = _member("NumberOfClusters", int, "Number of classes.")
    breaks = _member("Breaks", _doubles, "Break points.")
    goodness_of_variance_fit = _member("GoodnessOfVarianceFit", float,
                                       "Goodness of variance fit, 1 - SDCM / SDAM; closer to 1 is better.")

    @property
    def clusters(self) -> list[JenksCluster]:
        """The classes, in ascending order."""
        clusters = self.net.Clusters
        return [] if clusters is None else [JenksCluster._wrap(cluster) for cluster in clusters]


class PrincipalComponentAnalysis(_NetModel):
    """Principal component analysis (``Numerics.MachineLearning.PrincipalComponentAnalysis``).

    Parameters
    ----------
    x : array_like
        Data, shape (n, p), observations in rows.
    basis : PrincipalComponentBasis or str, optional
        Covariance (default, centered) or correlation (centered and scaled) basis.
    number_of_components : int, optional
        Components to retain; zero or negative (default) retains min(n, p).

    Attributes
    ----------
    net : Numerics.MachineLearning.PrincipalComponentAnalysis
        The wrapped .NET object.

    Notes
    -----
    Each principal axis is oriented so that its element of largest magnitude is positive, the
    scikit-learn convention. Variances use the n - 1 divisor.
    """

    def __init__(self, x: Any, basis: PrincipalComponentBasis | str = PrincipalComponentBasis.COVARIANCE,
                 number_of_components: int = 0) -> None:
        net_type = _net_type("PrincipalComponentAnalysis")
        ml = namespace("Numerics.MachineLearning")
        self.net = net_type(_net_columns(x), net_enum(ml.PrincipalComponentBasis, basis), int(number_of_components))

    def train(self) -> None:
        """Compute the principal components (``Train``)."""
        self.net.Train()

    def transform(self, x: Any) -> np.ndarray:
        """Project data onto the retained components (``Transform``).

        Parameters
        ----------
        x : array_like
            Data, shape (m, p), or a 1-D array holding the p values of one observation.

        Returns
        -------
        numpy.ndarray
            Scores, shape (m, k), or length k for one observation.
        """
        array = np.asarray(x, dtype=np.float64)
        if array.ndim == 1:
            return to_numpy_doubles(self.net.Transform(to_net_doubles(array)))
        return to_numpy_matrix(self.net.Transform(to_net_matrix(array)))

    def inverse_transform(self, scores: Any) -> np.ndarray:
        """Map component scores back to the data space (``InverseTransform``).

        Parameters
        ----------
        scores : array_like
            Scores, shape (m, k), or a 1-D array holding the k scores of one observation.

        Returns
        -------
        numpy.ndarray
            Reconstructed data, shape (m, p), or length p for one observation.
        """
        array = np.asarray(scores, dtype=np.float64)
        if array.ndim == 1:
            return to_numpy_doubles(self.net.InverseTransform(to_net_doubles(array)))
        return to_numpy_matrix(self.net.InverseTransform(to_net_matrix(array)))

    @property
    def basis(self) -> PrincipalComponentBasis:
        """Covariance or correlation basis."""
        return PrincipalComponentBasis(str(self.net.Basis))

    number_of_samples = _member("NumberOfSamples", int, "Number of observations n.")
    number_of_features = _member("NumberOfFeatures", int, "Number of variables p.")
    number_of_components = _member("NumberOfComponents", int, "Number of retained components k.")
    is_trained = _member("IsTrained", bool, "Whether the analysis has been trained.")
    means = _member("Means", _doubles, "Column means of the training data.")
    scales = _member("Scales", _doubles, "Column scales: standard deviations (correlation basis) or ones.")
    components = _member("Components", _matrix, "Principal axes, shape (k, p); row i is component i.")
    loadings = _member("Loadings", _matrix, "Loadings, shape (p, k): axes times the square root of the variance.")
    explained_variance = _member("ExplainedVariance", _doubles, "Variance explained by each component.")
    explained_variance_ratio = _member("ExplainedVarianceRatio", _doubles, "Fraction of the total variance.")
    cumulative_explained_variance_ratio = _member("CumulativeExplainedVarianceRatio", _doubles,
                                                  "Cumulative fraction of the total variance.")
    singular_values = _member("SingularValues", _doubles, "Singular values of the centered (scaled) data.")
    total_variance = _member("TotalVariance", float, "Trace of the covariance or correlation matrix.")
    scores = _member("Scores", _matrix, "Component scores of the training data, shape (n, k).")


__all__ = [
    "DecisionNode", "DecisionTree", "ElasticNetRegression", "GaussianMixtureModel", "GaussianProcessKernel",
    "GaussianProcessMean", "GaussianProcessRegression", "GeneralizedLinearModel", "JenksCluster",
    "JenksNaturalBreaks", "KMeans", "KNearestNeighbors", "Kriging", "KrigingType", "NaiveBayes",
    "PrincipalComponentAnalysis", "PrincipalComponentBasis", "RandomForest", "VariogramModel",
]

translate_public(globals())
