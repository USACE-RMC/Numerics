"""Cross-cutting behavior: Python exceptions, value semantics, the family registry and enum lookup."""

from __future__ import annotations

import copy
import importlib
import inspect
import math

import numpy as np
import pytest

from pynumerics import (
    bootstrap,
    copulas,
    data,
    distributions,
    enums,
    errors,
    linalg,
    multivariate,
    roots,
    sampling,
    statistics,
)

TRANSLATED_MODULES = (
    "pynumerics.distributions", "pynumerics.multivariate", "pynumerics.copulas", "pynumerics.data",
    "pynumerics.statistics", "pynumerics.sampling", "pynumerics.mcmc", "pynumerics.mathematics", "pynumerics.bootstrap",
    "pynumerics.callbacks", "pynumerics.integration", "pynumerics.differentiation", "pynumerics.roots",
    "pynumerics.ode", "pynumerics.machine_learning", "pynumerics.hypothesis_tests", "pynumerics.goodness_of_fit",
    "pynumerics.probability", "pynumerics.optimization", "pynumerics.graphs", "pynumerics.special", "pynumerics.linalg",
    "pynumerics.functions", "pynumerics.timeseries", "pynumerics.download",
)


def _members() -> list[str]:
    """Return every public function and class member of the wrapper modules, as dotted names."""
    names = []
    for module_name in TRANSLATED_MODULES:
        module = importlib.import_module(module_name)
        for name in module.__all__:
            obj = getattr(module, name)
            if getattr(obj, "__module__", None) != module_name:
                continue
            if inspect.isfunction(obj):
                names.append(f"{module_name}:{name}")
            elif inspect.isclass(obj) and not issubclass(obj, (BaseException,)) and not issubclass(obj, enums.NetEnum):
                for member, value in vars(obj).items():
                    if member.startswith("_") and member != "__init__":
                        continue
                    if inspect.isfunction(value) or isinstance(value, (property, classmethod, staticmethod)):
                        names.append(f"{module_name}:{name}.{member}")
    return names


@pytest.mark.parametrize("module_name", TRANSLATED_MODULES)
def test_every_exported_function_and_class_is_defined_in_its_module(module_name: str) -> None:
    # Arrange
    module = importlib.import_module(module_name)
    # Act
    foreign = [name for name in module.__all__
               if (inspect.isfunction(getattr(module, name)) or inspect.isclass(getattr(module, name)))
               and getattr(module, name).__module__ != module_name]
    # Assert: one home per member, the module that mirrors its Numerics namespace or class
    assert foreign == []


@pytest.mark.parametrize(("module_name", "name"), [
    ("pynumerics.mathematics", "beta_incomplete"), ("pynumerics.mathematics", "beta_incomplete_inverse"),
    ("pynumerics.mathematics", "LinkController"), ("pynumerics.mathematics", "make_symmetric_positive_definite"),
    ("pynumerics.mathematics", "cholesky"), ("pynumerics.mathematics", "matrix_product"),
    ("pynumerics.mathematics", "brent_solve"), ("pynumerics.mathematics", "brent_bracket"),
    ("pynumerics.mathematics", "brent_search"), ("pynumerics.mathematics", "BrentSearchResult"),
    ("pynumerics.mathematics", "integrate_adaptive_gauss_kronrod"), ("pynumerics.graphs", "knapsack"),
    ("pynumerics.graphs", "DynamicProgramming"), ("pynumerics.graphs", "DynamicProgrammingPath"),
    *(("pynumerics.distributions", name) for name in
      ("EstimationMethod", "Extrapolation", "KernelType", "PlottingPosition", "SortOrder", "Transform")),
])
def test_members_are_reachable_only_from_the_module_of_their_namespace(module_name: str, name: str) -> None:
    # Act / Assert
    assert not hasattr(importlib.import_module(module_name), name)


def test_gamma_max_is_a_special_function_constant() -> None:
    # Act / Assert: Gamma.GammaMax lives with the gamma function, not with the Tools constants
    from pynumerics import mathematics, special  # noqa: PLC0415

    assert not hasattr(mathematics.constants, "gamma_max") and special.gamma_max() > 171.0


@pytest.mark.parametrize("name", _members())
def test_every_public_member_translates_net_exceptions(name: str) -> None:
    # Arrange
    module_name, path = name.split(":")
    owner, _, member = path.partition(".")
    obj = getattr(importlib.import_module(module_name), owner)
    # Act
    target = vars(obj)[member] if member else obj
    # Assert
    assert errors.is_translated(target), f"{name} does not translate .NET exceptions"


@pytest.mark.parametrize(("call", "kind", "message"), [
    (lambda: distributions.Normal().inverse_cdf(1.5), ValueError, "Probability must be between 0 and 1."),
    (lambda: distributions.Normal().inverse_cdf([0.5, 1.5]), ValueError, "Probability must be between 0 and 1."),
    (lambda: distributions.Normal(0.0, -1.0).cdf(0.0), ValueError, "Standard deviation must be positive."),
    (lambda: statistics.percentile([1.0, 2.0, 3.0], 1.5), ValueError, r"k must be in \[0,1\]"),
    (lambda: statistics.wilson_score_interval(5, 3, 0.9), ValueError, "number of successes must be between"),
    (lambda: roots.brent_solve(lambda x: x * x + 1.0, 0.0, 1.0), ValueError, "root is not bracketed"),
    (lambda: linalg.CholeskyDecomposition([[1.0, 2.0], [2.0, 1.0]]), RuntimeError, "not positive-definite"),
    (lambda: copulas.ClaytonCopula(-5.0).cdf(0.5, 0.5), ValueError, "must be greater than or equal to -1"),
    (lambda: multivariate.MultivariateNormal([0.0, 0.0], [[1.0, 2.0], [2.0, 1.0]]), ValueError,
     "Covariance matrix is not positive-definite."),
    (lambda: multivariate.Dirichlet([1.0, 0.0]), ValueError, "must be positive"),
    (lambda: sampling.latin_hypercube(-1), OverflowError, "overflow"),
    (lambda: bootstrap.BootstrapAnalysis(distributions.Normal(), "MethodOfMoments", 0, 10), ValueError,
     "The sample size must be at least 10."),
    (lambda: data.Linear([1.0, 2.0, 3.0], [1.0, 2.0]), ValueError, None),
])
def test_net_exceptions_surface_as_python_exceptions_with_the_numerics_message(call, kind, message) -> None:
    # Act / Assert
    with pytest.raises(kind, match=message) as raised:
        call()
    import System  # noqa: PLC0415 - the call started the runtime

    assert isinstance(raised.value.__cause__, System.Exception)
    assert raised.type.__module__ == "builtins"


def test_python_exceptions_from_callbacks_propagate_unchanged() -> None:
    # Arrange
    def failing(x: float) -> float:
        raise KeyError("from the callback")

    # Act / Assert
    with pytest.raises(KeyError, match="from the callback"):
        roots.brent_solve(failing, 0.0, 1.0)


def test_translate_maps_each_net_family() -> None:
    # Arrange
    distributions.Normal()
    import System  # noqa: PLC0415

    cases = [
        (System.ArgumentOutOfRangeException("x", "bad x"), ValueError),
        (System.OverflowException("big"), OverflowError),
        (System.DivideByZeroException("zero"), ZeroDivisionError),
        (System.ArithmeticException("arith"), ArithmeticError),
        (System.NotSupportedException("unsupported"), NotImplementedError),
        (System.NotImplementedException("missing"), NotImplementedError),
        (System.InvalidOperationException("state"), RuntimeError),
        (System.Collections.Generic.KeyNotFoundException("key"), RuntimeError),
        (System.Reflection.TargetInvocationException(System.ArgumentException("inner")), ValueError),
    ]
    # Act / Assert
    for error, kind in cases:
        mapped = errors.translate(error)
        assert type(mapped) is kind
        assert str(mapped) == str(error.InnerException.Message if error.InnerException else error.Message)
    assert errors.translate(KeyError("python")) is None


def test_equality_compares_family_and_parameter_bits() -> None:
    # Arrange
    gev = distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1)
    same = distributions.create_distribution("GeneralizedExtremeValue", [100.0, 20.0, -0.1])
    # Act / Assert
    assert gev == same and type(same) is distributions.GeneralizedExtremeValue
    assert gev != distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1 + 1e-16)
    assert gev != distributions.GeneralizedPareto(100.0, 20.0, -0.1)
    assert distributions.Normal(0.0, 1.0) != distributions.Normal(-0.0, 1.0)
    assert gev != "GeneralizedExtremeValue"
    assert distributions.UnivariateDistribution(gev.net) == gev


def test_equality_includes_base_and_table_state() -> None:
    # Arrange
    base10 = distributions.LogNormal(2.0, 0.3)
    base_e = distributions.LogNormal(2.0, 0.3)
    base_e.net.Base = math.e
    table = distributions.Empirical([1.0, 2.0, 3.0], [0.1, 0.5, 0.9])
    # Act / Assert
    assert base10 == distributions.LogNormal(2.0, 0.3)
    assert base10 != base_e
    assert table == distributions.Empirical([1.0, 2.0, 3.0], [0.1, 0.5, 0.9])
    assert table != distributions.Empirical([1.0, 2.0, 3.0], [0.1, 0.6, 0.9])


def test_mutable_distributions_are_unhashable_and_frozen_ones_hash_by_value() -> None:
    # Arrange
    gev = distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1)
    frozen = gev.frozen()
    # Act
    gev.set_parameters([1.0, 2.0, 0.0])
    # Assert
    with pytest.raises(TypeError):
        hash(gev)
    assert isinstance(frozen, distributions.FrozenDistribution)
    assert isinstance(frozen, distributions.GeneralizedExtremeValue)
    assert frozen.parameters.tolist() == [100.0, 20.0, -0.1]
    assert hash(frozen) == hash(distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1).frozen())
    assert {frozen: 1}[distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1).frozen()] == 1
    assert frozen == distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1)
    assert frozen.frozen() is frozen and copy.deepcopy(frozen) is frozen and copy.copy(frozen) is frozen
    assert frozen.cdf(120.0) == distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1).cdf(120.0)


def test_frozen_distributions_refuse_changes() -> None:
    # Arrange
    frozen = distributions.Empirical([1.0, 2.0, 3.0], [0.1, 0.5, 0.9]).frozen()
    # Act / Assert
    with pytest.raises(TypeError, match="frozen"):
        frozen.set_parameters([])
    with pytest.raises(TypeError, match="frozen"):
        frozen.estimate([1.0, 2.0, 3.0])
    with pytest.raises(AttributeError, match="frozen"):
        frozen.x_transform = "logarithmic"
    with pytest.raises(AttributeError, match="frozen"):
        frozen.net = None
    thawed = frozen.clone()
    thawed.x_transform = "logarithmic"
    assert type(thawed) is distributions.Empirical and frozen.x_transform == "None"


def test_generic_wrapper_freezes_to_the_frozen_base() -> None:
    # Act
    frozen = distributions.UnivariateDistribution(distributions.Weibull(2.0, 3.0).net).frozen()
    # Assert
    assert type(frozen) is distributions.FrozenDistribution
    assert frozen == distributions.create_distribution("Weibull", [2.0, 3.0])


def test_families_registry_and_class_metadata() -> None:
    # Act
    families = distributions.FAMILIES
    # Assert
    assert families["Beta"] is distributions.Beta and families["Empirical"] is distributions.Empirical
    with pytest.raises(TypeError):
        families["Normal"] = distributions.Normal  # type: ignore[index]
    for name, cls in families.items():
        assert cls.type_name == name
    assert distributions.GeneralizedExtremeValue.parameter_names == ("Xi", "Alpha", "Kappa")
    assert distributions.Normal.parameter_names == distributions.Normal().parameter_names
    assert distributions.UnivariateDistribution.type_name is None
    assert distributions.Mixture.parameter_names is None
    assert distributions.create_distribution("Weibull", [2.0, 3.0]).type_name == "Weibull"


def test_validate_parameters_on_the_class_checks_arity_then_numerics() -> None:
    # Act / Assert
    assert distributions.GeneralizedExtremeValue.validate_parameters([1.0, 1.0, 0.0]) is None
    assert "must be positive" in distributions.GeneralizedExtremeValue.validate_parameters([1.0, -1.0, 0.0])
    with pytest.raises(ValueError, match="GeneralizedExtremeValue takes 3 parameters"):
        distributions.GeneralizedExtremeValue.validate_parameters([1.0, 1.0])
    with pytest.raises(ValueError, match="takes 2 parameters"):
        distributions.Gumbel(0.0, 1.0).validate_parameters([1.0])
    with pytest.raises(TypeError):
        distributions.UnivariateDistribution.validate_parameters([1.0])
    with pytest.raises(NotImplementedError, match="requires external components"):
        distributions.Mixture.validate_parameters([1.0])


def test_create_distribution_names_the_expected_parameter_count() -> None:
    # Act / Assert
    with pytest.raises(ValueError, match=r"Gumbel takes 2 parameters \(Xi, Alpha\); got 3\."):
        distributions.create_distribution("Gumbel", [1.0, 2.0, 3.0])
    assert distributions.create_distribution("Gumbel", np.array([1.0, 2.0])).parameters.tolist() == [1.0, 2.0]


@pytest.mark.parametrize(("enum_type", "spelling", "member"), [
    (enums.Transform, "normal_z", enums.Transform.NORMAL_Z),
    (enums.Transform, "NORMAL_Z", enums.Transform.NORMAL_Z),
    (enums.Transform, "NormalZ", enums.Transform.NORMAL_Z),
    (enums.ImportanceDistribution, "log10_uniform", enums.ImportanceDistribution.LOG10_UNIFORM),
    (enums.SortOrder, "descending", enums.SortOrder.DESCENDING),
    (enums.PlottingPosition, "kaplan_meier", enums.PlottingPosition.KAPLAN_MEIER),
    (enums.EstimationMethod, "method_of_linear_moments", enums.EstimationMethod.METHOD_OF_LINEAR_MOMENTS),
    (enums.WeightType, "Reliability", enums.WeightType.RELIABILITY),
    (enums.KernelType, "gaussian", enums.KernelType.GAUSSIAN),
    (enums.AutocorrelationType, "partial", enums.AutocorrelationType.PARTIAL),
    (enums.Extrapolation, "both", enums.Extrapolation.BOTH),
])
def test_enums_resolve_folded_spellings(enum_type, spelling, member) -> None:
    # Act / Assert
    assert enum_type(spelling) is member


def test_every_package_enum_derives_from_net_enum_and_refuses_unknown_names() -> None:
    # Arrange
    found = [obj for obj in vars(enums).values()
             if inspect.isclass(obj) and issubclass(obj, enums.enum.Enum) and obj is not enums.NetEnum]
    # Act / Assert
    assert found and all(issubclass(obj, enums.NetEnum) for obj in found)
    with pytest.raises(ValueError):
        enums.Transform("normal_y")
    with pytest.raises(ValueError):
        enums.Transform(3)
