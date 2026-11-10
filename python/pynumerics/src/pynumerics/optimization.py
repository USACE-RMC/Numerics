"""Local, global, constrained and multi-objective optimization from ``Numerics.Mathematics.Optimization``.

Every optimizer class wraps the Numerics optimizer of the same name and shares the interface of
:class:`Optimizer`: settings are properties (and keyword arguments of the constructor), and
:meth:`Optimizer.minimize` and :meth:`Optimizer.maximize` run the Numerics search and return an
:class:`OptimizationResult` with the best parameters, the objective value, the counts, the status,
the Hessian and the best-so-far trace.

Objectives, gradients and constraint functions are Python callables of a float64 array, C# source
compiled by :func:`pynumerics.callbacks.compile_csharp`, or .NET delegates. A compiled or .NET
objective runs with no Python crossing per evaluation. With a Python callable the population
optimizers (:class:`DifferentialEvolution`, :class:`ParticleSwarm`) always evaluate serially, because
a Python callback holds the GIL (:func:`pynumerics.callbacks.apply_parallel_policy`); with a compiled
objective ``use_parallel=True`` evaluates each generation on several threads.

Settings left unset keep the Numerics defaults. Numerics minimizes internally: for a maximization it
stores the negated objective as the fitness, and :class:`OptimizationResult` reports the objective
value itself.

:class:`ParetoFront` (dominance, non-dominated sorting, crowding distance, hypervolume and the
generational distances) takes one row of objective values per point and minimizes every objective.

:func:`knapsack` solves the exact 0/1 knapsack, single or multi-constraint.
:class:`DynamicProgramming` solves finite-horizon and stationary discrete dynamic programs with
deterministic or stochastic transitions; its stage-value and transition functions are Python
callables, C# source or .NET delegates, and with a Python callable the solve runs serially
(:func:`pynumerics.callbacks.apply_parallel_policy`). The shortest-path solvers of the same
namespace are in :mod:`pynumerics.graphs`.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import _opt, net_enum
from pynumerics.enums import NetEnum
from pynumerics.errors import translate_class, translate_public
from pynumerics.marshal import (
    _int32_array,
    _net_rectangular,
    _numpy_rectangular,
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_ints,
    to_numpy_matrix,
)


class OptimizationStatus(NetEnum):
    """Optimizer outcomes (``Numerics.Mathematics.Optimization.OptimizationStatus``)."""

    NONE = "None"
    SUCCESS = "Success"
    MAXIMUM_ITERATIONS_REACHED = "MaximumIterationsReached"
    MAXIMUM_FUNCTION_EVALUATIONS_REACHED = "MaximumFunctionEvaluationsReached"
    FAILURE = "Failure"
    LINE_SEARCH_FAILED = "LineSearchFailed"
    CANCELLED = "Cancelled"


class LocalMethod(NetEnum):
    """Local search methods of the multi-start and polishing steps (``LocalMethod``)."""

    ADAM = "ADAM"
    BFGS = "BFGS"
    GRADIENT_DESCENT = "GradientDescent"
    NELDER_MEAD = "NelderMead"
    POWELL = "Powell"
    LBFGSB = "LBFGSB"


class ConstraintType(NetEnum):
    """Constraint relations (``Numerics.Mathematics.Optimization.ConstraintType``)."""

    EQUAL_TO = "EqualTo"
    GREATER_THAN_OR_EQUAL_TO = "GreaterThanOrEqualTo"
    LESSER_THAN_OR_EQUAL_TO = "LesserThanOrEqualTo"


@dataclass(frozen=True, eq=False, slots=True)
class ParameterSet:
    """A parameter vector with its fitness and weight (``Numerics.Mathematics.Optimization.ParameterSet``).

    Parameters
    ----------
    values : numpy.ndarray
        The parameter values.
    fitness : float
        The fitness (the minimized objective; the negated objective of a maximization).
    weight : float, optional
        The weight; default 0.
    """

    values: np.ndarray
    fitness: float
    weight: float = 0.0

    @classmethod
    def from_net(cls, net: Any) -> ParameterSet:
        """Copy a .NET ``ParameterSet`` into a :class:`ParameterSet`.

        Parameters
        ----------
        net : Numerics.Mathematics.Optimization.ParameterSet
            The .NET value.

        Returns
        -------
        ParameterSet
            The copy; ``values`` is empty when the .NET values are null.
        """
        values = net.Values
        array = np.empty(0) if values is None else to_numpy_doubles(values)
        return cls(array, float(net.Fitness), float(net.Weight))

    def to_net(self) -> Any:
        """Return the .NET ``ParameterSet`` with these values, fitness and weight.

        Returns
        -------
        Numerics.Mathematics.Optimization.ParameterSet
            A new .NET value.
        """
        return _opt().ParameterSet(to_net_doubles(self.values), float(self.fitness), float(self.weight))

    def clone(self) -> ParameterSet:
        """Return a copy with its own values array (``ParameterSet.Clone``).

        Returns
        -------
        ParameterSet
            The copy.
        """
        return ParameterSet(np.array(self.values, dtype=np.float64), self.fitness, self.weight)

    def to_xml(self) -> str:
        """Return the Numerics XML form (``ParameterSet.ToXElement``).

        Returns
        -------
        str
            The ``ParameterSet`` element with its ``Values``, ``Fitness`` and ``Weight`` attributes.
        """
        return str(self.to_net().ToXElement().ToString())

    @classmethod
    def from_xml(cls, xml: str) -> ParameterSet:
        """Read the Numerics XML form (the ``ParameterSet(XElement)`` constructor).

        Parameters
        ----------
        xml : str
            A ``ParameterSet`` element, as written by :meth:`to_xml`.

        Returns
        -------
        ParameterSet
            The parameter set.
        """
        opt = _opt()
        from System.Xml.Linq import XElement  # noqa: PLC0415 - the runtime is running

        return cls.from_net(opt.ParameterSet(XElement.Parse(str(xml))))


@dataclass(frozen=True, eq=False, slots=True)
class OptimizationProgress:
    """A progress report of a running optimizer (``OptimizationProgress``).

    Parameters
    ----------
    iterations : int
        Iterations completed.
    function_evaluations : int
        Objective evaluations completed.
    best_fitness : float
        The best fitness so far (the negated objective during a maximization).
    best_values : numpy.ndarray
        The parameters of the best fitness.
    """

    iterations: int
    function_evaluations: int
    best_fitness: float
    best_values: np.ndarray

    @classmethod
    def from_net(cls, net: Any) -> OptimizationProgress:
        """Copy a .NET ``OptimizationProgress``.

        Parameters
        ----------
        net : Numerics.Mathematics.Optimization.OptimizationProgress
            The .NET report.

        Returns
        -------
        OptimizationProgress
            The copy.
        """
        return cls(int(net.Iterations), int(net.FunctionEvaluations), float(net.BestFitness),
                   to_numpy_doubles(net.BestValues))

    def to_net(self) -> Any:
        """Return the .NET ``OptimizationProgress`` with these values.

        Returns
        -------
        Numerics.Mathematics.Optimization.OptimizationProgress
            A new .NET report.
        """
        return _opt().OptimizationProgress(int(self.iterations), int(self.function_evaluations),
                                           float(self.best_fitness), to_net_doubles(self.best_values))


_TRACE_FLATTEN = """
static double[] Flatten(System.Collections.Generic.List<Numerics.Mathematics.Optimization.ParameterSet> trace)
{
    int n = trace.Count, p = n == 0 ? 0 : trace[0].Values.Length;
    var flat = new double[n * (p + 1)];
    for (int i = 0; i < n; i++)
    {
        var item = trace[i];
        System.Array.Copy(item.Values, 0, flat, i * (p + 1), p);
        flat[i * (p + 1) + p] = item.Fitness;
    }
    return flat;
}
"""

BULK_TRACE_LENGTH = 2000
"""Trace length from which a trace crosses in one copy through a compiled flattener, when Roslyn is staged."""


def _trace_arrays(net_trace: Any, scale: float, bulk: bool | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return the parameters and objective values of a .NET ``List<ParameterSet>`` trace.

    Long traces cross in one copy through a compiled C# flattener when Roslyn is staged; short ones, or
    any trace without Roslyn, cross one entry at a time. Both paths copy the same bits.
    """
    count = int(net_trace.Count)
    if count == 0:
        return np.empty((0, 0)), np.empty(0)
    if bulk is None:
        bulk = count >= BULK_TRACE_LENGTH and callbacks.roslyn_available()
    if bulk:
        import System  # noqa: PLC0415 - the runtime is running
        from System.Collections.Generic import List  # noqa: PLC0415

        dtype = System.Func[List[_opt().ParameterSet], System.Array[System.Double]]
        flat = to_numpy_doubles(callbacks.compile_csharp(_TRACE_FLATTEN, dtype)(net_trace))
        table = flat.reshape(count, -1)
        return np.ascontiguousarray(table[:, :-1]), scale * table[:, -1]
    parameters, fitness = [], np.empty(count)
    for i, item in enumerate(net_trace):
        parameters.append(to_numpy_doubles(item.Values))
        fitness[i] = float(item.Fitness)
    return np.vstack(parameters), scale * fitness


@dataclass(frozen=True, eq=False)
class OptimizationResult:
    """The outcome of an optimizer run, shared by every optimizer.

    Parameters
    ----------
    parameters : numpy.ndarray
        The best parameters (``BestParameterSet.Values``).
    value : float
        The objective at ``parameters``: ``BestParameterSet.Fitness``, with the sign restored for a
        maximization.
    iterations : int
        Iterations used (``Iterations``).
    function_evaluations : int
        Objective evaluations used (``FunctionEvaluations``).
    status : OptimizationStatus
        The outcome (``Status``).
    maximized : bool, optional
        Whether the run maximized; default False.
    hessian : numpy.ndarray or None, optional
        The Hessian of the objective at ``parameters`` when Numerics computed it (``ComputeHessian``).
    net_trace : object, optional
        The .NET ``ParameterSetTrace`` list that :attr:`trace_parameters` and :attr:`trace_values` copy on
        first use.
    """

    parameters: np.ndarray
    value: float
    iterations: int
    function_evaluations: int
    status: OptimizationStatus
    maximized: bool = False
    hessian: np.ndarray | None = None
    net_trace: Any = field(default=None, repr=False)
    _cache: dict = field(default_factory=dict, repr=False, init=False)

    @property
    def success(self) -> bool:
        """Whether the status is ``SUCCESS``."""
        return self.status is OptimizationStatus.SUCCESS

    def _trace(self) -> tuple[np.ndarray, np.ndarray]:
        """Copy the .NET trace once and cache it."""
        if "trace" not in self._cache:
            if self.net_trace is None:
                self._cache["trace"] = (np.empty((0, 0)), np.empty(0))
            else:
                self._cache["trace"] = _trace_arrays(self.net_trace, -1.0 if self.maximized else 1.0)
        return self._cache["trace"]

    @property
    def trace_parameters(self) -> np.ndarray:
        """The best parameters after each evaluation, one row per evaluation (``ParameterSetTrace``)."""
        return self._trace()[0]

    @property
    def trace_values(self) -> np.ndarray:
        """The best objective value after each evaluation, with the sign restored for a maximization."""
        return self._trace()[1]


def _vector(fn: Any) -> Any:
    """Return an objective or constraint function as a ``Func<double[], double>``."""
    return callbacks.vector_function(fn)


def _gradient(fn: Any) -> Any:
    """Return a gradient as a ``Func<double[], double[]>``, or None."""
    return None if fn is None else callbacks.to_delegate(fn, "Func<double[],double[]>")


def _setting(name: str, kind: type, doc: str) -> property:
    """Return a read-write property mirroring the .NET property ``name``, converted with ``kind``."""

    def fget(self: Any) -> Any:
        return kind(getattr(self.net, name))

    def fset(self: Any, value: Any) -> None:
        setattr(self.net, name, kind(value))

    return property(fget, fset, doc=doc)


def _array(name: str, doc: str) -> property:
    """Return a read-only property copying the .NET ``double[]`` property ``name`` (null gives None)."""

    def fget(self: Any) -> np.ndarray | None:
        values = getattr(self.net, name)
        return None if values is None else to_numpy_doubles(values)

    return property(fget, doc=doc)


def _method_setting(name: str, doc: str) -> property:
    """Return a read-write :class:`LocalMethod` property mirroring the .NET property ``name``."""

    def fget(self: Any) -> LocalMethod:
        return LocalMethod(str(getattr(self.net, name)))

    def fset(self: Any, value: Any) -> None:
        setattr(self.net, name, net_enum(_opt().LocalMethod, LocalMethod(value).value))

    return property(fget, fset, doc=doc)


_REPORTER: dict[str, Any] = {}


def _progress_reporter(fn: Callable[[OptimizationProgress], Any]) -> Any:
    """Return an ``IProgress<OptimizationProgress>`` that calls ``fn`` with each report, synchronously."""
    reporter = _REPORTER.get("type")
    if reporter is None:
        opt = _opt()
        from System import IProgress  # noqa: PLC0415 - the runtime is running

        class _Reporter(IProgress[opt.OptimizationProgress]):
            __namespace__ = "pynumerics.optimization"

            def Report(self, value: Any) -> None:  # noqa: N802 - the .NET interface member
                self.callback(OptimizationProgress.from_net(value))

        reporter = _REPORTER["type"] = _Reporter
    instance = reporter()
    instance.callback = fn
    return instance


class Optimizer:
    """The settings, run methods and results every optimizer shares (``Optimizer``).

    Construct one of the concrete optimizers; this base class only wraps an existing .NET optimizer.

    Parameters
    ----------
    net : Numerics.Mathematics.Optimization.Optimizer
        The .NET optimizer to wrap.
    functions : sequence, optional
        The callables given to it, for the parallel policy.

    Attributes
    ----------
    net : Numerics.Mathematics.Optimization.Optimizer
        The wrapped .NET optimizer.
    """

    def __init__(self, net: Any, functions: Sequence[Any] = ()) -> None:
        self.net = net
        self._functions = tuple(functions)
        self._progress: Callable[[OptimizationProgress], Any] | None = None
        self._maximized = False
        callbacks.apply_parallel_policy(net, *self._functions)

    def _configure(self, settings: dict[str, Any]) -> None:
        """Set each keyword setting through its property (None keeps the default), then re-apply the policy."""
        for name, value in settings.items():
            member = getattr(type(self), name, None)
            if not isinstance(member, property) or member.fset is None:
                raise TypeError(f"{type(self).__name__} has no setting {name!r}.")
            if value is not None:
                setattr(self, name, value)
        callbacks.apply_parallel_policy(self.net, *self._functions)

    max_iterations = _setting("MaxIterations", int, "Iteration limit (``MaxIterations``).")
    max_function_evaluations = _setting("MaxFunctionEvaluations", int,
                                        "Objective evaluation limit (``MaxFunctionEvaluations``).")
    absolute_tolerance = _setting("AbsoluteTolerance", float, "Absolute convergence tolerance (``AbsoluteTolerance``).")
    relative_tolerance = _setting("RelativeTolerance", float, "Relative convergence tolerance (``RelativeTolerance``).")
    report_failure = _setting("ReportFailure", bool, "Whether a failed search raises (``ReportFailure``).")
    record_traces = _setting("RecordTraces", bool, "Whether the best-so-far trace is recorded (``RecordTraces``).")
    compute_hessian = _setting("ComputeHessian", bool,
                               "Whether the Hessian is computed after convergence (``ComputeHessian``).")
    progress_interval = _setting("ProgressInterval", int,
                                 "Evaluations between progress reports (``ProgressInterval``).")

    @property
    def cancellation_token(self) -> Any:
        """The ``System.Threading.CancellationToken`` that stops the search (``CancellationToken``)."""
        return self.net.CancellationToken

    @cancellation_token.setter
    def cancellation_token(self, token: Any) -> None:
        self.net.CancellationToken = token

    @property
    def progress(self) -> Callable[[OptimizationProgress], Any] | None:
        """A callable given an :class:`OptimizationProgress` every ``progress_interval`` evaluations (``Progress``).

        It runs synchronously on the optimizer thread. None, the default, turns reporting off.
        """
        return self._progress

    @progress.setter
    def progress(self, fn: Callable[[OptimizationProgress], Any] | None) -> None:
        self.net.Progress = None if fn is None else _progress_reporter(fn)
        self._progress = fn

    @property
    def objective_function(self) -> Any:
        """The objective as the .NET ``Func<double[], double>`` Numerics calls (``ObjectiveFunction``).

        Setting it accepts a Python callable, C# source or a delegate, and re-applies the parallel policy.
        """
        return self.net.ObjectiveFunction

    @objective_function.setter
    def objective_function(self, fn: Any) -> None:
        self.net.ObjectiveFunction = _vector(fn)
        self._functions = (fn, *self._functions[1:])
        callbacks.apply_parallel_policy(self.net, *self._functions)

    @property
    def number_of_parameters(self) -> int:
        """The number of parameters (``NumberOfParameters``)."""
        return int(self.net.NumberOfParameters)

    @property
    def iterations(self) -> int:
        """Iterations used by the last run (``Iterations``)."""
        return int(self.net.Iterations)

    @property
    def function_evaluations(self) -> int:
        """Objective evaluations used by the last run (``FunctionEvaluations``)."""
        return int(self.net.FunctionEvaluations)

    @property
    def status(self) -> OptimizationStatus:
        """The outcome of the last run (``Status``)."""
        return OptimizationStatus(str(self.net.Status))

    @property
    def best_parameter_set(self) -> ParameterSet:
        """The best parameters and fitness (``BestParameterSet``; the fitness is negated after a maximization)."""
        return ParameterSet.from_net(self.net.BestParameterSet)

    @property
    def parameter_set_trace(self) -> list[ParameterSet]:
        """The best parameter set after each evaluation, fitness as Numerics stores it (``ParameterSetTrace``)."""
        return [ParameterSet.from_net(item) for item in self.net.ParameterSetTrace]

    @property
    def hessian(self) -> np.ndarray | None:
        """The Hessian at the best parameters, or None when it was not computed (``Hessian``)."""
        matrix = self.net.Hessian
        return None if matrix is None else to_numpy_matrix(matrix.ToArray())

    @property
    def result(self) -> OptimizationResult:
        """The :class:`OptimizationResult` of the last run."""
        best = self.net.BestParameterSet
        values = best.Values
        scale = -1.0 if self._maximized else 1.0
        return OptimizationResult(np.empty(0) if values is None else to_numpy_doubles(values),
                                  scale * float(best.Fitness), self.iterations, self.function_evaluations,
                                  self.status, self._maximized, self.hessian, self.net.ParameterSetTrace)

    def minimize(self) -> OptimizationResult:
        """Minimize the objective (``Minimize``).

        Returns
        -------
        OptimizationResult
            The best parameters, the objective value, the counts and the status.

        Raises
        ------
        RuntimeError
            Whatever Numerics throws, translated by ``pynumerics.errors``, for example when the search
            fails with ``report_failure`` set.
        """
        self._maximized = False
        self.net.Minimize()
        return self.result

    def maximize(self) -> OptimizationResult:
        """Maximize the objective (``Maximize``).

        Returns
        -------
        OptimizationResult
            The best parameters, the objective value (not negated), the counts and the status.

        Raises
        ------
        RuntimeError
            Whatever Numerics throws, translated by ``pynumerics.errors``.
        """
        self._maximized = True
        self.net.Maximize()
        return self.result

    def clear_results(self) -> None:
        """Reset the counts, best parameter set, trace, status and Hessian (``ClearResults``)."""
        self.net.ClearResults()


def _bounds(*arrays: Any) -> list[Any]:
    """Copy parameter vectors to ``System.Double[]``."""
    return [to_net_doubles(a) for a in arrays]


class _BoundedSearch(Optimizer):
    """Bounds shared by the multidimensional optimizers."""

    lower_bounds = _array("LowerBounds", "The parameter lower bounds (``LowerBounds``).")
    upper_bounds = _array("UpperBounds", "The parameter upper bounds (``UpperBounds``).")


class _LocalSearch(_BoundedSearch):
    """Initial values shared by the local optimizers."""

    initial_values = _array("InitialValues", "The starting parameters (``InitialValues``).")


class _GradientSearch(_LocalSearch):
    """The optional analytic gradient of the gradient-based local optimizers."""

    @property
    def gradient(self) -> Any:
        """The analytic gradient as a .NET ``Func<double[], double[]>``, or None for finite differences (``Gradient``).

        Setting it accepts a Python callable returning an array, C# source or a delegate.
        """
        return self.net.Gradient

    @gradient.setter
    def gradient(self, fn: Any) -> None:
        self.net.Gradient = _gradient(fn)
        self._functions = (self._functions[0], fn)
        callbacks.apply_parallel_policy(self.net, *self._functions)


class _ScalarSearch(Optimizer):
    """The interval of the one-dimensional optimizers."""

    @property
    def lower_bound(self) -> float:
        """The lower end of the search interval (``LowerBound``)."""
        return float(self.net.LowerBound)

    @property
    def upper_bound(self) -> float:
        """The upper end of the search interval (``UpperBound``)."""
        return float(self.net.UpperBound)


class BrentSearch(_ScalarSearch):
    """Brent's one-dimensional minimizer on a bounded interval (``BrentSearch``).

    The search is a port of ``scipy.optimize.minimize_scalar(method="bounded")``. Its step tolerance is
    ``relative_tolerance * |x| + absolute_tolerance / 3``; the constructor sets ``absolute_tolerance`` to about
    3.3e-19, so the default is ``1e-8 * |x| + 1.1e-19``. With ``relative_tolerance=math.sqrt(2.2e-16)``,
    ``absolute_tolerance=1e-5`` and ``max_function_evaluations=500`` (SciPy's ``xatol`` and ``maxiter``
    defaults) the minimizer is bitwise equal to SciPy's.

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of one float (``Func<double, double>``).
    lower_bound, upper_bound : float
        The search interval.
    **settings
        Initial values of any settings property, for example ``relative_tolerance``.
    """

    def __init__(self, objective: Any, lower_bound: float, upper_bound: float, **settings: Any) -> None:
        super().__init__(_opt().BrentSearch(callbacks.scalar_function(objective), float(lower_bound),
                                            float(upper_bound)), (objective,))
        self._configure(settings)

    def bracket(self, s: float = 1e-2, k: float = 2.0) -> tuple[float, float]:
        """Move the interval downhill until it brackets a minimum (``Bracket``).

        Parameters
        ----------
        s : float, optional
            Initial step; default 1e-2, as in Numerics.
        k : float, optional
            Step growth factor; default 2, as in Numerics.

        Returns
        -------
        tuple of float
            The new ``(lower_bound, upper_bound)``.
        """
        self.net.Bracket(float(s), float(k))
        return self.lower_bound, self.upper_bound


class GoldenSection(_ScalarSearch):
    """Golden-section one-dimensional minimizer on a bounded interval (``GoldenSection``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of one float (``Func<double, double>``).
    lower_bound, upper_bound : float
        The search interval.
    **settings
        Initial values of any settings property.
    """

    def __init__(self, objective: Any, lower_bound: float, upper_bound: float, **settings: Any) -> None:
        super().__init__(_opt().GoldenSection(callbacks.scalar_function(objective), float(lower_bound),
                                              float(upper_bound)), (objective,))
        self._configure(settings)


class NelderMead(_LocalSearch):
    """Bounded Nelder-Mead downhill simplex (``NelderMead``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array (``Func<double[], double>``).
    initial_values, lower_bounds, upper_bounds : array_like
        Start point and bounds, one entry per parameter.
    **settings
        Initial values of any settings property, for example ``simplex_scale``.
    """

    def __init__(self, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 **settings: Any) -> None:
        arrays = _bounds(initial_values, lower_bounds, upper_bounds)
        super().__init__(_opt().NelderMead(_vector(objective), int(arrays[0].Length), *arrays), (objective,))
        self._configure(settings)

    simplex_scale = _setting("SimplexScale", float,
                             "Initial simplex size as a fraction of the range (``SimplexScale``).")
    enable_start_point_probe = _setting("EnableStartPointProbe", bool,
                                        "Whether probes near the start may replace it (``EnableStartPointProbe``).")
    start_point_probe_count = _setting("StartPointProbeCount", int,
                                       "Number of start probes (``StartPointProbeCount``).")
    start_point_probe_radius = _setting("StartPointProbeRadius", float,
                                        "Start probe radius as a fraction of the range (``StartPointProbeRadius``).")


class Powell(_LocalSearch):
    """Bounded Powell conjugate-direction search (``Powell``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    initial_values, lower_bounds, upper_bounds : array_like
        Start point and bounds.
    **settings
        Initial values of any settings property.
    """

    def __init__(self, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 **settings: Any) -> None:
        arrays = _bounds(initial_values, lower_bounds, upper_bounds)
        super().__init__(_opt().Powell(_vector(objective), int(arrays[0].Length), *arrays), (objective,))
        self._configure(settings)


class BFGS(_GradientSearch):
    """Bounded BFGS quasi-Newton search (``BFGS``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    initial_values, lower_bounds, upper_bounds : array_like
        Start point and bounds.
    gradient : callable, str or System.Delegate, optional
        Analytic gradient returning an array; None, the default, uses finite differences.
    **settings
        Initial values of any settings property, for example ``scale_initial_inverse_hessian``.
    """

    def __init__(self, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 gradient: Any = None, **settings: Any) -> None:
        arrays = _bounds(initial_values, lower_bounds, upper_bounds)
        super().__init__(_opt().BFGS(_vector(objective), int(arrays[0].Length), *arrays, _gradient(gradient)),
                         (objective, gradient))
        self._configure(settings)

    scale_initial_inverse_hessian = _setting(
        "ScaleInitialInverseHessian", bool,
        "Whether the first inverse Hessian is scaled (``ScaleInitialInverseHessian``).")


class LBFGSB(_GradientSearch):
    """Limited-memory BFGS with bound constraints (``LBFGSB``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    initial_values, lower_bounds, upper_bounds : array_like
        Start point and bounds.
    gradient : callable, str or System.Delegate, optional
        Analytic gradient; None, the default, uses finite differences.
    **settings
        Initial values of any settings property, for example ``memory_size``.
    """

    def __init__(self, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 gradient: Any = None, **settings: Any) -> None:
        arrays = _bounds(initial_values, lower_bounds, upper_bounds)
        super().__init__(_opt().LBFGSB(_vector(objective), int(arrays[0].Length), *arrays, _gradient(gradient)),
                         (objective, gradient))
        self._configure(settings)

    memory_size = _setting("MemorySize", int, "Number of correction pairs kept (``MemorySize``).")
    max_line_search_steps = _setting("MaxLineSearchSteps", int, "Line-search step limit (``MaxLineSearchSteps``).")


class GradientDescent(_GradientSearch):
    """Bounded gradient descent (``GradientDescent``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    initial_values, lower_bounds, upper_bounds : array_like
        Start point and bounds.
    alpha : float, optional
        Learning rate; None keeps the Numerics default.
    gradient : callable, str or System.Delegate, optional
        Analytic gradient; None, the default, uses finite differences.
    **settings
        Initial values of any settings property.
    """

    def __init__(self, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 alpha: float | None = None, gradient: Any = None, **settings: Any) -> None:
        arrays = _bounds(initial_values, lower_bounds, upper_bounds)
        net = _opt().GradientDescent(_vector(objective), int(arrays[0].Length), *arrays)
        net.Gradient = _gradient(gradient)
        super().__init__(net, (objective, gradient))
        self._configure({"alpha": alpha, **settings})

    alpha = _setting("Alpha", float, "Learning rate (``Alpha``).")


class ADAM(_GradientSearch):
    """Bounded Adam optimizer (``ADAM``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    initial_values, lower_bounds, upper_bounds : array_like
        Start point and bounds.
    alpha : float, optional
        Learning rate; None keeps the Numerics default.
    gradient : callable, str or System.Delegate, optional
        Analytic gradient; None, the default, uses finite differences.
    **settings
        Initial values of any settings property, for example ``beta1``.
    """

    def __init__(self, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 alpha: float | None = None, gradient: Any = None, **settings: Any) -> None:
        arrays = _bounds(initial_values, lower_bounds, upper_bounds)
        net = _opt().ADAM(_vector(objective), int(arrays[0].Length), *arrays)
        net.Gradient = _gradient(gradient)
        super().__init__(net, (objective, gradient))
        self._configure({"alpha": alpha, **settings})

    alpha = _setting("Alpha", float, "Learning rate (``Alpha``).")
    beta1 = _setting("Beta1", float, "First-moment decay rate (``Beta1``).")
    beta2 = _setting("Beta2", float, "Second-moment decay rate (``Beta2``).")


class _GlobalSearch(_BoundedSearch):
    """The seed and polishing settings of the global optimizers."""

    prng_seed = _setting("PRNGSeed", int, "Random seed (``PRNGSeed``).")
    polish = _setting("Polish", bool, "Whether a local search polishes the best point (``Polish``).")
    polish_method = _method_setting("PolishMethod", "The local method of the polishing step (``PolishMethod``).")


class _PopulationSearch(_GlobalSearch):
    """Population size and parallel evaluation of the population optimizers."""

    population_size = _setting("PopulationSize", int, "Population size (``PopulationSize``).")
    max_degree_of_parallelism = _setting(
        "MaxDegreeOfParallelism", int,
        "Concurrent evaluations when parallel, or -1 for no limit (``MaxDegreeOfParallelism``).")

    @property
    def use_parallel(self) -> bool:
        """Whether each generation is evaluated on several threads (``UseParallel``; Numerics default False).

        Parallel evaluation runs the deferred (synchronous) variant of the search, whose result differs
        from the serial search but not with the thread count. It stays off with a Python objective, which
        serializes on the GIL.
        """
        return bool(self.net.UseParallel)

    @use_parallel.setter
    def use_parallel(self, value: bool) -> None:
        self.net.UseParallel = bool(value)
        callbacks.apply_parallel_policy(self.net, *self._functions)


class DifferentialEvolution(_PopulationSearch):
    """Differential evolution (``DifferentialEvolution``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    lower_bounds, upper_bounds : array_like
        Parameter bounds.
    **settings
        Initial values of any settings property, for example ``population_size`` or ``use_parallel``.
    """

    def __init__(self, objective: Any, lower_bounds: Any, upper_bounds: Any, **settings: Any) -> None:
        arrays = _bounds(lower_bounds, upper_bounds)
        super().__init__(_opt().DifferentialEvolution(_vector(objective), int(arrays[0].Length), *arrays),
                         (objective,))
        self._configure(settings)

    mutation = _setting("Mutation", float, "Mutation (differential weight) factor (``Mutation``).")
    dither_rate = _setting("DitherRate", float, "Mutation dithering rate (``DitherRate``).")
    crossover_probability = _setting("CrossoverProbability", float,
                                     "Crossover probability (``CrossoverProbability``).")

    @property
    def initial_population(self) -> np.ndarray | None:
        """Optional starting population, one row per member (``InitialPopulation``); None draws a Latin hypercube."""
        population = self.net.InitialPopulation
        return None if population is None else np.vstack([to_numpy_doubles(row) for row in population])

    @initial_population.setter
    def initial_population(self, rows: Any) -> None:
        if rows is None:
            self.net.InitialPopulation = None
            return
        import System  # noqa: PLC0415 - the runtime is running
        from System.Collections.Generic import List  # noqa: PLC0415

        population = List[System.Array[System.Double]]()
        for row in np.atleast_2d(np.asarray(rows, dtype=np.float64)):
            population.Add(to_net_doubles(row))
        self.net.InitialPopulation = population

    @property
    def integer_parameters(self) -> np.ndarray | None:
        """Optional flags marking integer-valued parameters (``IntegerParameters``)."""
        flags = self.net.IntegerParameters
        return None if flags is None else np.array([bool(f) for f in flags])

    @integer_parameters.setter
    def integer_parameters(self, flags: Any) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        self.net.IntegerParameters = None if flags is None else System.Array[System.Boolean](
            [bool(f) for f in flags])


class ParticleSwarm(_PopulationSearch):
    """Particle swarm optimization (``ParticleSwarm``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    lower_bounds, upper_bounds : array_like
        Parameter bounds.
    **settings
        Initial values of any settings property, for example ``population_size``.
    """

    def __init__(self, objective: Any, lower_bounds: Any, upper_bounds: Any, **settings: Any) -> None:
        arrays = _bounds(lower_bounds, upper_bounds)
        super().__init__(_opt().ParticleSwarm(_vector(objective), int(arrays[0].Length), *arrays), (objective,))
        self._configure(settings)


class _StartedGlobalSearch(_GlobalSearch):
    """The optional start point of CMA-ES and DDS."""

    initial_values = _array("InitialValues", "The optional starting parameters (``InitialValues``).")


def _started(factory: Any, objective: Any, lower_bounds: Any, upper_bounds: Any, initial_values: Any) -> Any:
    """Build a global optimizer with or without its optional start point."""
    lower, upper = _bounds(lower_bounds, upper_bounds)
    if initial_values is None:
        return factory(_vector(objective), int(lower.Length), lower, upper)
    return factory(_vector(objective), int(lower.Length), to_net_doubles(initial_values), lower, upper)


class CMAES(_StartedGlobalSearch):
    """Covariance matrix adaptation evolution strategy with restarts (``CMAES``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    lower_bounds, upper_bounds : array_like
        Parameter bounds.
    initial_values : array_like, optional
        Start point; None lets Numerics choose it.
    **settings
        Initial values of any settings property, for example ``restarts``.
    """

    def __init__(self, objective: Any, lower_bounds: Any, upper_bounds: Any, initial_values: Any = None,
                 **settings: Any) -> None:
        super().__init__(_started(_opt().CMAES, objective, lower_bounds, upper_bounds, initial_values), (objective,))
        self._configure(settings)

    population_size = _setting("PopulationSize", int, "Offspring per generation (``PopulationSize``).")
    initial_step_size = _setting("InitialStepSize", float,
                                 "Initial step size as a fraction of the range (``InitialStepSize``).")
    restarts = _setting("Restarts", int, "Maximum number of restarts (``Restarts``).")

    @property
    def restarts_performed(self) -> int:
        """Restarts performed by the last run (``RestartsPerformed``)."""
        return int(self.net.RestartsPerformed)


class DynamicallyDimensionedSearch(_StartedGlobalSearch):
    """Dynamically dimensioned search (``DynamicallyDimensionedSearch``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    lower_bounds, upper_bounds : array_like
        Parameter bounds.
    initial_values : array_like, optional
        Start point; None lets Numerics choose it.
    **settings
        Initial values of any settings property, for example ``evaluation_budget``.
    """

    def __init__(self, objective: Any, lower_bounds: Any, upper_bounds: Any, initial_values: Any = None,
                 **settings: Any) -> None:
        super().__init__(_started(_opt().DynamicallyDimensionedSearch, objective, lower_bounds, upper_bounds,
                                  initial_values), (objective,))
        self._configure(settings)

    evaluation_budget = _setting("EvaluationBudget", int, "Objective evaluation budget (``EvaluationBudget``).")
    perturbation_size = _setting("PerturbationSize", float,
                                 "Perturbation size as a fraction of the range (``PerturbationSize``).")


class ShuffledComplexEvolution(_GlobalSearch):
    """Shuffled complex evolution, SCE-UA (``ShuffledComplexEvolution``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    lower_bounds, upper_bounds : array_like
        Parameter bounds.
    **settings
        Initial values of any settings property, for example ``complexes``.
    """

    def __init__(self, objective: Any, lower_bounds: Any, upper_bounds: Any, **settings: Any) -> None:
        arrays = _bounds(lower_bounds, upper_bounds)
        super().__init__(_opt().ShuffledComplexEvolution(_vector(objective), int(arrays[0].Length), *arrays),
                         (objective,))
        self._configure(settings)

    complexes = _setting("Complexes", int, "Number of complexes (``Complexes``).")
    cce_iterations = _setting("CCEIterations", int,
                              "Competitive complex evolution steps per shuffle (``CCEIterations``).")
    tolerance_steps = _setting("ToleranceSteps", int,
                               "Shuffles without improvement before stopping (``ToleranceSteps``).")


class SimulatedAnnealing(_GlobalSearch):
    """Simulated annealing (``SimulatedAnnealing``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    lower_bounds, upper_bounds : array_like
        Parameter bounds.
    **settings
        Initial values of any settings property, for example ``initial_temperature``.
    """

    def __init__(self, objective: Any, lower_bounds: Any, upper_bounds: Any, **settings: Any) -> None:
        arrays = _bounds(lower_bounds, upper_bounds)
        super().__init__(_opt().SimulatedAnnealing(_vector(objective), int(arrays[0].Length), *arrays), (objective,))
        self._configure(settings)

    initial_temperature = _setting("InitialTemperature", float, "Starting temperature (``InitialTemperature``).")
    min_temperature = _setting("MinTemperature", float, "Final temperature (``MinTemperature``).")
    cooling_rate = _setting("CoolingRate", float, "Temperature reduction factor (``CoolingRate``).")
    update_cycles = _setting("UpdateCycles", int, "Step-size updates per temperature (``UpdateCycles``).")
    temperature_cycles = _setting("TemperatureCycles", int, "Moves per step-size update (``TemperatureCycles``).")
    tolerance_steps = _setting("ToleranceSteps", int,
                               "Temperatures without improvement before stopping (``ToleranceSteps``).")


class _MultiStartSearch(_LocalSearch):
    """The local method and tolerances of the multi-start optimizers."""

    prng_seed = _setting("PRNGSeed", int, "Random seed (``PRNGSeed``).")
    method = _method_setting("Method", "The local search method (``Method``).")
    local_absolute_tolerance = _setting("LocalAbsoluteTolerance", float,
                                        "Absolute tolerance of each local search (``LocalAbsoluteTolerance``).")
    local_relative_tolerance = _setting("LocalRelativeTolerance", float,
                                        "Relative tolerance of each local search (``LocalRelativeTolerance``).")
    polish = _setting("Polish", bool, "Whether the best point is polished at the end (``Polish``).")


def _multi_start(factory: Any, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 method: Any) -> Any:
    """Build MultiStart or MLSL, keeping the Numerics default method when ``method`` is None."""
    arrays = _bounds(initial_values, lower_bounds, upper_bounds)
    if method is None:
        return factory(_vector(objective), int(arrays[0].Length), *arrays)
    return factory(_vector(objective), int(arrays[0].Length), *arrays,
                   net_enum(_opt().LocalMethod, LocalMethod(method).value))


class MultiStart(_MultiStartSearch):
    """Multi-start local search (``MultiStart``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    initial_values, lower_bounds, upper_bounds : array_like
        First start point and bounds.
    method : LocalMethod or str, optional
        The local method; None keeps the Numerics default.
    **settings
        Initial values of any settings property; ``max_iterations`` is the number of starts.
    """

    def __init__(self, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 method: Any = None, **settings: Any) -> None:
        super().__init__(_multi_start(_opt().MultiStart, objective, initial_values, lower_bounds, upper_bounds,
                                      method), (objective,))
        self._configure(settings)


@dataclass(frozen=True, eq=False, slots=True)
class SamplePoint:
    """A sampled point of MLSL (``MLSL.SamplePoint``).

    Parameters
    ----------
    parameter_set : ParameterSet
        The point and its fitness.
    minimized : bool, optional
        Whether a local search started from it; default False.
    """

    parameter_set: ParameterSet
    minimized: bool = False

    @classmethod
    def from_net(cls, net: Any) -> SamplePoint:
        """Copy a .NET ``MLSL.SamplePoint``.

        Parameters
        ----------
        net : Numerics.Mathematics.Optimization.MLSL.SamplePoint
            The .NET point.

        Returns
        -------
        SamplePoint
            The copy.
        """
        return cls(ParameterSet.from_net(net.ParameterSet), bool(net.Minimized))


class MLSL(_MultiStartSearch):
    """Multi-level single linkage global search (``MLSL``).

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    initial_values, lower_bounds, upper_bounds : array_like
        First start point and bounds.
    method : LocalMethod or str, optional
        The local method; None keeps the Numerics default.
    **settings
        Initial values of any settings property, for example ``sample_size``.
    """

    def __init__(self, objective: Any, initial_values: Any, lower_bounds: Any, upper_bounds: Any,
                 method: Any = None, **settings: Any) -> None:
        super().__init__(_multi_start(_opt().MLSL, objective, initial_values, lower_bounds, upper_bounds, method),
                         (objective,))
        self._configure(settings)

    sample_size = _setting("SampleSize", int, "Points sampled per iteration (``SampleSize``).")
    gamma = _setting("Gamma", float, "Fraction of the best samples that may start local searches (``Gamma``).")
    sigma = _setting("Sigma", float, "Critical-distance scale (``Sigma``).")
    min_no_improvement = _setting("MinNoImprovement", int,
                                  "Iterations without improvement before stopping may begin (``MinNoImprovement``).")
    max_no_improvement = _setting("MaxNoImprovement", int,
                                  "Iterations without improvement that stop the search (``MaxNoImprovement``).")

    @property
    def sampled_points(self) -> list[SamplePoint]:
        """The points sampled by the last run (``SampledPoints``)."""
        points = self.net.SampledPoints
        return [] if points is None else [SamplePoint.from_net(p) for p in points]

    @property
    def local_minimums(self) -> list[SamplePoint]:
        """The local minima found by the last run (``LocalMinimums``)."""
        points = self.net.LocalMinimums
        return [] if points is None else [SamplePoint.from_net(p) for p in points]


class Constraint:
    """A constraint ``g(x) relation value`` for :class:`AugmentedLagrange` (``Constraint``).

    Parameters
    ----------
    function : callable, str or System.Delegate
        The constraint function g of a float64 array.
    number_of_parameters : int
        The number of parameters.
    value : float
        The right-hand side.
    type : ConstraintType or str
        The relation: g(x) = value, g(x) >= value or g(x) <= value.
    tolerance : float, optional
        The feasibility tolerance; None keeps the Numerics default.

    Attributes
    ----------
    net : Numerics.Mathematics.Optimization.Constraint
        The wrapped .NET constraint.
    """

    def __init__(self, function: Any, number_of_parameters: int, value: float, type: Any,  # noqa: A002
                 tolerance: float | None = None) -> None:
        opt = _opt()
        relation = net_enum(opt.ConstraintType, ConstraintType(type).value)
        arguments = [_vector(function), int(number_of_parameters), float(value), relation]
        if tolerance is not None:
            arguments.append(float(tolerance))
        self.net = opt.Constraint(*arguments)
        self._python_function = function

    @property
    def type(self) -> ConstraintType:
        """The relation (``Type``)."""
        return ConstraintType(str(self.net.Type))

    @property
    def number_of_parameters(self) -> int:
        """The number of parameters (``NumberOfParameters``)."""
        return int(self.net.NumberOfParameters)

    @property
    def function(self) -> Any:
        """The constraint function as the .NET ``Func<double[], double>`` (``Function``)."""
        return self.net.Function

    @property
    def value(self) -> float:
        """The right-hand side (``Value``)."""
        return float(self.net.Value)

    @property
    def tolerance(self) -> float:
        """The feasibility tolerance (``Tolerance``)."""
        return float(self.net.Tolerance)


class AugmentedLagrange(Optimizer):
    """Augmented Lagrangian constrained optimization around an inner optimizer (``AugmentedLagrange``).

    Numerics replaces the inner optimizer's objective with the augmented Lagrangian of ``objective``
    and the constraints; its bounds, start point and settings are kept. The parallel policy is applied
    to the inner optimizer for ``objective`` and every constraint function.

    Parameters
    ----------
    objective : callable, str or System.Delegate
        Objective of a float64 array.
    optimizer : Optimizer
        The inner optimizer, for example :class:`BFGS` or :class:`DifferentialEvolution`.
    constraints : sequence of Constraint
        At least one constraint.
    **settings
        Initial values of any settings property.
    """

    def __init__(self, objective: Any, optimizer: Optimizer, constraints: Sequence[Constraint],
                 **settings: Any) -> None:
        opt = _opt()
        from System.Collections.Generic import List  # noqa: PLC0415 - the runtime is running

        items = List[opt.IConstraint]()
        for constraint in constraints:
            items.Add(constraint.net)
        functions = (objective, *(c._python_function for c in constraints))
        net = opt.AugmentedLagrange(_vector(objective), optimizer.net, items)
        optimizer._functions = functions
        callbacks.apply_parallel_policy(optimizer.net, *functions)
        self._inner = optimizer
        self._constraints = tuple(constraints)
        super().__init__(net, functions)
        self._configure(settings)

    @property
    def optimizer(self) -> Optimizer:
        """The inner optimizer (``Optimizer``)."""
        return self._inner

    @property
    def constraints(self) -> tuple[Constraint, ...]:
        """The constraints (``Constraints``)."""
        return self._constraints

    @property
    def lambda_(self) -> np.ndarray:
        """The equality multipliers of the last run (``Lambda``)."""
        return np.array([float(v) for v in self.net.Lambda])

    @property
    def mu(self) -> np.ndarray:
        """The less-than inequality multipliers of the last run (``Mu``)."""
        return np.array([float(v) for v in self.net.Mu])

    @property
    def nu(self) -> np.ndarray:
        """The greater-than inequality multipliers of the last run (``Nu``)."""
        return np.array([float(v) for v in self.net.Nu])


class ParetoFront:
    """Pareto-front utilities for minimization (``Numerics.Mathematics.Optimization.ParetoFront``).

    Every function takes points as a two-dimensional array with one row of objective values per
    point; negate an objective to maximize it.
    """

    @staticmethod
    def dominates(a: Any, b: Any) -> bool:
        """Return whether point ``a`` Pareto dominates point ``b`` (``ParetoFront.Dominates``).

        Parameters
        ----------
        a, b : array_like
            Objective values of the two points.

        Returns
        -------
        bool
            True when ``a`` is no worse in every objective and better in one.
        """
        return bool(_opt().ParetoFront.Dominates(to_net_doubles(a), to_net_doubles(b)))

    @staticmethod
    def non_dominated_sort(points: Any) -> list[np.ndarray]:
        """Sort points into non-dominated fronts (``ParetoFront.NonDominatedSort``).

        Parameters
        ----------
        points : array_like
            One row of objective values per point.

        Returns
        -------
        list of numpy.ndarray
            The fronts, front 0 first, each the ascending int32 point indices.
        """
        fronts = _opt().ParetoFront.NonDominatedSort(to_net_matrix(points))
        return [to_numpy_ints(front.ToArray()) for front in fronts]

    @staticmethod
    def non_dominated_indices(points: Any) -> np.ndarray:
        """Return the indices of the non-dominated points (``ParetoFront.NonDominatedIndices``).

        Parameters
        ----------
        points : array_like
            One row of objective values per point.

        Returns
        -------
        numpy.ndarray
            The ascending int32 indices of front 0.
        """
        return to_numpy_ints(_opt().ParetoFront.NonDominatedIndices(to_net_matrix(points)))

    @staticmethod
    def crowding_distance(points: Any, front_indices: Any = None) -> np.ndarray:
        """Return the NSGA-II crowding distance of the points of one front (``ParetoFront.CrowdingDistance``).

        Parameters
        ----------
        points : array_like
            The front, or with ``front_indices`` the whole population; one row per point.
        front_indices : array_like of int, optional
            Indices of the front members in ``points``.

        Returns
        -------
        numpy.ndarray
            One distance per front member; boundary points are infinite.
        """
        pareto = _opt().ParetoFront
        if front_indices is None:
            return to_numpy_doubles(pareto.CrowdingDistance(to_net_matrix(points)))
        import System  # noqa: PLC0415 - the runtime is running

        rows = np.atleast_2d(np.asarray(points, dtype=np.float64))
        jagged = System.Array[System.Array[System.Double]]([to_net_doubles(r) for r in rows])
        return to_numpy_doubles(pareto.CrowdingDistance(jagged, to_net_ints(front_indices)))

    @staticmethod
    def hypervolume(points: Any, reference_point: Any) -> float:
        """Return the exact hypervolume the points dominate up to a reference point (``ParetoFront.Hypervolume``).

        Parameters
        ----------
        points : array_like
            One row of objective values per point.
        reference_point : array_like
            The reference point.

        Returns
        -------
        float
            The hypervolume.
        """
        return float(_opt().ParetoFront.Hypervolume(to_net_matrix(points), to_net_doubles(reference_point)))

    @staticmethod
    def generational_distance(approximation: Any, reference_front: Any) -> float:
        """Return the mean distance from the approximation to the reference front (``GenerationalDistance``).

        Parameters
        ----------
        approximation, reference_front : array_like
            One row of objective values per point.

        Returns
        -------
        float
            The generational distance (pymoo ``GD`` without normalization).
        """
        return float(_opt().ParetoFront.GenerationalDistance(to_net_matrix(approximation),
                                                              to_net_matrix(reference_front)))

    @staticmethod
    def inverted_generational_distance(approximation: Any, reference_front: Any) -> float:
        """Return the mean distance from the reference front to the approximation (``InvertedGenerationalDistance``).

        Parameters
        ----------
        approximation, reference_front : array_like
            One row of objective values per point.

        Returns
        -------
        float
            The inverted generational distance (pymoo ``IGD`` without normalization).
        """
        return float(_opt().ParetoFront.InvertedGenerationalDistance(to_net_matrix(approximation),
                                                                      to_net_matrix(reference_front)))


def _ints_matrix_to_numpy(table: Any) -> np.ndarray:
    """Copy a .NET ``int[,]`` into a two-dimensional int32 array in one copy."""
    return _numpy_rectangular(table, np.int32)


def _ints_matrix_to_net(values: Any) -> Any:
    """Copy a two-dimensional integer array into a new .NET ``int[,]``, refusing values outside the int32 range."""
    array = np.asarray(values)
    if array.ndim != 2:
        raise ValueError(f"Expected a 2-D weight array; got shape {array.shape}.")
    return _net_rectangular(_int32_array(array))


def knapsack(values: Any, weights: Any, capacity: Any) -> tuple[float, np.ndarray]:
    """Solve the exact 0/1 knapsack by dynamic programming (``Knapsack.Solve``).

    Parameters
    ----------
    values : array_like
        The value of each item.
    weights : array_like of int
        The non-negative integer weight of each item, or a ``(items, constraints)`` array for the
        multi-constraint problem.
    capacity : int or array_like of int
        The capacity, or one capacity per constraint.

    Returns
    -------
    tuple
        ``(total_value, selected_items)``, the selected indices ascending as int32.

    Raises
    ------
    ValueError
        If a weight or capacity is a floating-point value that is not a finite integer, or if Numerics
        refuses the problem (for example a negative weight).
    """
    for name, given in (("weights", weights), ("capacity", capacity)):
        array = np.asarray(given)
        if array.dtype.kind == "f" and not np.all(np.isfinite(array) & (array == np.floor(array))):
            raise ValueError(f"The knapsack {name} must be integers.")
    solver = _opt().Knapsack
    if np.ndim(weights) == 2:
        total, selected = solver.Solve(to_net_doubles(values), _ints_matrix_to_net(weights),
                                       to_net_ints(np.atleast_1d(capacity)), None)
    else:
        total, selected = solver.Solve(to_net_doubles(values), to_net_ints(weights), int(capacity), None)
    return float(total), to_numpy_ints(selected)


@dataclass(frozen=True, eq=False, slots=True)
class DynamicProgrammingPath:
    """A forward simulation of a dynamic programming policy (``DynamicProgrammingPath``).

    Parameters
    ----------
    states : numpy.ndarray
        The state at each stage, T + 1 entries.
    decisions : numpy.ndarray
        The decision at each stage, T entries.
    stage_values : numpy.ndarray
        The stage cost or reward at each stage.
    total_value : float
        The total value, terminal value included.
    """

    states: np.ndarray
    decisions: np.ndarray
    stage_values: np.ndarray
    total_value: float

    @classmethod
    def from_net(cls, net: Any) -> DynamicProgrammingPath:
        """Copy a .NET ``DynamicProgrammingPath``.

        Parameters
        ----------
        net : Numerics.Mathematics.Optimization.DynamicProgrammingPath
            The .NET path.

        Returns
        -------
        DynamicProgrammingPath
            The copy.
        """
        return cls(to_numpy_ints(net.States), to_numpy_ints(net.Decisions), to_numpy_doubles(net.StageValues),
                   float(net.TotalValue))


class DynamicProgramming:
    """Discrete dynamic programming over integer stages, states and decisions (``DynamicProgramming``).

    Parameters
    ----------
    stages, states, decisions : int
        The numbers of stages T, states S and decisions A.
    stage_value : callable, str or System.Delegate
        The stage cost or reward ``c(t, s, a)`` (``Func<int, int, int, double>``); a non-finite value
        marks the decision infeasible.
    transition : callable, str or System.Delegate
        Without ``max_outcomes``, the next state ``f(t, s, a)`` (``Func<int, int, int, int>``; -1 is
        infeasible). With ``max_outcomes``, a ``TransitionDistribution``
        ``(t, s, a, next_states, probabilities) -> count`` that fills the two .NET arrays and returns
        the number of outcomes written.
    max_outcomes : int, optional
        The largest number of outcomes of a stochastic transition; None builds a deterministic problem.

    Attributes
    ----------
    net : Numerics.Mathematics.Optimization.DynamicProgramming
        The wrapped .NET problem.
    """

    def __init__(self, stages: int, states: int, decisions: int, stage_value: Any, transition: Any,
                 max_outcomes: int | None = None) -> None:
        opt = _opt()
        value = callbacks.to_delegate(stage_value, "Func<int,int,int,double>")
        if max_outcomes is None:
            step = callbacks.to_delegate(transition, "Func<int,int,int,int>")
            self.net = opt.DynamicProgramming(int(stages), int(states), int(decisions), value, step)
        else:
            step = callbacks.to_delegate(transition, "TransitionDistribution")
            self.net = opt.DynamicProgramming(int(stages), int(states), int(decisions), value, step, int(max_outcomes))
        self._functions = (stage_value, transition)
        callbacks.apply_parallel_policy(self.net, *self._functions)

    @property
    def stages(self) -> int:
        """The number of stages T (``Stages``)."""
        return int(self.net.Stages)

    @property
    def states(self) -> int:
        """The number of states S (``States``)."""
        return int(self.net.States)

    @property
    def decisions(self) -> int:
        """The number of decisions A (``Decisions``)."""
        return int(self.net.Decisions)

    @property
    def max_outcomes(self) -> int:
        """The largest number of outcomes of a transition (``MaxOutcomes``)."""
        return int(self.net.MaxOutcomes)

    @property
    def is_stochastic(self) -> bool:
        """Whether the transitions are stochastic (``IsStochastic``)."""
        return bool(self.net.IsStochastic)

    @property
    def maximize(self) -> bool:
        """Whether the stage values are rewards to maximize (``Maximize``; default False)."""
        return bool(self.net.Maximize)

    @maximize.setter
    def maximize(self, value: bool) -> None:
        self.net.Maximize = bool(value)

    @property
    def discount_factor(self) -> float:
        """The discount factor in (0, 1] (``DiscountFactor``; default 1)."""
        return float(self.net.DiscountFactor)

    @discount_factor.setter
    def discount_factor(self, value: float) -> None:
        self.net.DiscountFactor = float(value)

    @property
    def terminal_value(self) -> np.ndarray | None:
        """The terminal value of each state, or None for zero (``TerminalValue``)."""
        values = self.net.TerminalValue
        return None if values is None else to_numpy_doubles(values)

    @terminal_value.setter
    def terminal_value(self, values: Any) -> None:
        self.net.TerminalValue = None if values is None else to_net_doubles(values)

    @property
    def use_parallel(self) -> bool:
        """Whether the states of each stage are solved in parallel (``UseParallel``; default False).

        The result equals the serial solve. It stays off with a Python callback.
        """
        return bool(self.net.UseParallel)

    @use_parallel.setter
    def use_parallel(self, value: bool) -> None:
        self.net.UseParallel = bool(value)
        callbacks.apply_parallel_policy(self.net, *self._functions)

    @property
    def is_solved(self) -> bool:
        """Whether :meth:`solve` has completed (``IsSolved``)."""
        return bool(self.net.IsSolved)

    @property
    def value(self) -> np.ndarray:
        """The optimal value table ``V[t, s]`` for t = 0..T; row T is the terminal value (``Value``)."""
        return to_numpy_matrix(self.net.Value)

    @property
    def policy(self) -> np.ndarray:
        """The optimal decision table ``[t, s]`` for t = 0..T-1, -1 where infeasible (``Policy``)."""
        return _ints_matrix_to_numpy(self.net.Policy)

    @property
    def steady_state_value(self) -> np.ndarray:
        """The stationary value of each state (``SteadyStateValue``)."""
        return to_numpy_doubles(self.net.SteadyStateValue)

    @property
    def steady_state_policy(self) -> np.ndarray:
        """The stationary decision of each state, -1 where infeasible (``SteadyStatePolicy``)."""
        return to_numpy_ints(self.net.SteadyStatePolicy)

    @property
    def steady_state_iterations(self) -> int:
        """Value-iteration sweeps of the last steady-state solve (``SteadyStateIterations``)."""
        return int(self.net.SteadyStateIterations)

    @property
    def steady_state_converged(self) -> bool:
        """Whether the last steady-state solve converged (``SteadyStateConverged``)."""
        return bool(self.net.SteadyStateConverged)

    def solve(self) -> None:
        """Solve the finite-horizon problem by backward recursion (``Solve``)."""
        self.net.Solve()

    def solve_steady_state(self, tolerance: float = 1e-8, max_iterations: int = 10000) -> bool:
        """Solve the stationary infinite-horizon problem by value iteration (``SolveSteadyState``).

        Parameters
        ----------
        tolerance : float, optional
            Convergence tolerance on the largest value change; default 1e-8, as in Numerics.
        max_iterations : int, optional
            Sweep limit; default 10,000, as in Numerics.

        Returns
        -------
        bool
            True when the iteration converged.
        """
        return bool(self.net.SolveSteadyState(float(tolerance), int(max_iterations)))

    def simulate(self, initial_state: int, outcomes: Sequence[int] | None = None,
                 seed: int | None = None) -> DynamicProgrammingPath:
        """Simulate the optimal policy forward from a state (``Simulate``).

        Parameters
        ----------
        initial_state : int
            The state at stage 0.
        outcomes : sequence of int, optional
            The outcome index at each stage (for example an observed inflow class sequence).
        seed : int, optional
            Seed of the ``System.Random`` that samples each stochastic outcome.

        Returns
        -------
        DynamicProgrammingPath
            The state and decision path. Without ``outcomes`` or ``seed`` the problem must be
            deterministic.

        Raises
        ------
        ValueError
            If both ``outcomes`` and ``seed`` are given.
        """
        if outcomes is not None and seed is not None:
            raise ValueError("Give outcomes or seed, not both.")
        if outcomes is not None:
            path = self.net.Simulate(int(initial_state), to_net_ints(outcomes))
        elif seed is not None:
            import System  # noqa: PLC0415 - the runtime is running

            path = self.net.Simulate(int(initial_state), System.Random(int(seed)))
        else:
            path = self.net.Simulate(int(initial_state))
        return DynamicProgrammingPath.from_net(path)


__all__ = [
    "ADAM", "BFGS", "BULK_TRACE_LENGTH", "CMAES", "LBFGSB", "MLSL", "AugmentedLagrange", "BrentSearch",
    "Constraint", "ConstraintType", "DifferentialEvolution", "DynamicProgramming", "DynamicProgrammingPath",
    "DynamicallyDimensionedSearch", "GoldenSection", "GradientDescent", "LocalMethod", "MultiStart", "NelderMead",
    "OptimizationProgress", "OptimizationResult", "OptimizationStatus", "Optimizer", "ParameterSet", "ParetoFront",
    "ParticleSwarm", "Powell", "SamplePoint", "ShuffledComplexEvolution", "SimulatedAnnealing", "knapsack",
]

translate_public(globals())
for _base in (_BoundedSearch, _LocalSearch, _GradientSearch, _ScalarSearch, _GlobalSearch, _PopulationSearch,
              _StartedGlobalSearch, _MultiStartSearch):
    translate_class(_base)
