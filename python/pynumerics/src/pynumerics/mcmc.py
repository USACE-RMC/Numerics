"""Bayesian sampling from ``Numerics.Sampling.MCMC``: the MCMC samplers, their results and diagnostics.

Every sampler takes a list of prior distributions and a log-likelihood. The priors bound the
parameters and seed the chains; the log-likelihood is the full log-posterior up to a constant
(Numerics adds no prior density), so uniform priors need no prior term and other priors need
their log-densities added (:func:`distribution_log_likelihood` does this with
``add_log_prior=True``). The log-likelihood is any of:

- a Python callable of the parameter array, called once per evaluation (about 1 microsecond of
  crossing each); the wrapper then sets ``ParallelizeChains = False``, because the chains would
  only contend for the GIL;
- C# source compiled with Roslyn (:func:`pynumerics.callbacks.compile_csharp`), or the in-.NET
  model of :func:`distribution_log_likelihood`: no Python runs per evaluation and the chains
  keep running in parallel;
- a ``Numerics.Sampling.MCMC.LogLikelihood`` delegate.

Runs are reproducible: every sampler draws from Mersenne Twister streams seeded by ``seed``
(``PRNGSeed``, default 12345), one per chain, so a seeded run returns the same chains whether
or not the chains run in parallel.

Chains and outputs come back as NumPy arrays: ``markov_chains`` has shape
``(chains, iterations, parameters)``. The chains are read from .NET in one pass per run (about
3 microseconds per recorded state) and cached, so later reads cost only a NumPy copy.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import (
    _opt,
    bind_list,
    cached_overload,
    namespace,
    net_enum,
    unwrap,
)
from pynumerics.enums import NetEnum
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _linear_algebra,
    _numerics_matrix,
    _numerics_vector,
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_ints,
    to_numpy_matrix,
)
from pynumerics.sampling import _parameter_set


def _mcmc() -> Any:
    """Return the ``Numerics.Sampling.MCMC`` namespace after starting the runtime."""
    return namespace("Numerics.Sampling.MCMC")


class InitializationType(NetEnum):
    """How the chains start (``MCMCSampler.InitializationType``)."""

    MAP = "MAP"
    RANDOMIZE = "Randomize"
    USER_DEFINED = "UserDefined"


def _priors(priors: Sequence[Any]) -> Any:
    """Build the ``List<IUnivariateDistribution>`` of priors from pynumerics or .NET distributions."""
    nd = namespace("Numerics.Distributions")
    from System.Collections.Generic import List  # noqa: PLC0415 - the runtime is running

    if isinstance(priors, List[nd.IUnivariateDistribution]):
        return priors
    items = List[nd.IUnivariateDistribution]()
    for prior in priors:
        items.Add(getattr(prior, "net", prior))
    return items


def _point(parameter_set: Any) -> tuple[np.ndarray, float]:
    """Return a ``ParameterSet`` as ``(values, fitness)``; an empty set gives an empty array."""
    values = parameter_set.Values
    return (np.empty(0) if values is None else to_numpy_doubles(values)), float(parameter_set.Fitness)


def _read_sets(sets: Any, width: int, *, fitness: bool = True,
               weights: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a .NET list of ``ParameterSet`` into ``(values (n, p), fitness (n,), weight (n,))``.

    The list is enumerated once and the value arrays are joined by one ``numpy.concatenate`` over
    their buffers, so each set costs one crossing for the set and one for its values. The fitness
    or weight array is empty, and not read, when ``fitness`` or ``weights`` is False.
    """
    items = list(sets)
    count = len(items)
    values = np.concatenate([item.Values for item in items]) if count else np.empty(0)
    values = values.astype(np.float64, copy=False).reshape(count, width)
    scores = (np.fromiter((item.Fitness for item in items), dtype=np.float64, count=count) if fitness
              else np.empty(0))
    weight = np.fromiter((item.Weight for item in items), dtype=np.float64, count=count) if weights else np.empty(0)
    return values, scores, weight


def _chain_values(chains: Any, width: int) -> tuple[np.ndarray, np.ndarray]:
    """Read ``List<ParameterSet>[]`` into ``(values (c, n, p), fitness (c, n))`` over the common length."""
    read = [_read_sets(chain, width, weights=False) for chain in chains]
    if not read:
        return np.empty((0, 0, width)), np.empty((0, 0))
    length = min(r[0].shape[0] for r in read)
    return np.stack([r[0][:length] for r in read]), np.stack([r[1][:length] for r in read])


def _cached_chain_values(owner: Any, chains: Any, width: int) -> tuple[np.ndarray, np.ndarray]:
    """Return :func:`_chain_values` of ``chains``, read once and kept on ``owner`` until the chains change.

    The cache holds while ``owner.net`` is the same object, ``chains`` is the same .NET array and
    every chain has the same count; :meth:`MCMCSampler.sample` also clears it. Callers return copies,
    so the cached arrays are never handed out.
    """
    import System  # noqa: PLC0415 - the runtime is running

    counts = tuple(int(chain.Count) for chain in chains)
    cache = getattr(owner, "_chain_cache", None)
    if (cache is not None and cache[0] is owner.net and System.Object.ReferenceEquals(cache[1], chains)
            and cache[2] == counts):
        return cache[3], cache[4]
    values, fitness = _chain_values(chains, width)
    owner._chain_cache = (owner.net, chains, counts, values, fitness)
    return values, fitness


def _chains_to_net(chains: Any) -> Any:
    """Return Markov chains as ``List<ParameterSet>[]``.

    Accepts a sampler or results wrapper, a .NET chain array, or an array of shape
    ``(chains, iterations, parameters)`` (a 2-D array is one parameter per chain).
    """
    net = getattr(chains, "net", None)
    if net is not None:
        return net.MarkovChains
    ps = _opt().ParameterSet
    import System  # noqa: PLC0415 - the runtime is running
    from System.Collections.Generic import List  # noqa: PLC0415

    if isinstance(chains, System.Array):
        return chains
    array = np.asarray(chains, dtype=np.float64)
    if array.ndim == 2:
        array = array[:, :, None]
    if array.ndim != 3:
        raise ValueError(f"Expected chains of shape (chains, iterations, parameters); got {array.shape}.")
    result = System.Array.CreateInstance(List[ps], array.shape[0])
    nan = float("nan")
    for c in range(array.shape[0]):
        chain = List[ps](array.shape[1])
        for row in array[c]:
            chain.Add(ps(to_net_doubles(row), nan))
        result[c] = chain
    return result


_LIKELIHOOD_FACTORY = """
static LogLikelihood Create(UnivariateDistributionBase template, double[] data, IUnivariateDistribution[] priors,
                            int[] indices)
{
    var baseline = template.GetParameters;
    return parameters =>
    {
        double logPrior = 0d;
        if (priors != null)
        {
            for (int i = 0; i < priors.Length; i++)
                logPrior += priors[i].LogPDF(parameters[i]);
            if (double.IsNaN(logPrior) || double.IsNegativeInfinity(logPrior))
                return double.NegativeInfinity;
        }
        var values = (double[])baseline.Clone();
        for (int i = 0; i < indices.Length; i++)
            values[indices[i]] = parameters[i];
        var model = template.Clone();
        model.SetParameters(values);
        if (!model.ParametersValid)
            return double.NegativeInfinity;
        return logPrior + model.LogLikelihood(data);
    };
}
"""
"""C# factory of the in-.NET log-likelihood: a fresh clone per evaluation, so parallel chains share nothing."""


def distribution_log_likelihood(distribution: Any, data: Any, priors: Sequence[Any] | None = None,
                                parameter_indices: Sequence[int] | None = None, add_log_prior: bool = False) -> Any:
    """Return the log-likelihood of data under a Numerics distribution as a delegate that runs in .NET.

    The delegate sets the sampled parameters on a clone of ``distribution`` and returns
    ``LogLikelihood(data)`` (``UnivariateDistributionBase.LogLikelihood``), plus the priors'
    log-densities when ``add_log_prior`` is True. No Python runs per evaluation, so samplers
    given it keep ``ParallelizeChains`` on. The factory is compiled once per process with Roslyn
    (:func:`pynumerics.callbacks.compile_csharp`).

    Parameters
    ----------
    distribution : pynumerics.distributions.UnivariateDistribution or UnivariateDistributionBase
        The model; its current parameters are the values of the parameters not sampled.
    data : array_like
        The observations.
    priors : sequence of distributions, optional
        The priors of the sampled parameters, in order; used only with ``add_log_prior``.
    parameter_indices : sequence of int, optional
        The positions in ``GetParameters`` of the sampled parameters, in sampling order; default
        all parameters in their order. Parameters not listed keep the template's values (for
        example a known scale).
    add_log_prior : bool, optional
        Add ``sum(prior.LogPDF(theta_i))`` so the delegate is the full log-posterior; default
        False, which is right for uniform priors.

    Returns
    -------
    Numerics.Sampling.MCMC.LogLikelihood
        The delegate. Parameters the distribution rejects (``ParametersValid`` false) and draws of
        zero prior density return negative infinity.

    Raises
    ------
    ValueError
        If ``add_log_prior`` is True without priors, or the priors and indices differ in length.
    pynumerics.RuntimeLoadError
        If the bundled Roslyn assemblies are missing.
    """
    mc = _mcmc()
    nd = namespace("Numerics.Distributions")
    import System  # noqa: PLC0415 - the runtime is running

    template = unwrap(getattr(distribution, "net", distribution))
    count = int(template.NumberOfParameters)
    indices = list(range(count)) if parameter_indices is None else [int(i) for i in parameter_indices]
    prior_array = None
    if add_log_prior:
        if priors is None:
            raise ValueError("add_log_prior=True needs the priors.")
        prior_list = [getattr(p, "net", p) for p in priors]
        if len(prior_list) != len(indices):
            raise ValueError(f"{len(prior_list)} priors for {len(indices)} sampled parameters.")
        prior_array = System.Array[nd.IUnivariateDistribution](prior_list)
    factory_type = System.Func[nd.UnivariateDistributionBase, System.Array[System.Double],
                               System.Array[nd.IUnivariateDistribution], System.Array[System.Int32], mc.LogLikelihood]
    factory = callbacks.compile_csharp(_LIKELIHOOD_FACTORY, factory_type, usings=("Numerics.Sampling.MCMC",))
    flat = np.ravel(np.asarray(data, dtype=np.float64))
    return factory(template.Clone(), to_net_doubles(flat), prior_array, to_net_ints(indices))


def _likelihood(log_likelihood: Any) -> Any:
    """Return the ``LogLikelihood`` delegate of a callable, C# source or delegate."""
    return callbacks.log_likelihood(log_likelihood)


class MCMCSampler:
    """Base of the Numerics MCMC samplers (``Numerics.Sampling.MCMC.MCMCSampler``).

    Use a concrete sampler (:class:`RWMH`, :class:`ARWMH`, :class:`DEMCz`, :class:`DEMCzs`,
    :class:`HMC`, :class:`NUTS`, :class:`Gibbs`, :class:`SNIS`). The keyword settings are applied
    after the .NET sampler is built; None keeps the Numerics default. Setting the seed or a run
    length resets the sampler, as in Numerics.

    Parameters
    ----------
    net : Numerics.Sampling.MCMC.MCMCSampler
        The .NET sampler.
    functions : sequence
        The functions given to it, for the parallel rule of :mod:`pynumerics.callbacks`.
    seed : int, optional
        ``PRNGSeed``; Numerics default 12345.
    number_of_chains : int, optional
        ``NumberOfChains``.
    iterations : int, optional
        ``Iterations``, the recorded iterations per chain including warmup.
    warmup_iterations : int, optional
        ``WarmupIterations``.
    thinning_interval : int, optional
        ``ThinningInterval``.
    output_length : int, optional
        ``OutputLength``, the posterior draws output after the iterations.
    initial_iterations : int, optional
        ``InitialIterations``, the initial population size.
    initialize : InitializationType or str, optional
        ``Initialize``; Numerics default ``Randomize``.
    parallelize_chains : bool, optional
        ``ParallelizeChains``; Numerics default True. A Python callback always turns it off.

    Attributes
    ----------
    net : Numerics.Sampling.MCMC.MCMCSampler
        The wrapped .NET sampler.
    """

    def __init__(self, net: Any, functions: Sequence[Any] = (), *, seed: int | None = None,
                 number_of_chains: int | None = None, iterations: int | None = None,
                 warmup_iterations: int | None = None, thinning_interval: int | None = None,
                 output_length: int | None = None, initial_iterations: int | None = None,
                 initialize: InitializationType | str | None = None, parallelize_chains: bool | None = None) -> None:
        self.net = net
        self._handlers: list[Any] = []
        settings = (("PRNGSeed", seed), ("NumberOfChains", number_of_chains), ("Iterations", iterations),
                    ("WarmupIterations", warmup_iterations), ("ThinningInterval", thinning_interval),
                    ("OutputLength", output_length), ("InitialIterations", initial_iterations))
        for name, value in settings:
            if value is not None:
                setattr(net, name, int(value))
        if initialize is not None:
            self.initialize = initialize
        if parallelize_chains is not None:
            net.ParallelizeChains = bool(parallelize_chains)
        callbacks.apply_parallel_policy(net, *functions)

    @property
    def seed(self) -> int:
        """The PRNG seed (``PRNGSeed``); setting it resets the sampler."""
        return int(self.net.PRNGSeed)

    @seed.setter
    def seed(self, value: int) -> None:
        self.net.PRNGSeed = int(value)

    @property
    def initial_iterations(self) -> int:
        """The size of the initial population (``InitialIterations``); setting it resets the sampler."""
        return int(self.net.InitialIterations)

    @initial_iterations.setter
    def initial_iterations(self, value: int) -> None:
        self.net.InitialIterations = int(value)

    @property
    def warmup_iterations(self) -> int:
        """The warmup iterations, part of ``iterations`` (``WarmupIterations``); setting it resets the sampler."""
        return int(self.net.WarmupIterations)

    @warmup_iterations.setter
    def warmup_iterations(self, value: int) -> None:
        self.net.WarmupIterations = int(value)

    @property
    def iterations(self) -> int:
        """The recorded iterations per chain (``Iterations``); setting it resets the sampler."""
        return int(self.net.Iterations)

    @iterations.setter
    def iterations(self, value: int) -> None:
        self.net.Iterations = int(value)

    @property
    def number_of_chains(self) -> int:
        """The number of chains (``NumberOfChains``); setting it resets the sampler."""
        return int(self.net.NumberOfChains)

    @number_of_chains.setter
    def number_of_chains(self, value: int) -> None:
        self.net.NumberOfChains = int(value)

    @property
    def thinning_interval(self) -> int:
        """Transitions per recorded iteration (``ThinningInterval``); setting it resets the sampler."""
        return int(self.net.ThinningInterval)

    @thinning_interval.setter
    def thinning_interval(self, value: int) -> None:
        self.net.ThinningInterval = int(value)

    @property
    def output_length(self) -> int:
        """The number of posterior draws output after the iterations, over all chains (``OutputLength``)."""
        return int(self.net.OutputLength)

    @output_length.setter
    def output_length(self, value: int) -> None:
        self.net.OutputLength = int(value)

    @property
    def initialize(self) -> InitializationType:
        """How the chains start (``Initialize``)."""
        return InitializationType(str(self.net.Initialize))

    @initialize.setter
    def initialize(self, value: InitializationType | str) -> None:
        self.net.Initialize = net_enum(_mcmc().MCMCSampler.InitializationType, value)

    @property
    def parallelize_chains(self) -> bool:
        """Whether chains advance on several threads (``ParallelizeChains``)."""
        return bool(self.net.ParallelizeChains)

    @parallelize_chains.setter
    def parallelize_chains(self, value: bool) -> None:
        self.net.ParallelizeChains = bool(value)

    @property
    def resume_simulation(self) -> bool:
        """Whether another :meth:`sample` continues the chains instead of restarting (``ResumeSimulation``)."""
        return bool(self.net.ResumeSimulation)

    @resume_simulation.setter
    def resume_simulation(self, value: bool) -> None:
        self.net.ResumeSimulation = bool(value)

    @property
    def progress_changed_rate(self) -> float:
        """The fraction of the run between progress reports (``ProgressChangedRate``)."""
        return float(self.net.ProgressChangedRate)

    @progress_changed_rate.setter
    def progress_changed_rate(self, value: float) -> None:
        self.net.ProgressChangedRate = float(value)

    @property
    def cancellation_token_source(self) -> Any:
        """The ``System.Threading.CancellationTokenSource`` of the current run (``CancellationTokenSource``)."""
        return self.net.CancellationTokenSource

    @cancellation_token_source.setter
    def cancellation_token_source(self, value: Any) -> None:
        self.net.CancellationTokenSource = value

    @property
    def cancellation_token(self) -> Any:
        """An external ``System.Threading.CancellationToken`` that stops the run (``CancellationToken``)."""
        return self.net.CancellationToken

    @cancellation_token.setter
    def cancellation_token(self, value: Any) -> None:
        self.net.CancellationToken = value

    @property
    def log_likelihood_function(self) -> Any:
        """The ``LogLikelihood`` delegate the sampler evaluates (``LogLikelihoodFunction``)."""
        return self.net.LogLikelihoodFunction

    @property
    def prior_distributions(self) -> list[Any]:
        """The prior distributions as pynumerics wrappers (``PriorDistributions``)."""
        from pynumerics.distributions import wrap  # noqa: PLC0415 - keeps the import of this module light

        return [wrap(p) for p in self.net.PriorDistributions]

    @property
    def number_of_parameters(self) -> int:
        """The number of model parameters (``NumberOfParameters``)."""
        return int(self.net.NumberOfParameters)

    @property
    def is_population_sampler(self) -> bool:
        """Whether proposals draw on the population of past states (``IsPopulationSampler``)."""
        return bool(self.net.IsPopulationSampler)

    @property
    def transition_count(self) -> int:
        """Chain transitions one chain performs per :meth:`sample` (``TransitionCount``)."""
        return int(self.net.TransitionCount)

    @property
    def total_transition_count(self) -> int:
        """Chain transitions over all chains per :meth:`sample` (``TotalTransitionCount``)."""
        return int(self.net.TotalTransitionCount)

    @property
    def markov_chains(self) -> np.ndarray:
        """The recorded chains, shape ``(chains, iterations, parameters)`` (``MarkovChains``).

        The chains are read from .NET once per run and cached; each access returns a new copy.
        """
        return _cached_chain_values(self, self.net.MarkovChains, self.number_of_parameters)[0].copy()

    @property
    def markov_chain_log_likelihoods(self) -> np.ndarray:
        """The log-likelihood of each recorded state, shape ``(chains, iterations)`` (``MarkovChains``).

        Read with :attr:`markov_chains` and cached the same way; each access returns a new copy.
        """
        return _cached_chain_values(self, self.net.MarkovChains, self.number_of_parameters)[1].copy()

    @property
    def output(self) -> list[np.ndarray]:
        """The posterior draws of each chain, one ``(draws, parameters)`` array per chain (``Output``)."""
        width = self.number_of_parameters
        return [_read_sets(chain, width, fitness=False, weights=False)[0] for chain in self.net.Output]

    @property
    def population_matrix(self) -> np.ndarray:
        """The population of a population sampler, shape ``(states, parameters)`` (``PopulationMatrix``)."""
        return _read_sets(self.net.PopulationMatrix, self.number_of_parameters, fitness=False, weights=False)[0]

    @property
    def accept_count(self) -> np.ndarray:
        """Accepted proposals per chain (``AcceptCount``)."""
        return to_numpy_ints(self.net.AcceptCount).astype(np.int64)

    @property
    def sample_count(self) -> np.ndarray:
        """Proposals per chain (``SampleCount``)."""
        return to_numpy_ints(self.net.SampleCount).astype(np.int64)

    @property
    def acceptance_rates(self) -> np.ndarray:
        """Accepted over proposed per chain (``AcceptanceRates``)."""
        return to_numpy_doubles(self.net.AcceptanceRates)

    @property
    def mean_log_likelihood(self) -> np.ndarray:
        """The log-likelihood averaged over chains at each iteration (``MeanLogLikelihood``)."""
        return to_numpy_doubles(self.net.MeanLogLikelihood)

    @property
    def map(self) -> tuple[np.ndarray, float]:
        """The output draw of largest log-likelihood as ``(values, log_likelihood)`` (``MAP``)."""
        return _point(self.net.MAP)

    @property
    def map_initialization_failed(self) -> bool:
        """Whether MAP initialization failed and the run started from the priors (``MAPInitializationFailed``)."""
        return bool(self.net.MAPInitializationFailed)

    def sample(self) -> MCMCResults:
        """Run the sampler (``Sample``) and post-process the run (``MCMCResults(sampler)``).

        Returns
        -------
        MCMCResults
            The results at the default 90 percent intervals; :meth:`results` gives another level.
        """
        self._chain_cache = None
        self.net.Sample()
        return MCMCResults(self)

    def results(self, alpha: float = 0.1) -> MCMCResults:
        """Return the post-processed results of the last run (``MCMCResults(sampler, alpha)``).

        Parameters
        ----------
        alpha : float, optional
            Two-sided interval level; default 0.1 (90 percent intervals).

        Returns
        -------
        MCMCResults
            The results.
        """
        return MCMCResults(self, alpha)

    def reset(self) -> None:
        """Clear the chains and results and reseed the generators (``Reset``)."""
        self.net.Reset()

    def cancel_simulation(self) -> None:
        """Stop the running :meth:`sample` after the current iteration (``CancelSimulation``)."""
        self.net.CancelSimulation()

    def report_progress(self, percent_complete: float) -> None:
        """Raise the progress event (``ReportProgress``).

        Parameters
        ----------
        percent_complete : float
            Fraction complete in [0, 1].
        """
        self.net.ReportProgress(float(percent_complete))

    def add_progress_handler(self, handler: Any) -> Any:
        """Subscribe to progress reports (the ``ProgressChanged`` event).

        Parameters
        ----------
        handler : callable, str or System.Delegate
            ``handler(percent_complete, progress_text)``; called on the thread running :meth:`sample`.

        Returns
        -------
        Numerics.Sampling.MCMC.MCMCSampler.ProgressChangedEventHandler
            The subscribed delegate, for :meth:`remove_progress_handler`.
        """
        delegate = callbacks.to_delegate(handler, "MCMCSampler.ProgressChangedEventHandler")
        self.net.ProgressChanged += delegate
        self._handlers.append(delegate)
        return delegate

    def remove_progress_handler(self, delegate: Any) -> None:
        """Unsubscribe a delegate returned by :meth:`add_progress_handler`.

        Parameters
        ----------
        delegate : Numerics.Sampling.MCMC.MCMCSampler.ProgressChangedEventHandler
            The delegate.
        """
        self.net.ProgressChanged -= delegate
        if delegate in self._handlers:
            self._handlers.remove(delegate)

    @staticmethod
    def current_chain_index() -> int:
        """Return the chain the calling thread is advancing, or -1 outside one (``CurrentChainIndex``).

        Returns
        -------
        int
            The zero-based chain index.
        """
        return int(_mcmc().MCMCSampler.CurrentChainIndex)


class RWMH(MCMCSampler):
    """Random-walk Metropolis-Hastings (``Numerics.Sampling.MCMC.RWMH``).

    Parameters
    ----------
    priors : sequence of distributions
        Prior of each parameter (pynumerics or .NET distributions).
    log_likelihood : callable, str or System.Delegate
        Log-posterior of the parameter array.
    proposal_sigma : array_like, optional
        Proposal covariance, shape ``(p, p)``; default a zero matrix, which suits only
        ``initialize="MAP"`` (the MAP covariance replaces it).
    **settings
        Keyword settings of :class:`MCMCSampler`.
    """

    def __init__(self, priors: Sequence[Any], log_likelihood: Any, proposal_sigma: Any = None, **settings: Any) -> None:
        net_priors = _priors(priors)
        sigma = (_linear_algebra().Matrix(int(net_priors.Count)) if proposal_sigma is None
                 else _numerics_matrix(proposal_sigma, pass_through=True))
        super().__init__(_mcmc().RWMH(net_priors, _likelihood(log_likelihood), sigma), (log_likelihood,), **settings)

    @property
    def proposal_sigma(self) -> np.ndarray:
        """The proposal covariance (``ProposalSigma``)."""
        return to_numpy_matrix(self.net.ProposalSigma.ToArray())


class ARWMH(MCMCSampler):
    """Adaptive random-walk Metropolis-Hastings (``Numerics.Sampling.MCMC.ARWMH``).

    Parameters
    ----------
    priors : sequence of distributions
        Prior of each parameter.
    log_likelihood : callable, str or System.Delegate
        Log-posterior of the parameter array.
    **settings
        Keyword settings of :class:`MCMCSampler`.
    """

    def __init__(self, priors: Sequence[Any], log_likelihood: Any, **settings: Any) -> None:
        super().__init__(_mcmc().ARWMH(_priors(priors), _likelihood(log_likelihood)), (log_likelihood,), **settings)

    @property
    def scale(self) -> float:
        """The scale of the adapted covariance (``Scale``)."""
        return float(self.net.Scale)

    @scale.setter
    def scale(self, value: float) -> None:
        self.net.Scale = float(value)

    @property
    def beta(self) -> float:
        """The fraction of proposals drawn from the small identity covariance (``Beta``)."""
        return float(self.net.Beta)

    @beta.setter
    def beta(self, value: float) -> None:
        self.net.Beta = float(value)

    @property
    def proposal_sigma(self) -> list[np.ndarray]:
        """The current proposal covariance of each chain (``ProposalSigma``)."""
        return [to_numpy_matrix(m.ToArray()) for m in self.net.ProposalSigma]


class DEMCz(MCMCSampler):
    """Differential evolution Markov chain with past-state sampling (``Numerics.Sampling.MCMC.DEMCz``).

    Parameters
    ----------
    priors : sequence of distributions
        Prior of each parameter.
    log_likelihood : callable, str or System.Delegate
        Log-posterior of the parameter array.
    **settings
        Keyword settings of :class:`MCMCSampler`; at least three chains are needed.
    """

    _NET_TYPE = "DEMCz"

    def __init__(self, priors: Sequence[Any], log_likelihood: Any, **settings: Any) -> None:
        net_type = getattr(_mcmc(), self._NET_TYPE)
        super().__init__(net_type(_priors(priors), _likelihood(log_likelihood)), (log_likelihood,), **settings)

    @property
    def jump(self) -> float:
        """The jump scale gamma (``Jump``)."""
        return float(self.net.Jump)

    @jump.setter
    def jump(self, value: float) -> None:
        self.net.Jump = float(value)

    @property
    def jump_threshold(self) -> float:
        """The fraction of iterations with a jump of 1 (``JumpThreshold``)."""
        return float(self.net.JumpThreshold)

    @jump_threshold.setter
    def jump_threshold(self, value: float) -> None:
        self.net.JumpThreshold = float(value)

    @property
    def noise(self) -> float:
        """The standard deviation b of the proposal noise (``Noise``)."""
        return float(self.net.Noise)

    @noise.setter
    def noise(self, value: float) -> None:
        self.net.Noise = float(value)


class DEMCzs(DEMCz):
    """DEMCz with snooker updates (``Numerics.Sampling.MCMC.DEMCzs``).

    Parameters
    ----------
    priors : sequence of distributions
        Prior of each parameter.
    log_likelihood : callable, str or System.Delegate
        Log-posterior of the parameter array.
    **settings
        Keyword settings of :class:`MCMCSampler`; at least three chains are needed.
    """

    _NET_TYPE = "DEMCzs"

    @property
    def snooker_threshold(self) -> float:
        """The fraction of snooker updates (``SnookerThreshold``)."""
        return float(self.net.SnookerThreshold)

    @snooker_threshold.setter
    def snooker_threshold(self, value: float) -> None:
        self.net.SnookerThreshold = float(value)


class HMC(MCMCSampler):
    """Hamiltonian Monte Carlo (``Numerics.Sampling.MCMC.HMC``).

    Parameters
    ----------
    priors : sequence of distributions
        Prior of each parameter.
    log_likelihood : callable, str or System.Delegate
        Log-posterior of the parameter array.
    mass : array_like, optional
        Mass of each parameter; default ones.
    step_size : float, optional
        Leapfrog step size; default 0.1.
    steps : int, optional
        Leapfrog steps; default 50.
    gradient : callable, str or System.Delegate, optional
        Gradient of the log-posterior (``HMC.Gradient``), returning an array; default Numerics
        finite differences.
    **settings
        Keyword settings of :class:`MCMCSampler`.
    """

    def __init__(self, priors: Sequence[Any], log_likelihood: Any, mass: Any = None, step_size: float = 0.1,
                 steps: int = 50, gradient: Any = None, **settings: Any) -> None:
        grad = None if gradient is None else callbacks.to_delegate(gradient, "HMC.Gradient")
        net = _mcmc().HMC(_priors(priors), _likelihood(log_likelihood), _numerics_vector(mass, pass_through=True),
                          float(step_size), int(steps), grad)
        super().__init__(net, (log_likelihood, gradient), **settings)

    @property
    def mass(self) -> np.ndarray:
        """The mass vector (``Mass``)."""
        return to_numpy_doubles(self.net.Mass.Array)

    @property
    def step_size(self) -> float:
        """The leapfrog step size (``StepSize``)."""
        return float(self.net.StepSize)

    @step_size.setter
    def step_size(self, value: float) -> None:
        self.net.StepSize = float(value)

    @property
    def steps(self) -> int:
        """The number of leapfrog steps (``Steps``)."""
        return int(self.net.Steps)

    @steps.setter
    def steps(self, value: int) -> None:
        self.net.Steps = int(value)

    @property
    def gradient_function(self) -> Any:
        """The ``HMC.Gradient`` delegate in use (``GradientFunction``)."""
        return self.net.GradientFunction

    @staticmethod
    def quadratic_form(phi: Any, inverse_mass: Any) -> float:
        """Return the kinetic quadratic form ``sum(inverse_mass * phi**2)`` (``HMC.QuadraticForm``).

        Parameters
        ----------
        phi : array_like
            Momentum.
        inverse_mass : array_like
            Inverse mass of each parameter.

        Returns
        -------
        float
            The quadratic form.
        """
        return float(_mcmc().HMC.QuadraticForm(_numerics_vector(phi, pass_through=True),
                                               _numerics_vector(inverse_mass, pass_through=True)))


class NUTS(MCMCSampler):
    """The No-U-Turn sampler (``Numerics.Sampling.MCMC.NUTS``).

    Parameters
    ----------
    priors : sequence of distributions
        Prior of each parameter.
    log_likelihood : callable, str or System.Delegate
        Log-posterior of the parameter array.
    mass : array_like, optional
        Initial mass vector, adapted during warmup; default ones.
    step_size : float, optional
        Initial step size, adapted during warmup; default 0.1.
    max_tree_depth : int, optional
        Maximum tree depth; default 10.
    gradient : callable, str or System.Delegate, optional
        Gradient of the log-posterior (``HMC.Gradient``); default finite differences.
    log_likelihood_and_gradient : callable, str or System.Delegate, optional
        ``(value, gradient)`` in one call (``NUTS.LogLikelihoodAndGradient``).
    **settings
        Keyword settings of :class:`MCMCSampler`.
    """

    def __init__(self, priors: Sequence[Any], log_likelihood: Any, mass: Any = None, step_size: float = 0.1,
                 max_tree_depth: int = 10, gradient: Any = None, log_likelihood_and_gradient: Any = None,
                 **settings: Any) -> None:
        grad = None if gradient is None else callbacks.to_delegate(gradient, "HMC.Gradient")
        net = _mcmc().NUTS(_priors(priors), _likelihood(log_likelihood), _numerics_vector(mass, pass_through=True),
                           float(step_size), int(max_tree_depth), grad)
        super().__init__(net, (log_likelihood, gradient, log_likelihood_and_gradient), **settings)
        if log_likelihood_and_gradient is not None:
            self.log_likelihood_and_gradient_function = log_likelihood_and_gradient

    @property
    def mass(self) -> np.ndarray:
        """The initial mass vector (``Mass``)."""
        return to_numpy_doubles(self.net.Mass.Array)

    @property
    def max_tree_depth(self) -> int:
        """The maximum tree depth (``MaxTreeDepth``)."""
        return int(self.net.MaxTreeDepth)

    @max_tree_depth.setter
    def max_tree_depth(self, value: int) -> None:
        self.net.MaxTreeDepth = int(value)

    @property
    def gradient_function(self) -> Any:
        """The ``HMC.Gradient`` delegate in use (``GradientFunction``)."""
        return self.net.GradientFunction

    @property
    def log_likelihood_and_gradient_function(self) -> Any:
        """The joint value-and-gradient delegate, or None (``LogLikelihoodAndGradientFunction``)."""
        return self.net.LogLikelihoodAndGradientFunction

    @log_likelihood_and_gradient_function.setter
    def log_likelihood_and_gradient_function(self, value: Any) -> None:
        delegate = None if value is None else callbacks.to_delegate(value, "NUTS.LogLikelihoodAndGradient")
        self.net.LogLikelihoodAndGradientFunction = delegate
        callbacks.apply_parallel_policy(self.net, value)

    @property
    def target_acceptance_rate(self) -> float:
        """The acceptance probability step-size adaptation targets (``TargetAcceptanceRate``)."""
        return float(self.net.TargetAcceptanceRate)

    @target_acceptance_rate.setter
    def target_acceptance_rate(self, value: float) -> None:
        self.net.TargetAcceptanceRate = float(value)

    @property
    def adapt_mass_matrix(self) -> bool:
        """Whether warmup adapts the diagonal mass (``AdaptMassMatrix``)."""
        return bool(self.net.AdaptMassMatrix)

    @adapt_mass_matrix.setter
    def adapt_mass_matrix(self, value: bool) -> None:
        self.net.AdaptMassMatrix = bool(value)

    @property
    def hamiltonian_acceptance_rates(self) -> np.ndarray:
        """Mean Hamiltonian acceptance statistic per chain (``HamiltonianAcceptanceRates``)."""
        return to_numpy_doubles(self.net.HamiltonianAcceptanceRates)

    @property
    def diagnostic_sample_counts(self) -> np.ndarray:
        """Post-warmup transitions per chain behind the diagnostics (``DiagnosticSampleCounts``)."""
        return to_numpy_ints(self.net.DiagnosticSampleCounts).astype(np.int64)

    @property
    def divergence_counts(self) -> np.ndarray:
        """Divergent transitions per chain (``DivergenceCounts``)."""
        return to_numpy_ints(self.net.DivergenceCounts).astype(np.int64)

    @property
    def max_tree_depth_hit_counts(self) -> np.ndarray:
        """Transitions per chain that reached the maximum depth (``MaxTreeDepthHitCounts``)."""
        return to_numpy_ints(self.net.MaxTreeDepthHitCounts).astype(np.int64)

    @property
    def mean_tree_depths(self) -> np.ndarray:
        """Mean tree depth per chain (``MeanTreeDepths``)."""
        return to_numpy_doubles(self.net.MeanTreeDepths)

    @property
    def mean_leapfrog_steps(self) -> np.ndarray:
        """Mean leapfrog steps per transition per chain (``MeanLeapfrogSteps``)."""
        return to_numpy_doubles(self.net.MeanLeapfrogSteps)

    @property
    def step_sizes(self) -> np.ndarray:
        """The adapted step size of each chain (``StepSizes``)."""
        return to_numpy_doubles(self.net.StepSizes)

    @property
    def energy_bayesian_fraction_of_missing_information(self) -> np.ndarray:
        """E-BFMI per chain (``EnergyBayesianFractionOfMissingInformation``)."""
        return to_numpy_doubles(self.net.EnergyBayesianFractionOfMissingInformation)


class Gibbs(MCMCSampler):
    """Gibbs sampling from a user proposal of the full conditionals (``Numerics.Sampling.MCMC.Gibbs``).

    Every proposal is accepted, so the proposal must draw from the full conditional distributions.

    Parameters
    ----------
    priors : sequence of distributions
        Prior of each parameter.
    log_likelihood : callable, str or System.Delegate
        Log-posterior of the parameter array (recorded with each state).
    proposal : callable, str or System.Delegate
        ``proposal(values, prng)`` returning the next state (``Gibbs.Proposal``); ``prng`` is the
        chain's ``System.Random`` (a Numerics ``MersenneTwister``), so draws from it are reproducible.
    **settings
        Keyword settings of :class:`MCMCSampler`.
    """

    def __init__(self, priors: Sequence[Any], log_likelihood: Any, proposal: Any, **settings: Any) -> None:
        net = _mcmc().Gibbs(_priors(priors), _likelihood(log_likelihood),
                            callbacks.to_delegate(proposal, "Gibbs.Proposal"))
        super().__init__(net, (log_likelihood, proposal), **settings)

    @property
    def proposal_function(self) -> Any:
        """The ``Gibbs.Proposal`` delegate (``ProposalFunction``)."""
        return self.net.ProposalFunction


class SNIS(MCMCSampler):
    """Self-normalized importance sampling (``Numerics.Sampling.MCMC.SNIS``).

    Draws ``iterations`` independent points from the priors, or from ``importance`` when given
    (or from the inflated MAP normal with ``initialize="MAP"``), weights them and resamples
    ``output_length`` draws. Numerics evaluates the draws with ``Parallel.For`` while
    ``parallelize_chains`` is True; a Python log-likelihood turns it off, so the draws are then
    evaluated in order on the calling thread, with the same bits. C# source or
    :func:`distribution_log_likelihood` keep the draws parallel.

    Parameters
    ----------
    priors : sequence of distributions
        Prior of each parameter.
    log_likelihood : callable, str or System.Delegate
        Log-posterior of the parameter array.
    importance : pynumerics.multivariate.MultivariateNormal or Numerics MultivariateNormal, optional
        Importance distribution; default the priors (naive Monte Carlo).
    **settings
        Keyword settings of :class:`MCMCSampler`.
    """

    def __init__(self, priors: Sequence[Any], log_likelihood: Any, importance: Any = None, **settings: Any) -> None:
        mvn = None if importance is None else getattr(importance, "net", importance)
        super().__init__(_mcmc().SNIS(_priors(priors), _likelihood(log_likelihood), mvn), (log_likelihood,),
                         **settings)

    @property
    def weights(self) -> np.ndarray:
        """Normalized importance weight of each draw, in the order of ``markov_chains`` (``Weight``)."""
        return _read_sets(self.net.MarkovChains[0], self.number_of_parameters)[2]


class ParameterStatistics:
    """Posterior summary of one parameter (``Numerics.Sampling.MCMC.ParameterStatistics``).

    Parameters
    ----------
    net : Numerics.Sampling.MCMC.ParameterStatistics, optional
        The .NET object to wrap; default a new empty one.

    Attributes
    ----------
    net : Numerics.Sampling.MCMC.ParameterStatistics
        The wrapped .NET object.
    """

    def __init__(self, net: Any = None) -> None:
        self.net = _mcmc().ParameterStatistics() if net is None else net

    def __repr__(self) -> str:
        return (f"ParameterStatistics(mean={self.mean!r}, standard_deviation={self.standard_deviation!r}, "
                f"lower_ci={self.lower_ci!r}, upper_ci={self.upper_ci!r}, rhat={self.rhat!r}, ess={self.ess!r})")

    @property
    def rhat(self) -> float:
        """Rank-normalized split R-hat (``Rhat``)."""
        return float(self.net.Rhat)

    @rhat.setter
    def rhat(self, value: float) -> None:
        self.net.Rhat = float(value)

    @property
    def ess(self) -> float:
        """Conservative effective sample size (``ESS``)."""
        return float(self.net.ESS)

    @ess.setter
    def ess(self, value: float) -> None:
        self.net.ESS = float(value)

    @property
    def n(self) -> int:
        """Number of draws (``N``)."""
        return int(self.net.N)

    @n.setter
    def n(self, value: int) -> None:
        self.net.N = int(value)

    @property
    def mean(self) -> float:
        """Posterior mean (``Mean``)."""
        return float(self.net.Mean)

    @mean.setter
    def mean(self, value: float) -> None:
        self.net.Mean = float(value)

    @property
    def median(self) -> float:
        """Posterior median (``Median``)."""
        return float(self.net.Median)

    @median.setter
    def median(self, value: float) -> None:
        self.net.Median = float(value)

    @property
    def standard_deviation(self) -> float:
        """Posterior standard deviation (``StandardDeviation``)."""
        return float(self.net.StandardDeviation)

    @standard_deviation.setter
    def standard_deviation(self, value: float) -> None:
        self.net.StandardDeviation = float(value)

    @property
    def lower_ci(self) -> float:
        """Lower credible bound, the alpha/2 percentile (``LowerCI``)."""
        return float(self.net.LowerCI)

    @lower_ci.setter
    def lower_ci(self, value: float) -> None:
        self.net.LowerCI = float(value)

    @property
    def upper_ci(self) -> float:
        """Upper credible bound, the 1 - alpha/2 percentile (``UpperCI``)."""
        return float(self.net.UpperCI)

    @upper_ci.setter
    def upper_ci(self, value: float) -> None:
        self.net.UpperCI = float(value)


class ParameterResults:
    """Posterior summary, density and histogram of one parameter (``Numerics.Sampling.MCMC.ParameterResults``).

    Parameters
    ----------
    values : array_like or Numerics.Sampling.MCMC.ParameterResults
        Posterior draws of the parameter, or a .NET object to wrap.
    alpha : float, optional
        Two-sided interval level; default 0.1.
    is_sorted : bool, optional
        Whether ``values`` is already ascending (``sorted``); default False.

    Attributes
    ----------
    net : Numerics.Sampling.MCMC.ParameterResults
        The wrapped .NET object.
    """

    def __init__(self, values: Any, alpha: float = 0.1, is_sorted: bool = False) -> None:
        mc = _mcmc()
        if isinstance(values, mc.ParameterResults):
            self.net = values
        else:
            draws = np.ravel(np.asarray(values, dtype=np.float64))
            self.net = mc.ParameterResults(to_net_doubles(draws), float(alpha), bool(is_sorted))

    @property
    def summary_statistics(self) -> ParameterStatistics:
        """The summary statistics (``SummaryStatistics``)."""
        return ParameterStatistics(self.net.SummaryStatistics)

    @property
    def kernel_density(self) -> np.ndarray:
        """The kernel density graph, columns x and density (``KernelDensity``)."""
        return to_numpy_matrix(self.net.KernelDensity)

    @property
    def histogram(self) -> Any:
        """The ``Numerics.Data.Statistics.Histogram`` of the draws (``Histogram``)."""
        return self.net.Histogram

    @property
    def autocorrelation(self) -> np.ndarray:
        """The chain-averaged autocorrelation at lags 0 to 50 in column 1 (``Autocorrelation``); column 0 is zero."""
        return to_numpy_matrix(self.net.Autocorrelation)

    @autocorrelation.setter
    def autocorrelation(self, value: Any) -> None:
        self.net.Autocorrelation = to_net_matrix(value)


class MCMCResults:
    """Post-processed MCMC output (``Numerics.Sampling.MCMC.MCMCResults``).

    Built from a sampler after a run: R-hat over the post-warmup chains, ESS and autocorrelation
    over the output, and per-parameter summaries.

    Parameters
    ----------
    sampler : MCMCSampler or Numerics.Sampling.MCMC.MCMCSampler or Numerics.Sampling.MCMC.MCMCResults
        The sampler to post-process, or .NET results to wrap.
    alpha : float, optional
        Two-sided interval level; default 0.1 (90 percent intervals).

    Attributes
    ----------
    net : Numerics.Sampling.MCMC.MCMCResults
        The wrapped .NET object.
    """

    def __init__(self, sampler: Any, alpha: float = 0.1) -> None:
        mc = _mcmc()
        net = getattr(sampler, "net", sampler)
        self.net = net if isinstance(net, mc.MCMCResults) else mc.MCMCResults(net, float(alpha))

    @classmethod
    def from_output(cls, values: Any, map_values: Any, map_log_likelihood: float = float("nan"),
                    alpha: float = 0.1) -> MCMCResults:
        """Summarize posterior draws from elsewhere (``MCMCResults(ParameterSet, IList<ParameterSet>, alpha)``).

        Parameters
        ----------
        values : array_like
            Draws, shape ``(draws, parameters)``.
        map_values : array_like
            The MAP parameters.
        map_log_likelihood : float, optional
            The MAP log-likelihood; default NaN.
        alpha : float, optional
            Two-sided interval level; default 0.1.

        Returns
        -------
        MCMCResults
            Results with output, MAP and parameter summaries (no chains or diagnostics).
        """
        ps = _opt().ParameterSet
        from System.Collections.Generic import List  # noqa: PLC0415 - the runtime is running

        draws = np.atleast_2d(np.asarray(values, dtype=np.float64))
        sets = List[ps](draws.shape[0])
        nan = float("nan")
        for row in draws:
            sets.Add(ps(to_net_doubles(row), nan))
        return cls(_mcmc().MCMCResults(_parameter_set(map_values, map_log_likelihood), sets, float(alpha)))

    @classmethod
    def from_bytes(cls, data: bytes) -> MCMCResults:
        """Read results written by :meth:`to_bytes` (``MCMCResults.FromByteArray``).

        Parameters
        ----------
        data : bytes
            The serialized results.

        Returns
        -------
        MCMCResults
            The results.

        Raises
        ------
        ValueError
            If the bytes hold no results.
        """
        mc = _mcmc()
        import System  # noqa: PLC0415 - the runtime is running

        net = mc.MCMCResults.FromByteArray(System.Array[System.Byte](bytes(data)))
        if net is None:
            raise ValueError("The bytes hold no MCMC results.")
        return cls(net)

    def to_bytes(self) -> bytes:
        """Serialize the results (``MCMCResults.ToByteArray``).

        Returns
        -------
        bytes
            The serialized results.
        """
        return bytes(_mcmc().MCMCResults.ToByteArray(self.net))

    def _width(self) -> int:
        """Return the number of parameters."""
        return int(self.net.ParameterResults.Length)

    @property
    def markov_chains(self) -> np.ndarray:
        """The recorded chains, shape ``(chains, iterations, parameters)`` (``MarkovChains``).

        The chains are read from .NET once and cached; each access returns a new copy.
        """
        chains = self.net.MarkovChains
        if chains is None:
            return np.empty((0, 0, self._width()))
        return _cached_chain_values(self, chains, self._width())[0].copy()

    @property
    def output(self) -> np.ndarray:
        """All posterior draws, chain after chain, shape ``(draws, parameters)`` (``Output``)."""
        return _read_sets(self.net.Output, self._width(), fitness=False, weights=False)[0]

    @property
    def mean_log_likelihood(self) -> np.ndarray:
        """The log-likelihood averaged over chains at each iteration (``MeanLogLikelihood``)."""
        values = self.net.MeanLogLikelihood
        return np.empty(0) if values is None else to_numpy_doubles(values)

    @property
    def acceptance_rates(self) -> np.ndarray:
        """Acceptance rate per chain; the mean Hamiltonian acceptance statistic for NUTS (``AcceptanceRates``)."""
        values = self.net.AcceptanceRates
        return np.empty(0) if values is None else to_numpy_doubles(values)

    @property
    def parameter_results(self) -> list[ParameterResults]:
        """The per-parameter results (``ParameterResults``)."""
        return [ParameterResults(r) for r in self.net.ParameterResults]

    @property
    def map(self) -> tuple[np.ndarray, float]:
        """The MAP draw as ``(values, log_likelihood)`` (``MAP``)."""
        return _point(self.net.MAP)

    @property
    def posterior_mean(self) -> tuple[np.ndarray, float]:
        """The posterior mean and its log-likelihood as ``(values, log_likelihood)`` (``PosteriorMean``)."""
        return _point(self.net.PosteriorMean)

    def recompute_parameter_results(self, alpha: float) -> None:
        """Recompute the summaries at another interval level, keeping R-hat, ESS and autocorrelation.

        Calls ``RecomputeParameterResults``.

        Parameters
        ----------
        alpha : float
            The new two-sided level.
        """
        self.net.RecomputeParameterResults(float(alpha))


def _diagnostics() -> Any:
    """Return ``MCMCDiagnostics``."""
    return _mcmc().MCMCDiagnostics


def effective_sample_size(series: Any) -> float:
    """Return the conservative rank-normalized ESS of one series (``MCMCDiagnostics.EffectiveSampleSize``).

    The minimum of bulk and tail ESS of the series split in two halves (Vehtari et al. 2021). This is
    the effective size of a correlated series; for the Kish effective size of a set of weights use
    :func:`pynumerics.statistics.effective_sample_size`.

    Parameters
    ----------
    series : array_like
        The draws.

    Returns
    -------
    float
        The ESS, or NaN for an insufficient or constant series.
    """
    member = cached_overload("MCMCDiagnostics.EffectiveSampleSize(IList<double>)",
                             lambda: bind_list(_diagnostics().EffectiveSampleSize))
    return float(member(to_net_doubles(np.ravel(np.asarray(series, dtype=np.float64)))))


def effective_sample_sizes(chains: Any) -> tuple[np.ndarray, list[np.ndarray]]:
    """Return the conservative ESS of each parameter over several chains, and the averaged autocorrelation.

    Calls ``MCMCDiagnostics.EffectiveSampleSize(IList<List<ParameterSet>>, out double[][,])``.

    Parameters
    ----------
    chains : array_like or MCMCSampler or MCMCResults
        Chains of shape ``(chains, iterations, parameters)``, or a sampler or results whose
        ``MarkovChains`` are used.

    Returns
    -------
    ess : numpy.ndarray
        ESS per parameter.
    autocorrelation : list of numpy.ndarray
        The chain-averaged autocorrelation of each parameter at lags 0 to 50, shape ``(51, 2)``:
        column 1 holds the values; Numerics leaves column 0 at zero.
    """
    ess, acf = _diagnostics().EffectiveSampleSize(_chains_to_net(chains), None)
    return to_numpy_doubles(ess), [to_numpy_matrix(a) for a in acf]


def gelman_rubin(chains: Any, warmup_iterations: int = 0) -> np.ndarray:
    """Return rank-normalized split R-hat of each parameter (``MCMCDiagnostics.GelmanRubin``).

    The maximum of the rank-normalized and the folded rank-normalized split R-hat (Vehtari et al. 2021).

    Parameters
    ----------
    chains : array_like or MCMCSampler or MCMCResults
        Equal-length chains of shape ``(chains, iterations, parameters)``, or a sampler or results.
    warmup_iterations : int, optional
        Leading iterations of each chain to discard; default 0.

    Returns
    -------
    numpy.ndarray
        R-hat per parameter; NaN for a single chain.
    """
    return to_numpy_doubles(_diagnostics().GelmanRubin(_chains_to_net(chains), int(warmup_iterations)))


def minimum_sample_size(quantile: float, tolerance: float, probability: float) -> int:
    """Return the Raftery-Lewis minimum sample size, rounded to hundreds (``MCMCDiagnostics.MinimumSampleSize``).

    Parameters
    ----------
    quantile : float
        Posterior quantile of interest, for example 0.975.
    tolerance : float
        Acceptable error of that quantile's probability, for example 0.005.
    probability : float
        Probability of being within the tolerance, for example 0.95.

    Returns
    -------
    int
        ``q (1 - q) z^2 / r^2`` rounded to hundreds, with z the standard normal ``(1 + s) / 2`` quantile.
    """
    return int(_diagnostics().MinimumSampleSize(float(quantile), float(tolerance), float(probability)))


__all__ = [
    "ARWMH", "DEMCz", "DEMCzs", "Gibbs", "HMC", "InitializationType", "MCMCResults", "MCMCSampler", "NUTS",
    "ParameterResults", "ParameterStatistics", "RWMH", "SNIS", "distribution_log_likelihood",
    "effective_sample_size", "effective_sample_sizes", "gelman_rubin", "minimum_sample_size",
]

translate_public(globals())
