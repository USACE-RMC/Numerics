"""The MCMC samplers, results and diagnostics of ``pynumerics.mcmc``.

Posterior moments are checked on conjugate models against their closed forms within Monte Carlo
error (a few posterior standard deviations over the square root of the effective sample size),
the diagnostics against ArviZ (rank-normalized split R-hat and bulk and tail ESS, Vehtari et al.
2021) and the Raftery-Lewis sample size against its textbook formula.
"""

from __future__ import annotations

import math
import threading
import time

import numpy as np
import pytest
from scipy import stats

from pynumerics import callbacks, distributions, mcmc, sampling
from pynumerics._util import namespace
from pynumerics.marshal import to_net_doubles

needs_roslyn = pytest.mark.skipif(not callbacks.roslyn_available(), reason="Roslyn assemblies not staged")

RUN = {"iterations": 2000, "warmup_iterations": 500, "thinning_interval": 1, "output_length": 4000}
"""Short runs: 2000 recorded iterations per chain and 4000 output draws."""

X = stats.norm(3.0, 2.0).ppf(sampling.MersenneTwister(7).next_doubles(30))
"""Thirty observations from N(3, 2), drawn from the Numerics Mersenne Twister so every platform sees them."""


def _within_mc_error(results: mcmc.MCMCResults, index: int, mean: float, sd: float, k: float = 5.0) -> None:
    """Assert a posterior mean and standard deviation match the exact ones within Monte Carlo error."""
    summary = results.parameter_results[index].summary_statistics
    ess = summary.ess if math.isfinite(summary.ess) and summary.ess > 0 else summary.n
    assert abs(summary.mean - mean) < k * sd / math.sqrt(ess), (summary, mean, sd)
    assert summary.standard_deviation == pytest.approx(sd, rel=k / math.sqrt(2.0 * ess) + 0.02)


def _normal_normal() -> tuple[float, float]:
    """Return the exact posterior mean and sd of mu for N(mu, 2) data under a N(0, 5) prior."""
    precision = 1.0 / 25.0 + X.size / 4.0
    return float((X.sum() / 4.0) / precision), float(1.0 / math.sqrt(precision))


def _normal_normal_python(theta: np.ndarray) -> float:
    """Log-posterior of mu: N(mu, 2) likelihood plus the N(0, 5) prior, up to a constant."""
    mu = theta[0]
    return float(-0.5 * np.sum((X - mu) ** 2) / 4.0 - 0.5 * mu * mu / 25.0)


BETA_POSTERIOR = stats.beta(2 + 7, 3 + 13)
"""Beta(2, 3) prior and 7 successes in 20 trials: the posterior is Beta(9, 16)."""


def _beta_binomial_python(theta: np.ndarray) -> float:
    """Log-posterior of p up to a constant: 8 log p + 15 log(1 - p)."""
    p = theta[0]
    if not 0.0 < p < 1.0:
        return -math.inf
    return 8.0 * math.log(p) + 15.0 * math.log1p(-p)


BETA_BINOMIAL_CSHARP = ("p => (p[0] <= 0 || p[0] >= 1) ? double.NegativeInfinity"
                        " : 8 * Math.Log(p[0]) + 15 * Math.Log(1 - p[0])")


# Conjugate oracles


def test_rwmh_python_callback_recovers_the_normal_normal_posterior() -> None:
    # Arrange
    mean, sd = _normal_normal()
    sampler = mcmc.RWMH([distributions.Normal(0.0, 5.0)], _normal_normal_python, [[0.3]], **RUN)
    # Act
    results = sampler.sample()
    # Assert
    assert sampler.parallelize_chains is False
    _within_mc_error(results, 0, mean, sd)
    assert sampler.markov_chains.shape == (4, 2000, 1)
    assert results.output.shape == (4000, 1)


@needs_roslyn
def test_distribution_log_likelihood_with_prior_recovers_the_normal_normal_posterior_in_net() -> None:
    # Arrange: sample mu only; sigma = 2 stays fixed in the template
    mean, sd = _normal_normal()
    prior = distributions.Normal(0.0, 5.0)
    model = mcmc.distribution_log_likelihood(distributions.Normal(0.0, 2.0), X, priors=[prior], parameter_indices=[0],
                                             add_log_prior=True)
    sampler = mcmc.ARWMH([prior], model, **RUN)
    # Act
    results = sampler.sample()
    # Assert: no Python per evaluation, so the chains stay parallel
    assert sampler.parallelize_chains is True
    _within_mc_error(results, 0, mean, sd)
    expected = stats.norm(1.5, 2.0).logpdf(X).sum() + stats.norm(0.0, 5.0).logpdf(1.5)
    assert model(to_net_doubles([1.5])) == pytest.approx(expected, rel=1e-12)
    assert model(to_net_doubles([math.inf])) == -math.inf


@needs_roslyn
@pytest.mark.parametrize("sampler_type", [mcmc.DEMCz, mcmc.DEMCzs])
def test_population_samplers_with_csharp_source_recover_the_beta_binomial_posterior(sampler_type: type) -> None:
    # Arrange
    sampler = sampler_type([distributions.Uniform(0.0, 1.0)], BETA_BINOMIAL_CSHARP, **RUN)
    # Act
    results = sampler.sample()
    # Assert
    assert sampler.is_population_sampler and sampler.parallelize_chains
    _within_mc_error(results, 0, BETA_POSTERIOR.mean(), BETA_POSTERIOR.std())
    assert sampler.population_matrix.shape[1] == 1


def test_hmc_with_python_gradient_recovers_the_beta_binomial_posterior() -> None:
    # Arrange
    def gradient(theta: np.ndarray) -> np.ndarray:
        p = theta[0]
        return np.array([8.0 / p - 15.0 / (1.0 - p)])

    sampler = mcmc.HMC([distributions.Uniform(0.0, 1.0)], _beta_binomial_python, step_size=0.05, steps=5,
                       gradient=gradient, iterations=1000, warmup_iterations=200, thinning_interval=1,
                       output_length=2000)
    # Act
    results = sampler.sample()
    # Assert
    assert sampler.parallelize_chains is False
    assert sampler.steps == 5 and sampler.step_size == 0.05 and sampler.mass.tolist() == [1.0]
    _within_mc_error(results, 0, BETA_POSTERIOR.mean(), BETA_POSTERIOR.std())


def test_nuts_with_value_and_gradient_recovers_the_beta_binomial_posterior() -> None:
    # Arrange
    def value_and_gradient(theta: np.ndarray) -> tuple[float, np.ndarray]:
        p = theta[0]
        return _beta_binomial_python(theta), np.array([8.0 / p - 15.0 / (1.0 - p)])

    sampler = mcmc.NUTS([distributions.Uniform(0.0, 1.0)], _beta_binomial_python,
                        log_likelihood_and_gradient=value_and_gradient, iterations=600, warmup_iterations=200,
                        thinning_interval=1, output_length=1200, number_of_chains=2)
    # Act
    results = sampler.sample()
    # Assert
    assert sampler.parallelize_chains is False
    _within_mc_error(results, 0, BETA_POSTERIOR.mean(), BETA_POSTERIOR.std())
    assert sampler.divergence_counts.tolist() == [0, 0]
    assert np.all(sampler.step_sizes > 0) and np.all(sampler.mean_tree_depths >= 1)
    assert np.all((sampler.hamiltonian_acceptance_rates > 0) & (sampler.hamiltonian_acceptance_rates <= 1))
    assert sampler.diagnostic_sample_counts.tolist() == [600 - 200 + 1200 // 2] * 2
    assert sampler.max_tree_depth_hit_counts.tolist() == [0, 0]
    assert np.all(sampler.mean_leapfrog_steps >= 1)
    assert np.all(np.isfinite(sampler.energy_bayesian_fraction_of_missing_information))
    assert results.acceptance_rates.tolist() == sampler.hamiltonian_acceptance_rates.tolist()


def test_gibbs_with_exact_conditional_draws_recovers_the_beta_binomial_posterior() -> None:
    # Arrange: the full conditional of the one-parameter model is the posterior itself
    posterior = namespace("Numerics.Distributions").BetaDistribution(9.0, 16.0)

    def proposal(theta: np.ndarray, prng: object) -> np.ndarray:
        return np.array([posterior.InverseCDF(prng.NextDouble())])

    sampler = mcmc.Gibbs([distributions.Uniform(0.0, 1.0)], _beta_binomial_python, proposal, iterations=2000,
                         output_length=2000)
    # Act
    results = sampler.sample()
    # Assert
    _within_mc_error(results, 0, BETA_POSTERIOR.mean(), BETA_POSTERIOR.std())
    assert sampler.proposal_function is not None


def test_snis_recovers_the_beta_binomial_posterior() -> None:
    # Arrange
    sampler = mcmc.SNIS([distributions.Uniform(0.0, 1.0)], _beta_binomial_python, iterations=20000,
                        output_length=4000)
    # Act
    results = sampler.sample()
    weights = sampler.weights
    draws = sampler.markov_chains[0, :, 0]
    # Assert: the importance-weighted moments are the posterior moments
    assert weights.sum() == pytest.approx(1.0, abs=1e-12)
    weighted_mean = float(np.sum(weights * draws))
    assert weighted_mean == pytest.approx(BETA_POSTERIOR.mean(), abs=4 * BETA_POSTERIOR.std() / math.sqrt(2000))
    _within_mc_error(results, 0, BETA_POSTERIOR.mean(), BETA_POSTERIOR.std(), k=6.0)


def test_snis_with_a_python_likelihood_runs_serially_with_the_parallel_bits() -> None:
    # Arrange: the Python likelihood turns ParallelizeChains off; SNIS honours it
    threads: set[int] = set()

    def log_likelihood(theta: np.ndarray) -> float:
        threads.add(threading.get_ident())
        return _beta_binomial_python(theta)

    serial = mcmc.SNIS([distributions.Uniform(0.0, 1.0)], log_likelihood, iterations=3000, output_length=500)
    parallel = mcmc.SNIS([distributions.Uniform(0.0, 1.0)], _beta_binomial_python, iterations=3000,
                         output_length=500)
    parallel.parallelize_chains = True
    # Act
    serial.sample()
    parallel.sample()
    # Assert
    assert not serial.parallelize_chains and parallel.parallelize_chains
    assert threads == {threading.get_ident()}
    np.testing.assert_array_equal(serial.markov_chains.view(np.uint64), parallel.markov_chains.view(np.uint64))
    np.testing.assert_array_equal(serial.weights.view(np.uint64), parallel.weights.view(np.uint64))


@needs_roslyn
def test_flat_prior_normal_model_gives_the_student_t_marginal_of_the_mean() -> None:
    # Arrange: p(mu, sigma) flat, so mu | x = xbar + sqrt(S / (n (n - 2))) t_{n - 2}
    n, xbar, ss = X.size, float(X.mean()), float(np.sum((X - X.mean()) ** 2))
    sd = math.sqrt(ss / (n * (n - 4)))
    model = mcmc.distribution_log_likelihood(distributions.Normal(0.0, 1.0), X)
    sampler = mcmc.DEMCzs([distributions.Uniform(-20.0, 20.0), distributions.Uniform(0.01, 20.0)], model, **RUN)
    # Act
    results = sampler.sample()
    # Assert
    _within_mc_error(results, 0, xbar, sd)


# Blueprint shapes (Numerics-Python-Examples notebooks 03 to 06 and examples/bayesian_regression.py)


def test_rwmh_with_map_initialization_replaces_the_blank_proposal() -> None:
    # Arrange: notebook 03 passes Matrix(n) and initializes at the MAP
    sampler = mcmc.RWMH([distributions.Uniform(-20, 20), distributions.Uniform(0.1, 20)],
                        lambda t: float(np.sum(stats.norm(t[0], t[1]).logpdf(X))), initialize="MAP", **RUN)
    # Act
    results = sampler.sample()
    # Assert
    assert sampler.initialize is mcmc.InitializationType.MAP and not sampler.map_initialization_failed
    assert np.all(np.diag(sampler.proposal_sigma) > 0)
    values, log_likelihood = results.map
    assert values == pytest.approx([X.mean(), X.std()], rel=0.05)
    assert log_likelihood == pytest.approx(stats.norm(X.mean(), X.std()).logpdf(X).sum(), abs=0.1)


def test_default_settings_and_derived_counts_match_numerics() -> None:
    # Arrange
    sampler = mcmc.DEMCzs([distributions.Uniform(0, 1)], _beta_binomial_python)
    # Act / Assert: Numerics defaults and TransitionCount = (Iterations + OutputLength / chains) * thinning
    assert (sampler.seed, sampler.iterations, sampler.warmup_iterations, sampler.number_of_chains,
            sampler.thinning_interval, sampler.output_length) == (12345, 3500, 1750, 4, 20, 10000)
    assert sampler.transition_count == (3500 + 2500) * 20
    assert sampler.total_transition_count == 4 * sampler.transition_count
    assert sampler.number_of_parameters == 1
    assert isinstance(sampler.prior_distributions[0], distributions.Uniform)
    assert sampler.log_likelihood_function is not None
    assert sampler.resume_simulation is True and sampler.progress_changed_rate == 0.01
    sampler.jump, sampler.jump_threshold, sampler.noise, sampler.snooker_threshold = 2.0, 0.2, 1e-4, 0.2
    assert (sampler.jump, sampler.jump_threshold, sampler.noise, sampler.snooker_threshold) == (2.0, 0.2, 1e-4, 0.2)


def test_arwmh_settings_and_acceptance_bookkeeping() -> None:
    # Arrange
    sampler = mcmc.ARWMH([distributions.Uniform(0, 1)], _beta_binomial_python, **RUN)
    sampler.scale, sampler.beta = sampler.scale, sampler.beta
    # Act
    sampler.sample()
    # Assert
    assert np.all(sampler.accept_count <= sampler.sample_count)
    assert sampler.acceptance_rates == pytest.approx(sampler.accept_count / sampler.sample_count)
    assert sampler.mean_log_likelihood.shape == (2000,)
    assert len(sampler.proposal_sigma) == 4 and sampler.proposal_sigma[0].shape == (1, 1)
    assert [o.shape for o in sampler.output] == [(1000, 1)] * 4
    assert sampler.markov_chain_log_likelihoods.shape == (4, 2000)


def test_blueprint_hmc_and_nuts_constructor_arguments() -> None:
    # Arrange / Act: notebook 05 builds HMC(stepSize=0.03, steps=3) and NUTS(stepSize=0.1, maxTreeDepth=7)
    priors = [distributions.Uniform(0, 1)]
    hmc = mcmc.HMC(priors, _beta_binomial_python, step_size=0.03, steps=3)
    nuts = mcmc.NUTS(priors, _beta_binomial_python, step_size=0.1, max_tree_depth=7, mass=[2.0])
    nuts.target_acceptance_rate, nuts.adapt_mass_matrix = 0.9, False
    # Assert
    assert (hmc.step_size, hmc.steps, nuts.max_tree_depth, nuts.mass.tolist()) == (0.03, 3, 7, [2.0])
    assert nuts.target_acceptance_rate == 0.9 and nuts.adapt_mass_matrix is False
    assert hmc.gradient_function is not None and nuts.gradient_function is not None
    assert nuts.log_likelihood_and_gradient_function is None


def test_hmc_quadratic_form_is_the_mass_weighted_sum_of_squares() -> None:
    # Act / Assert
    assert mcmc.HMC.quadratic_form([1.0, -2.0, 3.0], [0.5, 1.0, 2.0]) == pytest.approx(0.5 + 4.0 + 18.0)
    with pytest.raises(ValueError, match="same length"):
        mcmc.HMC.quadratic_form([1.0], [1.0, 2.0])


# Reproducibility and the parallel rule


def test_a_seed_reproduces_the_chains_and_another_seed_changes_them() -> None:
    # Arrange
    def run(seed: int) -> np.ndarray:
        return mcmc.ARWMH([distributions.Uniform(0, 1)], _beta_binomial_python, seed=seed, **RUN).sample().output

    # Act
    first, again, other = run(99), run(99), run(100)
    # Assert
    assert np.array_equal(first, again)
    assert not np.array_equal(first, other)


@needs_roslyn
def test_parallel_compiled_chains_equal_serial_python_chains_bit_for_bit() -> None:
    # Arrange: the same model as a Python callback (serial) and as compiled C# (parallel)
    nd = namespace("Numerics.Distributions")
    data = to_net_doubles(X)
    priors = [distributions.Uniform(-20, 20), distributions.Uniform(0.1, 20)]
    python = mcmc.ARWMH(priors, lambda t: nd.Normal(t[0], t[1]).LogLikelihood(data), **RUN)
    compiled = mcmc.ARWMH(priors, mcmc.distribution_log_likelihood(distributions.Normal(0, 1), X), **RUN)
    # Act
    python.sample()
    compiled.sample()
    # Assert
    assert python.parallelize_chains is False and compiled.parallelize_chains is True
    assert np.array_equal(python.markov_chains, compiled.markov_chains)


def test_an_explicit_parallel_setting_is_kept_for_net_delegates_and_overridden_for_python() -> None:
    # Arrange / Act
    python = mcmc.DEMCz([distributions.Uniform(0, 1)], _beta_binomial_python, parallelize_chains=True)
    native = mcmc.DEMCz([distributions.Uniform(0, 1)], python.log_likelihood_function, parallelize_chains=False)
    # Assert: a delegate built from a Python callable still calls Python
    assert python.parallelize_chains is False
    assert native.parallelize_chains is False
    native.parallelize_chains = True
    assert native.parallelize_chains is True


@needs_roslyn
def test_compiled_likelihood_runs_rwmh_faster_than_a_python_likelihood() -> None:
    # Arrange: the notebook 03 model, mu and sigma of a normal sample
    nd = namespace("Numerics.Distributions")
    data = to_net_doubles(X)
    priors = [distributions.Uniform(-20, 20), distributions.Uniform(0.1, 20)]
    proposal = [[0.2, 0.0], [0.0, 0.1]]
    compiled = mcmc.distribution_log_likelihood(distributions.Normal(0, 1), X)

    def timed(log_likelihood: object, parallel: bool | None) -> float:
        best = math.inf
        for _ in range(3):
            sampler = mcmc.RWMH(priors, log_likelihood, proposal, parallelize_chains=parallel, **RUN)
            start = time.perf_counter()
            sampler.sample()
            best = min(best, time.perf_counter() - start)
        return best

    # Act
    python_seconds = timed(lambda t: nd.Normal(t[0], t[1]).LogLikelihood(data), None)
    compiled_serial = timed(compiled, False)
    compiled_parallel = timed(compiled, True)
    # Assert
    timings = (f"RWMH, 12000 evaluations: Python {python_seconds:.3f} s, compiled serial {compiled_serial:.3f} s, "
               f"compiled parallel {compiled_parallel:.3f} s")
    assert compiled_serial < python_seconds / 2, timings
    assert compiled_parallel < python_seconds / 2, timings


# Results


def test_results_serialize_and_recompute_intervals_without_touching_the_diagnostics() -> None:
    # Arrange
    results = mcmc.ARWMH([distributions.Uniform(0, 1)], _beta_binomial_python, **RUN).sample()
    before = results.parameter_results[0].summary_statistics
    rhat, ess, lower90, upper90 = before.rhat, before.ess, before.lower_ci, before.upper_ci
    # Act
    copy = mcmc.MCMCResults.from_bytes(results.to_bytes())
    results.recompute_parameter_results(0.05)
    after = results.parameter_results[0].summary_statistics
    # Assert
    assert np.array_equal(copy.output, results.output)
    assert np.array_equal(copy.markov_chains, results.markov_chains)
    assert copy.map[1] == results.map[1]
    assert (after.rhat, after.ess) == (rhat, ess)
    assert after.lower_ci < lower90 and after.upper_ci > upper90
    assert results.mean_log_likelihood.shape == (2000,)
    assert results.posterior_mean[0] == pytest.approx([after.mean])
    with pytest.raises(ValueError, match="no MCMC results"):
        mcmc.MCMCResults.from_bytes(b"null")


def test_parameter_results_summarize_draws_like_numpy() -> None:
    # Arrange
    draws = stats.norm(1.0, 2.0).ppf(sampling.MersenneTwister(3).next_doubles(2001))
    # Act
    result = mcmc.ParameterResults(draws, alpha=0.1)
    summary = result.summary_statistics
    # Assert: Type-7 percentiles and the sample moments
    assert summary.n == 2001
    assert summary.median == pytest.approx(np.percentile(draws, 50), rel=1e-12)
    assert summary.lower_ci == pytest.approx(np.percentile(draws, 5), rel=1e-12)
    assert summary.upper_ci == pytest.approx(np.percentile(draws, 95), rel=1e-12)
    assert summary.mean == pytest.approx(draws.mean(), rel=1e-9)
    assert summary.standard_deviation == pytest.approx(draws.std(ddof=1), rel=0.02)
    assert result.kernel_density.shape[1] == 2 and result.histogram is not None
    result.autocorrelation = np.array([[0.0, 1.0], [1.0, 0.5]])
    assert result.autocorrelation.tolist() == [[0.0, 1.0], [1.0, 0.5]]
    assert math.isnan(summary.rhat) and math.isnan(summary.ess)


def test_results_from_output_and_empty_statistics() -> None:
    # Arrange
    draws = np.column_stack([np.linspace(0, 1, 101), np.linspace(1, 3, 101)])
    # Act
    results = mcmc.MCMCResults.from_output(draws, [0.5, 2.0], -1.5, alpha=0.2)
    stats_obj = mcmc.ParameterStatistics()
    stats_obj.mean, stats_obj.median, stats_obj.n = 1.0, 2.0, 3
    stats_obj.standard_deviation, stats_obj.lower_ci, stats_obj.upper_ci = 0.5, -1.0, 4.0
    stats_obj.rhat, stats_obj.ess = 1.01, 250.0
    # Assert
    assert results.map[0].tolist() == [0.5, 2.0] and results.map[1] == -1.5
    assert results.output.shape == (101, 2)
    assert results.parameter_results[1].summary_statistics.median == pytest.approx(2.0)
    assert results.markov_chains.shape == (0, 0, 2)
    assert results.mean_log_likelihood.size == 0 and results.acceptance_rates.size == 0
    assert (stats_obj.mean, stats_obj.median, stats_obj.n, stats_obj.standard_deviation) == (1.0, 2.0, 3, 0.5)
    assert (stats_obj.lower_ci, stats_obj.upper_ci, stats_obj.rhat, stats_obj.ess) == (-1.0, 4.0, 1.01, 250.0)
    assert "mean=1.0" in repr(stats_obj)


# Diagnostics


def _ar1_chains(seed: int, shift: float = 0.0) -> np.ndarray:
    """Four AR(1) chains (phi 0.6) of 400 draws of two parameters from the Numerics Mersenne Twister."""
    e = stats.norm.ppf(sampling.MersenneTwister(seed).next_doubles(4 * 400 * 2).reshape(4, 400, 2))
    x = np.empty_like(e)
    x[:, 0] = e[:, 0]
    for t in range(1, 400):
        x[:, t] = 0.6 * x[:, t - 1] + e[:, t]
    x[0] += shift
    return x


ARVIZ_REFERENCE = {
    0.0: ([1.0163622345152468, 1.004331087025769], [356.5131362059822, 517.3949251670716]),
    0.5: ([1.028876833214592, 1.019040378202041], [212.92645062938462, 481.02246099413975]),
}
"""ArviZ 0.23.4 ``rhat(method="rank")`` and ``min(ess bulk, ess tail)`` of :func:`_ar1_chains(2024, shift)`."""


@pytest.mark.parametrize("shift", [0.0, 0.5])
def test_gelman_rubin_and_ess_match_arviz(shift: float) -> None:
    # Arrange
    chains = _ar1_chains(2024, shift)
    rhat_expected, ess_expected = ARVIZ_REFERENCE[shift]
    # Act
    rhat = mcmc.gelman_rubin(chains)
    ess, acf = mcmc.effective_sample_sizes(chains)
    # Assert
    assert rhat == pytest.approx(rhat_expected, rel=1e-10)
    assert ess == pytest.approx(ess_expected, rel=1e-10)
    assert len(acf) == 2 and acf[0].shape == (51, 2)
    assert mcmc.effective_sample_size(chains[0, :, 0]) == pytest.approx(54.56986461700465, rel=1e-10)


def test_gelman_rubin_and_ess_match_live_arviz() -> None:
    # Arrange
    az = pytest.importorskip("arviz")
    chains = _ar1_chains(77, 0.3)
    data = az.convert_to_dataset({"a": chains[..., 0], "b": chains[..., 1]})
    # Act
    rhat = mcmc.gelman_rubin(chains)
    ess, _ = mcmc.effective_sample_sizes(chains)
    # Assert
    expected_rhat = [float(az.rhat(data, method="rank")[k]) for k in "ab"]
    expected_ess = [min(float(az.ess(data, method="bulk")[k]), float(az.ess(data, method="tail")[k])) for k in "ab"]
    assert rhat == pytest.approx(expected_rhat, rel=1e-10)
    assert ess == pytest.approx(expected_ess, rel=1e-10)


def test_diagnostics_accept_samplers_results_and_respect_warmup() -> None:
    # Arrange
    sampler = mcmc.ARWMH([distributions.Uniform(0, 1)], _beta_binomial_python, **RUN)
    results = sampler.sample()
    # Act
    from_sampler = mcmc.gelman_rubin(sampler, warmup_iterations=sampler.warmup_iterations)
    from_array = mcmc.gelman_rubin(sampler.markov_chains, warmup_iterations=sampler.warmup_iterations)
    # Assert: MCMCResults uses the post-warmup R-hat
    assert from_sampler.tolist() == from_array.tolist()
    assert from_sampler[0] == results.parameter_results[0].summary_statistics.rhat
    assert mcmc.gelman_rubin(results).shape == (1,)
    assert mcmc.gelman_rubin(sampler.markov_chains[:1]).tolist() == [pytest.approx(math.nan, nan_ok=True)]


def test_diagnostics_reject_bad_input() -> None:
    # Act / Assert
    with pytest.raises(ValueError, match="chains, iterations, parameters"):
        mcmc.gelman_rubin(np.zeros(5))
    with pytest.raises(ValueError):
        mcmc.gelman_rubin(np.zeros((0, 10, 1)))
    assert math.isnan(mcmc.effective_sample_size(np.ones(50)))


@pytest.mark.parametrize(("q", "r", "s"), [(0.025, 0.005, 0.95), (0.975, 0.005, 0.95), (0.5, 0.01, 0.99)])
def test_minimum_sample_size_is_the_raftery_lewis_formula(q: float, r: float, s: float) -> None:
    # Act
    n = mcmc.minimum_sample_size(q, r, s)
    # Assert: q (1 - q) z^2 / r^2, z the (1 + s) / 2 normal quantile, rounded to hundreds
    expected = round(q * (1 - q) * stats.norm.ppf(0.5 * (s + 1)) ** 2 / r ** 2 / 100) * 100
    assert n == expected


# Errors and events


def test_invalid_settings_raise_value_errors() -> None:
    # Arrange
    demcz = mcmc.DEMCz([distributions.Uniform(0, 1)], _beta_binomial_python, number_of_chains=1)
    rwmh = mcmc.RWMH([distributions.Uniform(0, 1)], _beta_binomial_python, [[1.0, 0.0], [0.0, 1.0]])
    # Act / Assert
    with pytest.raises(ValueError):
        demcz.sample()
    with pytest.raises(ValueError, match="same number of rows"):
        rwmh.sample()
    with pytest.raises(ValueError, match="is not a member"):
        rwmh.initialize = "warm"
    with pytest.raises(ValueError, match="needs the priors"):
        mcmc.distribution_log_likelihood(distributions.Normal(0, 1), X, add_log_prior=True)


def test_progress_handler_receives_reports_and_unsubscribes() -> None:
    # Arrange
    seen: list[float] = []
    sampler = mcmc.ARWMH([distributions.Uniform(0, 1)], _beta_binomial_python, **RUN)
    sampler.progress_changed_rate = 0.25
    handler = sampler.add_progress_handler(lambda fraction, text: seen.append(fraction))
    # Act
    sampler.sample()
    sampler.report_progress(0.5)
    sampler.remove_progress_handler(handler)
    sampler.report_progress(0.75)
    # Assert
    assert seen[:4] == [0.25, 0.5, 0.75, 1.0] and seen[-1] == 0.5


def test_reset_cancel_and_chain_index_outside_a_run() -> None:
    # Arrange
    sampler = mcmc.ARWMH([distributions.Uniform(0, 1)], _beta_binomial_python, **RUN)
    sampler.sample()
    # Act
    sampler.reset()
    sampler.cancel_simulation()
    # Assert
    assert sampler.markov_chains.shape == (4, 0, 1)
    assert mcmc.MCMCSampler.current_chain_index() == -1
    assert sampler.cancellation_token_source is not None
    token = sampler.cancellation_token
    sampler.cancellation_token = token
    sampler.cancellation_token_source = None
    assert sampler.cancellation_token_source is None


def _chains_one_set_at_a_time(chains: object) -> tuple[np.ndarray, np.ndarray]:
    """Read ``List<ParameterSet>[]`` element by element, the reference for the bulk read."""
    from pynumerics.marshal import to_numpy_doubles  # noqa: PLC0415

    values = [[to_numpy_doubles(chains[c][i].Values) for i in range(chains[c].Count)] for c in range(chains.Length)]
    fitness = [[float(chains[c][i].Fitness) for i in range(chains[c].Count)] for c in range(chains.Length)]
    return np.array(values, dtype=np.float64), np.array(fitness, dtype=np.float64)


def _two_parameter_rwmh(seed: int) -> mcmc.RWMH:
    def target(p: np.ndarray) -> float:
        return float(-0.5 * ((p[0] - 1.0) / 0.7) ** 2 - 0.5 * ((p[1] + 2.0) / 0.5) ** 2)

    priors = [distributions.Uniform(-10.0, 10.0), distributions.Uniform(-10.0, 10.0)]
    return mcmc.RWMH(priors, target, [[0.5, 0.0], [0.0, 0.3]], seed=seed, iterations=1500, warmup_iterations=200)


def test_markov_chains_read_in_bulk_are_bitwise_the_per_set_values_and_are_cached() -> None:
    # Arrange
    sampler = _two_parameter_rwmh(12345)
    results = sampler.sample()
    expected_values, expected_fitness = _chains_one_set_at_a_time(sampler.net.MarkovChains)
    # Act
    chains = sampler.markov_chains
    log_likelihoods = sampler.markov_chain_log_likelihoods
    results_chains = results.markov_chains
    chains[0, 0, 0] = 1e300
    again = sampler.markov_chains
    # Assert
    assert chains.shape == expected_values.shape == (4, 1500, 2)
    np.testing.assert_array_equal(again.view(np.uint64), expected_values.view(np.uint64))
    np.testing.assert_array_equal(results_chains.view(np.uint64), expected_values.view(np.uint64))
    np.testing.assert_array_equal(log_likelihoods.view(np.uint64), expected_fitness.view(np.uint64))
    assert again is not chains and again.flags.owndata and again.flags.writeable
    assert sampler._chain_cache is not None and results._chain_cache is not None


def test_markov_chains_cache_is_cleared_when_the_sampler_runs_again() -> None:
    # Arrange
    sampler = _two_parameter_rwmh(1)
    sampler.sample()
    first = sampler.markov_chains
    # Act
    sampler.seed = 2
    sampler.sample()
    second = sampler.markov_chains
    # Assert
    expected, _ = _chains_one_set_at_a_time(sampler.net.MarkovChains)
    assert not np.array_equal(first, second)
    np.testing.assert_array_equal(second.view(np.uint64), expected.view(np.uint64))


def test_markov_chains_cache_follows_a_run_started_on_the_net_sampler() -> None:
    # Arrange
    sampler = _two_parameter_rwmh(3)
    sampler.sample()
    first = sampler.markov_chains
    # Act: a run through the .NET object bypasses sample(), so the chain identity check must notice
    sampler.net.PRNGSeed = 4
    sampler.net.Sample()
    second = sampler.markov_chains
    # Assert
    expected, _ = _chains_one_set_at_a_time(sampler.net.MarkovChains)
    assert not np.array_equal(first, second)
    np.testing.assert_array_equal(second.view(np.uint64), expected.view(np.uint64))


@pytest.mark.benchmark
def test_markov_chains_cost_one_bulk_read_then_a_copy() -> None:
    # Arrange: measured on a desktop, the first read costs about 3 us per set (22 us before the
    # bulk read) and a cached read about 10 us in all; the bounds leave a wide margin.
    sampler = _two_parameter_rwmh(5)
    sampler.sample()
    sets = 4 * 1500
    # Act
    start = time.perf_counter()
    _ = sampler.markov_chains
    first = time.perf_counter() - start
    start = time.perf_counter()
    for _ in range(100):
        _ = sampler.markov_chains
    cached = (time.perf_counter() - start) / 100
    # Assert
    assert first / sets < 12e-6, first / sets
    assert cached < first / 20, (cached, first)
