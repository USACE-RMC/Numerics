# Bootstrap

[← Back to Index](../index.md)

`Bootstrap<TData>` in `Numerics.Sampling` is a general-purpose bootstrap engine. You supply delegates
that resample the data, fit the model and compute the statistics of interest, and the class runs the
replicates in parallel and builds confidence intervals from them. It supports the percentile,
bias-corrected, BCa, normal and bootstrap-t intervals [[1]](#1), a double bootstrap and a
covariance-aware pivotal bootstrap.

## A Parametric Bootstrap

The example below bootstraps the 1% and 99% quantiles of a Normal distribution fitted by the method of
moments.

```cs
using Numerics.Distributions;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

var original = new Normal(100, 15);
double[] data = original.GenerateRandomValues(50, 12345);
var fitted = new Normal();
((IEstimation)fitted).Estimate(data, ParameterEstimationMethod.MethodOfMoments);

var boot = new Bootstrap<double[]>(data, new ParameterSet(fitted.GetParameters, double.NaN))
{
    Replicates = 2000,
    PRNGSeed = 12345,
    // Draw a new sample of the same size from the fitted distribution
    ResampleFunction = (sample, parameters, rng) =>
        new Normal(parameters.Values[0], parameters.Values[1]).GenerateRandomValues(sample.Length, rng.Next()),
    // Refit the model to a resampled data set
    FitFunction = sample =>
    {
        var d = new Normal();
        ((IEstimation)d).Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
        return new ParameterSet(d.GetParameters, double.NaN);
    },
    // The statistics to bootstrap
    StatisticFunction = parameters =>
    {
        var d = new Normal(parameters.Values[0], parameters.Values[1]);
        return new[] { d.InverseCDF(0.01), d.InverseCDF(0.99) };
    }
};

boot.Run();
var results = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile, alpha: 0.1);
foreach (var r in results.StatisticResults)
    Console.WriteLine($"{r.PopulationEstimate:F2}  90% CI [{r.LowerCI:F2}, {r.UpperCI:F2}]");
```

Each replicate draws from its own `MersenneTwister`, seeded from `PRNGSeed`, so the results are
reproducible regardless of how the replicates are scheduled across threads. A replicate whose resample,
fit or statistic throws, or returns non-finite values, is retried up to `MaxRetries` times, continuing
the replicate's own generator, and then counted in `FailedReplicates`. `FailedAttempts` counts every attempt that threw an exception, including
the ones a retry recovered from and the skipped inner replicates of the double and studentized bootstraps.

## Serial Runs

`UseParallel` (default `true`) runs the outer replicates with `Parallel.For`, so the delegates must be thread safe. Set it to `false` to run them one after another, in index order, on the calling thread, for delegates that are not thread safe or that gain nothing from threads (for example callbacks into Python, which hold the interpreter lock). Each replicate keeps its own seed, so a serial run gives bit-identical results. The BCa leave-one-out jackknife follows the same setting; its moments are merged in a fixed chunk order, so its acceleration constants are also bit-identical in both modes.

```cs
boot.UseParallel = false;   // same results, one thread
boot.Run();
```

## Progress and Cancellation

A long bootstrap can report progress and be stopped early:

- `Progress` (an `IProgress<int>`) receives the number of outer replicates completed so far, once per
  replicate, after its retries finish and whether or not it succeeded. The last value of a completed run
  is `Replicates`. With `UseParallel` on, `Report` is called concurrently from worker threads
  and the values can arrive out of order. Inner replicates of the double and studentized bootstraps and
  the BCa jackknife are not reported.
- `CancellationToken` stops `Run`, `RunDoubleBootstrap`, `RunWithStudentizedBootstrap` and
  `RunPivotalBootstrap` from starting new outer replicates. The replicates already running finish, and
  then the method throws `OperationCanceledException`. A cancelled run leaves no usable results, so run
  the bootstrap again before requesting confidence intervals.

Neither hook changes the random streams, so a run with them gives exactly the same results as a run
without them.

```cs
using System.Threading;

using var cancel = new CancellationTokenSource(TimeSpan.FromMinutes(5));
boot.Progress = new ConsoleProgress();
boot.CancellationToken = cancel.Token;
try
{
    boot.Run();
}
catch (OperationCanceledException)
{
    Console.WriteLine("The bootstrap was cancelled.");
}

// Called from worker threads; keep it thread-safe and cheap.
sealed class ConsoleProgress : IProgress<int>
{
    public void Report(int completed)
    {
        if (completed % 500 == 0) Console.WriteLine($"{completed} replicates completed");
    }
}
```

`System.Progress<T>` also works, but it posts each report to the captured synchronization context or the
thread pool, so in a console application some reports can print after `Run` has returned.

---

## References

<a id="1">[1]</a> Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the Bootstrap*. Chapman & Hall/CRC.

---

[← Previous: Convergence Diagnostics](convergence-diagnostics.md) | [Back to Index](../index.md)
