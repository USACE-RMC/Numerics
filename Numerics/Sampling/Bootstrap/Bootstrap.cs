/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Numerics.Sampling
{
    public class Bootstrap<TData>
    {

        /// <summary>
        /// Delegate function for resampling the original data and model fit.
        /// </summary>
        public Func<TData, ParameterSet, Random, TData> ResampleFunction { get; set; } = null!;

        /// <summary>
        /// Delegate function for fitting a model.
        /// </summary>
        public Func<TData, ParameterSet> FitFunction { get; set; } = null!;

        /// <summary>
        /// Delegate function for extracting a statistic from the fit result.
        /// </summary>
        public Func<ParameterSet, double[]> StatisticFunction { get; set; } = null!;

        /// <summary>
        /// Number of bootstrap replicates.
        /// </summary>
        public int Replicates { get; set; } = 10000;

        /// <summary>
        /// Gets and sets the PRNG seed for reproducibility.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// Constructs a new bootstrap class.
        /// </summary>
        /// <param name="originalData">The original data.</param>
        /// <param name="originalParameters">The original fitted parameter set.</param>
        public Bootstrap(TData originalData, ParameterSet originalParameters) 
        { 
            _originalData = originalData;
            _originalParameters = originalParameters;
        }

        private TData _originalData;
        private ParameterSet _originalParameters;
        private ParameterSet[] _bootstrapParameterSets = null!;
        private double[][] _bootstrapStatistics = null!;

        /// <summary>
        /// Gets the bootstrapped model parameter sets.
        /// </summary>
        public ParameterSet[] BootstrapParameterSets => _bootstrapParameterSets;

        /// <summary>
        /// Gets the bootstrapped statistics.
        /// </summary>
        public IReadOnlyList<double[]> BootstrapStatistics => _bootstrapStatistics;

        /// <summary>
        /// Runs the basic bootstrap procedure.
        /// </summary>
        public void Run()
        {
            if (ResampleFunction == null)
                throw new InvalidOperationException("Bootstrap Sample Function must be set.");
            if (FitFunction == null)
                throw new InvalidOperationException("Fit Function must be set.");
            if (StatisticFunction == null)
                throw new InvalidOperationException("Statistic Function must be set.");

            _bootstrapParameterSets = new ParameterSet[Replicates];
            _bootstrapStatistics = new double[Replicates][];

            var prng = new MersenneTwister(PRNGSeed);
            var seeds = prng.NextIntegers(Replicates);

            Parallel.For(0, Replicates, idx =>
            {
                var rng = new MersenneTwister(seeds[idx]);
                var sample = ResampleFunction(_originalData, _originalParameters, rng);
                var fit = FitFunction(sample);
                var stat = StatisticFunction(fit);
                _bootstrapParameterSets[idx] = fit;
                _bootstrapStatistics[idx] = stat;
            });
        }

        public void RunDoubleBootstrap(int innerReplicates = 300)
        {
            if (ResampleFunction == null)
                throw new InvalidOperationException("Bootstrap Sample Function must be set.");
            if (FitFunction == null)
                throw new InvalidOperationException("Fit Function must be set.");
            if (StatisticFunction == null)
                throw new InvalidOperationException("Statistic Function must be set.");

            _bootstrapParameterSets = new ParameterSet[Replicates];
            _bootstrapStatistics = new double[Replicates][];

            var prng = new MersenneTwister(PRNGSeed);
            var seeds = prng.NextIntegers(Replicates);

            Parallel.For(0, Replicates, idx =>
            {
                var rng = new MersenneTwister(seeds[idx]);

                // Step 1: outer bootstrap
                var outerSample = ResampleFunction(_originalData, _originalParameters, rng);
                var outerFit = FitFunction(outerSample);
                var outerStat = StatisticFunction(outerFit);

                // Step 2: inner bootstrap
                int p = outerFit.Values.Length;
                var parmsInnerSum = new double[p];
                int s = outerStat.Length;
                var statsInnerSum = new double[s];

                for (int k = 0; k < innerReplicates; k++)
                {
                    var innerSample = ResampleFunction(outerSample, outerFit, rng);
                    var innerFit = FitFunction(innerSample);
                    var innerStat = StatisticFunction(innerFit);

                    for (int i = 0; i < p; i++)
                        parmsInnerSum[i] += innerFit.Values[i];
                    for (int i = 0; i < s; i++)
                        statsInnerSum[i] += innerStat[i];
                }

                // bias correct the parameters
                var biasCorrectedParms = new double[p];
                for (int i = 0; i < p; i++)
                {
                    double innerMean = parmsInnerSum[i] / innerReplicates;
                    biasCorrectedParms[i] = outerFit.Values[i] - (innerMean - outerFit.Values[i]);
                }

                // bias correct the statistics
                var biasCorrectedStats = new double[s];
                for (int i = 0; i < s; i++)
                {
                    double innerMean = statsInnerSum[i] / innerReplicates;
                    biasCorrectedStats[i] = outerStat[i] - (innerMean - outerStat[i]);
                }

                _bootstrapParameterSets[idx].Values = biasCorrectedParms;
                _bootstrapStatistics[idx] = biasCorrectedStats;
            });
        }

        /// <summary>
        /// Returns percentile confidence intervals for each parameter.
        /// </summary>
        /// <param name="confidenceLevel">Confidence level. Default = 0.9.</param>
        public (double[] Lower, double[] Upper) GetPercentileConfidenceIntervals(double confidenceLevel = 0.90)
        {
            if (_bootstrapStatistics is null)
                throw new InvalidOperationException("Run() must be called before requesting CIs.");

            int nParams = _bootstrapStatistics[0].Length;
            int n = _bootstrapStatistics.Length;
            var lower = new double[nParams];
            var upper = new double[nParams];

            double alpha = (1.0 - confidenceLevel) / 2.0;
            int iLo = Math.Max(0, (int)Math.Floor(alpha * n));
            int iHi = Math.Min(n - 1, (int)Math.Ceiling((1.0 - alpha) * n) - 1);

            for (int j = 0; j < nParams; j++)
            {
                double[] values = _bootstrapStatistics.Select(s => s[j]).ToArray();
                Array.Sort(values);
                lower[j] = values[iLo];
                upper[j] = values[iHi];
            }

            return (lower, upper);
        }
    }
}
