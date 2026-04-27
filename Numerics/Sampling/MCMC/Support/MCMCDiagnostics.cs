using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Sampling.MCMC
{
    /// <summary>
    /// A class for assessing Bayesian MCMC convergence diagnostics.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class MCMCDiagnostics
    {

        /// <summary>
        /// Compute the effective sample size.
        /// </summary>
        /// <param name="series">The series of posterior samples to evaluate.</param>
        public static double EffectiveSampleSize(IList<double> series)
        {
            //https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.4/topics/ESS
            int N = series.Count;
            var acf = Fourier.Autocorrelation(series, (int)Math.Ceiling((double)N / 2));
            if (acf == null) return N;
            double rho = 0;
            for (int i = 1; i < acf.GetLength(0); i++)
            {
                if (acf[i, 1] < 0.0) break;
                rho += acf[i, 1];
            }
            return Math.Min(N / (1d + 2d * rho), N);
        }

        /// <summary>
        /// Computes the effective samples size for each model parameter.
        /// </summary>
        /// <param name="markovChains">The list of Markov Chains to be evaluated. The chains must be of equal length.</param>
        /// <param name="averageACF">Output. A jagged array of averaged autocorrelation functions, one for each parameter.</param>
        public static double[] EffectiveSampleSize(IList<List<ParameterSet>> markovChains, out double[][,] averageACF)
        {
            // Get number of chains
            int M = markovChains.Count;
            if (M == 0) throw new ArgumentException(nameof(markovChains), "No chains provided.");
            // Get number of parameters
            int P = markovChains[0][0].Values.Length;
            // Get the minimum number of iterations
            int N = markovChains.Min(chain => chain.Count);

            // Validation checks
            if (N < 2) throw new ArgumentOutOfRangeException(nameof(markovChains), "There must be at least two iterations to evaluate.");
            if (P < 1) throw new ArgumentOutOfRangeException(nameof(markovChains), "There must be at least one parameter to evaluate.");

            // Create result arrays. 
            var ESS= new double[P];
            averageACF = new double[P][,];

            for (int p = 0; p < P; p++)
            {
                // Compute the Autocorrelation Function (ACF) and Effective Sample Size (ESS)
                // Average the ACF and sum the ESS
                averageACF[p] = new double[51, 2];
                double meanRho = 0;

                for (int i = 0; i < M; i++)
                {
                    // Get values for this parameter within this chain
                    var values = markovChains[i].Select(set => set.Values[p]).ToArray();

                    // Get ACF for this chain
                    var acf = Fourier.Autocorrelation(values, (int)Math.Ceiling((double)N / 2));
                    if (acf == null) continue;
                    // Update the average ACF across all chains
                    for (int j = 0; j < acf.GetLength(0); j++)
                    {
                        if (j > 50) break;
                        averageACF[p][j, 1] += acf[j, 1] / M;
                    }
                    // https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.4/topics/ESS
                    // Get rho for the current chain
                    double rho = 0;
                    for (int j = 1; j < acf.GetLength(0); j++)
                    {
                        if (acf[j, 1] < 0.0) break;
                        rho += acf[j, 1];
                    }
                    meanRho += rho / M;
                    
                }
                ESS[p] += Math.Min(N * M / (1d + 2d * meanRho), N * M);
            }

            return ESS;
        }

        /// <summary>
        /// The Gelman-Rubin diagnostic.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     The Gelman-Rubin diagnostic tests for lack of convergence by comparing the variance between multiple chains
        ///     to the variance within each chain. If convergence has been achieved, the between-chain and within-chain
        ///     variances should be identical. To be most effective in detecting evidence for non convergence, each chain should
        ///     have been initialized to starting values that are dispersed relative to the target distribution.
        /// </para>
        /// </remarks>
        /// <param name="markovChains">The list of Markov Chains to be evaluated. The chains must be of equal length.</param>
        /// <param name="warmupIterations">The number of warm up MCMC iterations to discard at the beginning of the chains.</param>
        public static double[] GelmanRubin(IList<List<ParameterSet>> markovChains, int warmupIterations = 0)
        {

            // Get number of chains
            int M = markovChains.Count;
            if (M == 0) throw new ArgumentException(nameof(markovChains), "No chains provided.");
            // Get number of parameters
            int P = markovChains[0][0].Values.Length;
            // Get the number of iterations
            int N = markovChains[0].Count;
            foreach (var chain in markovChains)
            {
                if (chain.Count != N)
                    throw new ArgumentException("All chains must have the same length.");
            }

            // Create result array. 
            var Rhat = new double[P];
            Rhat.Fill(double.NaN);


            // Validation checks
            if (M < 2) return Rhat;           
            if (N < 2) throw new ArgumentOutOfRangeException(nameof(markovChains), "There must be at least two iterations to evaluate.");
            if (P < 1) throw new ArgumentOutOfRangeException(nameof(markovChains), "There must be at least one parameter to evaluate.");
            if (warmupIterations < 0) throw new ArgumentOutOfRangeException(nameof(warmupIterations), "The warm up iterations must be non-negative.");
            int startIndex = Math.Max(0, warmupIterations);

            // Compute R-hat or each parameter
            for (int p = 0; p < P; p++)
            {
                // Step 1. Compute between- and within-chain mean
                var chainMeans = new double[M];
                double overallMean = 0;         
                for (int i = 0; i < M; i++)
                {
                    for (int j = startIndex; j < N; j++)
                    {
                        chainMeans[i] += markovChains[i][j].Values[p];
                    }
                    // Get within-chain mean
                    chainMeans[i] /= (N - startIndex);
                    overallMean += chainMeans[i];
                }
                // Get between-chain mean
                overallMean /= M;

                // Step 2. Compute between- and within-chain variance
                int n = N - startIndex;
                double B = 0, W = 0;
                for (int i = 0; i < M; i++)
                {
                    double sum = 0;
                    for (int j = startIndex; j < N; j++)
                    {
                        sum += Tools.Sqr(markovChains[i][j].Values[p] - chainMeans[i]);
                    }
                    // within-chain variance
                    W += sum / (n - 1);
                    // between-chain variance
                    B += Tools.Sqr(chainMeans[i] - overallMean);
                }
                // Set between- and within-chain variance
                W /= M;
                B *= n / (double)(M - 1);

                // Step 3. Compute the pooled variance
                double V = ((n - 1d) / n) * W + (1d / n) * B;

                // Step 4. Compute R-hat
                Rhat[p] = Math.Sqrt(V / W);
            }

            return Rhat;
        }

        /// <summary>
        /// Computes the minimum sample size rounded to the nearest 100 based on the Raftery-Lewis method.
        /// </summary>
        /// <param name="quantile">The posterior quantile of interest; e.g., 0.975.</param>
        /// <param name="tolerance">The acceptable tolerance for this quantile; e.g., ±0.005.</param>
        /// <param name="probability">Probability of being within the range of tolerance; e.g., 0.95.</param>
        /// <returns>The minimum sample size as an integer, rounded to the nearest thousands place.</returns>
        public static int MinimumSampleSize(double quantile, double tolerance, double probability)
        {
            double q = quantile;
            double r = tolerance;
            double s = probability;
            double N = (q * (1d - q) * Math.Pow(Normal.StandardZ(0.5d * (s + 1d)), 2d)) / Math.Pow(r, 2d);
            int Nmin = (int)Math.Round(decimal.Parse(N.ToString()) / 100M, 0) * 100;
            return Nmin;
        }

    }

   
}
