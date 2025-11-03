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

using System;
using System.Collections.Generic;
using System.Linq;
using Numerics.Distributions;

namespace Numerics.Sampling
{
    /// <summary>
    /// A class for stratifying probabilities for sampling, or values for numerical integration. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <b> References:</b>
    /// <list type="bullet">
    /// <item>
    /// <see href = "https://en.wikipedia.org/wiki/Stratified_sampling" />
    /// </item>
    /// </list>
    /// </remarks>
    public class Stratify
    {

        /// <summary>
        /// Returns a list of stratified x value bins.
        /// </summary>
        /// <param name="options">The stratification options.</param>
        /// <param name="isLogarithmic">Optional. Determines if the x values should be stratified with a linear or logarithmic (base 10) scale. Default = False.</param>
        public static List<StratificationBin> XValues(StratificationOptions options, bool isLogarithmic = false)
        {
            var bins = new List<StratificationBin>();
            if (!options.IsValid || options.IsProbability)
                return bins;

            double delta, offset, min, max, xl, xu;
            
            if (isLogarithmic)
            {
                // first determine the offset
                offset = options.LowerBound <= 0d ? 1.0d - options.LowerBound : 0.0d;
                // get delta
                min = Math.Log10(options.LowerBound + offset);
                max = Math.Log10(options.UpperBound + offset);
                delta = (max - min) / options.NumberOfBins;
                // stratify first bin
                xl = Math.Pow(10d, min) - offset;
                xu = Math.Pow(10d, Math.Log10(xl + offset) + delta) - offset;
                bins.Add(new StratificationBin(xl, xu));
                // stratify remaining bins
                for (int i = 1; i < options.NumberOfBins; i++)
                {
                    xl = xu;
                    xu = Math.Pow(10d, Math.Log10(xl + offset) + delta) - offset;
                    bins.Add(new StratificationBin(xl, xu));
                }
            }
            else
            {
                // get delta
                delta = (options.UpperBound - options.LowerBound) / options.NumberOfBins;
                // stratify first bin
                xl = options.LowerBound;
                xu = xl + delta;
                bins.Add(new StratificationBin(xl, xu));
                // stratify remaining bins
                for (int i = 1; i < options.NumberOfBins; i++)
                {
                    xl = xu;
                    xu = xl + delta;
                    bins.Add(new StratificationBin(xl, xu));
                }
            }
             
            return bins;
        }

        /// <summary>
        /// Returns a list of stratified x value bins.
        /// </summary>
        /// <param name="options">a list of stratification options.</param>
        /// <param name="isLogarithmic">Optional. Determines if the x values should be stratified with a linear or logarithmic (base 10) scale. Default = False.</param>
        public static List<StratificationBin> XValues(List<StratificationOptions> options, bool isLogarithmic = false)
        {
            if (options == null || options.Count == 0)
                return new List<StratificationBin>();

            var sortedOptions = options.ToArray();
            Array.Sort(sortedOptions, (x, y) => x.UpperBound.CompareTo(y.UpperBound));
            var bins = new List<StratificationBin>();
            foreach (var opt in sortedOptions)
                bins.AddRange(XValues(opt, isLogarithmic));
            return bins;
        }

        /// <summary>
        /// Returns a list of stratified exceedance probability bins given a list of stratified x value bins and a function to transform x values to exceedance probabilities.
        /// The transform function must convert x values in ascending order to exceedance probabilities in descending order.
        /// </summary>
        /// <param name="xValues">List of stratified x value bins.</param>
        /// <param name="transformFunction">Function to transform x values to exceedance probabilities.
        /// The transform function must convert x values in ascending order to exceedance probabilities in descending order.</param>
        /// <param name="isExhaustive">Determines if the probability bin weights should be collectively exhaustive (i.e., sum to 1).</param>
        public static List<StratificationBin> XToExceedanceProbability(List<StratificationBin> xValues, Func<double, double> transformFunction, bool isExhaustive = true)
        {
            var bins = new List<StratificationBin>();
            if (xValues == null || xValues.Count == 0)
                return bins;

            double xl, xu;
            xu = transformFunction(xValues.First().LowerBound);
            xl = transformFunction(xValues.First().UpperBound);
            bins.Add(new StratificationBin(xl, xu));

            // transform the bins
            for (int i = 1; i < xValues.Count; i++)
            {
                xu = xl;
                xl = transformFunction(xValues[i].UpperBound);
                bins.Add(new StratificationBin(xl, xu));
            }

            if (isExhaustive)
            {
                // adjust end bins to ensure exhaustivity
                bins.First().Weight = 1.0d - bins.First().LowerBound;
                bins.Last().Weight = bins.Last().UpperBound;
            }

            return bins;
        }

        /// <summary>
        /// Returns a list of stratified non-exceedance probability bins given a list of stratified x value bins and a function to transform x values to non-exceedance probabilities.
        /// The transform function must convert x values in ascending order to non-exceedance probabilities in ascending order.
        /// </summary>
        /// <param name="xValues">List of stratified x value bins.</param>
        /// <param name="transformFunction">Function to transform x values to non-exceedance probabilities.
        /// The transform function must convert x values in ascending order to non-exceedance probabilities in ascending order.</param>
        /// <param name="isExhaustive">Determines if the probability bin weights should be collectively exhaustive (i.e., sum to 1).</param>
        public static List<StratificationBin> XToProbability(List<StratificationBin> xValues, Func<double, double> transformFunction, bool isExhaustive = true)
        {
            var bins = new List<StratificationBin>();
            if (xValues == null || xValues.Count == 0)
                return bins;

            double xl, xu;
            // transform first bin
            xl = transformFunction(xValues.First().LowerBound);
            xu = transformFunction(xValues.First().UpperBound);
            // determine if the first bin weight should ensure exhaustivity
            bins.Add(new StratificationBin(xl, xu, isExhaustive ? xu : xu - xl));

            // transform the remaining inner bins
            for (int i = 1; i < xValues.Count - 1; i++)
            {
                xl = xu;
                xu = transformFunction(xValues[i].UpperBound);
                bins.Add(new StratificationBin(xl, xu));
            }
            // transform the last bin
            xl = xu;
            xu = transformFunction(xValues.Last().UpperBound);
            // determine if the last bin weight should ensure exhaustivity
            bins.Add(new StratificationBin(xl, xu, isExhaustive ? 1.0d - xl : xu - xl));
            
            return bins;
        }

        /// <summary>
        /// Returns a list of stratified probability bins.
        /// </summary>
        /// <param name="options">The stratification options.</param>
        /// <param name="distributionType">Optional. The importance distribution type to stratify with. Default = Uniform.</param>
        /// <param name="isExhaustive">Determines if the probability bin weights should be collectively exhaustive (i.e., sum to 1).</param>
        public static List<StratificationBin> Probabilities(StratificationOptions options, ImportanceDistribution distributionType = ImportanceDistribution.Uniform, bool isExhaustive = true)
        {
            var bins = new List<StratificationBin>();
            if (options == null || !options.IsValid || !options.IsProbability)
                return bins;

            double delta, offset = 0, min, max, pl, pu, xl, xu, w;

            UnivariateDistributionBase distribution = distributionType switch
            {
                ImportanceDistribution.Gumbel => new Gumbel(0d, 1d),
                ImportanceDistribution.Normal => new Normal(0d, 1d),
                ImportanceDistribution.Uniform => new Uniform(0d, 1d),
                ImportanceDistribution.Log10Uniform => new Uniform(0d, 1d), // same distribution, stratified differently
                _ => throw new ArgumentOutOfRangeException(nameof(distributionType), "Invalid importance distribution type.")
            };

            if (distributionType == ImportanceDistribution.Log10Uniform)
            {
                // first determine the offset
                if (options.LowerBound == 0d)
                    offset = 1.0d;
                // get delta
                min = Math.Log10(distribution.InverseCDF(options.LowerBound) + offset);
                max = Math.Log10(distribution.InverseCDF(options.UpperBound) + offset);
                delta = (max - min) / options.NumberOfBins;
                // stratify first bin
                xl = Math.Pow(10d, min) - offset;
                xu = Math.Pow(10d, Math.Log10(xl + offset) + delta) - offset;
                pl = distribution.CDF(xl);
                pu = distribution.CDF(xu);
                // determine if the first bin weight should ensure exhaustivity
                bins.Add(new StratificationBin(pl, pu, isExhaustive ? pu : pu - pl));
                
                // stratify inner bins
                for (int i = 1; i < options.NumberOfBins - 1; i++)
                {
                    xl = xu;
                    xu = Math.Pow(10d, Math.Log10(xl + offset) + delta) - offset;
                    pl = distribution.CDF(xl);
                    pu = distribution.CDF(xu);
                    bins.Add(new StratificationBin(pl, pu));
                    // The weight of the inner intervals is simply the 
                    // width of the probability interval under consideration.
                }
                // stratify last bin
                xl = xu;
                xu = Math.Pow(10d, Math.Log10(xl + offset) + delta) - offset;
                pl = distribution.CDF(xl);
                pu = distribution.CDF(xu);
                w = isExhaustive ? 1.0d - pl : pu - pl;
                bins.Add(new StratificationBin(pl, pu, w));
            }
            else
            {
                // get delta
                min = distribution.InverseCDF(options.LowerBound);
                max = distribution.InverseCDF(options.UpperBound);
                delta = (max - min) / options.NumberOfBins;
                // stratify first bin
                xl = min;
                xu = xl + delta;
                pl = distribution.CDF(xl);
                pu = distribution.CDF(xu);
                // determine if the first bin weight should ensure exhaustivity
                bins.Add(new StratificationBin(pl, pu, isExhaustive ? pu : pu - pl));
                // stratify inner bins
                for (int i = 1; i < options.NumberOfBins - 1; i++)
                {
                    xl = xu;
                    xu = xl + delta;
                    pl = distribution.CDF(xl);
                    pu = distribution.CDF(xu);
                    bins.Add(new StratificationBin(pl, pu));
                    // The weight of the inner intervals is simply the 
                    // width of the probability interval under consideration.
                }
                // stratify last bin
                xl = xu;
                xu = xl + delta;
                pl = distribution.CDF(xl);
                pu = distribution.CDF(xu);
                w = isExhaustive ? 1.0d - pl : pu - pl;
                bins.Add(new StratificationBin(pl, pu, w));
            }
            
            return bins;
        }

        /// <summary>
        /// Returns a multivariate list of stratified probability bins.
        /// </summary>
        /// <param name="options">The stratification options.</param>
        /// <param name="distributionType">Optional. The importance distribution type to stratify with. Default = Uniform.</param>
        /// <param name="isExhaustive">Determines if the probability bin weights should be collectively exhaustive (i.e., sum to 1).</param>
        /// <param name="dimension">The number of dimensions to stratify.</param>
        /// <param name="seed"> Seed for random number generator. </param>
        /// <param name="correlation">The correlation matrix. If null, independence is assumed.</param>
        public static List<List<StratificationBin>> MultivariateProbabilities(StratificationOptions options, ImportanceDistribution distributionType = ImportanceDistribution.Uniform, bool isExhaustive = true, int dimension = 1, int seed = -1, double[,] correlation = null)
        {
            // Validate inputs
            var output = new List<List<StratificationBin>>();
            if (options == null || !options.IsValid || !options.IsProbability || dimension < 1)
                return output;

            // Generate univariate bins for each dimension
            var binPool = new List<List<StratificationBin>>();
            for (int i = 0; i < dimension; i++)
            {
                var bins = Probabilities(options, distributionType, isExhaustive);
                if (bins.Count < options.NumberOfBins)
                    throw new InvalidOperationException($"Insufficient bins ({bins.Count}) for dimension {i} to sample {options.NumberOfBins} points without replacement.");
                binPool.Add(bins);
                output.Add(new List<StratificationBin>());
            }

            // Generate correlation matrix if null (identity = independence)
            if (correlation == null)
            {
                correlation = new double[dimension, dimension];
                for (int i = 0; i < dimension; i++)
                    for (int j = 0; j < dimension; j++)
                        correlation[i, j] = i == j ? 1.0 : 0.0;
            }

            // Create multivariate normal distribution
            var mean = new double[dimension]; // all zero
            var mvn = new MultivariateNormal(dimension);
            mvn.SetParameters(mean, correlation);

            // Latin Hypercube Samples (LHS) in standard normal space
            var lhsSamples = mvn.LatinHypercubeRandomValues(options.NumberOfBins, seed);

            // Sample bins without replacement in each dimension
            var rng = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            for (int row = 0; row < options.NumberOfBins; row++)
            {
                for (int col = 0; col < dimension; col++)
                {
                    var u = Normal.StandardCDF(lhsSamples[row, col]);
                    var availableBins = binPool[col];

                    int r = (int)Math.Floor(u * availableBins.Count);
                    r = Math.Min(r, availableBins.Count - 1); // Clamp in case u == 1.0

                    output[col].Add((StratificationBin)availableBins[r].Clone());
                    availableBins.RemoveAt(r);
                }
            }

            return output;
        }

        /// <summary>
        /// Returns a list of stratified probability bins.
        /// </summary>
        /// <param name="options">The list of stratification options.</param>
        /// <param name="distributionType">Optional. The importance distribution type to stratify with. Default = Uniform.</param>
        /// <param name="isExhaustive">Determines if the probability bin weights should be collectively exhaustive (i.e., sum to 1).</param>
        public static List<StratificationBin> Probabilities(List<StratificationOptions> options, ImportanceDistribution distributionType = ImportanceDistribution.Uniform, bool isExhaustive = true)
        {
            if (options == null || options.Count == 0)
                return new List<StratificationBin>();

            var sortedOptions = options.ToArray();
            Array.Sort(sortedOptions, (x, y) => x.UpperBound.CompareTo(y.UpperBound));

            var bins = new List<StratificationBin>();
            for (int i = 0; i < sortedOptions.Length; i++)
                bins.AddRange(Probabilities(sortedOptions[i], distributionType, false));

            if (isExhaustive)
            {
                // adjust end bins to ensure exhaustivity
                bins.First().Weight = bins.First().UpperBound;
                bins.Last().Weight = 1.0d - bins.Last().LowerBound;
            }

            return bins;
        }

        /// <summary>
        /// Returns a list of stratified exceedance probability bins.
        /// </summary>
        /// <param name="options">The stratification options.</param>
        /// <param name="distributionType">Optional. The importance distribution type to stratify with. Default = Uniform.</param>
        /// <param name="isExhaustive">Determines if the probability bin weights should be collectively exhaustive (i.e., sum to 1).</param>
        public static List<StratificationBin> ExceedanceProbabilities(StratificationOptions options, ImportanceDistribution distributionType = ImportanceDistribution.Uniform, bool isExhaustive = true)
        {
            var bins = new List<StratificationBin>();
            if (options == null || !options.IsValid || !options.IsProbability)
                return bins;

            double delta, offset = 0, min, max, pl, pu, xl, xu;

            UnivariateDistributionBase distribution = distributionType switch
            {
                ImportanceDistribution.Gumbel => new Gumbel(0d, 1d),
                ImportanceDistribution.Normal => new Normal(0d, 1d),
                ImportanceDistribution.Uniform => new Uniform(0d, 1d),
                ImportanceDistribution.Log10Uniform => new Uniform(0d, 1d),
                _ => throw new ArgumentOutOfRangeException(nameof(distributionType), "Invalid importance distribution type.")
            };

            if (distributionType == ImportanceDistribution.Log10Uniform)
            {
                // first determine the offset
                if (options.LowerBound == 0d)
                    offset = 1.0d;
                // get delta
                min = Math.Log10(distribution.InverseCDF(options.LowerBound) + offset);
                max = Math.Log10(distribution.InverseCDF(options.UpperBound) + offset);
                delta = (max - min) / options.NumberOfBins;
                // stratify first bin
                xu = Math.Pow(10d, max) - offset;
                xl = Math.Pow(10d, Math.Log10(xu + offset) - delta) - offset;
                pl = distribution.CDF(xl);
                pu = distribution.CDF(xu);
                // determine if the first bin weight should ensure exhaustivity
                bins.Add(new StratificationBin(pl, pu, isExhaustive ? 1.0d - pl : pu - pl));
                // stratify inner bins
                for (int i = 1; i < options.NumberOfBins - 1; i++)
                {
                    xu = xl;
                    xl = Math.Pow(10d, Math.Log10(xu + offset) - delta) - offset;
                    pl = distribution.CDF(xl);
                    pu = distribution.CDF(xu);
                    bins.Add(new StratificationBin(pl, pu));
                    // The weight of the inner intervals is simply the 
                    // width of the probability interval under consideration.
                }
                // stratify last bin
                xu = xl;
                xl = Math.Pow(10d, Math.Log10(xu + offset) - delta) - offset;
                pl = distribution.CDF(xl);
                pu = distribution.CDF(xu);
                bins.Add(new StratificationBin(pl, pu, isExhaustive ? pu : pu - pl));
            }
            else
            {
                // get delta
                min = distribution.InverseCDF(options.LowerBound);
                max = distribution.InverseCDF(options.UpperBound);
                delta = (max - min) / options.NumberOfBins;
                // stratify first bin
                xu = max;
                xl = xu - delta;
                pl = distribution.CDF(xl);
                pu = distribution.CDF(xu);
                // determine if the first bin weight should ensure exhaustivity
                bins.Add(new StratificationBin(pl, pu, isExhaustive ? 1.0d - pl : pu - pl));
                // stratify inner bins
                for (int i = 1; i < options.NumberOfBins - 1; i++)
                {
                    xu = xl;
                    xl = xu - delta;
                    pl = distribution.CDF(xl);
                    pu = distribution.CDF(xu);
                    bins.Add(new StratificationBin(pl, pu));
                    // The weight of the inner intervals is simply the 
                    // width of the probability interval under consideration.
                }
                // stratify last bin
                xu = xl;
                xl = xu - delta;
                pl = distribution.CDF(xl);
                pu = distribution.CDF(xu);
                bins.Add(new StratificationBin(pl, pu, isExhaustive ? pu : pu - pl));
            }
            
            return bins;
        }

        /// <summary>
        /// Returns a list of stratified exceedance probability bins.
        /// </summary>
        /// <param name="options">The list of stratification options.</param>
        /// <param name="distributionType">Optional. The importance distribution type to stratify with. Default = Uniform.</param>
        /// <param name="isExhaustive">Determines if the probability bin weights should be collectively exhaustive (i.e., sum to 1).</param>
        public static List<StratificationBin> ExceedanceProbabilities(List<StratificationOptions> options, ImportanceDistribution distributionType = ImportanceDistribution.Uniform, bool isExhaustive = true)
        {
            if (options == null || options.Count == 0)
                return new List<StratificationBin>();

            var sortedOptions = options.ToArray();
            Array.Sort(sortedOptions, (x, y) => -x.UpperBound.CompareTo(y.UpperBound));

            var bins = new List<StratificationBin>();
            for (int i = 0; i < sortedOptions.Length; i++)
                bins.AddRange(ExceedanceProbabilities(sortedOptions[i], distributionType, false));

            if (isExhaustive)
            {
                // adjust end bins to ensure exhaustivity
                bins.First().Weight = 1.0d - bins.First().LowerBound;
                bins.Last().Weight = bins.Last().UpperBound;
            }

            return bins;
        }

        /// <summary>
        /// Returns a list of stratified x value bins given a list of stratified non-exceedance probability bins and a function to transform non-exceedance probabilities to x values.
        /// </summary>
        /// <param name="probabilities">List of stratified probabilities.</param>
        /// <param name="transformFunction">Function to transform non-exceedance probabilities to x values.</param>
        public static List<StratificationBin> ProbabilityToX(List<StratificationBin> probabilities, Func<double, double> transformFunction)
        {
            var bins = new List<StratificationBin>();
            if (probabilities == null || probabilities.Count == 0)
                return bins;

            double xl, xu;
            xl = transformFunction(probabilities.First().LowerBound);
            xu = transformFunction(probabilities.First().UpperBound);
            bins.Add(new StratificationBin(xl, xu));

            // transform the bins
            for (int i = 1; i < probabilities.Count; i++)
            {
                xl = xu;
                xu = transformFunction(probabilities[i].UpperBound);
                bins.Add(new StratificationBin(xl, xu));
            }

            return bins;
        }

        /// <summary>
        /// Returns a list of stratified x value bins given a list of stratified exceedance probability bins and a function to transform exceedance probabilities to x values.
        /// </summary>
        /// <param name="probabilities">List of stratified exceedance probabilities.</param>
        /// <param name="transformFunction">Function to transform exceedance probabilities to x values.</param>
        public static List<StratificationBin> ExceedanceProbabilityToX(List<StratificationBin> probabilities, Func<double, double> transformFunction)
        {
            var bins = new List<StratificationBin>();
            if (probabilities == null || probabilities.Count == 0)
                return bins;

            double xl, xu;
            xl = transformFunction(probabilities.First().UpperBound);
            xu = transformFunction(probabilities.First().LowerBound);
            bins.Add(new StratificationBin(xl, xu));

            // transform the bins
            for (int i = 1; i < probabilities.Count; i++)
            {
                xl = xu;
                xu = transformFunction(probabilities[i].LowerBound);
                bins.Add(new StratificationBin(xl, xu));
            }

            return bins;
        }

        /// <summary>
        /// Enumeration of importance distributions for stratified probabilities.
        /// </summary>
        public enum ImportanceDistribution
        {
            /// <summary>
            /// The Gumbel (Extreme Value Type I) distribution.
            /// </summary>
            Gumbel,
            /// <summary>
            /// The normal (Gaussian) distribution.
            /// </summary>
            Normal,
            /// <summary>
            /// The uniform distribution.
            /// </summary>
            Uniform,
            /// <summary>
            /// The log (base 10) uniform distribution.
            /// </summary>
            Log10Uniform
        }
    }
}