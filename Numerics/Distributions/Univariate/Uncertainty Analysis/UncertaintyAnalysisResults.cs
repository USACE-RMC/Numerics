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

using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling.MCMC;
using Numerics.Utilities;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Serialization;
using System.Threading.Tasks;
using System.Xml.Linq;

namespace Numerics.Distributions
{

    /// <summary>
    /// A class for storing distribution uncertainty analysis results.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class UncertaintyAnalysisResults
    {

        /// <summary>
        /// Construct an instance of the UncertaintyAnalysisResults class.
        /// </summary>
        public UncertaintyAnalysisResults() { }

        /// <summary>
        /// Constructs a new instance with computed uncertainty metrics.
        /// </summary>
        /// <param name="parentDistribution">The parent (mode/point estimate) distribution.</param>
        /// <param name="sampledDistributions">Array of sampled distributions from posterior or bootstrap.</param>
        /// <param name="probabilities">Array of non-exceedance probabilities for quantile estimation.</param>
        /// <param name="alpha">The confidence level (default = 0.1 for 90% CI).</param>
        /// <param name="minProbability">Minimum probability for mean curve computation (default = 0.001).</param>
        /// <param name="maxProbability">Maximum probability for mean curve computation (default = 1 - 1e-9).</param>
        /// <param name="recordParameterSets">If true, stores all parameter sets from sampled distributions.</param>
        public UncertaintyAnalysisResults(UnivariateDistributionBase parentDistribution,
                                          UnivariateDistributionBase[] sampledDistributions,
                                          double[] probabilities,
                                          double alpha = 0.1,
                                          double minProbability = 0.001,
                                          double maxProbability = 1 - 1e-9,
                                          bool recordParameterSets = false)
        {
            if (parentDistribution == null!)
                throw new ArgumentNullException(nameof(parentDistribution));
            if (sampledDistributions == null || sampledDistributions.Length == 0)
                throw new ArgumentException("Sampled distributions cannot be null or empty.", nameof(sampledDistributions));
            if (probabilities == null || probabilities.Length == 0)
                throw new ArgumentException("Probabilities cannot be null or empty.", nameof(probabilities));

            ProcessModeCurve(parentDistribution, probabilities);
            ProcessConfidenceIntervals(sampledDistributions, probabilities, alpha);
            ProcessMeanCurve(sampledDistributions, probabilities, minProbability, maxProbability);

            if (recordParameterSets)
                ProcessParameterSets(sampledDistributions);

            // Set default values
            AIC = double.NaN;
            BIC = double.NaN;
            DIC = double.NaN;
            RMSE = double.NaN;
            ERL = double.NaN;
        }

        /// <summary>
        /// The parent probability distribution.
        /// </summary>
        public UnivariateDistributionBase ParentDistribution { get; set; } = null!;

        /// <summary>
        /// The array of parameter sets.
        /// </summary>
        public ParameterSet[] ParameterSets { get; set; } = null!;

        /// <summary>
        /// The confidence intervals. 
        /// </summary>
        public double[,] ConfidenceIntervals { get; set; } = null!;

        /// <summary>
        /// The mode (or computed) curve from the parent distribution. 
        /// </summary>
        public double[] ModeCurve { get; set; } = null!;

        /// <summary>
        /// The mean (or predictive) curve. 
        /// </summary>
        public double[] MeanCurve { get; set; } = null!;

        /// <summary>
        /// Gets or sets the Akaike information criteria (AIC) of the fit.
        /// </summary>
        public double AIC { get; set; }

        /// <summary>
        /// Gets or sets the Bayesian information criteria (BIC) of the fit.
        /// </summary>
        public double BIC { get; set; }

        /// <summary>
        /// Gets or sets the Deviance Information Criterion (DIC) of the fit. 
        /// </summary>
        public double DIC { get; set; }

        /// <summary>
        /// Gets or sets the Root Mean Square Error (RMSE) of the fit. 
        /// </summary>
        public double RMSE { get; set; }

        /// <summary>
        /// Gets or sets the Effective Record Length (ERL).
        /// </summary>
        public double ERL { get; set; }

        /// <summary>
        /// Returns the class as a byte array. 
        /// </summary>
        /// <param name="results">The uncertainty analysis results.</param>
        public static byte[] ToByteArray(UncertaintyAnalysisResults results)
        {
            var options = new JsonSerializerOptions
            {
                WriteIndented = false,
                DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
                IncludeFields = true
            };
            // Add custom converters for unsupported types
            options.Converters.Add(new Double2DArrayConverter());
            options.Converters.Add(new String2DArrayConverter());
            options.Converters.Add(new UnivariateDistributionConverter());
            return JsonSerializer.SerializeToUtf8Bytes(results, options);
        }

        /// <summary>
        /// Returns the class from a byte array. 
        /// </summary>
        /// <param name="bytes">Byte array.</param>
        public static UncertaintyAnalysisResults FromByteArray(byte[] bytes)
        {
            try
            {
                var options = new JsonSerializerOptions
                {
                    DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
                    IncludeFields = true
                };
                // Add custom converters for unsupported types
                options.Converters.Add(new Double2DArrayConverter());
                options.Converters.Add(new String2DArrayConverter());
                options.Converters.Add(new UnivariateDistributionConverter());
                var result = JsonSerializer.Deserialize<UncertaintyAnalysisResults>(bytes, options);
                return result ?? FromByteArrayLegacy(bytes);
            }
            catch (Exception)
            {
                // An error can occur because we're trying to deserialize a blob written with binary formatter, 
                //as a blob of json bytes. If that happens, fall back to the old.  
                return FromByteArrayLegacy(bytes);
            }          
        }

         /// <summary>
        /// Returns the class from a byte array. 
        /// </summary>
        /// <param name="bytes">Byte array.</param>
        private static UncertaintyAnalysisResults FromByteArrayLegacy(byte[] bytes)
        {
            try
            {
                using (var memStream = new MemoryStream())
                {
                    #pragma warning disable SYSLIB0011 // Suppress obsolete BinaryFormatter warning for legacy support
                    var binForm = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
                    memStream.Write(bytes, 0, bytes.Length);
                    memStream.Seek(0, SeekOrigin.Begin);
                    var obj = binForm.Deserialize(memStream);
                    #pragma warning disable SYSLIB0011 // Suppress obsolete BinaryFormatter warning for legacy support
                    return (UncertaintyAnalysisResults)obj;
                }


            }
            catch (Exception)
            {
                // An error can occur because of differences in versions. 
                // If there is an error, just catch it and force the user to rerun the
                // uncertainty analysis. 
            }
            return null!;
        }

        /// <summary>
        /// Returns the class as XElement. The parameter sets are not included. 
        /// </summary>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(UncertaintyAnalysisResults));         
            if (ParentDistribution != null!) result.Add(ParentDistribution.ToXElement());
            result.SetAttributeValue(nameof(AIC), AIC.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(BIC), BIC.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(DIC), DIC.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(RMSE), RMSE.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(ERL), ERL.ToString("G17", CultureInfo.InvariantCulture));
            if (ModeCurve != null) result.SetAttributeValue(nameof(ModeCurve), String.Join("|", ModeCurve.Select(d => d.ToString("G17", CultureInfo.InvariantCulture))));
            if (MeanCurve != null) result.SetAttributeValue(nameof(MeanCurve), String.Join("|", MeanCurve.Select(d => d.ToString("G17", CultureInfo.InvariantCulture))));
            // CIs
            if (ConfidenceIntervals != null)
            {
                string CIstring = "";
                for (int i = 0; i < ConfidenceIntervals.GetLength(0); i++)
                {
                    for (int j = 0; j < ConfidenceIntervals.GetLength(1); j++)
                    {
                        CIstring += ConfidenceIntervals[i, j].ToString("G17", CultureInfo.InvariantCulture) + (j < ConfidenceIntervals.GetLength(1) - 1 ? "," : "");
                    }
                    CIstring += (i < ConfidenceIntervals.GetLength(0) - 1 ? "|" : "");
                }
                result.SetAttributeValue(nameof(ConfidenceIntervals), CIstring);
            }
            return result;
        }

        /// <summary>
        /// Returns the class from an XElement. 
        /// </summary>
        /// <param name="xElement">XElement to deserialize.</param>
        public static UncertaintyAnalysisResults FromXElement(XElement xElement)
        {
            var ua = new UncertaintyAnalysisResults();    
            // Parent distribution
            var distElement = xElement.Element("Distribution");
            if (distElement != null)
            {
                var parentDist = UnivariateDistributionFactory.CreateDistribution(distElement);
                if (parentDist is not null)
                {
                    ua.ParentDistribution = parentDist;
                }
                else
                {
                    throw new InvalidDataException("Unable to deserialize parent distribution from XElement.");
                }    
            }
                

            // AIC
            var aicElement = xElement.Attribute(nameof(AIC));
            if (aicElement != null)
            {
                double.TryParse(aicElement.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var aic);
                ua.AIC = aic;
            }
            // BIC
            var bicElement = xElement.Attribute(nameof(BIC));
            if (bicElement != null)
            {
                double.TryParse(bicElement.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var bic);
                ua.BIC = bic;
            }
            // DIC
            var dicElement = xElement.Attribute(nameof(DIC));
            if (dicElement != null)
            {
                double.TryParse(dicElement.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var dic);
                ua.DIC = dic;
            }
            // RMSE
            var rmseElement = xElement.Attribute(nameof(RMSE));
            if (rmseElement != null)
            {
                double.TryParse(rmseElement.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var rmse);
                ua.RMSE = rmse;
            }
            // ERL
            var erlElement = xElement.Attribute(nameof(ERL));
            if (erlElement != null)
            {
                double.TryParse(erlElement.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var erl);
                ua.ERL = erl;
            }

            // Mode Curve
            var modeAttr = xElement.Attribute(nameof(ua.ModeCurve));
            if (modeAttr != null)
            {
                var vals = modeAttr.Value.Split('|');
                if (vals.Length > 0)
                {
                    ua.ModeCurve = new double[vals.Length];
                    for (int i = 0; i < vals.Length; i++)
                    {
                        double.TryParse(vals[i], NumberStyles.Any, CultureInfo.InvariantCulture, out ua.ModeCurve[i]);
                    }
                }
            }
            // Mean Curve
            var meanAttr = xElement.Attribute(nameof(ua.MeanCurve));
            if (meanAttr != null)
            {
                var vals = meanAttr.Value.Split('|');
                if (vals.Length > 0)
                {
                    ua.MeanCurve = new double[vals.Length];
                    for (int i = 0; i < vals.Length; i++)
                    {
                        double.TryParse(vals[i], NumberStyles.Any, CultureInfo.InvariantCulture, out ua.MeanCurve[i]);
                    }
                }
            }
            // Confidence Intervals
            var ciAttr = xElement.Attribute(nameof(ua.ConfidenceIntervals));
            if (ciAttr != null)
            {
                var vals = ciAttr.Value.Split('|');
                if (vals.Length > 0)
                {
                    ua.ConfidenceIntervals = new double[vals.Length, vals[0].Split(',').Length];
                    for (int i = 0; i < vals.Length; i++)
                    {
                        var jVals = vals[i].Split(',');
                        for (int j = 0; j < jVals.Length; j++)
                        {
                            double.TryParse(jVals[j], NumberStyles.Any, CultureInfo.InvariantCulture, out ua.ConfidenceIntervals[i, j]);
                        }
                    }
                }
            }
            return ua;
        }

        /// <summary>
        /// Process and set the parent distribution and computed curve (mode, plug-in, point estimate).
        /// </summary>
        /// <param name="parentDistribution">The parent distribution.</param>
        /// <param name="probabilities">Array of non-exceedance probabilities.</param>
        public void ProcessModeCurve(UnivariateDistributionBase parentDistribution, double[] probabilities)
        {
            if (parentDistribution == null!)
                throw new ArgumentNullException(nameof(parentDistribution));
            if (probabilities == null || probabilities.Length == 0)
                throw new ArgumentException("Probabilities cannot be null or empty.", nameof(probabilities));

            ParentDistribution = parentDistribution;
            ModeCurve = ParentDistribution.InverseCDF(probabilities);
        }

        /// <summary>
        /// Process and set the confidence intervals from a list of sampled distributions.
        /// </summary>
        /// <param name="sampledDistributions">The list of sampled distributions to process.</param>
        /// <param name="probabilities">Array of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public void ProcessConfidenceIntervals(UnivariateDistributionBase[] sampledDistributions, double[] probabilities, double alpha = 0.1)
        {
            if (sampledDistributions == null || sampledDistributions.Length == 0)
                throw new ArgumentException("Sampled distributions cannot be null or empty.", nameof(sampledDistributions));
            if (probabilities == null || probabilities.Length == 0)
                throw new ArgumentException("Probabilities cannot be null or empty.", nameof(probabilities));
            if (alpha <= 0 || alpha >= 1)
                throw new ArgumentOutOfRangeException(nameof(alpha), "Alpha must be between 0 and 1.");

            int B = sampledDistributions.Length;
            int p = probabilities.Length;
            double lowerCI = alpha / 2d;
            double upperCI = 1d - alpha / 2d;
            ConfidenceIntervals = new double[p, 2];

            // Loop over probabilities and record percentiles
            for (int i = 0; i < p; i++)
            {
                var XValues = new double[B];

                // Compute quantiles in parallel
                Parallel.For(0, B, idx => {
                    XValues[idx] = sampledDistributions[idx]?.InverseCDF(probabilities[i]) ?? double.NaN;
                });

                // Filter valid values and sort
                int validCount = 0;
                for (int j = 0; j < B; j++)
                {
                    if (!double.IsNaN(XValues[j])) validCount++;
                }

                var validValues = new double[validCount];
                int writeIdx = 0;
                for (int j = 0; j < B; j++)
                {
                    if (!double.IsNaN(XValues[j]))
                        validValues[writeIdx++] = XValues[j];
                }

                Array.Sort(validValues);

                // Record percentiles for CIs
                ConfidenceIntervals[i, 0] = Statistics.Percentile(validValues, lowerCI, true);
                ConfidenceIntervals[i, 1] = Statistics.Percentile(validValues, upperCI, true);
            }
        }

        /// <summary>
        /// Computes the mean (predictive) curve by averaging CDFs across all sampled distributions.
        /// Uses log-spaced quantiles for efficient computation across wide ranges.
        /// </summary>
        /// <param name="sampledDistributions">Array of sampled distributions from posterior or bootstrap.</param>
        /// <param name="probabilities">Array of non-exceedance probabilities for interpolation.</param>
        /// <param name="minProbability">Minimum probability for range determination (default = 0.001).</param>
        /// <param name="maxProbability">Maximum probability for range determination (default = 1 - 1e-9).</param>
        public void ProcessMeanCurve(UnivariateDistributionBase[] sampledDistributions, double[] probabilities, double minProbability = 0.001, double maxProbability = 1 - 1e-9)
        {
            if (sampledDistributions == null || sampledDistributions.Length == 0)
                throw new ArgumentException("Sampled distributions cannot be null or empty.", nameof(sampledDistributions));
            if (probabilities == null || probabilities.Length == 0)
                throw new ArgumentException("Probabilities cannot be null or empty.", nameof(probabilities));

            int B = sampledDistributions.Length;

            // Compute min and max X values across all distributions
            double minX = double.MaxValue;
            double maxX = double.MinValue;
            object lockObject = new object();

            Parallel.For(0, B, j =>
            {
                if (sampledDistributions[j] != null!)
                {
                    var innerMin = sampledDistributions[j].InverseCDF(minProbability);
                    var innerMax = sampledDistributions[j].InverseCDF(maxProbability);

                    lock (lockObject)
                    {
                        if (innerMin < minX) minX = innerMin;
                        if (innerMax > maxX) maxX = innerMax;
                    }
                }
            });

            // Create log-spaced quantiles for efficient coverage
            double shift = minX <= 0 ? Math.Abs(minX) + 1d : 0;
            double min = minX + shift;
            double max = maxX + shift;
            int order = (int)Math.Floor(Math.Log10(max) - Math.Log10(min));
            int bins = Math.Max(200, Math.Min(1000, 100 * order));

            var quantiles = new double[bins];
            double delta = (Math.Log10(max) - Math.Log10(min)) / (bins - 1);

            for (int i = 0; i < bins; i++)
            {
                double logX = Math.Log10(min) + i * delta;
                quantiles[i] = Math.Pow(10, logX) - shift;
            }

            // Compute expected probability for each quantile
            var expected = new double[bins];
            for (int i = 0; i < bins; i++)
            {
                double total = 0d;
                Parallel.For(0, B, () => 0d, (j, loop, sum) =>
                {
                    if (sampledDistributions[j] != null!)
                    {
                        sum += sampledDistributions[j].CDF(quantiles[i]);
                    }
                    return sum;
                }, z => Tools.ParallelAdd(ref total, z));
                expected[i] = total / B;
            }

            // Build monotonic interpolation points
            var yVals = new List<double> { quantiles[0] };
            var xVals = new List<double> { expected[0] };
            double minY = quantiles[0];
            double maxY = quantiles[0];

            for (int i = 1; i < bins; i++)
            {
                if (expected[i] > xVals[xVals.Count - 1])
                {
                    minY = Math.Min(minY, quantiles[i]);
                    maxY = Math.Max(maxY, quantiles[i]);
                    yVals.Add(quantiles[i]);
                    xVals.Add(expected[i]);
                }
            }

            // Determine if log transform is appropriate
            bool useLogTransform = minY > 0 && (Math.Log10(maxY) - Math.Log10(minY)) > 1;

            // Interpolate mean curve at requested probabilities
            var linint = new Linear(xVals, yVals)
            {
                XTransform = Transform.NormalZ,
                YTransform = useLogTransform ? Transform.Logarithmic : Transform.None
            };
            MeanCurve = linint.Interpolate(probabilities);
        }

        /// <summary>
        /// Processes and stores the parameter sets from all sampled distributions.
        /// </summary>
        /// <param name="sampledDistributions">Array of sampled distributions to extract parameters from.</param>
        public void ProcessParameterSets(UnivariateDistributionBase[] sampledDistributions)
        {
            if (sampledDistributions == null || sampledDistributions.Length == 0)
                throw new ArgumentException("Sampled distributions cannot be null or empty.", nameof(sampledDistributions));

            int B = sampledDistributions.Length;
            ParameterSets = new ParameterSet[B];

            Parallel.For(0, B, idx =>
            {
                if (sampledDistributions[idx] != null!)
                {
                    ParameterSets[idx] = new ParameterSet(sampledDistributions[idx].GetParameters, double.NaN);
                }
            });
        }

    }
}
