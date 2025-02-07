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
using Numerics.Sampling.MCMC;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
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
    [Serializable]
    public class UncertaintyAnalysisResults
    {

        /// <summary>
        /// The parent probability distribution.
        /// </summary>
        public UnivariateDistributionBase ParentDistribution { get; set; }

        /// <summary>
        /// The array of parameter sets.
        /// </summary>
        public ParameterSet[] ParameterSets { get; set; }

        /// <summary>
        /// The confidence intervals. 
        /// </summary>
        public double[,] ConfidenceIntervals { get; set; }

        /// <summary>
        /// The mode (or computed) curve from the parent distribution. 
        /// </summary>
        public double[] ModeCurve { get; set; }

        /// <summary>
        /// The mean (or predictive) curve. 
        /// </summary>
        public double[] MeanCurve { get; set; }

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
            BinaryFormatter bf = new BinaryFormatter();
            using (var ms = new MemoryStream())
            {
                bf.Serialize(ms, results);
                return ms.ToArray();
            }
        }

        /// <summary>
        /// Returns the class from a byte array. 
        /// </summary>
        /// <param name="bytes">Byte array.</param>
        public static UncertaintyAnalysisResults FromByteArray(byte[] bytes)
        {
            try
            {
                using (var memStream = new MemoryStream())
                {
                    var binForm = new BinaryFormatter();
                    memStream.Write(bytes, 0, bytes.Length);
                    memStream.Seek(0, SeekOrigin.Begin);
                    var obj = binForm.Deserialize(memStream);
                    return (UncertaintyAnalysisResults)obj;
                }
            }
            catch (Exception)
            {
                // An error can occur because of differences in versions. 
                // If there is an error, just catch it and force the user to rerun the
                // uncertainty analysis. 
            }
            return null;
        }

        /// <summary>
        /// Returns the class as XElement. The parameter sets are not included. 
        /// </summary>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(UncertaintyAnalysisResults));
            if (ParentDistribution != null) result.Add(ParentDistribution.ToXElement());
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
            if (xElement.Element("Distribution") != null)
                ua.ParentDistribution = UnivariateDistributionFactory.CreateDistribution(xElement.Element("Distribution"));

            // AIC
            if (xElement.Attribute(nameof(AIC)) != null)
            {
                double.TryParse(xElement.Attribute(nameof(AIC)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var aic);
                ua.AIC = aic;
            }
            // BIC
            if (xElement.Attribute(nameof(BIC)) != null)
            {
                double.TryParse(xElement.Attribute(nameof(BIC)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var bic);
                ua.BIC = bic;
            }
            // DIC
            if (xElement.Attribute(nameof(DIC)) != null)
            {
                double.TryParse(xElement.Attribute(nameof(DIC)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var dic);
                ua.DIC = dic;
            }
            // RMSE
            if (xElement.Attribute(nameof(RMSE)) != null)
            {
                double.TryParse(xElement.Attribute(nameof(RMSE)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var rmse);
                ua.RMSE = rmse;
            }
            // ERL
            if (xElement.Attribute(nameof(ERL)) != null)
            {
                double.TryParse(xElement.Attribute(nameof(ERL)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var erl);
                ua.ERL = erl;
            }

            // Mode Curve
            if (xElement.Attribute(nameof(ua.ModeCurve)) != null)
            {
                var vals = xElement.Attribute(nameof(ua.ModeCurve)).Value.Split('|');
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
            if (xElement.Attribute(nameof(ua.MeanCurve)) != null)
            {
                var vals = xElement.Attribute(nameof(ua.MeanCurve)).Value.Split('|');
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
            if (xElement.Attribute(nameof(ua.ConfidenceIntervals)) != null)
            {
                var vals = xElement.Attribute(nameof(ua.ConfidenceIntervals)).Value.Split('|');
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
        /// Returns uncertainty analysis results for the distribution estimated from MCMC.
        /// </summary>
        /// <param name="results">The MCMC results.</param>
        /// <param name="distribution">The parent distribution.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public static UncertaintyAnalysisResults FromMCMCResults(MCMCResults results, UnivariateDistributionBase distribution, IList<double> probabilities, double alpha = 0.1)
        {
            if (results == null|| results.Output == null ||results.Output.Count == 0) return null;
            if (results.Output[0].Values.Length != distribution.NumberOfParameters) return null;

            int realz = results.Output.Count;
            var dists = new UnivariateDistributionBase[results.Output.Count];

            Parallel.For(0, realz, idx => 
            {
                var d = distribution.Clone();
                d.SetParameters(results.Output[idx].Values);
                dists[idx] = d;
            });

            // Set up dummy bootstrap analysis
            var boot = new BootstrapAnalysis(distribution, ParameterEstimationMethod.MaximumLikelihood, 100, realz);
            var uar = boot.Estimate(probabilities, alpha, dists);
            uar.ParameterSets = null;
            return uar;
        }

    }
}
