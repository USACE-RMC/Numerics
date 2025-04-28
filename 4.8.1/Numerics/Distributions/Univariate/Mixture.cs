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
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Xml.Linq;

namespace Numerics.Distributions
{
    /// <summary>
    /// A Mixture distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Mixture_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Mixture : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {

        /// <summary>
        /// Construct new mixture distribution.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public Mixture(double[] weights, UnivariateDistributionBase[] distributions)
        {
            SetParameters(weights, distributions);
        }

        /// <summary>
        /// Construct new mixture distribution.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public Mixture(double[] weights, IUnivariateDistribution[] distributions)
        {
            SetParameters(weights, distributions);
        }

        private double[] _weights;
        private UnivariateDistributionBase[] _distributions;
        private EmpiricalDistribution _empiricalCDF;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;
        private bool _empiricalCDFCreated = false;

        /// <summary>
        /// Returns the array of distribution weights.
        /// </summary>
        public double[] Weights => _weights;

        /// <summary>
        /// Returns the array of univariate probability distributions.
        /// </summary>
        public UnivariateDistributionBase[] Distributions => _distributions;

        /// <summary>
        /// Determines whether to use assign a weight to all data points less than or equal to zero.
        /// </summary>
        public bool IsZeroInflated { get; set; } = false;

        /// <summary>
        /// The zero-value weight used if the mixture is zero-inflated.
        /// </summary>
        public double ZeroWeight { get; set; } = 0.0;

        /// <summary>
        /// Determines the interpolation transform for the X-values.
        /// </summary>
        public Transform XTransform { get; set; } = Transform.None;

        /// <summary>
        /// Determines the interpolation transform for the Probability-values.
        /// </summary>
        public Transform ProbabilityTransform { get; set; } = Transform.NormalZ;

        /// <summary>
        /// The maximum iterations in the Expectation Maximization algorithm. Default = 1,000. 
        /// </summary>
        public int MaxIterations { get; set; } = 1000;

        /// <summary>
        /// The relative tolerance for convergence. Default = 1E-8.
        /// </summary>
        public double Tolerance { get; set; } = 1E-8;

        /// <summary>
        /// The total number of iterations required to find the MLE.
        /// </summary>
        public int Iterations { get; private set; }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get
            {
                int sum = 0;
                sum += Distributions.Count();
                for (int i = 0; i < Distributions.Count(); i++)
                    sum += Distributions[i].NumberOfParameters;
                return sum;
            }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type => UnivariateDistributionType.Mixture;

        /// <inheritdoc/>
        public override string DisplayName => "Mixture";

        /// <inheritdoc/>
        public override string ShortDisplayName => "MIX";

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                string Wstring = "{";
                string Dstring = "{";
                for (int i = 1; i < Weights.Count() - 1; i++)
                {
                    Wstring += Weights[i].ToString();
                    Dstring += Distributions[i].DisplayName;
                    if (i < Weights.Count() - 2)
                    {
                        Wstring += ",";
                        Dstring += ",";
                    }
                }
                Wstring += "}";
                Dstring += "}";
                parmString[0, 0] = "Weights";
                parmString[1, 0] = "Distributions";
                parmString[0, 1] = Wstring;
                parmString[1, 1] = Dstring;
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNames
        {
            get
            {
                var result = new List<string>();
                for (int i = 1; i <= Distributions.Count(); i++)
                {
                    result.Add("Weight " + i.ToString());
                }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    for (int j = 0; j < Distributions[i].ParameterNames.Length; j++)
                    {
                        result.Add("D" + (i + 1).ToString() + " " + Distributions[i].ParameterNames[j]);
                    }
                }
                return result.ToArray();
            }
        }


        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get
            {
                var result = new List<string>();
                for (int i = 1; i <= Distributions.Count(); i++)
                {
                    result.Add("W" + i.ToString());
                }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    for (int j = 0; j < Distributions[i].ParameterNamesShortForm.Length; j++)
                    {
                        result.Add("D" + (i + 1).ToString() + " " + Distributions[i].ParameterNamesShortForm[j]);
                    }
                }
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get
            {
                var result = new List<double>();
                result.AddRange(Weights);
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.AddRange(Distributions[i].GetParameters);
                }                  
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Weights), nameof(Distributions)]; }
        }

        /// <summary>
        /// Compute central moments of the distribution.
        /// </summary>
        private void ComputeMoments()
        {
            var mom = CentralMoments(1000);
            u1 = mom[0];
            u2 = mom[1];
            u3 = mom[2];
            u4 = mom[3];
            _momentsComputed = true;
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (!_momentsComputed) 
                    ComputeMoments();
                return u1;
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                var brent = new BrentSearch(PDF, InverseCDF(0.001), InverseCDF(0.999));
                brent.Maximize();
                return brent.BestParameterSet.Values[0];
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (!_momentsComputed) 
                    ComputeMoments();
                return u2;
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (!_momentsComputed) 
                    ComputeMoments();
                return u3;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (!_momentsComputed) 
                    ComputeMoments();
                return u4;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return Distributions.Min(p => p.Minimum); }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Distributions.Max(p => p.Maximum); }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        { 
            get 
            {
                var result = new List<double>();
                if (IsZeroInflated) { result.Add(0.0); }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.Add(0.0);
                }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.AddRange(Distributions[i].MinimumOfParameters);
                }
                return result.ToArray(); 
            }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get
            {
                var result = new List<double>();
                if (IsZeroInflated) { result.Add(1.0); }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.Add(1.0);
                }
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.AddRange(Distributions[i].MaximumOfParameters);
                }
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                SetParameters(MLE(sample));
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <inheritdoc/>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var newDistribution = (Mixture)Clone();
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public void SetParameters(double[] weights, UnivariateDistributionBase[] distributions)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            if (weights.Length != distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));

            _weights = weights;
            _distributions = distributions;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public void SetParameters(double[] weights, IUnivariateDistribution[] distributions)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            if (weights.Length != distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));

            _weights = weights;
            _distributions = new UnivariateDistributionBase[distributions.Length];
            for (int i = 0; i < distributions.Length; i++)
            {
                _distributions[i] = (UnivariateDistributionBase)distributions[i];
            }
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="parameters">The mixture distribution parameters.</param>
        public void SetParameters(double[] weights, double[] parameters)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (weights.Length != Distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));
            if (parameters.Length != Distributions.Sum(x => x.NumberOfParameters))
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }

            // Set weights
            _weights = weights.ToArray();
            // Set distribution parameters
            int t = 0;
            for (int i = 0; i < Distributions.Count(); i++)
            {
                var parms = new List<double>();
                for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                {
                    parms.Add(parameters[j]);
                }
                Distributions[i].SetParameters(parms);
                t += Distributions[i].NumberOfParameters;
            }
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            if (parameters.Count != NumberOfParameters)
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }

            // Set the weights
            int t = 0;
            for (int i = 0; i < Distributions.Count(); i++)
            {
                Weights[i] = parameters[i];
                t++;
            }

            // Set the distribution parameters
            for (int i = 0; i < Distributions.Count(); i++)
            {
                var parms = new List<double>();
                for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                {
                    parms.Add(parameters[j]);
                }
                Distributions[i].SetParameters(parms);
                t += Distributions[i].NumberOfParameters;
            }
            // Validate parameters
            _parametersValid = ValidateParameters(parameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters from a referenced array. Weights are normalized to sum to 1.
        /// </summary>
        /// <param name="parameters">The array of parameters.</param>
        public void SetParameters(ref double[] parameters)
        {
            if (parameters == null) return;
            if (Weights == null || Weights.Length == 0) return;
            if (Distributions == null || Distributions.Count() == 0) return;
            if (Distributions.Count() == 1 && parameters.Length == Distributions[0].NumberOfParameters)
            {
                if (IsZeroInflated)
                {
                    Weights[0] = 1 - ZeroWeight;
                }
                else
                {
                    Weights[0] = 1;
                }
                Distributions[0].SetParameters(parameters);
            }
            else
            {
                   
                // Get the weights
                int K = Distributions.Count();
                int t = 0; // keep track of parameter index

                // Get weights
                double sum = 0.0;
                for (int i = 0; i < K; i++)
                {
                    Weights[i] = parameters[i];
                    sum += Weights[i];
                    t++;
                }

                // Check if weights need to be normalized
                if (sum <= 0.0)
                {
                    // If weights sum to 0, reset to be uniformly distributed
                    double w = IsZeroInflated ? (1d - ZeroWeight) / K : 1d / K;
                    for (int i = 0; i < K; i++)
                    {
                        Weights[i] = w;
                        parameters[i] = Weights[i];
                    }
                }                  
                else
                {
                    // Normalize weights to sum to 1.
                    var c = IsZeroInflated ? (1d - ZeroWeight) / sum : 1d / sum;
                    // Normalize weights
                    for (int i = 0; i < K; i++)
                    {
                        Weights[i] *= c;
                        parameters[i] = Weights[i];
                    }
                }
                

                // Set distribution parameters
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    var parms = new List<double>();
                    for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                    {
                        parms.Add(parameters[j]);
                    }
                    Distributions[i].SetParameters(parms);
                    t += Distributions[i].NumberOfParameters;
                }
            }

            // Validate parameters
            _parametersValid = ValidateParameters(parameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            // Check if weights are between 0 and 1.
            if (IsZeroInflated && (ZeroWeight < 0.0 || ZeroWeight > 1.0))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(ZeroWeight), "The zero value weight must be between 0 and 1.");
                return new ArgumentOutOfRangeException(nameof(ZeroWeight), "The zero value weight must be between 0 and 1.");
            }
            for (int i = 0; i < Distributions.Count(); i++)
            {
                if (Weights[i] < 0.0 || Weights[i] > 1.0)
                {
                    if (throwException)
                        throw new ArgumentOutOfRangeException(nameof(Weights), "The weights must be between 0 and 1.");
                    return new ArgumentOutOfRangeException(nameof(Weights), "The weights must be between 0 and 1.");
                }
            }
            // Check if weights sum to 1.
            double sum = IsZeroInflated ? ZeroWeight : 0.0;
            for (int i = 0; i < Distributions.Count(); i++)
                sum += Weights[i];
            if (sum.AlmostEquals(1d, 1E-8) == false)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Weights), "The weights must sum to 1.0.");
                return new ArgumentOutOfRangeException(nameof(Weights), "The weights must sum to 1.0.");
            }
            // Check if distributions are valid
            for (int i = 0; i < Distributions.Count(); i++)
            {
                if (Distributions[i].ParametersValid == false)
                {
                    if (throwException)
                        throw new ArgumentOutOfRangeException(nameof(Distributions), "Distribution " + (i + 1).ToString() + " has invalid parameters.");
                    return new ArgumentOutOfRangeException(nameof(Distributions), "Distribution " + (i + 1).ToString() + " has invalid parameters.");
                }
            }
            return null;
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];

            // Weights are first
            int t = 0;
            for (int i = 0; i < Distributions.Count(); i++)
            {
                initialVals[i] = IsZeroInflated ? (1d - ZeroWeight) / Distributions.Count() : 1d / Distributions.Count();
                lowerVals[i] = 0.0;
                upperVals[i] = 1.0;
                t += 1;
            }

            for (int i = 0; i < Distributions.Count(); i++)
            {
                var tuple = ((IMaximumLikelihoodEstimation)Distributions[i]).GetParameterConstraints(sample);
                var initials = tuple.Item1;
                var lowers = tuple.Item2;
                var uppers = tuple.Item3;

                for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                {
                    initialVals[j] = initials[j - t];
                    lowerVals[j] = lowers[j - t];
                    upperVals[j] = uppers[j - t];
                }
                t += Distributions[i].NumberOfParameters;
            }
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        public double[] MLE(IList<double> sample)
        {

            int N = sample.Count;
            int Np = Distributions.Sum(x => x.NumberOfParameters);
            int K = Distributions.Count();

            // Set constraints
            var tuple = GetParameterConstraints(sample);
            var Initials = tuple.Item1.Subset(K);
            var Lowers = tuple.Item2.Subset(K);
            var Uppers = tuple.Item3.Subset(K);

            // Set up EM parameters
            var mleWeights = tuple.Item1.Subset(0 , K - 1);
            var mleParameters = Initials;
            var likelihood = new double[N, K];
            double oldLogLH = double.MinValue, newLogLH = double.MinValue;

            // The Expectation step. 
            double EStep(double[] x)
            {
                var dist = (Mixture)Clone();
                dist.SetParameters(mleWeights, x);
                // Outer loop for computing the likelihoods
                for (int k = 0; k < K; k++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        if (IsZeroInflated && sample[i] <= 0.0)
                        {
                            likelihood[i, k] = Math.Log(ZeroWeight);
                        }
                        else
                        {
                            likelihood[i, k] = Math.Log(mleWeights[k]) + dist.Distributions[k].LogPDF(sample[i]);
                        }
                    }
                }
                // At this point we have unnormalized log likelihoods.
                // We need to normalize using log-sum-exp and compute the true log-likelihoods.
                double logLH = 0;
                for (int i = 0; i < N; i++)
                {
                    // Get max likelihood
                    double max = double.MinValue;
                    for (int k = 0; k < K; k++)
                    {
                        if (likelihood[i, k] > max)
                        {
                            max = likelihood[i, k];
                        }
                    }
                    // log-sum-exp trick begins here
                    double sum = 0;
                    for (int k = 0; k < K; k++)
                        sum += Math.Exp(likelihood[i, k] - max);
                    double tmp = max + Math.Log(sum);
                    for (int k = 0; k < K; k++)
                        likelihood[i, k] = Math.Exp(likelihood[i, k] - tmp);
                    logLH += tmp;
                }
                return logLH;
            }

            // The Maximization step
            double[] MStep(double[] x)
            {
                // Get updated weights
                for (int k = 0; k < K; k++)
                {
                    double wgt = 0d;
                    for (int i = 0; i < N; i++)
                    {
                        if (!IsZeroInflated || sample[i] > 0.0 )
                        {
                            wgt += likelihood[i, k];
                        }
                    }  
                    mleWeights[k] = wgt / N;
                }
                // MLE
                var solver = new NelderMead(logLH, Np, x, Lowers, Uppers);
                solver.Maximize();
                return solver.BestParameterSet.Values;
            }

            // The log-likelihood to maximize in the M-Step
            // Weights are held fixed, only the distribution parameters are solved.
            double logLH(double[] x)
            {
                var dist = (Mixture)Clone();
                dist.SetParameters(mleWeights, x);
                double lh = dist.LogLikelihood(sample);
                if (double.IsNaN(lh) || double.IsInfinity(lh)) return double.MinValue;
                return lh;
            }

            // Estimate using the EM Algorithm
            for (Iterations = 1; Iterations <= MaxIterations; Iterations++)
            {
                // Perform the expectation step
                newLogLH = EStep(mleParameters);

                // Check convergence
                if (Math.Abs((oldLogLH - newLogLH) / oldLogLH) < Tolerance)
                    break;

                // Perform the maximization step
                mleParameters = MStep(mleParameters);

                // Update log-likelihood state
                oldLogLH = newLogLH;

            }

            // Return the full list of distribution parameters
            var result = new List<double>();
            result.AddRange(mleWeights);
            result.AddRange(mleParameters);
            return result.ToArray();
        }


        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            double f = 0.0;
            if (IsZeroInflated && x <= 0.0)
            {
                f = ZeroWeight;
            }
            else
            {
                for (int i = 0; i < Distributions.Count(); i++)
                    f += Weights[i] * Distributions[i].PDF(x);
            }
            return f < 0d ? 0d : f;
        }

        /// <inheritdoc/>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            var lnf = new List<double>();
            if (IsZeroInflated && x <= 0.0)
            {
                lnf.Add(Math.Log(ZeroWeight));
            }
            else
            {
                for (int i = 0; i < Distributions.Count(); i++)
                    lnf.Add(Math.Log(Weights[i]) + Distributions[i].LogPDF(x));
            }
            var f = Tools.LogSumExp(lnf);
            return f;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            double F = 0.0;
            if (IsZeroInflated)
            {
                F = ZeroWeight;
                if (x > 0.0)
                {
                    for (int i = 0; i < Distributions.Count(); i++)
                        F += Weights[i] * Distributions[i].CDF(x);
                }
            }
            else
            {
                for (int i = 0; i < Distributions.Count(); i++)
                    F += Weights[i] * Distributions[i].CDF(x);
            }
            return F < 0d ? 0d : F > 1d ? 1d : F;
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            var lnF = new List<double>();
            if (IsZeroInflated)
            {
                lnF.Add(Math.Log(ZeroWeight));
                if (x > 0.0)
                {
                    for (int i = 0; i < Distributions.Count(); i++)
                        lnF.Add(Math.Log(Weights[i]) + Distributions[i].LogCDF(x));
                }             
            }
            else
            {
                for (int i = 0; i < Distributions.Count(); i++)
                    lnF.Add(Math.Log(Weights[i]) + Distributions[i].LogCDF(x));
            }
            var F = Tools.LogSumExp(lnF);
            return F;
        }

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            var lnF = new List<double>();
            if (IsZeroInflated)
            {
                if (x > 0.0)
                {
                    for (int i = 0; i < Distributions.Count(); i++)
                        lnF.Add(Math.Log(Weights[i]) + Distributions[i].LogCCDF(x));
                }
            }
            else
            {
                for (int i = 0; i < Distributions.Count(); i++)
                    lnF.Add(Math.Log(Weights[i]) + Distributions[i].LogCCDF(x));
            }
            var F = Tools.LogSumExp(lnF);
            return F;
        }


        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            if (IsZeroInflated && probability <= ZeroWeight) return 0;

            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            // If there is only one distribution, return its inverse CDF
            if (Distributions.Count() == 1)
            {
                return Distributions[0].InverseCDF(probability);
            }

            double x = 0;
            if (_empiricalCDFCreated == true)
            {
                x = _empiricalCDF.InverseCDF(probability);
            }
            else
            {
                // Use a root finder to solve the inverse CDF
                var xVals = Distributions.Select(d => d.InverseCDF(probability));
                double minX = xVals.Min();
                double maxX = xVals.Max();
                try
                {
                    if (IsZeroInflated)
                    {
                        Brent.Bracket((y) => { return probability - CDF(y); }, ref minX, ref maxX, out var f1, out var f2);
                    }
                    x = Brent.Solve((y) => { return probability - CDF(y); }, minX, maxX, 1E-6, 100, true);
                }
                catch (Exception)
                {
                    // If the root finder fails, create an empirical CDF
                    if (_empiricalCDFCreated == false)
                        CreateEmpiricalCDF();
                    x = _empiricalCDF.InverseCDF(probability);
                }
            }
            double min = Minimum;
            double max = Maximum;
            return x < min ? min : x > max ? max : x;
        }
    
        /// <summary>
        /// Create empirical distribution for the CDF.
        /// </summary>
        public void CreateEmpiricalCDF()
        {
            // Get min & max
            double minP = 1E-16;
            double maxP = 1 - 1E-16;
            double minX = Distributions.Min(d => d.InverseCDF(minP));
            double maxX = Distributions.Max(d => d.InverseCDF(maxP));
            // Get number of bins
            double shift = 0;
            if (minX <= 0) shift = Math.Abs(minX) + 1d;
            double min = minX + shift;
            double max = maxX + shift;
            int order = (int)Math.Floor(Math.Log10(max) - Math.Log10(min));
            int binN = Math.Max(200, 100 * order) - 1;
            // Create bins
            var bins = Stratify.XValues(new StratificationOptions(minX, maxX, binN, false), true);
            var xValues = new List<double>();
            var pValues = new List<double>();
            var x = bins.First().LowerBound;
            var p = CDF(bins.First().LowerBound);
            xValues.Add(x);
            pValues.Add(p);
            for (int i = 1; i < bins.Count; i++)
            {
                x = bins[i].LowerBound;
                p = CDF(x);
                if (x > xValues.Last() && p > pValues.Last())
                {
                    xValues.Add(x);
                    pValues.Add(p);
                }
            }
            x = maxX;
            p = CDF(x);
            if (x > xValues.Last() && p > pValues.Last())
            {
                xValues.Add(x);
                pValues.Add(p);
            }
            _empiricalCDF = new EmpiricalDistribution(xValues, pValues) { XTransform = XTransform, ProbabilityTransform = ProbabilityTransform };
            _empiricalCDFCreated = true;
            _momentsComputed = false;
        }

        /// <inheritdoc/>
        public override double[] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            // Create PRNG for generating random numbers
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var weights = new List<double>();
            var distributions = new List<UnivariateDistributionBase>();
            if (IsZeroInflated)
            {
                weights.Add(ZeroWeight);
                distributions.Add(new Deterministic(0.0));
            }
            weights.AddRange(Weights);
            distributions.AddRange(Distributions);

            var sample = new double[sampleSize];
            // Generate values
            for (int i = 0; i < sampleSize; i++)
            {
                var u = rnd.NextDouble();
                var cdfW = new double[distributions.Count()];
                for (int j = 0; j < distributions.Count(); j++)
                {
                    cdfW[j] = j == 0 ? weights[j] : cdfW[j - 1] + weights[j];
                    if (u <= cdfW[j])
                    {
                        sample[i] = distributions[j].InverseCDF(rnd.NextDouble());
                        break;
                    }
                }
            }
            // Return array of random values
            return sample;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            var dists = new UnivariateDistributionBase[Distributions.Count()];
            for (int i = 0; i < Distributions.Count(); i++)
                dists[i] = Distributions[i].Clone();

            return new Mixture(Weights.ToArray(), dists)
            {
                IsZeroInflated = IsZeroInflated,
                ZeroWeight = ZeroWeight,
                XTransform = XTransform,
                ProbabilityTransform = ProbabilityTransform
            };
        }

        /// <inheritdoc/>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue(nameof(IsZeroInflated), IsZeroInflated.ToString());
            result.SetAttributeValue(nameof(ZeroWeight), ZeroWeight.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(XTransform), XTransform.ToString());
            result.SetAttributeValue(nameof(ProbabilityTransform), ProbabilityTransform.ToString());
            result.SetAttributeValue(nameof(Distributions), String.Join("|", Distributions.Select(x => x.Type)));
            // Weights
            var weights = Weights;
            var weightStrings = new string[Weights.Length];
            for (int i = 0; i < Weights.Length; i++)
            {
                weightStrings[i] = weights[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue(nameof(Weights), String.Join("|", weightStrings));
            // Parameters
            var parms = GetParameters;
            var parmStrings = new string[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("Parameters", String.Join("|", parmStrings));
            return result;
        }

        /// <summary>
        /// Create a mixture distribution from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new mixture distribution.</returns>
        public static Mixture FromXElement(XElement xElement)
        {
            UnivariateDistributionType type = UnivariateDistributionType.Deterministic;
            if (xElement.Attribute(nameof(UnivariateDistributionBase.Type)) != null)
            {
                Enum.TryParse(xElement.Attribute(nameof(UnivariateDistributionBase.Type)).Value, out type);

            }
            if (type == UnivariateDistributionType.Mixture)
            {
                var weights = new List<double>();
                var distributions = new List<UnivariateDistributionBase>();
                if (xElement.Attribute(nameof(Weights)) != null)
                {
                    var w = xElement.Attribute(nameof(Weights)).Value.Split('|');
                    for (int i = 0; i < w.Length; i++)
                    {
                        double.TryParse(w[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var weight);
                        weights.Add(weight);
                    }
                }
                if (xElement.Attribute(nameof(Distributions)) != null)
                {
                    var types = xElement.Attribute(nameof(Distributions)).Value.Split('|');
                    for (int i = 0; i < types.Length; i++)
                    {
                        Enum.TryParse(types[i], out UnivariateDistributionType distType);
                        distributions.Add(UnivariateDistributionFactory.CreateDistribution(distType));
                    }
                }
                var mixture = new Mixture(weights.ToArray(), distributions.ToArray());

                if (xElement.Attribute(nameof(IsZeroInflated)) != null)
                {
                    bool.TryParse(xElement.Attribute(nameof(IsZeroInflated)).Value, out var isZeroInflated);
                    mixture.IsZeroInflated = isZeroInflated;
                }
                if (xElement.Attribute(nameof(ZeroWeight)) != null)
                {
                    double.TryParse(xElement.Attribute(nameof(ZeroWeight)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var zeroWeight);
                    mixture.ZeroWeight = zeroWeight;
                }
                if (xElement.Attribute(nameof(XTransform)) != null)
                {
                    Enum.TryParse(xElement.Attribute(nameof(XTransform)).Value, out Transform xTransform);
                    mixture.XTransform = xTransform;
                }
                if (xElement.Attribute(nameof(ProbabilityTransform)) != null)
                {
                    Enum.TryParse(xElement.Attribute(nameof(ProbabilityTransform)).Value, out Transform probabilityTransform);
                    mixture.ProbabilityTransform = probabilityTransform;
                }
                if (xElement.Attribute("Parameters") != null)
                {
                    var vals = xElement.Attribute("Parameters").Value.Split('|');
                    var parameters = new List<double>();
                    for (int i = 0; i < vals.Length; i++)
                    {
                        double.TryParse(vals[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var parm);
                        parameters.Add(parm);
                    }
                    mixture.SetParameters(parameters);
                }

                return mixture;
            }
            else
            {
                return null;
            }
        }

    }
}
