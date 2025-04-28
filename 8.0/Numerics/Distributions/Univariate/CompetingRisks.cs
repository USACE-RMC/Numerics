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
using Numerics.Mathematics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Xml.Linq;

namespace Numerics.Distributions
{
    /// <summary>
    /// A competing risks distribution.
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
    /// <see href = "https://reliability.readthedocs.io/en/latest/Competing%20risk%20models.html" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class CompetingRisks : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {
        /// <summary>
        /// Construct new competing risks distribution.
        /// </summary>
        /// <param name="distributions">The competing distributions.</param>
        public CompetingRisks(UnivariateDistributionBase[] distributions)
        {
            SetParameters(distributions);
        }

        /// <summary>
        /// Construct new competing risks distribution.
        /// </summary>
        /// <param name="distributions">The competing distributions.</param>
        public CompetingRisks(IUnivariateDistribution[] distributions)
        {
            SetParameters(distributions);
        }

        private UnivariateDistributionBase[] _distributions;
        private EmpiricalDistribution _empiricalCDF;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;
        private bool _empiricalCDFCreated = false;
        private double[,] _correlationMatrix;
        private bool _mvnCreated = false;
        private MultivariateNormal _mvn;

        /// <summary>
        /// Returns the array of univariate probability distributions.
        /// </summary>
        public ReadOnlyCollection<UnivariateDistributionBase> Distributions => new(_distributions);

        /// <summary>
        /// Determines the interpolation transform for the X-values.
        /// </summary>
        public Transform XTransform { get; set; } = Transform.None;

        /// <summary>
        /// Determines the interpolation transform for the Probability-values.
        /// </summary>
        public Transform ProbabilityTransform { get; set; } = Transform.NormalZ;

        /// <summary>
        /// If true, the competing risks model computes the minimum of the random variables. If false, it computes the maximum of random variables. 
        /// </summary>
        public bool MinimumOfRandomVariables { get; set; } = true;

        /// <summary>
        /// The dependency between random variables. 
        /// </summary>
        public Probability.DependencyType Dependency { get; set; } = Probability.DependencyType.Independent;

        /// <summary>
        /// The correlation matrix used for modeling dependency between the marginal distributions.
        /// This is only used when the Dependency Type = CorrelationMatrix.
        /// </summary>
        public double[,] CorrelationMatrix
        {
            get { return _correlationMatrix; }
            set
            {
                _correlationMatrix = value;
                _mvnCreated = false;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get
            {
                int sum = 0;
                for (int i = 0; i < Distributions.Count; i++)
                    sum += Distributions[i].NumberOfParameters;
                return sum;
            }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type => UnivariateDistributionType.CompetingRisks;

        /// <inheritdoc/>
        public override string DisplayName => "Competing Risks";

        /// <inheritdoc/>
        public override string ShortDisplayName => "CR";

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[1, 2];
                string Dstring = "{";
                for (int i = 1; i < Distributions.Count - 1; i++)
                {
                    Dstring += Distributions[i].DisplayName;
                    if (i < Distributions.Count - 2)
                    {
                        Dstring += ",";
                    }
                }
                Dstring += "}";
                parmString[0, 0] = "Distributions";
                parmString[0, 1] = Dstring;
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNames
        {
            get
            {
                var result = new List<string>();
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
                for (int i = 0; i < Distributions.Count; i++)
                    result.AddRange(Distributions[i].GetParameters);
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Distributions)]; }
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
                if (!_momentsComputed) ComputeMoments();
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
                if (!_momentsComputed) ComputeMoments();
                return u2;
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
                return u3;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (!_momentsComputed) ComputeMoments();
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
            var newDistribution = (CompetingRisks)Clone();
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="distributions">The competing distributions.</param>
        public void SetParameters(UnivariateDistributionBase[] distributions)
        {
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            _distributions = distributions;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
            _mvnCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="distributions">The competing distributions.</param>
        public void SetParameters(IUnivariateDistribution[] distributions)
        {
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            _distributions = new UnivariateDistributionBase[distributions.Length];
            for (int i = 0; i < distributions.Length; i++)
            {
                _distributions[i] = (UnivariateDistributionBase)distributions[i];
            }
            _momentsComputed = false;
            _empiricalCDFCreated = false;
            _mvnCreated = false;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            if (Distributions == null || Distributions.Count == 0) return;

            int t = 0;
            for (int i = 0; i < Distributions.Count; i++)
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
            _mvnCreated = false;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (Distributions.Count == 0) return new ArgumentOutOfRangeException(nameof(Distributions), "There must be at least 1 distribution");
            for (int i = 0; i < Distributions.Count; i++)
            {
                if (Distributions[i].ParametersValid == false)
                {
                    if (throwException)
                        throw new ArgumentOutOfRangeException(nameof(Distributions), "One of the distributions have invalid parameters.");
                    return new ArgumentOutOfRangeException(nameof(Distributions), "One of the distributions have invalid parameters.");
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

            int t = 0;
            for (int i = 0; i < Distributions.Count; i++)
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
            // Set constraints
            var tuple = GetParameterConstraints(sample);
            var Initials = tuple.Item1;
            var Lowers = tuple.Item2;
            var Uppers = tuple.Item3;

            // Solve using Nelder-Mead (Downhill Simplex)
            double logLH(double[] x)
            {
                var dist = (CompetingRisks)Clone();
                dist.SetParameters(x);
                return dist.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            if (Distributions.Count == 1)
            {
                return Distributions[0].PDF(x);
            }

            double f = double.NaN;

            // Only compute the exact PDF for independent random variables
            if (Dependency == Probability.DependencyType.Independent)
            {
                double hf = 0d;
                double sf = 1d;
                for (int i = 0; i < Distributions.Count; i++)
                {

                    if (MinimumOfRandomVariables == true)
                    {
                        hf += Distributions[i].HF(x);
                        sf *= Distributions[i].CCDF(x);
                    }
                    else
                    {
                        hf += Distributions[i].PDF(x) / Distributions[i].CDF(x);
                        sf *= Distributions[i].CDF(x);
                    }
                }
                f = hf * sf;
            }
            else
            {
                // Compute the PDF using numerical differentiation
                f = NumericalDerivative.Derivative(CDF, x);
            }

            return f < 0d ? 0d : f;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            if (Distributions.Count == 1)
            {
                return Distributions[0].CDF(x);
            }

            double p = double.NaN;
            var ind = new int[Distributions.Count];
            var cdf = new double[Distributions.Count];
            for (int i = 0; i < Distributions.Count; i++)
            {
                ind[i] = 1;
                cdf[i] = Distributions[i].CDF(x);
            }

            if (MinimumOfRandomVariables == true)
            {

                if (Dependency == Probability.DependencyType.PerfectlyNegative || Dependency == Probability.DependencyType.CorrelationMatrix)
                {
                    if (_mvnCreated == false)
                        CreateMultivariateNormal();
                    p = Probability.UnionPCM(cdf, _mvn.Covariance);
                }
                else
                {
                    p = Probability.Union(cdf, Dependency);
                }
            }
            else
            {
                if (Dependency == Probability.DependencyType.PerfectlyNegative || Dependency == Probability.DependencyType.CorrelationMatrix)
                {
                    if (_mvnCreated == false)
                        CreateMultivariateNormal();
                    p = Probability.JointProbability(cdf, ind, _mvn.Covariance);
                }
                else
                {
                    p = Probability.JointProbability(cdf, Dependency);
                }

            }
            return p < 0d ? 0d : p > 1d ? 1d : p;
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([0], true);

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
                // use a root finder to solve the inverse CDF
                var xVals = Distributions.Select(d => d.InverseCDF(probability));
                double minX = xVals.Min();
                double maxX = xVals.Max();
                try
                {
                    Brent.Bracket((y) => { return probability - CDF(y); }, ref minX, ref maxX, out var f1, out var f2);
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
        /// Returns a list of cumulative incidence functions. 
        /// </summary>
        /// <param name="bins">Optional. The stratification bins to integrate over. Default is 200 bins.</param>
        public List<EmpiricalDistribution> CumulativeIncidenceFunctions(List<StratificationBin> bins = null)
        {
            // Get stratification bins
            if (bins == null)
            {
                double minP = 1E-16;
                double maxP = 1 - 1E-16;
                double minX = Distributions.Min(d => d.InverseCDF(minP));
                double maxX = Distributions.Max(d => d.InverseCDF(maxP));
                bins = Stratify.XValues(new StratificationOptions(minX, maxX, 200, false), XTransform == Transform.Logarithmic ? true : false);
            }

            var D = Distributions.Count();
            var x = new List<double[]>();
            var p = new List<double[]>();
            var dF = new List<double[]>();

            if (Dependency == Probability.DependencyType.PerfectlyNegative || Dependency == Probability.DependencyType.CorrelationMatrix)
            {
                /* 
                * For perfect negative dependency or a custom correlation matrix,
                * use the Genz Method. This method is slow but accurate.
                */

                if (_mvnCreated == false)
                    CreateMultivariateNormal();

                var lower = new double[D];
                var upper = new double[D];
                for (int i = 0; i < D; i++)
                {
                    x.Add(new double[bins.Count + 1]);
                    p.Add(new double[bins.Count + 1]);
                    dF.Add(new double[bins.Count + 1]);

                    // Record the first bin
                    x[i][0] = bins[0].LowerBound;
                    for (int k = 0; k < D; k++)
                    {
                        if (MinimumOfRandomVariables == true)
                        {
                            lower[k] = k == i ? Normal.StandardZ(1E-16) : Normal.StandardZ(Distributions[k].CDF(bins[0].LowerBound));
                            upper[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[0].LowerBound)) : Normal.StandardZ(1 - 1E-16);
                        }
                        else
                        {
                            lower[k] = Normal.StandardZ(1E-16);
                            upper[k] = Normal.StandardZ(Distributions[i].CDF(bins[0].LowerBound));
                        }
                    }
                    dF[i][0] = _mvn.Interval(lower, upper);
                    if (double.IsNaN(dF[i][0])) dF[i][0] = 0;
                    dF[i][0] = Math.Max(0, Math.Min(1, dF[i][0]));
                    p[i][0] = dF[i][0];

                    // Record the remaining bins
                    for (int j = 0; j < bins.Count; j++)
                    {
                        x[i][j + 1] = bins[j].UpperBound;
                        for (int k = 0; k < D; k++)
                        {
                            if (MinimumOfRandomVariables == true)
                            {
                                lower[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[j].LowerBound)) : Normal.StandardZ(Distributions[k].CDF(bins[j].Midpoint));
                                upper[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[j].UpperBound)) : Normal.StandardZ(1 - 1E-16);
                            }
                            else
                            {
                                lower[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[j].LowerBound)) : Normal.StandardZ(1E-16);
                                upper[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[j].UpperBound)) : Normal.StandardZ(Distributions[k].CDF(bins[j].Midpoint));
                            }
                        }
                        dF[i][j + 1] = _mvn.Interval(lower, upper);
                        if (double.IsNaN(dF[i][j + 1])) dF[i][j + 1] = 0;
                        dF[i][j + 1] = Math.Max(0, Math.Min(1, dF[i][j + 1]));
                    }
                }
            }
            else if (Dependency == Probability.DependencyType.Independent || Dependency == Probability.DependencyType.PerfectlyPositive)
            {
                /* 
                * For perfect independent or perfectly positive,
                * use the "Delta Method" developed by Haden Smith and Dave Margo.
                * It is fast and accurate.
                */

                double F1 = 0, F2 = 0;
                var pm = new double[D];
                var pl = new double[D];
                var pu = new double[D];
                var ind = new int[D];
                ind.Fill(1);

                for (int i = 0; i < D; i++)
                {
                    x.Add(new double[bins.Count + 1]);
                    p.Add(new double[bins.Count + 1]);
                    dF.Add(new double[bins.Count + 1]);

                    // Record first bin
                    x[i][0] = bins[0].LowerBound;
                    for (int k = 0; k < D; k++)
                    {
                        if (MinimumOfRandomVariables == true)
                        {
                            pl[k] = k == i ? Distributions[k].CCDF(bins[0].LowerBound) : Distributions[k].CCDF(bins[0].LowerBound);
                            pu[k] = k == i ? 1.0 : Distributions[k].CCDF(bins[0].LowerBound);
                        }
                        else
                        {
                            pu[k] = Distributions[k].CDF(bins[0].LowerBound);
                        }
                    }
                    F1 = Probability.JointProbability(pl, Dependency);
                    F2 = Probability.JointProbability(pu, Dependency);
                    dF[i][0] = F2 - F1;
                    if (double.IsNaN(dF[i][0]))
                        dF[i][0] = 0;
                    dF[i][0] = Math.Max(0, Math.Min(1, dF[i][0]));
                    p[i][0] = dF[i][0];

                    // Record remaining bins
                    for (int j = 0; j < bins.Count; j++)
                    {
                        x[i][j + 1] = bins[j].UpperBound;
                        for (int k = 0; k < D; k++)
                        {
                            if (MinimumOfRandomVariables == true)
                            {
                                pl[k] = k == i ? Distributions[k].CCDF(bins[j].UpperBound) : Distributions[k].CCDF(bins[j].Midpoint);
                                pu[k] = k == i ? Distributions[k].CCDF(bins[j].LowerBound) : Distributions[k].CCDF(bins[j].Midpoint);
                            }
                            else
                            {
                                pl[k] = k == i ? Distributions[k].CDF(bins[j].LowerBound) : Distributions[k].CDF(bins[j].Midpoint);
                                pu[k] = k == i ? Distributions[k].CDF(bins[j].UpperBound) : Distributions[k].CDF(bins[j].Midpoint);
                            }

                        }
                        F1 = Probability.JointProbability(pl, Dependency);
                        F2 = Probability.JointProbability(pu, Dependency);
                        dF[i][j + 1] = F2 - F1;
                        if (double.IsNaN(dF[i][j + 1]))
                            dF[i][j + 1] = 0;
                        dF[i][j + 1] = Math.Max(0, Math.Min(1, dF[i][j + 1]));
                    }
                }
            }

            // Get cumulative probabilities and make sure they sum <= 1 across D
            bool fixDF = false;
            var sum = new double[bins.Count + 1];
            for (int j = 1; j <= bins.Count; j++)
            {
                for (int i = 0; i < D; i++)
                {
                    sum[j] += p[i][j - 1] + dF[i][j];
                    p[i][j] = Math.Max(0, Math.Min(1, p[i][j - 1] + dF[i][j]));
                }
                if (sum[j] > 1 && sum[j] != sum[j - 1] && fixDF == false)
                {
                    double s = 0;
                    for (int i = 0; i < D; i++)
                    {
                        dF[i][j] *= (1 - sum[j - 1]) / (sum[j] - sum[j - 1]);
                        s += p[i][j - 1] + dF[i][j];

                        p[i][j] = Math.Max(0, Math.Min(1, p[i][j - 1] + dF[i][j]));
                    }
                    sum[j] = s;
                    fixDF = true;
                }
                else if (fixDF == true)
                {
                    for (int i = 0; i < D; i++)
                    {
                        dF[i][j] = 0;
                        p[i][j] = Math.Max(0, Math.Min(1, p[i][j - 1] + dF[i][j]));
                    }
                }
            }

            // Return CIFs
            var CIFs = new List<EmpiricalDistribution>();
            for (int i = 0; i < D; i++)
                CIFs.Add(new EmpiricalDistribution(x[i], p[i]));

            return CIFs;

        }

        /// <summary>
        /// Create a Multivariate Normal distribution used for modeling dependency between the marginal distributions.
        /// </summary>
        private void CreateMultivariateNormal()
        {
            var D = Distributions.Count();
            var mu = new double[D];
            var sigma = new double[D, D];
            if (Dependency == Probability.DependencyType.PerfectlyNegative)
            {
                CorrelationMatrix = new double[D, D];
                double rho = -1d / (D - 1d) + Math.Sqrt(Tools.DoubleMachineEpsilon);
                for (int i = 0; i < D; i++)
                {
                    mu[i] = 0d;
                    for (int j = 0; j < D; j++)
                        sigma[i, j] = i == j ? 1d : rho;
                }
            }
            else
            {
                for (int i = 0; i < D; i++)
                {
                    mu[i] = 0d;
                    for (int j = 0; j < D; j++)
                        sigma[i, j] = CorrelationMatrix[i, j];
                }
            }
            _mvn = new MultivariateNormal(mu, sigma);
            _mvnCreated = true;
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
            var bins = Stratify.XValues(new StratificationOptions(minX, maxX, binN, false), XTransform == Transform.Logarithmic ? true : false);
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
            var sample = new double[sampleSize];
            // Generate values
            for (int i = 0; i < sampleSize; i++)
            {
                double xMin = double.MaxValue;
                double xMax = double.MinValue;
                for (int j = 0; j < Distributions.Count; j++)
                {
                    var x = Distributions[j].InverseCDF(rnd.NextDouble());
                    if (x < xMin) xMin = x;
                    if (x > xMax) xMax = x;
                }
                sample[i] = MinimumOfRandomVariables == true ? xMin : xMax;
            }
            // Return array of random values
            return sample;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            var dists = new UnivariateDistributionBase[Distributions.Count];
            for (int i = 0; i < Distributions.Count; i++)
                dists[i] = Distributions[i].Clone();

            var cr = new CompetingRisks(dists)
            {
                MinimumOfRandomVariables = MinimumOfRandomVariables,
                Dependency = Dependency,
                XTransform = XTransform,
                ProbabilityTransform = ProbabilityTransform
            };
            if (CorrelationMatrix != null)
                cr.CorrelationMatrix = CorrelationMatrix.Clone() as double[,];

            return cr;
        }

        /// <inheritdoc/>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue(nameof(XTransform), XTransform.ToString());
            result.SetAttributeValue(nameof(ProbabilityTransform), ProbabilityTransform.ToString());
            result.SetAttributeValue(nameof(Distributions), String.Join("|", Distributions.Select(x => x.Type)));
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
        /// Create a competing risks distribution from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new competing risks distribution.</returns>
        public static CompetingRisks FromXElement(XElement xElement)
        {
            UnivariateDistributionType type = UnivariateDistributionType.Deterministic;
            if (xElement.Attribute(nameof(UnivariateDistributionBase.Type)) != null)
            {
                Enum.TryParse(xElement.Attribute(nameof(UnivariateDistributionBase.Type)).Value, out type);

            }
            if (type == UnivariateDistributionType.CompetingRisks)
            {
                var distributions = new List<UnivariateDistributionBase>();
                if (xElement.Attribute(nameof(Distributions)) != null)
                {
                    var types = xElement.Attribute(nameof(Distributions)).Value.Split('|');
                    for (int i = 0; i < types.Length; i++)
                    {
                        Enum.TryParse(types[i], out UnivariateDistributionType distType);
                        distributions.Add(UnivariateDistributionFactory.CreateDistribution(distType));
                    }
                }
                var competingRisks = new CompetingRisks(distributions.ToArray());

                if (xElement.Attribute(nameof(XTransform)) != null)
                {
                    Enum.TryParse(xElement.Attribute(nameof(XTransform)).Value, out Transform xTransform);
                    competingRisks.XTransform = xTransform;
                }
                if (xElement.Attribute(nameof(ProbabilityTransform)) != null)
                {
                    Enum.TryParse(xElement.Attribute(nameof(ProbabilityTransform)).Value, out Transform probabilityTransform);
                    competingRisks.ProbabilityTransform = probabilityTransform;
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
                    competingRisks.SetParameters(parameters);
                }

                return competingRisks;
            }
            else
            {
                return null;
            }
        }
    }
}
