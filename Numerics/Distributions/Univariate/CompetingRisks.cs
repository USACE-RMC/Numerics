using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics;
using Numerics.Mathematics.LinearAlgebra;
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

        private UnivariateDistributionBase[] _distributions = null!;
        private EmpiricalDistribution _empiricalCDF = null!;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;
        private bool _empiricalCDFCreated = false;
        private double[,] _correlationMatrix = null!;
        private bool _mvnCreated = false;
        private Probability.DependencyType _dependency = Probability.DependencyType.Independent;
        private MultivariateNormal _mvn = null!;
        private int _prngSeed = MultivariateNormal.DefaultMVNUNISeed;

        // Soft finite floor used in tail arithmetic before returning the final log-density.
        private const double _logZero = -745.0;
        private const double _minDensity = 1E-300;

        /// <summary>
        /// Returns the array of univariate probability distributions.
        /// </summary>
        public ReadOnlyCollection<UnivariateDistributionBase> Distributions => new(_distributions);

        /// <summary>
        /// The seed for the multivariate normal's quadrature randomizer, used by the dependent
        /// (perfectly negative and correlation-matrix) branches.
        /// </summary>
        /// <remarks>
        /// Those branches evaluate a randomized-lattice rectangle integral drawing from
        /// <see cref="MultivariateNormal.MVNUNI"/>, so their results reproduce only when the seed
        /// is fixed. Applied when the multivariate normal is built — set it before the first
        /// dependent evaluation.
        /// </remarks>
        public int PRNGSeed
        {
            get { return _prngSeed; }
            set
            {
                if (_prngSeed == value) return;
                _prngSeed = value;
                _mvnCreated = false;
                _empiricalCDFCreated = false;
            }
        }

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
        public Probability.DependencyType Dependency
        {
            get { return _dependency; }
            set
            {
                if (_dependency != value)
                {
                    _dependency = value;
                    _mvnCreated = false;
                }
            }
        }

        /// <summary>
        /// The correlation matrix used for modeling dependency between the marginal distributions.
        /// This is only used when the Dependency Type = CorrelationMatrix.
        /// </summary>
        public double[,] CorrelationMatrix 
        { 
            get {  return _correlationMatrix; } 
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
                for (int i = 0; i < Distributions.Count; i++)
                {
                    Dstring += Distributions[i].DisplayName;
                    if (i < Distributions.Count - 1)
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
            _parametersValid = ValidateParameters(Array.Empty<double>(), false) is null;
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
            _parametersValid = ValidateParameters(Array.Empty<double>(), false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
            _mvnCreated = false;
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentException">Thrown when the flattened parameter count does not match the component distributions.</exception>
        public override void SetParameters(IList<double> parameters)
        {
            if (Distributions == null || Distributions.Count == 0)
            {
                _parametersValid = false;
                return;
            }
            if (parameters.Count != NumberOfParameters)
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }

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

            _parametersValid = ValidateParameters(parameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
            _mvnCreated = false;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (Distributions.Count == 0)
            {
                var exception = new ArgumentOutOfRangeException(nameof(Distributions), "There must be at least 1 distribution.");
                if (throwException) throw exception;
                return exception;
            }
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
                if (MinimumOfRandomVariables)
                {
                    f = PDFMinimumIndependent(x);
                }
                else
                {
                    f = PDFMaximumIndependent(x);
                }
            }
            else
            {
                // Compute the PDF using numerical differentiation
                f = NumericalDerivative.Derivative(CDF, x);
            }

            // Return minimum density instead of zero to prevent log-likelihood issues
            return f < _minDensity ? _minDensity : f;
        }

        /// <summary>
        /// Computes PDF for minimum of independent random variables.
        /// f(x) = h(x) * S(x) where h(x) = sum of hazard rates, S(x) = product of survival functions
        /// </summary>
        private double PDFMinimumIndependent(double x)
        {
            double sumHazard = 0.0;
            double productSurvival = 1.0;

            for (int i = 0; i < Distributions.Count; i++)
            {
                double ccdf = Distributions[i].CCDF(x);
                double pdf = Distributions[i].PDF(x);

                productSurvival *= ccdf;

                // Safe hazard calculation
                if (ccdf > _minDensity)
                {
                    sumHazard += pdf / ccdf;
                }
                else if (pdf > _minDensity)
                {
                    // CCDF ≈ 0 but PDF > 0: we're at the boundary
                    // The hazard is very large, but productSurvival will be ≈ 0
                    // so the contribution is negligible
                    sumHazard += pdf / _minDensity; // Cap the hazard
                }
            }

            return sumHazard * productSurvival;
        }

        /// <summary>
        /// Computes PDF for maximum of independent random variables.
        /// f(x) = sum_i [f_i(x) * prod_{j≠i} F_j(x)]
        ///      = [sum_i (f_i/F_i)] * [prod_j F_j]
        /// </summary>
        private double PDFMaximumIndependent(double x)
        {
            double sumRatio = 0.0;
            double productCdf = 1.0;

            for (int i = 0; i < Distributions.Count; i++)
            {
                double cdf = Distributions[i].CDF(x);
                double pdf = Distributions[i].PDF(x);

                productCdf *= cdf;

                // Safe ratio calculation
                if (cdf > _minDensity)
                {
                    sumRatio += pdf / cdf;
                }
                else if (pdf > _minDensity)
                {
                    // CDF ≈ 0 but PDF > 0: we're at the left boundary
                    // The ratio is very large, but productCdf will be ≈ 0
                    // so the contribution is negligible
                    sumRatio += pdf / _minDensity; // Cap the ratio
                }
            }

            return sumRatio * productCdf;
        }

        /// <inheritdoc/>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            if (Distributions.Count == 1)
            {
                return Distributions[0].LogPDF(x);
            }

            // Only compute the exact LogPDF for independent random variables
            if (Dependency == Probability.DependencyType.Independent)
            {
                if (MinimumOfRandomVariables)
                {
                    return LogPDFMinimumIndependent(x);
                }
                else
                {
                    return LogPDFMaximumIndependent(x);
                }
            }
            else
            {
                // For dependent cases, fall back to numerical differentiation
                // but use a more stable approach
                double pdf = NumericalDerivative.Derivative(CDF, x);
                return pdf > _minDensity ? Math.Log(pdf) : double.NegativeInfinity;
            }

        }

        /// <summary>
        /// Computes log-PDF for minimum of independent random variables.
        /// Uses the formula: log(f(x)) = log(sum of hazards) + sum of log(survival functions)
        /// 
        /// For minimum: f(x) = h(x) * S(x) where:
        ///   h(x) = sum_i h_i(x) = sum_i [f_i(x) / S_i(x)]
        ///   S(x) = prod_i S_i(x)
        /// 
        /// In log space: log(f) = log(h(x)) + sum_i log(S_i(x))
        /// </summary>
        private double LogPDFMinimumIndependent(double x)
        {
            int n = Distributions.Count;
            var logSurvival = new double[n];
            var logHazard = new double[n];

            double sumLogSurvival = 0.0;
            bool allSurvivalZero = true;

            for (int i = 0; i < n; i++)
            {
                double ccdf = Distributions[i].CCDF(x);
                double pdf = Distributions[i].PDF(x);

                if (ccdf > _minDensity)
                {
                    logSurvival[i] = Math.Log(ccdf);
                    allSurvivalZero = false;
                }
                else
                {
                    // Survival is essentially zero - we're far in the right tail
                    logSurvival[i] = _logZero;
                }

                sumLogSurvival += logSurvival[i];

                // Compute log-hazard: log(f_i / S_i) = log(f_i) - log(S_i)
                if (pdf > _minDensity && ccdf > _minDensity)
                {
                    logHazard[i] = Math.Log(pdf) - Math.Log(ccdf);
                }
                else if (pdf <= _minDensity)
                {
                    logHazard[i] = _logZero;
                }
                else
                {
                    // pdf > 0 but ccdf ≈ 0: hazard is very large
                    // This happens in the far right tail
                    logHazard[i] = Math.Log(pdf) - _logZero; // Will be very large
                }
            }

            // If all survival functions are zero, density is zero
            if (allSurvivalZero)
                return _logZero;

            // Compute log of sum of hazards using log-sum-exp trick
            double logSumHazard = Tools.LogSumExp(logHazard);

            // Final result: log(f) = log(sum h_i) + sum log(S_i)
            double logPdf = logSumHazard + sumLogSurvival;

            return double.IsNaN(logPdf) || double.IsInfinity(logPdf) ? double.NegativeInfinity : logPdf;
        }

        /// <summary>
        /// Computes log-PDF for maximum of independent random variables.
        /// Uses the formula: f(x) = sum_i [f_i(x) * prod_{j≠i} F_j(x)]
        /// 
        /// In log space, we use the identity:
        ///   f(x) = [sum_i (f_i/F_i)] * [prod_j F_j]
        /// 
        /// So: log(f) = log(sum_i f_i/F_i) + sum_j log(F_j)
        /// </summary>
        private double LogPDFMaximumIndependent(double x)
        {
            int n = Distributions.Count;
            var logCdf = new double[n];
            var logRatio = new double[n]; // log(f_i / F_i)

            double sumLogCdf = 0.0;
            bool allCdfZero = true;

            for (int i = 0; i < n; i++)
            {
                double cdf = Distributions[i].CDF(x);
                double pdf = Distributions[i].PDF(x);

                if (cdf > _minDensity)
                {
                    logCdf[i] = Math.Log(cdf);
                    allCdfZero = false;
                }
                else
                {
                    // CDF is essentially zero - we're far in the left tail
                    logCdf[i] = _logZero;
                }

                sumLogCdf += logCdf[i];

                // Compute log(f_i / F_i) = log(f_i) - log(F_i)
                if (pdf > _minDensity && cdf > _minDensity)
                {
                    logRatio[i] = Math.Log(pdf) - Math.Log(cdf);
                }
                else if (pdf <= _minDensity)
                {
                    logRatio[i] = _logZero;
                }
                else
                {
                    // pdf > 0 but cdf ≈ 0: ratio is very large
                    // This can happen but contributes negligibly when multiplied by near-zero CDF product
                    logRatio[i] = Math.Log(pdf) - _logZero;
                }
            }

            // If all CDFs are zero, density is zero
            if (allCdfZero)
                return _logZero;

            // Compute log of sum of ratios using log-sum-exp trick
            double logSumRatio = Tools.LogSumExp(logRatio);

            // Final result: log(f) = log(sum f_i/F_i) + sum log(F_i)
            double logPdf = logSumRatio + sumLogCdf;

            return double.IsNaN(logPdf) || double.IsInfinity(logPdf) ? double.NegativeInfinity : logPdf;
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
        public List<EmpiricalDistribution> CumulativeIncidenceFunctions(List<StratificationBin>? bins = null)
        {
            // Get stratification bins
            if (bins == null)
            {
                double minP = 1E-16;
                double maxP = 1 - 1E-16;
                double minX = Distributions.Min(d => d.InverseCDF(minP));
                double maxX = Distributions.Max(d => d.InverseCDF(maxP));
                bins = Stratify.XValues(new StratificationOptions(minX, maxX, 200, false), true);
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
        /// Validates the user-supplied correlation matrix used by the Gaussian copula.
        /// </summary>
        /// <exception cref="ArgumentException">
        /// Thrown when the matrix is missing, has the wrong dimensions, contains invalid
        /// entries, is not symmetric with unit diagonal, or is not positive definite.
        /// </exception>
        private void ValidateCorrelationMatrix()
        {
            int dimension = Distributions.Count;
            if (CorrelationMatrix == null)
                throw new ArgumentException("A correlation-matrix dependency requires a correlation matrix.", nameof(CorrelationMatrix));
            if (CorrelationMatrix.GetLength(0) != dimension || CorrelationMatrix.GetLength(1) != dimension)
                throw new ArgumentException("The correlation matrix dimensions must match the number of distributions.", nameof(CorrelationMatrix));

            for (int i = 0; i < dimension; i++)
            {
                if (Math.Abs(CorrelationMatrix[i, i] - 1d) > 1E-12)
                    throw new ArgumentException("The correlation matrix must have unit diagonal entries.", nameof(CorrelationMatrix));

                for (int j = 0; j < dimension; j++)
                {
                    double value = CorrelationMatrix[i, j];
                    if (!Tools.IsFinite(value) || value < -1d || value > 1d)
                        throw new ArgumentException("The correlation matrix must contain finite values between -1 and 1.", nameof(CorrelationMatrix));
                    if (j > i && Math.Abs(value - CorrelationMatrix[j, i]) > 1E-12)
                        throw new ArgumentException("The correlation matrix must be symmetric.", nameof(CorrelationMatrix));
                }
            }

            try
            {
                var cholesky = new CholeskyDecomposition(new Matrix(CorrelationMatrix));
                if (!cholesky.IsPositiveDefinite)
                    throw new ArgumentException("The correlation matrix must be positive definite.", nameof(CorrelationMatrix));
            }
            catch (ArgumentException)
            {
                throw;
            }
            catch (Exception exception)
            {
                throw new ArgumentException("The correlation matrix must be positive definite.", nameof(CorrelationMatrix), exception);
            }
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
                ValidateCorrelationMatrix();
                for (int i = 0; i < D; i++)
                {
                    mu[i] = 0d;
                    for (int j = 0; j < D; j++)
                        sigma[i, j] = CorrelationMatrix[i, j];
                }
            }
            _mvn = new MultivariateNormal(mu, sigma) { MVNUNI = new MersenneTwister(PRNGSeed) };
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
            return GenerateRandomValuesWithDependency(sampleSize, seed);
        }

        /// <summary>
        /// Generates random values accounting for dependency structure.
        /// </summary>
        /// <param name="sampleSize"> Size of random sample to generate. </param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <returns>Array of random samples.</returns>
        public double[] GenerateRandomValuesWithDependency(int sampleSize, int seed = -1)
        {
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize];

            if (Dependency == Probability.DependencyType.Independent)
            {
                // Original implementation is correct for independent case
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
                    sample[i] = MinimumOfRandomVariables ? xMin : xMax;
                }
            }
            else if (Dependency == Probability.DependencyType.PerfectlyPositive)
            {
                // For perfectly positive dependency, all variables share the same quantile
                for (int i = 0; i < sampleSize; i++)
                {
                    double u = rnd.NextDouble();
                    double xMin = double.MaxValue;
                    double xMax = double.MinValue;
                    for (int j = 0; j < Distributions.Count; j++)
                    {
                        var x = Distributions[j].InverseCDF(u); // Same u for all
                        if (x < xMin) xMin = x;
                        if (x > xMax) xMax = x;
                    }
                    sample[i] = MinimumOfRandomVariables ? xMin : xMax;
                }
            }
            else if (Dependency == Probability.DependencyType.PerfectlyNegative ||
                     Dependency == Probability.DependencyType.CorrelationMatrix)
            {
                // Use Gaussian copula for correlation structure
                if (_mvnCreated == false)
                    CreateMultivariateNormal();

                // Generate correlated standard normal samples
                var mvnSamples = _mvn.GenerateRandomValues(sampleSize, seed);

                for (int i = 0; i < sampleSize; i++)
                {
                    double xMin = double.MaxValue;
                    double xMax = double.MinValue;

                    for (int j = 0; j < Distributions.Count; j++)
                    {
                        // Transform standard normal to uniform via Phi, then to marginal via inverse CDF
                        double z = mvnSamples[i, j];
                        double u = Normal.StandardCDF(z);
                        var x = Distributions[j].InverseCDF(u);

                        if (x < xMin) xMin = x;
                        if (x > xMax) xMax = x;
                    }
                    sample[i] = MinimumOfRandomVariables ? xMin : xMax;
                }
            }

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
                ProbabilityTransform = ProbabilityTransform,
                PRNGSeed = PRNGSeed
            };
            if (CorrelationMatrix != null)
                cr.CorrelationMatrix = (double[,])CorrelationMatrix.Clone();

            return cr;
        }

        /// <inheritdoc/>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue(nameof(XTransform), XTransform.ToString());
            result.SetAttributeValue(nameof(ProbabilityTransform), ProbabilityTransform.ToString());
            result.SetAttributeValue(nameof(MinimumOfRandomVariables), MinimumOfRandomVariables.ToString());
            result.SetAttributeValue(nameof(Dependency), Dependency.ToString());
            result.SetAttributeValue(nameof(PRNGSeed), PRNGSeed.ToString(CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Distributions), String.Join("|", Distributions.Select(x => x.Type)));
            // Parameters
            var parms = GetParameters;
            var parmStrings = new string[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("Parameters", String.Join("|", parmStrings));

            // Correlation matrix
            var corrMatrixElement = new XElement(nameof(CorrelationMatrix));
            if (CorrelationMatrix != null 
                && CorrelationMatrix.GetLength(0) == Distributions.Count 
                && CorrelationMatrix.GetLength(1) == Distributions.Count)
            {
                int rows = Distributions.Count;
                int cols = Distributions.Count;
                var row = new double[cols];

                for (int i = 0; i < rows; i++)
                {
                    var corrRowElement = new XElement("Correlation_Row");

                    // collect one row of the 2D array
                    for (int j = 0; j < cols; j++)
                    {
                        row[j] = _correlationMatrix[i, j];
                    }

                    // format each double to "G17" with invariant culture
                    var formatted = row.Select(v => v.ToString("G17", CultureInfo.InvariantCulture));
                    // join with '|' and set as the element's text
                    corrRowElement.Value = string.Join("|", formatted);
                    corrMatrixElement.Add(corrRowElement);
                }
            }
            result.Add(corrMatrixElement);

            return result;
        }

        /// <summary>
        /// Creates a competing-risks distribution from its serialized representation.
        /// </summary>
        /// <param name="xElement">The element to deserialize.</param>
        /// <returns>A validated competing-risks distribution, or <see langword="null"/> when the element identifies another distribution type.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when serialized configuration, parameters, or correlation data is malformed.</exception>
        public static CompetingRisks? FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));

            var typeAttribute = xElement.Attribute(nameof(UnivariateDistributionBase.Type));
            if (typeAttribute == null
                || !Enum.TryParse(typeAttribute.Value, out UnivariateDistributionType type)
                || !Enum.IsDefined(typeof(UnivariateDistributionType), type))
                throw new ArgumentException("The serialized distribution type is missing or invalid.", nameof(xElement));
            if (type != UnivariateDistributionType.CompetingRisks) return null;

            var distributionsAttribute = xElement.Attribute(nameof(Distributions));
            if (distributionsAttribute == null || string.IsNullOrWhiteSpace(distributionsAttribute.Value))
                throw new ArgumentException("The serialized competing-risks distribution has no component distributions.", nameof(xElement));

            string[] typeTokens = distributionsAttribute.Value.Split('|');
            var distributions = new UnivariateDistributionBase[typeTokens.Length];
            for (int i = 0; i < typeTokens.Length; i++)
            {
                if (!Enum.TryParse(typeTokens[i], out UnivariateDistributionType componentType)
                    || !Enum.IsDefined(typeof(UnivariateDistributionType), componentType))
                    throw new ArgumentException("The serialized competing-risks distribution contains an invalid component type.", nameof(xElement));
                distributions[i] = UnivariateDistributionFactory.CreateDistribution(componentType);
            }

            var competingRisks = new CompetingRisks(distributions);

            var xTransformAttribute = xElement.Attribute(nameof(XTransform));
            if (xTransformAttribute != null)
            {
                if (!Enum.TryParse(xTransformAttribute.Value, out Transform xTransform)
                    || !Enum.IsDefined(typeof(Transform), xTransform))
                    throw new ArgumentException("The serialized X transform is invalid.", nameof(xElement));
                competingRisks.XTransform = xTransform;
            }

            var probabilityTransformAttribute = xElement.Attribute(nameof(ProbabilityTransform));
            if (probabilityTransformAttribute != null)
            {
                if (!Enum.TryParse(probabilityTransformAttribute.Value, out Transform probabilityTransform)
                    || !Enum.IsDefined(typeof(Transform), probabilityTransform))
                    throw new ArgumentException("The serialized probability transform is invalid.", nameof(xElement));
                competingRisks.ProbabilityTransform = probabilityTransform;
            }

            var minimumAttribute = xElement.Attribute(nameof(MinimumOfRandomVariables));
            if (minimumAttribute != null)
            {
                if (!bool.TryParse(minimumAttribute.Value, out bool minimumOfRandomVariables))
                    throw new ArgumentException("The serialized minimum-selection flag is invalid.", nameof(xElement));
                competingRisks.MinimumOfRandomVariables = minimumOfRandomVariables;
            }

            var dependencyAttribute = xElement.Attribute(nameof(Dependency));
            if (dependencyAttribute != null)
            {
                if (!Enum.TryParse(dependencyAttribute.Value, out Probability.DependencyType dependency)
                    || !Enum.IsDefined(typeof(Probability.DependencyType), dependency))
                    throw new ArgumentException("The serialized dependency type is invalid.", nameof(xElement));
                competingRisks.Dependency = dependency;
            }

            var seedAttribute = xElement.Attribute(nameof(PRNGSeed));
            if (seedAttribute != null)
            {
                if (!int.TryParse(seedAttribute.Value, NumberStyles.Integer, CultureInfo.InvariantCulture, out int seed))
                    throw new ArgumentException("The serialized competing-risks seed is invalid.", nameof(xElement));
                competingRisks.PRNGSeed = seed;
            }

            var parametersAttribute = xElement.Attribute("Parameters");
            if (parametersAttribute == null)
                throw new ArgumentException("The serialized competing-risks parameters are missing.", nameof(xElement));
            string[] parameterTokens = parametersAttribute.Value.Split('|');
            if (parameterTokens.Length != competingRisks.NumberOfParameters)
                throw new ArgumentException("The serialized competing-risks parameter count is invalid.", nameof(xElement));
            var parameters = new double[parameterTokens.Length];
            for (int i = 0; i < parameters.Length; i++)
            {
                if (!double.TryParse(parameterTokens[i], NumberStyles.Any, CultureInfo.InvariantCulture, out parameters[i])
                    || !Tools.IsFinite(parameters[i]))
                    throw new ArgumentException("The serialized competing-risks parameters contain an invalid value.", nameof(xElement));
            }

            int offset = 0;
            for (int i = 0; i < distributions.Length; i++)
            {
                int count = distributions[i].NumberOfParameters;
                var componentParameters = new double[count];
                Array.Copy(parameters, offset, componentParameters, 0, count);
                distributions[i].ValidateParameters(componentParameters, true);
                offset += count;
            }
            competingRisks.SetParameters(parameters);
            if (!competingRisks.ParametersValid)
                throw new ArgumentException("The serialized competing-risks parameters are invalid.", nameof(xElement));

            var correlationElement = xElement.Element(nameof(CorrelationMatrix));
            var correlationRows = correlationElement?.Elements("Correlation_Row").ToArray() ?? Array.Empty<XElement>();
            if (correlationRows.Length > 0)
            {
                int dimension = distributions.Length;
                if (correlationRows.Length != dimension)
                    throw new ArgumentException("The serialized correlation matrix has an invalid row count.", nameof(xElement));

                var correlation = new double[dimension, dimension];
                for (int i = 0; i < dimension; i++)
                {
                    string[] entries = correlationRows[i].Value.Split('|');
                    if (entries.Length != dimension)
                        throw new ArgumentException("The serialized correlation matrix has an invalid column count.", nameof(xElement));
                    for (int j = 0; j < dimension; j++)
                    {
                        if (!double.TryParse(entries[j], NumberStyles.Any, CultureInfo.InvariantCulture, out correlation[i, j])
                            || !Tools.IsFinite(correlation[i, j])
                            || correlation[i, j] < -1d
                            || correlation[i, j] > 1d)
                            throw new ArgumentException("The serialized correlation matrix contains an invalid value.", nameof(xElement));
                    }
                }

                for (int i = 0; i < dimension; i++)
                {
                    if (Math.Abs(correlation[i, i] - 1d) > 1E-12)
                        throw new ArgumentException("The serialized correlation matrix must have unit diagonal entries.", nameof(xElement));
                    for (int j = i + 1; j < dimension; j++)
                    {
                        if (Math.Abs(correlation[i, j] - correlation[j, i]) > 1E-12)
                            throw new ArgumentException("The serialized correlation matrix must be symmetric.", nameof(xElement));
                    }
                }
                competingRisks.CorrelationMatrix = correlation;
            }
            else if (competingRisks.Dependency == Probability.DependencyType.CorrelationMatrix)
            {
                throw new ArgumentException("A correlation-matrix dependency requires serialized correlation data.", nameof(xElement));
            }

            return competingRisks;
        }

    }
}
