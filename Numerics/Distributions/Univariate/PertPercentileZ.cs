using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Pert percentile z distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    ///     In probability and statistics, the PERT distribution is a family of continuous probability distributions
    ///     defined by the minimum (a), most likely (c) and maximum (b) values that a variable can take.
    ///     It is a transformation of the four-parameter Beta distribution.
    /// </para>
    /// <para>
    ///     This version of the PERT is parameterized using the 5th, 50th, and 95th percentiles in Normal Z space.
    ///     The variable is X = Φ(Z), where Φ is the standard Normal CDF and Z follows a PERT distribution fitted to the
    ///     z-scores of the three percentiles, so X lies in [Φ(a), Φ(b)] within (0, 1), where a and b are the PERT
    ///     min and max in z-space. The density of X is f(x) = g(z) / φ(z) with z = Φ⁻¹(x), where g is the PERT density
    ///     and φ the standard Normal density. The moments are those of X.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/PERT_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class PertPercentileZ : UnivariateDistributionBase
    {

        /// <summary>
        /// Constructs a PERT distribution with 5th = 0.05 (z = -1.64), 50th = 0.5 (z = 0.0), and 95th = 0.95 (z = 1.64).
        /// </summary>
        public PertPercentileZ()
        {
            SetParameters(0.05, 0.5d, 0.95d);
        }

        /// <summary>
        /// Constructs a PERT distribution with specified 5th, 50th, and 95th percentiles.
        /// </summary>
        /// <param name="fifth">The 5th percentile, a probability in (0, 1) whose standard Normal z-score is the 5th percentile of the underlying PERT.</param>
        /// <param name="fiftieth">The 50th percentile, a probability in (0, 1) whose standard Normal z-score is the 50th percentile of the underlying PERT.</param>
        /// <param name="ninetyFifth">The 95th percentile, a probability in (0, 1) whose standard Normal z-score is the 95th percentile of the underlying PERT.</param>
        public PertPercentileZ(double fifth, double fiftieth, double ninetyFifth)
        {
            SetParameters(fifth, fiftieth, ninetyFifth);
        }

        private bool _parametersSolved = false;
        private GeneralizedBeta _beta = new GeneralizedBeta();
        private double _5th;
        private double _50th;
        private double _95th;

        /// <summary>
        /// Gets and sets the 5th percentile.
        /// </summary>
        public double Percentile5th
        {
            get { return _5th; }
            set
            {
                if (_5th != value)
                {
                    _5th = value;
                    // validate parameters
                    _parametersValid = ValidateParameters(_5th, _50th, _95th, false) is null;
                    _parametersSolved = false;
                }
            }
        }

        /// <summary>
        /// Gets and sets the 50th percentile.
        /// </summary>
        public double Percentile50th
        {
            get { return _50th; }
            set
            {
                if (_50th != value)
                {
                    _50th = value;
                    // validate parameters
                    _parametersValid = ValidateParameters(_5th, _50th, _95th, false) is null;
                    _parametersSolved = false;
                }
            }
        }

        /// <summary>
        /// Gets and sets the 95th percentile.
        /// </summary>
        public double Percentile95th
        {
            get { return _95th; }
            set
            {
                if (_95th != value)
                {
                    _95th = value;
                    // validate parameters
                    _parametersValid = ValidateParameters(_5th, _50th, _95th, false) is null;
                    _parametersSolved = false;
                }
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.PertPercentileZ; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "PERT-Percentile Z"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "PERT-% Z"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[3, 2];
                parmString[0, 0] = "5%";
                parmString[1, 0] = "50%";
                parmString[2, 0] = "95%";
                parmString[0, 1] = Percentile5th.ToString();
                parmString[1, 1] = Percentile50th.ToString();
                parmString[2, 1] = Percentile95th.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["5%", "50%", "95%"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Percentile5th), nameof(Percentile50th), nameof(Percentile95th)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Percentile5th, Percentile50th, Percentile95th]; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The mean of X = Φ(Z), E[X] = ∫ Φ(z) g(z) dz over the z-space support, evaluated by adaptive Gauss-Kronrod
        /// integration. It differs from Φ(E[Z]) because Φ is nonlinear.
        /// </remarks>
        public override double Mean
        {
            get
            {
                if (_parametersSolved == false) SolveParameters();
                if (_beta.Min.AlmostEquals(_beta.Max)) return Normal.StandardCDF(_beta.Min);
                return CentralMomentsOfX()[0];
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Φ is increasing, so the median of X is Φ applied to the median of the underlying PERT.
        /// </remarks>
        public override double Median
        {
            get
            {
                if (_parametersSolved == false) SolveParameters();
                if (_beta.Min.AlmostEquals(_beta.Max)) return Normal.StandardCDF(_beta.Min);
                return Normal.StandardCDF(_beta.Median);
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The maximizer of the density f(x) = g(z) / φ(z) of X. Because Φ is increasing, it is Φ(z*), where z*
        /// maximizes log g(z) - log φ(z) over the z-space support.
        /// </remarks>
        public override double Mode
        {
            get
            {
                if (_parametersSolved == false) SolveParameters();
                if (_beta.Min.AlmostEquals(_beta.Max)) return Normal.StandardCDF(_beta.Min);
                return ModeOfX();
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The standard deviation of X = Φ(Z), evaluated by adaptive Gauss-Kronrod integration over the z-space support.
        /// A fixed-value distribution has a standard deviation of zero.
        /// </remarks>
        public override double StandardDeviation
        {
            get { return CentralMomentsOfX()[1]; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The skewness of X = Φ(Z), evaluated by adaptive Gauss-Kronrod integration over the z-space support.
        /// </remarks>
        public override double Skewness
        {
            get { return CentralMomentsOfX()[2]; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The kurtosis (not excess kurtosis) of X = Φ(Z), evaluated by adaptive Gauss-Kronrod integration over the
        /// z-space support.
        /// </remarks>
        public override double Kurtosis
        {
            get { return CentralMomentsOfX()[3]; }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get
            {
                if (_parametersSolved == false) SolveParameters();
                return Normal.StandardCDF(_beta.Minimum);
            }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get
            {
                if (_parametersSolved == false) SolveParameters();
                return Normal.StandardCDF(_beta.Maximum);
            }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [0d, Percentile5th, Percentile50th]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [Percentile50th, Percentile95th, 1d]; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="fifth">The 5th percentile value of the distribution.</param>
        /// <param name="fiftieth">The 50th percentile value of the distribution.</param>
        /// <param name="ninetyFifth">The 95th percentile value of the distribution.</param>
        public void SetParameters(double fifth, double fiftieth, double ninetyFifth)
        {
            Percentile5th = fifth;
            Percentile50th = fiftieth;
            Percentile95th = ninetyFifth;
            // The setters validate only on change, so validate the full set once more
            _parametersValid = ValidateParameters(_5th, _50th, _95th, false) is null;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="fifth">The 5th percentile value of the distribution.</param>
        /// <param name="fiftieth">The 50th percentile value of the distribution.</param>
        /// <param name="ninetyFifth">The 95th percentile value of the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        private ArgumentOutOfRangeException? ValidateParameters(double fifth, double fiftieth, double ninetyFifth, bool throwException)
        {
            if (double.IsNaN(fifth) || double.IsInfinity(fifth) ||
                double.IsNaN(ninetyFifth) || double.IsInfinity(ninetyFifth) || fifth > ninetyFifth)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Percentile5th), "The 5% cannot be greater than the 95%.");
                return new ArgumentOutOfRangeException(nameof(Percentile5th), "The 5% cannot be greater than the 95%.");
            }
            else if (double.IsNaN(fiftieth) || double.IsInfinity(fiftieth) || fiftieth < fifth || fiftieth > ninetyFifth)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Percentile50th), "The 50% must be between the 5% and 95%.");
                return new ArgumentOutOfRangeException(nameof(Percentile50th), "The 50% must be between the 5% and 95%.");
            }
            else if (fifth <= 0 || fifth >= 1)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Percentile5th), "The percentiles must be strictly between 0 and 1.");
                return new ArgumentOutOfRangeException(nameof(Percentile5th), "The percentiles must be strictly between 0 and 1.");
            }
            else if (fiftieth <= 0 || fiftieth >= 1)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Percentile50th), "The percentiles must be strictly between 0 and 1.");
                return new ArgumentOutOfRangeException(nameof(Percentile50th), "The percentiles must be strictly between 0 and 1.");
            }
            else if (ninetyFifth <= 0 || ninetyFifth >= 1)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Percentile95th), "The percentiles must be strictly between 0 and 1.");
                return new ArgumentOutOfRangeException(nameof(Percentile95th), "The percentiles must be strictly between 0 and 1.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }

        
        /// <summary>
        /// Solve the PERT parameters (min, mode, max) given the percentiles.
        /// </summary>
        public void SolveParameters()
        {
            if (_parametersValid == false) return;
            if (_parametersSolved == true) return;

            double fifth =  Percentile5th;
            double fiftieth = Percentile50th;
            double ninetyFifth = Percentile95th;

            if (fifth == fiftieth && fiftieth == ninetyFifth)
            {
                _beta = GeneralizedBeta.PERT(Normal.StandardZ(fifth), Normal.StandardZ(fiftieth), Normal.StandardZ(ninetyFifth));
                _parametersSolved = true;
                return;
            }
          
            fifth = Normal.StandardZ(Percentile5th);
            fiftieth = Normal.StandardZ(Percentile50th);
            ninetyFifth = Normal.StandardZ(Percentile95th);

            _beta = GeneralizedBeta.PERT(fifth, fiftieth, ninetyFifth);
            // Set wide bounds
            double min = Math.Max(Normal.StandardZ(1E-16), fifth - (ninetyFifth - fifth) * 2);
            double max = Math.Min(Normal.StandardZ(1 - 1E-16), ninetyFifth + (ninetyFifth - fifth) * 2);
            var Initials = new double[] { _beta.Min, _beta.Mode, _beta.Max };
            var Lowers = new double[] { min, min, min };
            var Uppers = new double[] { max, max, max };

            // Objective function is the sum of squared errors
            double sse(double[] x)
            {
                var dist = new Pert();
                try
                {
                    dist = new Pert(x[0], x[1], x[2]);
                    if (dist.ParametersValid == false) return double.MaxValue;

                    double SSE = 0d;
                    SSE += Math.Pow(fifth - dist.InverseCDF(0.05), 2d);
                    SSE += Math.Pow(fiftieth - dist.InverseCDF(0.5), 2d);
                    SSE += Math.Pow(ninetyFifth - dist.InverseCDF(0.95), 2d);
                    return SSE;
                }
                catch
                {
                    return double.MaxValue;
                }
            }

            var solver = new NelderMead(sse, 3, Initials, Lowers, Uppers);
            solver.MaxFunctionEvaluations = 200;
            solver.RelativeTolerance = 1E-8;
            solver.AbsoluteTolerance = 1E-8;
            solver.ComputeHessian = false;
            solver.ReportFailure = false;
            solver.Minimize();
            var solution = solver.BestParameterSet.Values;
            _beta = GeneralizedBeta.PERT(solution[0], solution[1], solution[2]);
            _parametersSolved = true;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Returns f(x) = g(z) / φ(z) with z = Φ⁻¹(x), where g is the PERT density in z-space and φ the standard Normal
        /// density; 1 / φ(z) is the Jacobian dz/dx. Returns zero outside the support, including x ≤ 0 and x ≥ 1.
        /// </remarks>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters(Percentile5th, Percentile50th, Percentile95th, true);
            if (_parametersSolved == false) SolveParameters();
            if (x <= 0.0d || x >= 1.0d) return 0.0d;
            // These checks are done specifically for an application where a
            // user inputs min = max = mode
            if (_beta.Min.AlmostEquals(_beta.Max) && _beta.Min.AlmostEquals(_beta.Mode)) return 0.0d;
            if (double.IsNaN(_beta.Mode)) return 0.0d;
            //
            double z = Normal.StandardZ(x);
            double density = _beta.PDF(z);
            if (density == 0.0d) return 0.0d;
            return density / Normal.StandardPDF(z);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Returns G(Φ⁻¹(x)), where G is the PERT CDF in z-space. Returns zero for x ≤ 0 and one for x ≥ 1.
        /// </remarks>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters(Percentile5th, Percentile50th, Percentile95th, true);
            if (_parametersSolved == false) SolveParameters();
            if (x <= 0.0d) return 0.0d;
            if (x >= 1.0d) return 1.0d;
            // These checks are done specifically for an application where a 
            // user inputs min = max = mode
            if (_beta.Min.AlmostEquals(_beta.Max) && _beta.Min.AlmostEquals(_beta.Mode)) return 1d;
            if (double.IsNaN(_beta.Mode)) return 1d;
            //
            return _beta.CDF(Normal.StandardZ(x));
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters(Percentile5th, Percentile50th, Percentile95th, true);
            if (_parametersSolved == false) SolveParameters();
            var x = Normal.StandardCDF(_beta.Min);
            // These checks are done specifically for an application where a 
            // user inputs min = max = mode
            if (_beta.Min.AlmostEquals(_beta.Max) && _beta.Min.AlmostEquals(_beta.Mode)) return x;
            if (double.IsNaN(_beta.Mode)) return x;
            //
            try
            {
                x = Normal.StandardCDF(_beta.InverseCDF(probability));
            }
            catch (ArithmeticException)
            {
                return x;
            }
            return x;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new PertPercentileZ(Percentile5th, Percentile50th, Percentile95th) { _parametersSolved = _parametersSolved, _parametersValid = _parametersValid, _beta = (GeneralizedBeta)_beta.Clone()};
        }

        /// <summary>
        /// Returns the PERT distribution of the standard Normal z-scores that underlies this distribution.
        /// </summary>
        /// <returns>The PERT distribution with the solved z-space min, mode and max. When the 5th and 95th percentiles
        /// are equal, this is the fixed-value PERT whose min, mode and max all equal their common z-score.</returns>
        public Pert ToPert()
        {
            if (_parametersSolved == false) SolveParameters();
            if (_parametersValid && Percentile5th == Percentile95th)
            {
                double z = Normal.StandardZ(Percentile5th);
                return new Pert(z, z, z);
            }
            return new Pert(_beta.Min, _beta.Mode, _beta.Max);
        }

        /// <summary>
        /// Computes the mean, standard deviation, skewness and kurtosis of X = Φ(Z).
        /// </summary>
        /// <returns>The mean, standard deviation, skewness and kurtosis, in that order.</returns>
        /// <remarks>
        /// Each expectation E[h(X)] = ∫ h(Φ(z)) g(z) dz is integrated over the finite z-space support [a, b] of the
        /// underlying PERT with adaptive Gauss-Kronrod quadrature at a relative tolerance of 1E-12. Integrating in z-space
        /// keeps the integrand bounded and smooth. A fixed-value distribution returns its value, a standard deviation of
        /// zero, and undefined skewness and kurtosis. A moment whose integration fails is returned as
        /// <see cref="double.NaN"/>.
        /// </remarks>
        private double[] CentralMomentsOfX()
        {
            if (_parametersSolved == false) SolveParameters();
            double a = _beta.Min, b = _beta.Max;
            if (a.AlmostEquals(b)) return [Normal.StandardCDF(a), 0d, double.NaN, double.NaN];
            double Expectation(Func<double, double> h)
            {
                var agk = new AdaptiveGaussKronrod(z => h(Normal.StandardCDF(z)) * _beta.PDF(z), a, b)
                {
                    RelativeTolerance = 1E-12, AbsoluteTolerance = 1E-15, ReportFailure = false
                };
                agk.Integrate();
                return agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            }
            double mean = Expectation(x => x);
            double sd = Math.Sqrt(Expectation(x => (x - mean) * (x - mean)));
            double skew = Expectation(x => Math.Pow((x - mean) / sd, 3d));
            double kurt = Expectation(x => Math.Pow((x - mean) / sd, 4d));
            return [mean, sd, skew, kurt];
        }

        /// <summary>
        /// Returns the mode of X = Φ(Z) by maximizing the log density of X over the z-space support.
        /// </summary>
        /// <returns>The value of X at which its density is largest.</returns>
        /// <remarks>
        /// The objective is log g(z) - log φ(z), whose maximizer z* gives the mode Φ(z*). The added z²/2 term is convex,
        /// so the objective can have more than one local maximum and can peak at an end of the support. The objective is
        /// therefore evaluated at 101 equally spaced points, including both ends, and the best point is refined with
        /// Brent's method on its neighbouring intervals. An interior maximum is then polished with Newton steps on the
        /// stationarity condition (α - 1) / (z - a) - (β - 1) / (b - z) + z = 0, where α and β are the PERT shapes.
        /// </remarks>
        private double ModeOfX()
        {
            double a = _beta.Min, b = _beta.Max;
            double LogDensity(double z) => _beta.LogPDF(z) - Normal.StandardLogPDF(z);
            const int intervals = 100;
            double h = (b - a) / intervals;
            int best = 0;
            double bestValue = double.NegativeInfinity;
            for (int i = 0; i <= intervals; i++)
            {
                double value = LogDensity(i == intervals ? b : a + i * h);
                if (value > bestValue)
                {
                    bestValue = value;
                    best = i;
                }
            }
            double zMode = best == intervals ? b : a + best * h;
            double lower = best == 0 ? a : a + (best - 1) * h;
            double upper = best >= intervals - 1 ? b : a + (best + 1) * h;
            var brent = new BrentSearch(LogDensity, lower, upper) { RelativeTolerance = 1E-12, ReportFailure = false };
            brent.Maximize();
            if (brent.Status == OptimizationStatus.Success && -brent.BestParameterSet.Fitness >= bestValue)
                zMode = brent.BestParameterSet.Values[0];
            // Polish an interior maximum with Newton steps on the stationarity condition
            double alpha1 = _beta.Alpha - 1d, beta1 = _beta.Beta - 1d;
            for (int k = 0; k < 3; k++)
            {
                double left = zMode - a, right = b - zMode;
                double d1 = alpha1 / left - beta1 / right + zMode;
                double d2 = 1d - alpha1 / (left * left) - beta1 / (right * right);
                if (!Tools.IsFinite(d1) || !(d2 < 0d)) break;
                double next = zMode - d1 / d2;
                if (!(next > lower && next < upper)) break;
                zMode = next;
            }
            return Normal.StandardCDF(zMode);
        }

        /// <inheritdoc/>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue("ParametersSolved", _parametersSolved.ToString());
            // Parameters
            var parms = GetParameters;
            var parmStrings = new string[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("Parameters", String.Join("|", parmStrings));
            // Beta Distribution parameters
            parms = _beta.GetParameters;
            parmStrings = new string[_beta.NumberOfParameters];
            for (int i = 0; i < _beta.NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("BetaParameters", String.Join("|", parmStrings));
            return result;
        }

        /// <summary>
        /// Create a Pert-Percentile distribution from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new PERT-PercentileZ distribution, or null if the element cannot be read.</returns>
        public static PertPercentileZ? FromXElement(XElement xElement)
        {
            UnivariateDistributionType type = UnivariateDistributionType.Deterministic;
            var typeAttr = xElement.Attribute(nameof(UnivariateDistributionBase.Type));
            if (typeAttr != null)
            {
                Enum.TryParse(typeAttr.Value, out type);

            }
            if (type == UnivariateDistributionType.PertPercentileZ)
            {
                bool parametersSolved = false;
                var paramsSolvedAttr = xElement.Attribute("ParametersSolved");
                if (paramsSolvedAttr != null)
                {
                    bool.TryParse(paramsSolvedAttr.Value, out parametersSolved);
                }
                else
                {
                    // This was saved from Base method
                    var dist = new PertPercentileZ();
                    var names = dist.GetParameterPropertyNames;
                    var parms = dist.GetParameters;
                    var vals = new double[dist.NumberOfParameters];
                    for (int i = 0; i < dist.NumberOfParameters; i++)
                    {
                        var nameAttr = xElement.Attribute(names[i]);
                        if (nameAttr != null)
                        {
                            double.TryParse(nameAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out vals[i]);
                        }
                    }
                    dist.SetParameters(vals);
                    return dist;

                }

                var beta = new GeneralizedBeta();
                var betaParamsAttr = xElement.Attribute("BetaParameters");
                if (betaParamsAttr != null)
                {
                    var vals = betaParamsAttr.Value.Split('|');
                    var parameters = new List<double>();
                    for (int i = 0; i < vals.Length; i++)
                    {
                        double.TryParse(vals[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var parm);
                        parameters.Add(parm);
                    }
                    beta.SetParameters(parameters);
                }
                double _5th = 0, _50th = 0, _95th = 0;
                var paramsAttr = xElement.Attribute("Parameters");
                if (paramsAttr != null)
                {
                    var vals = paramsAttr.Value.Split('|');
                    double.TryParse(vals[0], NumberStyles.Any, CultureInfo.InvariantCulture, out _5th);
                    double.TryParse(vals[1], NumberStyles.Any, CultureInfo.InvariantCulture, out _50th);
                    double.TryParse(vals[2], NumberStyles.Any, CultureInfo.InvariantCulture, out _95th);
                }

                var pert = new PertPercentileZ()
                {
                    _parametersSolved = parametersSolved,
                    _beta = beta,
                    _5th = _5th,
                    _50th = _50th,
                    _95th = _95th
                };

                return pert;
            }
            else
            {
                return null;
            }
        }

    }
}