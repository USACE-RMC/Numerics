using System;
using System.Collections.Generic;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Noncentral t probability distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Noncentral_t-distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class NoncentralT : UnivariateDistributionBase
    {
 
        /// <summary>
        /// Constructs a Noncentral t distribution with 10 degrees of freedom and noncentrality = 0.
        /// </summary>
        public NoncentralT()
        {
            SetParameters(10d, 0d);
        }

        /// <summary>
        /// Constructs a Noncentral t distribution with given degrees of freedom and noncentrality.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom.</param>
        /// <param name="noncentrality">The noncentrality parameter.</param>
        public NoncentralT(double degreesOfFreedom, double noncentrality)
        {
            SetParameters(degreesOfFreedom, noncentrality);
        }

        private double _degreesOfFreedom;
        private double _noncentrality;
       
        /// <summary>
        /// Gets and sets the degrees of freedom ν (nu) of the distribution.
        /// </summary>
        public double DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
            set
            {
                _parametersValid = ValidateParameters(value, Noncentrality, false) is null;
                _degreesOfFreedom = value;
            }
        }

        /// <summary>
        /// Gets and sets the noncentrality parameter μ (mu) of the distribution.
        /// </summary>
        public double Noncentrality
        {
            get { return _noncentrality; }
            set
            {
                _parametersValid = ValidateParameters(DegreesOfFreedom, value, false) is null;
                _noncentrality = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.NoncentralT; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Noncentral t"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "NCT"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Degrees of Freedom (ν)";
                parmString[1, 0] = "Noncentrality (μ)";
                parmString[0, 1] = DegreesOfFreedom.ToString();
                parmString[1, 1] = Noncentrality.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ν", "μ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(DegreesOfFreedom), nameof(Noncentrality)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [DegreesOfFreedom, Noncentrality]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (_degreesOfFreedom > 1)
                {
                    return Noncentrality * Math.Sqrt(DegreesOfFreedom / 2d) * Gamma.Function((DegreesOfFreedom - 1) / 2d) / Gamma.Function(DegreesOfFreedom / 2d);
                }
                else
                {
                    return double.NaN;
                }
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
                if (DegreesOfFreedom > 2)
                {
                    double a = DegreesOfFreedom * (1d + Math.Pow(Noncentrality, 2d)) / (DegreesOfFreedom - 2);
                    double b = DegreesOfFreedom * Math.Pow(Noncentrality, 2d) / 2d;
                    double c = Gamma.Function((DegreesOfFreedom - 1) / 2d) / Gamma.Function(DegreesOfFreedom / 2d);
                    return Math.Sqrt(a - b * c * c);
                }
                else
                {
                    return double.NaN;
                }
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {                
                return CentralMoments(1E-8)[2];
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                return CentralMoments(1E-8)[3];
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return double.NegativeInfinity; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [1.0d, double.NegativeInfinity]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="v">The degrees of freedom ν (nu). Range: ν > 0.</param>
        /// <param name="mu">The noncentrality parameter μ (mu).</param>
        public void SetParameters(double v, double mu)
        {
            _degreesOfFreedom = v;
            _noncentrality = mu;
            _parametersValid = ValidateParameters(v, mu, false) is null;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="v">The degrees of freedom ν (nu). Range: ν > 0.</param>
        /// <param name="mu">The noncentrality parameter μ (mu).</param>
        /// <param name="throwException">Whether to throw the validation exception.</param>
        /// <returns><see langword="null"/> when the parameters are valid; otherwise, the validation exception.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="throwException"/> is true and either parameter is invalid.</exception>
        public ArgumentOutOfRangeException? ValidateParameters(double v, double mu, bool throwException)
        {
            if (double.IsNaN(v) || double.IsInfinity(v) || v < 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
                return new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
            }
            if (double.IsNaN(mu) || double.IsInfinity(mu))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Noncentrality), "The noncentrality parameter μ (mu) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Noncentrality), "The noncentrality parameter μ (mu) must be a number.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            double u = Noncentrality; 
            double v = DegreesOfFreedom;
            if (x != 0)
            {
                double A = EvaluateCdfWithFallback(x * Math.Sqrt(1 + 2 / v), v + 2, u);
                double B = EvaluateCdfWithFallback(x, v, u);
                double C = v / x;
                return C * (A - B);
            }
            else
            {
                double A = Gamma.Function((v + 1) / 2);
                double B = Math.Sqrt(Math.PI * v) * Gamma.Function(v / 2);
                double C = Math.Exp(-(u * u) / 2);
                return (A / B) * C;
            }
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            return EvaluateCdfSeries(x, DegreesOfFreedom, Noncentrality);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            // 
            return InverseCdf(probability, DegreesOfFreedom, Noncentrality);
        }

        
        /// <summary>
        /// Cumulative probability at T of the non-central t-distribution
        /// with DF degrees of freedom (may be fractional) and non-centrality
        /// parameter DELTA.
        /// </summary>
        /// <param name="t">A single point in the distribution range.</param>
        /// <param name="df">The degrees of freedom.</param>
        /// <param name="delta">The noncentrality parameter.</param>
        /// <returns>The cumulative probability, using a normal approximation if the series does not converge.</returns>
        private static double EvaluateCdfWithFallback(double t, double df, double delta)
        {
            double result;
            try
            {
                result = EvaluateCdfSeries(t, df, delta);
            }
            catch (ArgumentException)
            {
                // If the series fails to converge, use normal approximation
                double z = (t * (1.0d - 1.0d / (4.0d * df)) - delta) / Math.Sqrt(1.0d + t * t / (2.0d * df));
                result = Normal.StandardCDF(z);
                if (result < 0d) result = 0d;
                if (result > 1d) result = 1d;
            }
            return result;
        }

        /// <summary>
        /// Cumulative probability at T of the non-central t-distribution
        /// with DF degrees of freedom (may be fractional) and non-centrality
        /// parameter DELTA.
        /// </summary>
        /// <param name="t">A single point in the distribution range.</param>
        /// <param name="df">The degrees of freedom.</param>
        /// <param name="delta">The noncentrality parameter.</param>
        /// <returns>The cumulative probability at <paramref name="t"/>.</returns>
        /// <exception cref="ArgumentException">Thrown when the AS 243 series exceeds its iteration limit.</exception>
        /// <remarks>
        /// The function is based on ALGORITHM AS 243  APPL. STATIST. (1989), VOL.38, NO. 1.
        /// Original FORTRAN code can be found at:
        /// http://people.sc.fsu.edu/~jburkardt/f77_src/asa243/asa243.html
        /// </remarks>
        private static double EvaluateCdfSeries(double t, double df, double delta)
        {
            double a;
            double logBeta;
            double b;
            double adjustedDelta;
            int iteration = 0;
            double errorBound;
            double gEven;
            double gOdd;
            double lambda;
            double pWeight;
            double qWeight;
            double oneMinusXPowerB;
            double remainingWeight;
            double x;
            double xEven;
            double xOdd;
            double result;
            bool reflected;
            const int maxIterations = 10000;
            const double maxError = 0.000000001d;
            const double zero = 0d;
            const double half = 0.5d;
            const double one = 1.0d;
            const double two = 2.0d;
            const double sqrtTwoOverPi = 0.797884560802865d;
            const double logSqrtPi = 0.5723649429247d;
            result = zero;
            adjustedDelta = delta;
            reflected = false;
            if (t < zero)
            {
                reflected = true;
                adjustedDelta = -adjustedDelta;
            }

            // Initialize the twin series of Guenther (1978).
            x = t * t / (t * t + df);
            if (x > zero)
            {
                lambda = adjustedDelta * adjustedDelta;
                pWeight = half * Math.Exp(-half * lambda);
                qWeight = sqrtTwoOverPi * pWeight * adjustedDelta;
                remainingWeight = half - pWeight;
                a = half;
                b = half * df;
                oneMinusXPowerB = Math.Pow(one - x, b);
                logBeta = logSqrtPi + Gamma.LogGamma(b) - Gamma.LogGamma(a + b);
                xOdd = Beta.IncompleteRatio(x, a, b, logBeta);
                gOdd = two * oneMinusXPowerB * Math.Exp(a * Math.Log(x) - logBeta);
                xEven = one - oneMinusXPowerB;
                gEven = b * x * oneMinusXPowerB;
                result = pWeight * xOdd + qWeight * xEven;
                iteration = 1;
                do
                {
                    a = a + one;
                    xOdd = xOdd - gOdd;
                    xEven = xEven - gEven;
                    gOdd = gOdd * x * (a + b - one) / a;
                    gEven = gEven * x * (a + b - half) / (a + half);
                    pWeight = pWeight * lambda / (two * iteration);
                    qWeight = qWeight * lambda / (two * iteration + one);
                    remainingWeight = remainingWeight - pWeight;
                    iteration = iteration + 1;
                    result = result + pWeight * xOdd + qWeight * xEven;
                    errorBound = two * remainingWeight * (xOdd - gOdd);
                }
                while (errorBound > maxError && iteration <= maxIterations);
            }
            if (iteration > maxIterations)
            {
                throw new ArgumentException("Max number of iterations were exceeded.");
            }

            if (reflected)
            {
                result = Normal.StandardCDF(adjustedDelta) - result;
            }
            else
            {
                result = result + (1.0d - Normal.StandardCDF(adjustedDelta));
            }
            return result;
        }

        /// <summary>
        /// The inverse of the non-central t distribution
        /// </summary>
        /// <param name="p">Probability between 0 and 1.</param>
        /// <param name="df">The degrees of freedom.</param>
        /// <param name="delta">The noncentrality parameter.</param>
        /// <returns>The quantile associated with <paramref name="p"/>.</returns>
        private static double InverseCdf(double p, double df, double delta)
        {
            double lower;
            double upper;
            double candidate;
            double lowerError;
            double upperError;
            double increment;
            double slope;
            int iteration;
            const double probabilityTolerance = 0.0000001d;
            const double quantileTolerance = 0.0000001d;
            const int maxIterations = 50;
            lower = InitialInverseCdfEstimate(p, df, delta);
            lowerError = EvaluateCdfSeries(lower, df, delta) - p;
            if (lowerError > probabilityTolerance)
            {
                increment = -1.0d;
            }
            else if (lowerError < -probabilityTolerance)
            {
                increment = 1.0d;
            }
            else
            {
                return lower;
            }

            // Find a bracket by secant extrapolation with a factor-of-two overshoot.
            upper = lower + increment;
            upperError = EvaluateCdfSeries(upper, df, delta) - p;
            iteration = 0;
            while ((lowerError < 0d) != (upperError > 0d) && Math.Abs(upper - lower) > quantileTolerance && iteration < maxIterations)
            {
                slope = (upperError - lowerError) / (upper - lower);
                if (slope == 0d)
                {
                    return (upper + lower) / 2.0d;
                }

                candidate = lower - 2.0d * lowerError / slope;
                lower = upper;
                upper = candidate;
                lowerError = upperError;
                upperError = EvaluateCdfSeries(upper, df, delta) - p;
                iteration = iteration + 1;
            }
            return Brent.Solve(x => EvaluateCdfSeries(x, df, delta) - p, Math.Min(lower, upper), Math.Max(lower, upper), quantileTolerance, reportFailure: false);
        }

        /// <summary>
        /// Computes the Johnson and Kotz starting estimate for noncentral-t inversion.
        /// </summary>
        /// <param name="probability">The cumulative probability.</param>
        /// <param name="degreesOfFreedom">The degrees of freedom.</param>
        /// <param name="noncentrality">The noncentrality parameter.</param>
        /// <returns>An initial quantile estimate for the root search.</returns>
        /// <remarks>Uses the Jennett and Welch approximation and falls back to an offset central-t quantile when necessary.</remarks>
        private static double InitialInverseCdfEstimate(double probability, double degreesOfFreedom, double noncentrality)
        {
            // Approximates percentage points of the non-central t distribution.
            // Source: Johnson & Kotz, Continuous Univariate Distributions, Volume 2.
            double z = Normal.StandardZ(probability);

            // Jennett & Welch approximation, formula (14.1).
            // Intended for the large noncentralities common in tolerance-interval calculations.
            double gammaRatio = Math.Exp(Gamma.LogGamma((degreesOfFreedom + 1d) / 2d) - Gamma.LogGamma(degreesOfFreedom / 2d)) * Math.Sqrt(2d / degreesOfFreedom);
            double zSquared = z * z;
            double gammaRatioSquared = gammaRatio * gammaRatio;
            double denominator = gammaRatioSquared - zSquared * (1d - gammaRatioSquared);
            double discriminant = gammaRatioSquared + (1d - gammaRatioSquared) * (noncentrality * noncentrality - zSquared);
            if (discriminant > 0d && Math.Abs(denominator) > 1e-12)
            {
                return (noncentrality * gammaRatio + z * Math.Sqrt(discriminant)) / denominator;
            }

            // Fall back to an offset central-t quantile.
            var studentT = new StudentT(degreesOfFreedom);
            return studentT.InverseCDF(probability) + noncentrality;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new NoncentralT(DegreesOfFreedom, Noncentrality);
        }

    }
}