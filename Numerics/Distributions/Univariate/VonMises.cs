using System;
using System.Collections.Generic;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// The von Mises distribution for circular data.
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
    /// The von Mises distribution (also known as the circular normal distribution) is a continuous
    /// probability distribution on the circle. It is the circular analogue of the normal distribution
    /// and is commonly used to model directional data such as flood seasonality (timing of annual
    /// maximum flows), wind directions, and other angular measurements.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Mardia, K.V. and Jupp, P.E. (2000). "Directional Statistics." John Wiley and Sons.
    /// </description></item>
    /// <item><description>
    /// Best, D.J. and Fisher, N.I. (1979). "Efficient simulation of the von Mises distribution."
    /// Applied Statistics, 28(2), 152-157.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Von_Mises_distribution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class VonMises : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {

        /// <summary>
        /// Constructs a von Mises distribution with μ = 0 and κ = 1.
        /// </summary>
        public VonMises()
        {
            SetParameters(0d, 1d);
        }

        /// <summary>
        /// Constructs a von Mises distribution with given μ and κ.
        /// </summary>
        /// <param name="mu">The mean direction μ (mu), in radians. Must be in [-π, π].</param>
        /// <param name="kappa">The concentration parameter κ (kappa). Must be ≥ 0.</param>
        public VonMises(double mu, double kappa)
        {
            SetParameters(mu, kappa);
        }

        private double _mu;    // mean direction
        private double _kappa; // concentration

        /// <summary>
        /// Gets and sets the mean direction parameter μ (mu), in radians.
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters([value, Kappa], false) is null;
                _mu = value;
            }
        }

        /// <summary>
        /// Gets and sets the concentration parameter κ (kappa).
        /// </summary>
        public double Kappa
        {
            get { return _kappa; }
            set
            {
                _parametersValid = ValidateParameters([Mu, value], false) is null;
                _kappa = value;
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
            get { return UnivariateDistributionType.VonMises; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Von Mises"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "VM"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Mean Direction (μ)";
                parmString[1, 0] = "Concentration (κ)";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Kappa.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["μ", "κ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mu), nameof(Kappa)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mu, Kappa]; }
        }

        /// <summary>
        /// Gets the circular mean direction μ (Mu).
        /// </summary>
        public override double Mean
        {
            get { return Mu; }
        }

        /// <summary>
        /// Gets the circular median, which equals the mean direction μ (Mu) for this symmetric unimodal distribution.
        /// </summary>
        public override double Median
        {
            get { return Mu; }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Mu; }
        }

        /// <summary>
        /// Gets the square root of the circular variance, √(1 - A(κ)), where A(κ) = I₁(κ)/I₀(κ).
        /// </summary>
        /// <remarks>
        /// The circular variance is V = 1 - R̄, where R̄ = A(κ) = I₁(κ)/I₀(κ) is the mean resultant length, and this
        /// property returns √V. It differs from the circular standard deviation √(-2 ln R̄) of Mardia and Jupp.
        /// For κ = 0 (uniform), V = 1. As κ → ∞, V → 0.
        /// Variance is computed as StandardDeviation² by the base class, so it equals V.
        /// A(κ) is evaluated as <see cref="Bessel.I1e(double)"/>/<see cref="Bessel.I0e(double)"/>, which equals
        /// I₁(κ)/I₀(κ) and stays finite for every κ.
        /// </remarks>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(1d - Bessel.I1e(_kappa) / Bessel.I0e(_kappa)); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 0d; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                // Circular kurtosis is not directly comparable to linear kurtosis
                return double.NaN;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return -Math.PI; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Math.PI; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [-Math.PI, 0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [Math.PI, double.PositiveInfinity]; }
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
            var newDistribution = new VonMises(Mu, Kappa);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="mu">The mean direction μ (mu), in radians.</param>
        /// <param name="kappa">The concentration parameter κ (kappa).</param>
        public void SetParameters(double mu, double kappa)
        {
            _parametersValid = ValidateParameters(new[] { mu, kappa }, false) is null;
            _mu = mu;
            _kappa = kappa;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "The mean direction parameter μ (mu) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Mu), "The mean direction parameter μ (mu) must be a number.");
            }
            if (parameters[0] < -Math.PI || parameters[0] > Math.PI)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "The mean direction parameter μ (mu) must be in [-π, π].");
                return new ArgumentOutOfRangeException(nameof(Mu), "The mean direction parameter μ (mu) must be in [-π, π].");
            }
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] < 0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Kappa), "The concentration parameter κ (kappa) must be non-negative.");
                return new ArgumentOutOfRangeException(nameof(Kappa), "The concentration parameter κ (kappa) must be non-negative.");
            }
            return null;
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = MLE(sample);
            var lowerVals = new double[] { -Math.PI, 0d };
            var upperVals = new double[] { Math.PI, Math.Max(initialVals[1] * 10d, 100d) };
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>
        /// Estimates the parameters by maximum likelihood.
        /// </summary>
        /// <param name="sample">The sample of angles, in radians.</param>
        /// <returns>The estimated parameters [μ, κ].</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than two observations or contains a NaN or
        /// infinite observation.</exception>
        /// <exception cref="ArgumentException">The sample has no angular dispersion: every observation is the same, or
        /// the mean resultant length R̄ is at least one, so the likelihood increases without bound in κ.</exception>
        /// <remarks>
        /// μ is the circular mean atan2(Σ sin xᵢ, Σ cos xᵢ). κ solves A(κ) = R̄, where R̄ is the mean resultant length,
        /// starting from the Mardia and Jupp (2000) approximation and refined by Newton-Raphson iterations. A sample
        /// with R̄ = 0 gives κ = 0.
        /// </remarks>
        public double[] MLE(IList<double> sample)
        {
            if (sample == null)
                throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2)
                throw new ArgumentOutOfRangeException(nameof(sample), "At least two observations are required to estimate the von Mises distribution.");

            // Compute mean direction
            double sumSin = 0d, sumCos = 0d;
            bool distinct = false;
            for (int i = 0; i < sample.Count; i++)
            {
                if (!Tools.IsFinite(sample[i]))
                    throw new ArgumentOutOfRangeException(nameof(sample), "Observations must be finite.");
                distinct |= sample[i] != sample[0];
                sumSin += Math.Sin(sample[i]);
                sumCos += Math.Cos(sample[i]);
            }
            double mu = Math.Atan2(sumSin, sumCos);

            // Compute mean resultant length R_bar
            double rBar = Math.Sqrt(sumSin * sumSin + sumCos * sumCos) / sample.Count;
            if (!distinct || rBar >= 1d)
                throw new ArgumentException("The sample has no angular dispersion (mean resultant length of one), so the maximum likelihood concentration κ (kappa) is unbounded.", nameof(sample));

            // Solve A(kappa) = R_bar for kappa, where A(kappa) = I1(kappa)/I0(kappa)
            // Use the approximation from Mardia & Jupp (2000) as initial estimate
            double kappa;
            if (rBar < 0.53)
            {
                kappa = 2d * rBar + rBar * rBar * rBar + 5d / 6d * Math.Pow(rBar, 5);
            }
            else if (rBar < 0.85)
            {
                kappa = -0.4 + 1.39 * rBar + 0.43 / (1d - rBar);
            }
            else
            {
                kappa = 1d / (rBar * rBar * rBar - 4d * rBar * rBar + 3d * rBar);
            }

            // Refine with Newton-Raphson iterations: A(kappa) = I1/I0, A'(kappa) = 1 - A(kappa)^2 - A(kappa)/kappa.
            // The scaled functions give the same ratio without overflow for large kappa.
            if (rBar > 0 && rBar < 1)
            {
                for (int i = 0; i < 20; i++)
                {
                    double a = Bessel.I1e(kappa) / Bessel.I0e(kappa);
                    double aPrime = 1d - a * a - a / kappa;
                    if (Math.Abs(aPrime) < 1e-30) break;
                    double delta = (a - rBar) / aPrime;
                    kappa -= delta;
                    if (kappa < 0) kappa = Tools.DoubleMachineEpsilon;
                    if (Math.Abs(delta) < 1e-12) break;
                }
            }
            else
            {
                kappa = 0d;
            }

            return [mu, kappa];
        }

        /// <summary>
        /// Gets the probability density function (PDF) of the distribution evaluated at a point x.
        /// </summary>
        /// <param name="x">A single point in the distribution range, in radians.</param>
        /// <returns>The density at x; zero outside [-π, π].</returns>
        /// <remarks>
        /// Evaluated as exp(-2κ sin²((x - μ)/2)) / (2π I₀e(κ)), where I₀e(κ) = exp(-κ) I₀(κ) is
        /// <see cref="Bessel.I0e(double)"/>. Since cos(x - μ) - 1 = -2 sin²((x - μ)/2), this equals
        /// exp(κ cos(x - μ)) / (2π I₀(κ)), but neither factor overflows and the exponent does not cancel, so the
        /// density is finite and accurate for every κ.
        /// </remarks>
        public override double PDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters([Mu, Kappa], true);
            if (x < -Math.PI || x > Math.PI)
                return 0d;
            double halfAngle = Math.Sin(0.5d * (x - _mu));
            return Math.Exp(-2d * _kappa * halfAngle * halfAngle) / (2d * Math.PI * Bessel.I0e(_kappa));
        }

        /// <summary>
        /// Gets the cumulative distribution function (CDF) of the distribution evaluated at a point x.
        /// </summary>
        /// <param name="x">A single point in the distribution range, in radians.</param>
        /// <returns>The probability that X ≤ x, measured from -π.</returns>
        /// <remarks>
        /// <para>
        /// The kernel exp(c - 2κ sin²((t - μ)/2)), with c = min(κ, 700), is integrated from -π to x by adaptive
        /// Gauss-Kronrod quadrature and divided by 2π I₀e(κ) exp(c), where I₀e(κ) = exp(-κ) I₀(κ) is
        /// <see cref="Bessel.I0e(double)"/>. Since cos(t - μ) - 1 = -2 sin²((t - μ)/2), the kernel equals
        /// exp(κ cos(t - μ)) exp(c - κ): it is the unnormalized density itself for κ ≤ 700, and above 700 its peak is
        /// capped at exp(700), below the double overflow threshold, so the CDF is finite for every κ. The half-angle
        /// form keeps the relative rounding error of the kernel near the peak at about 1E-13 instead of κ times the
        /// rounding error of cos, so the quadrature's error estimate converges for large κ. The integral is split at
        /// μ - w, μ and μ + w where they fall inside (-π, x), with w = 8/√max(κ, 1) (eight standard deviations of the
        /// normal limit), so the quadrature resolves the peak even when κ is large and the peak is narrow.
        /// </para>
        /// <para>
        /// The kernel stays large compared with the integrator's absolute tolerance wherever the density has
        /// appreciable mass, so the relative tolerance (1E-8) of each subinterval governs the accuracy. Against mpmath
        /// quadrature of the density the CDF agrees to within 2E-14 absolute for κ from 0.01 to 1E6; probabilities
        /// below about 1E-40 in the far tails carry relative errors up to about 1E-2.
        /// </para>
        /// </remarks>
        public override double CDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters([Mu, Kappa], true);
            if (x <= -Math.PI) return 0d;
            if (x >= Math.PI) return 1d;

            // Numerical integration from -π to x, split at μ - w, μ and μ + w, where w is eight
            // standard deviations of the normal limit, so narrow peaks are resolved
            double c = Math.Min(_kappa, 700d);
            Func<double, double> kernel = (t) =>
            {
                double s = Math.Sin(0.5 * (t - _mu));
                return Math.Exp(c - 2d * _kappa * s * s);
            };
            double w = 8d / Math.Sqrt(Math.Max(_kappa, 1d));
            double lower = -Math.PI;
            double integral = 0d;
            foreach (double breakpoint in new[] { _mu - w, _mu, _mu + w, x })
            {
                if (breakpoint <= lower || breakpoint > x) continue;
                var integrator = new AdaptiveGaussKronrod(kernel, lower, breakpoint);
                integrator.Integrate();
                integral += integrator.Result;
                lower = breakpoint;
            }
            return integral / (2d * Math.PI * Bessel.I0e(_kappa) * Math.Exp(c));
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            if (probability < 0d || probability > 1d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0d) return Minimum;
            if (probability == 1d) return Maximum;
            if (_parametersValid == false)
                ValidateParameters([Mu, Kappa], true);

            // Use Brent's method to solve CDF(x) = probability
            return Brent.Solve((x) => CDF(x) - probability, -Math.PI, Math.PI);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new VonMises(Mu, Kappa);
        }

        /// <summary>
        /// Generates random values from the von Mises distribution using Best's algorithm.
        /// </summary>
        /// <remarks>
        /// Best, D.J. and Fisher, N.I. (1979). "Efficient simulation of the von Mises distribution."
        /// </remarks>
        public override double[] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            if (_parametersValid == false)
                ValidateParameters([Mu, Kappa], true);

            var rng = seed < 0 ? new MersenneTwister() : new MersenneTwister(seed);
            var values = new double[sampleSize];

            if (_kappa < 1e-10)
            {
                // For kappa ≈ 0, the distribution is approximately uniform on [-π, π]
                for (int i = 0; i < sampleSize; i++)
                    values[i] = -Math.PI + 2d * Math.PI * rng.NextDouble();
                return values;
            }

            // Best's algorithm
            double tau = 1d + Math.Sqrt(1d + 4d * _kappa * _kappa);
            double rho = (tau - Math.Sqrt(2d * tau)) / (2d * _kappa);
            double r = (1d + rho * rho) / (2d * rho);

            for (int i = 0; i < sampleSize; i++)
            {
                double f, c;
                while (true)
                {
                    double u1 = rng.NextDouble();
                    double z = Math.Cos(Math.PI * u1);
                    f = (1d + r * z) / (r + z);
                    c = _kappa * (r - f);

                    double u2 = rng.NextDouble();
                    if (c * (2d - c) > u2 || Math.Log(c / u2) + 1d >= c)
                        break;
                }

                double u3 = rng.NextDouble();
                double theta = (u3 > 0.5 ? 1d : -1d) * Math.Acos(f) + _mu;

                // Wrap to [-π, π]
                theta = ((theta + Math.PI) % (2d * Math.PI) + 2d * Math.PI) % (2d * Math.PI) - Math.PI;
                values[i] = theta;
            }

            return values;
        }

    }
}
