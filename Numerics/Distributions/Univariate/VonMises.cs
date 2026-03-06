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
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

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

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Mu; }
        }

        /// <inheritdoc/>
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
        /// Gets the circular standard deviation, defined as √(1 - I₁(κ)/I₀(κ)).
        /// </summary>
        /// <remarks>
        /// The circular variance is defined as V = 1 - A(κ), where A(κ) = I₁(κ)/I₀(κ)
        /// is the mean resultant length. The standard deviation is √V.
        /// For κ = 0 (uniform), V = 1. As κ → ∞, V → 0.
        /// Note: Variance is computed as StandardDeviation² by the base class.
        /// </remarks>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(1d - Bessel.I1(_kappa) / Bessel.I0(_kappa)); }
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
                // Circular kurtosis: κ₂ = A₂(κ) / V² where A₂ = I₂/I₀
                // For now, return a numerical approximation
                double a1 = Bessel.I1(_kappa) / Bessel.I0(_kappa);
                double v = 1d - a1;
                if (v <= 0) return double.NaN;
                // The circular kurtosis is not directly comparable to linear kurtosis,
                // so we return the excess kurtosis of the wrapped distribution
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

        /// <inheritdoc/>
        public double[] MLE(IList<double> sample)
        {
            // Compute mean direction
            double sumSin = 0d, sumCos = 0d;
            for (int i = 0; i < sample.Count; i++)
            {
                sumSin += Math.Sin(sample[i]);
                sumCos += Math.Cos(sample[i]);
            }
            double mu = Math.Atan2(sumSin, sumCos);

            // Compute mean resultant length R_bar
            double rBar = Math.Sqrt(sumSin * sumSin + sumCos * sumCos) / sample.Count;

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

            // Refine with Newton-Raphson iterations: A(kappa) = I1/I0, A'(kappa) = 1 - A(kappa)^2 - A(kappa)/kappa
            if (rBar > 0 && rBar < 1)
            {
                for (int i = 0; i < 20; i++)
                {
                    double a = Bessel.I1(kappa) / Bessel.I0(kappa);
                    double aPrime = 1d - a * a - a / kappa;
                    double delta = (a - rBar) / aPrime;
                    kappa -= delta;
                    if (kappa < 0) kappa = Tools.DoubleMachineEpsilon;
                    if (Math.Abs(delta) < 1e-12) break;
                }
            }
            else if (rBar >= 1)
            {
                kappa = double.MaxValue;
            }
            else
            {
                kappa = 0d;
            }

            return [mu, kappa];
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters([Mu, Kappa], true);
            if (x < -Math.PI || x > Math.PI)
                return 0d;
            return Math.Exp(_kappa * Math.Cos(x - _mu)) / (2d * Math.PI * Bessel.I0(_kappa));
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters([Mu, Kappa], true);
            if (x <= -Math.PI) return 0d;
            if (x >= Math.PI) return 1d;

            // Numerical integration from -π to x
            var integrator = new AdaptiveGaussKronrod((t) => Math.Exp(_kappa * Math.Cos(t - _mu)), -Math.PI, x);
            integrator.Integrate();
            return integrator.Result / (2d * Math.PI * Bessel.I0(_kappa));
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

            var rng = seed < 0 ? new Random() : new Random(seed);
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
