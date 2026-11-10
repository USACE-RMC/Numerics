using System;
using Numerics.Mathematics.RootFinding;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// Declares common functionality of all Archimedean Copulas.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class ArchimedeanCopula : BivariateCopula, IArchimedeanCopula
    {

        /// <inheritdoc/>
        public override string[,] ParameterToString
        {
            get
            {
              var parmString = new string[2, 2];
              parmString[0, 0] = "Dependency (θ)";
              parmString[0, 1] = Theta.ToString();
              return parmString;
            }
        }

        /// <inheritdoc/>
        public override string ParameterNameShortForm
        {
            get { return "θ"; }
        }

        /// <inheritdoc/>
        public override int NumberOfCopulaParameters => 1;

        /// <inheritdoc/>
        public override double[] GetCopulaParameters => new double[] { Theta };

        /// <inheritdoc/>
        public override void SetCopulaParameters(double[] parameters)
        {
            Theta = parameters[0];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameter(double parameter, bool throwException)
        {
            if (double.IsNaN(parameter) || double.IsInfinity(parameter))
            {
                var exception = new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter must be finite.");
                if (throwException) throw exception;
                return exception;
            }
            if (parameter < ThetaMinimum)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter θ (theta) must be greater than or equal to " + ThetaMinimum.ToString() + ".");
                return new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter θ (theta) must be greater than or equal to " + ThetaMinimum.ToString() + ".");
            }

            if (parameter > ThetaMaximum)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter θ (theta) must be less than or equal to " + ThetaMaximum.ToString() + ".");
                return new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter θ (theta) must be less than or equal to " + ThetaMaximum.ToString() + ".");
            }
            // A valid parameter must return null: the Theta setter derives ParametersValid
            // from 'ValidateParameter(value, false) is null'.
            return null;
        }

        /// <summary>
        /// Solves g(z) = 0 for an increasing function g with a safeguarded Newton iteration.
        /// </summary>
        /// <param name="g">
        /// Returns the increasing function and its derivative, which is positive, at a point; the two are evaluated
        /// together so they can share their intermediate terms.
        /// </param>
        /// <param name="start">A starting point near the root.</param>
        /// <returns>The root.</returns>
        /// <remarks>
        /// Each evaluation narrows a bracket [lower, upper] by the sign of g. A Newton step that leaves the bracket is
        /// replaced by bisection once both ends are known, and by a step of 8 toward the root while one end is still
        /// open. The iteration stops when the step, or the bracket, is below 2 ulps of the iterate; a Newton step that
        /// small is accepted before the bracket test, so a step that rounds onto the end of the bracket at the root does
        /// not trigger the safeguard. No evaluations are spent
        /// on finding a bracket first, and from a start near the root Newton converges quadratically in a few steps. The
        /// generator-slope inverses of the Gumbel and Joe copulas use it in a log variable.
        /// </remarks>
        protected static double SolveIncreasing(Func<double, (double Value, double Slope)> g, double start)
        {
            double lower = double.NegativeInfinity, upper = double.PositiveInfinity;
            double z = start;
            for (int i = 0; i < 200; i++)
            {
                var (value, slope) = g(z);
                if (value == 0d) return z;
                if (value < 0d) lower = z; else upper = z;
                double next = z - value / slope;
                // A converged step ends the iteration even when rounding leaves it on the bracket end z; otherwise
                // the bracket test below would replace it with a step of 8 away from the root.
                if (Math.Abs(next - z) <= 4.4E-16 * Math.Max(1d, Math.Abs(next))) return next;
                if (!(next > lower && next < upper))
                {
                    if (double.IsInfinity(lower)) next = Math.Max(z - 8d, upper - 8d);
                    else if (double.IsInfinity(upper)) next = Math.Min(z + 8d, lower + 8d);
                    else next = 0.5d * (lower + upper);
                }
                double tolerance = 4.4E-16 * Math.Max(1d, Math.Abs(next));
                if (Math.Abs(next - z) <= tolerance || upper - lower <= tolerance) return next;
                z = next;
            }
            return z;
        }

        /// <summary>
        /// Returns ln h(v|u), the log of <see cref="BivariateCopula.ConditionalCDF(double, double)"/>.
        /// </summary>
        /// <param name="u">The conditioning variate.</param>
        /// <param name="v">The dependent variate.</param>
        /// <returns>The log of the conditional CDF.</returns>
        /// <remarks>
        /// The base implementation takes the log of <see cref="BivariateCopula.ConditionalCDF(double, double)"/>. Families
        /// with a closed form override it with the exponent itself, so 1 − h = −Expm1(ln h) keeps its relative accuracy
        /// where h is close to 1; <see cref="InvertConditionalCDF(double, double, double)"/> relies on it.
        /// </remarks>
        protected virtual double LogConditionalCDF(double u, double v)
        {
            return Math.Log(ConditionalCDF(u, v));
        }

        /// <summary>
        /// Returns ln h(v|u) and the density c(u,v) together, for the Newton iteration of
        /// <see cref="InvertConditionalCDF(double, double, double)"/>.
        /// </summary>
        /// <param name="u">The conditioning variate.</param>
        /// <param name="v">The dependent variate.</param>
        /// <param name="density">Receives the copula density at (u, v).</param>
        /// <returns>The log of the conditional CDF.</returns>
        /// <remarks>
        /// The base implementation calls <see cref="LogConditionalCDF(double, double)"/> and
        /// <see cref="BivariateCopula.PDF(double, double)"/>; families whose closed forms share terms override it to
        /// evaluate them once.
        /// </remarks>
        protected virtual double LogConditionalAndDensity(double u, double v, out double density)
        {
            density = PDF(u, v);
            return LogConditionalCDF(u, v);
        }

        /// <summary>
        /// Solves h(v|u) = t for v in (0, 1) with a safeguarded Newton iteration on the conditional CDF, whose derivative in
        /// v is the copula density.
        /// </summary>
        /// <param name="u">The conditioning variate, strictly inside (0, 1).</param>
        /// <param name="t">The conditional probability, strictly inside (0, 1).</param>
        /// <param name="start">The starting value of v, inside (0, 1).</param>
        /// <returns>The v with <see cref="BivariateCopula.ConditionalCDF(double, double)"/>(u, v) = t.</returns>
        /// <remarks>
        /// <para>
        /// Each evaluation of h narrows the bracket [lower, upper], which starts at [0, 1]. A Newton step v − (h(v) − t)/c(u,v)
        /// that leaves the bracket, or a density that is zero or not finite, is replaced by bisection. For t > 1/2 the
        /// residual is (1 − t) − (1 − h) with 1 − h from <see cref="LogConditionalCDF(double, double)"/>, so the root keeps
        /// its accuracy where t is close to 1. The bisection is geometric while the
        /// bracket spans more than a factor 4 above zero, so small roots are reached in few steps, and arithmetic
        /// otherwise. The iteration stops after a Newton step below 1E-9 of min(v, 1 − v), whose error is of the order of
        /// the step squared, or when the bracket is below one ulp of v; after 100 iterations
        /// the bracket, with the same residual, is handed to <see cref="Brent.Solve(Func{double, double}, double, double, double, double, int, bool)"/>
        /// at full precision. The family's ConditionalCDF and PDF must be accurate to the last bits for the inverse to be.
        /// </para>
        /// </remarks>
        protected double InvertConditionalCDF(double u, double t, double start)
        {
            double lower = 0d, upper = 1d, v = start;
            // Above t = 1/2 the residual is formed from the complements, (1 − t) − (1 − h), with 1 − h = −Expm1(ln h), so
            // a t close to 1 keeps its relative distance from 1.
            bool upperTail = t > 0.5d;
            double complement = 1d - t;
            for (int i = 0; i < 100; i++)
            {
                double logH = LogConditionalAndDensity(u, v, out double density);
                double residual = upperTail ? complement + Tools.Expm1(logH) : Math.Exp(logH) - t;
                if (residual == 0d) return v;
                if (double.IsNaN(residual)) return double.NaN;
                if (residual < 0d) lower = v; else upper = v;
                double next = v - residual / density;
                if (next > lower && next < upper)
                {
                    // Newton converges quadratically here: after a step below 1E-9 of the distance to the nearer end of
                    // (0, 1) the error is of the order of the step squared, far below an ulp.
                    if (Math.Abs(next - v) <= 1E-9 * Math.Min(next, 1d - next)) return next;
                }
                else
                {
                    next = lower > 0d && upper > 4d * lower ? Math.Sqrt(lower * upper) : lower == 0d ? 0.25d * upper : 0.5d * (lower + upper);
                }
                if (upper - lower <= 2.2E-16 * upper) return next;
                v = next;
            }
            Func<double, double> residualAt = upperTail
                ? x => complement + Tools.Expm1(LogConditionalCDF(u, x))
                : x => Math.Exp(LogConditionalCDF(u, x)) - t;
            return Brent.Solve(residualAt, lower, upper, 0d, Brent.MinimumRelativeTolerance, 1000, false);
        }

        /// <inheritdoc/>
        public abstract double Generator(double t);

        /// <inheritdoc/>
        public abstract double GeneratorInverse(double t);

        /// <inheritdoc/>
        public abstract double GeneratorPrime(double t);

        /// <inheritdoc/>
        public abstract double GeneratorPrime2(double t);

        /// <inheritdoc/>
        public abstract double GeneratorPrimeInverse(double t);

        /// <inheritdoc/>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double num = -GeneratorPrime2(CDF(u, v)) * GeneratorPrime(u) * GeneratorPrime(v);
            double den = Math.Pow(GeneratorPrime(CDF(u, v)), 3d);
            return num / den;
        }

        /// <inheritdoc/>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            return GeneratorInverse(Generator(u) + Generator(v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// For an Archimedean copula with generator φ, differentiating C(u,v) = φ⁻¹(φ(u) + φ(v))
        /// with respect to u gives the exact generic conditional h(v|u) = φ′(u) / φ′(C(u,v)).
        /// The ratio is insensitive to a family's internal sign convention for φ′.
        /// </para>
        /// <para>
        /// Per family the ratio reduces to the closed forms:
        /// Clayton h = u^(−θ−1)·(u^(−θ) + v^(−θ) − 1)^(−1−1/θ);
        /// Frank h = e^(−θu)(e^(−θv) − 1) / [(e^(−θ) − 1) + (e^(−θu) − 1)(e^(−θv) − 1)];
        /// Gumbel h = C(u,v)·A^(1/θ−1)·(−ln u)^(θ−1)/u with A = (−ln u)^θ + (−ln v)^θ;
        /// Joe h = (1−u)^(θ−1)·[1 − (1−v)^θ]·A^(1/θ−1) with A = (1−u)^θ + (1−v)^θ − (1−u)^θ(1−v)^θ;
        /// Ali-Mikhail-Haq h = v(1 − θ(1−v))/D² with D = 1 − θ(1−u)(1−v).
        /// </para>
        /// <para>
        /// The generic generator ratio has conditioning domain 0 &lt; u &lt; 1. Within that domain,
        /// the dependent-variable boundaries are exact: h(0|u) = 0 and h(1|u) = 1. Conditioning
        /// endpoint limits are family-specific, so this base implementation does not impose a
        /// family-independent value at u = 0 or u = 1.
        /// </para>
        /// </remarks>
        public override double ConditionalCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (v == 0d) return 0d;
            if (v == 1d) return 1d;
            return GeneratorPrime(u) / GeneratorPrime(CDF(u, v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// This method is based on Genest et al. 1986
        /// 1) Two independent uniformly distributed U(0,1) random variates, u and t, are generated.
        /// 2) Two new variables, s and w, are obtained as s = GeneratorPrime(u) / t and w = GeneratorPrimeInverse(s).
        /// 3) The dependent variate is obtained as v = GeneratorInverse(Generator(w) - Generator(u))
        /// 4) The pairs u and v are the simulated pair, preserving the dependence structure.
        /// 5) Both these u and v are in the range [0,1]. Simulated pairs of u and v are then
        /// back-transformed through their corresponding marginal distributions.
        /// This is the exact inverse of the generic conditional h(v|u) = φ′(u)/φ′(C(u,v)):
        /// setting h = t and solving gives C = φ′⁻¹(φ′(u)/t) and v = φ⁻¹(φ(C) − φ(u)).
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double s = GeneratorPrime(u) / t;
            double w = GeneratorPrimeInverse(s);
            return GeneratorInverse(Generator(w) - Generator(u));
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            return [u, InverseConditionalCDF(u, v)];
        }

    }
}
