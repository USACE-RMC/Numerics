using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{
    /// <summary>Expected fixed-observation score information for regular Kappa Four and fixed-hondo families.</summary>
    /// <remarks>
    /// Scores are integrated over both complete probability tails. The square-root Jacobian is
    /// combined with exponential scores before products are formed. Checked quadrature and
    /// diagonally scaled Cholesky inversion report numerical failure without regularization.
    /// </remarks>
    internal static class KappaExpectedInformation
    {
        private const double RelativeTolerance = 1E-10;
        private const double AbsoluteTolerance = 1E-12;
        private const double Roundoff = 2.2204460492503131E-16;
        private const int MaximumIntervals = 2048;
        private static readonly double[] Nodes =
        {
            .995657163025808080735527280689003, .973906528517171720077964012084452,
            .930157491355708226001207180059508, .865063366688984510732096688423493,
            .780817726586416897063717578345042, .679409568299024406234327365114874,
            .562757134668604683339000099272694, .433395394129247190799265943165784,
            .294392862701460198131126603103866, .148874338981631210884826001129720, 0
        };
        private static readonly double[] KronrodWeights =
        {
            .011694638867371874278064396062192, .032558162307964727478818972459390,
            .054755896574351996031381300244580, .075039674810919952767043140916190,
            .093125454583697605535065465083366, .109387158802297641899210590325805,
            .123491976262065851077958109831074, .134709217311473325928054001771707,
            .142775938577060080797094273138717, .147739104901338491374841515972068,
            .149445554002916905664936468389821
        };
        private static readonly double[] GaussWeights =
        {
            .066671344308688137593568809893332, .149451349150580593145776339657697,
            .219086362515982043995534934228163, .269266719309996355091226921569469,
            .295524224714752870173892994651338
        };

        /// <summary>Returns local asymptotic MLE covariance in xi, alpha, kappa[, hondo] coordinates.</summary>
        /// <param name="alpha">The finite positive scale.</param>
        /// <param name="k">The finite kappa shape.</param>
        /// <param name="h">The finite hondo shape, fixed when three parameters are requested.</param>
        /// <param name="sampleSize">The positive number of independent observations.</param>
        /// <param name="parameterCount">Three for a fixed-hondo family, or four for full Kappa Four.</param>
        /// <returns>The inverse expected information, transformed to public scale and divided by sample size.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Scale, sample size, dimension or information regularity is invalid.</exception>
        /// <exception cref="InvalidOperationException">Quadrature, positive definiteness, inversion or representability checks fail.</exception>
        internal static double[,] ParameterCovariance(double alpha, double k, double h, int sampleSize, int parameterCount)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (!DistributionNumerics.IsFinite(alpha) || alpha <= 0) throw new ArgumentOutOfRangeException(nameof(alpha));
            double[,] information = ExpectedInformation(k, h, parameterCount, out _, out _, out _);
            double[,] covariance = InvertInformation(information);
            for (int i = 0; i < parameterCount; i++)
            for (int j = i; j < parameterCount; j++)
            {
                double value = covariance[i, j];
                if (value != 0)
                {
                    double logarithm = Math.Log(Math.Abs(value)) - Math.Log(sampleSize)
                        + (i < 2 ? Math.Log(alpha) : 0) + (j < 2 ? Math.Log(alpha) : 0);
                    value = Math.Sign(value) * Math.Exp(logarithm);
                    if (!DistributionNumerics.IsFinite(value))
                        throw new InvalidOperationException("The local MLE covariance is outside the finite floating-point range.");
                }
                covariance[i, j] = covariance[j, i] = value;
            }
            return covariance;
        }

        /// <summary>Returns standardized information and complete-tail integration diagnostics.</summary>
        /// <param name="k">The kappa shape.</param>
        /// <param name="h">The hondo shape.</param>
        /// <param name="parameterCount">The estimated parameter dimension.</param>
        /// <param name="scoreMeans">Integrated fixed-observation score means.</param>
        /// <param name="scoreMeanErrors">Estimated absolute errors of the score means.</param>
        /// <param name="informationErrors">Estimated absolute errors of information entries.</param>
        /// <returns>The symmetric per-observation standardized information matrix.</returns>
        internal static double[,] ExpectedInformation(double k, double h, int parameterCount,
            out double[] scoreMeans, out double[] scoreMeanErrors, out double[,] informationErrors)
        {
            ValidateDomain(k, h, parameterCount);
            double boundary = Math.Max(k, Math.Max(h, k * h));
            // The stronger map regularizes integrable powers close to the information boundary.
            // Both maps cover the entire p interval; neither truncates a probability tail.
            int power = boundary > .4 ? 128 : 8;
            var integrals = Integrate(k, h, parameterCount, power);
            scoreMeans = new double[parameterCount];
            scoreMeanErrors = new double[parameterCount];
            var information = new double[parameterCount, parameterCount];
            informationErrors = new double[parameterCount, parameterCount];
            for (int i = 0; i < parameterCount; i++)
            {
                scoreMeans[i] = integrals.Values[i];
                scoreMeanErrors[i] = integrals.Errors[i];
                if (Math.Abs(scoreMeans[i]) > 4 * scoreMeanErrors[i] + 32 * Roundoff)
                    throw new InvalidOperationException("The integrated fixed-observation score does not have zero mean within numerical error.");
            }
            int index = parameterCount;
            for (int i = 0; i < parameterCount; i++)
            for (int j = i; j < parameterCount; j++, index++)
            {
                information[i, j] = information[j, i] = integrals.Values[index];
                informationErrors[i, j] = informationErrors[j, i] = integrals.Errors[index];
            }
            return information;
        }

        /// <summary>Checks information regularity without changing distribution or fitting validity.</summary>
        private static void ValidateDomain(double k, double h, int count)
        {
            if (count != 3 && count != 4) throw new ArgumentOutOfRangeException(nameof(count));
            if (!DistributionNumerics.IsFinite(k) || !DistributionNumerics.IsFinite(h))
                throw new ArgumentOutOfRangeException(nameof(k), "Information requires finite shape parameters.");
            if (k >= .5 || h >= .5 || k * h >= .5)
                throw new ArgumentOutOfRangeException(nameof(k), "Regular information requires kappa < 1/2, hondo < 1/2 and kappa*hondo < 1/2.");
        }

        /// <summary>Evaluates expm1(a*z)/a and its second divided difference without cancellation in weighted scores.</summary>
        private static double SecondExponentialRelative(double x)
        {
            if (Math.Abs(x) >= .1) return (Tools.Expm1(x) - x) / x / x;
            double sum = .5, term = .5;
            for (int n = 1; n < 24; n++)
            {
                term *= x / (n + 2);
                sum += term;
                if (Math.Abs(term) < Math.Abs(sum) * 1E-17) break;
            }
            return sum;
        }

        /// <summary>Returns log(abs(expm1(x))) without overflowing the exponential.</summary>
        private static double LogAbsoluteExpm1(double x) => x > 36 ? x + Tools.Log1p(-Math.Exp(-x)) : Math.Log(Math.Abs(Tools.Expm1(x)));

        /// <summary>Forms scores times the square-root integration Jacobian before multiplying exponential terms.</summary>
        private static double[] WeightedScores(double logp, double logq, double k, double h, double logroot)
        {
            double u = -logp;
            double logt, logA;
            double logAbsoluteH = Math.Log(Math.Abs(h));
            if (logq < -36 && logAbsoluteH + logq < -36)
            {
                // Corrections are below a binary64 rounding unit, including when q underflows.
                logt = logA = logq;
            }
            else if (h == 0) logt = logA = Math.Log(u);
            else
            {
                // Retain h*u through the log survival coordinate even if u is subnormal.
                double hu = logq < -36 ? Math.Sign(h) * Math.Exp(logAbsoluteH + logq) : h * u;
                logt = LogAbsoluteExpm1(-hu) - logAbsoluteH;
                logA = LogAbsoluteExpm1(hu) - logAbsoluteH;
            }
            double y = -logt;
            double root = Math.Exp(logroot);
            double logOneMinusK = Tools.Log1p(-k), logOneMinusH = Tools.Log1p(-h);
            double cRoot = Math.Exp(logOneMinusK + logroot) - Math.Exp(logOneMinusH + logA + logroot);
            double cExponentialRoot = Math.Exp(logOneMinusK + k * y + logroot)
                - Math.Exp(logOneMinusH + logA + k * y + logroot);
            double scaleScore, kScore;
            if (Math.Abs(k * y) < .1 || k == 0)
            {
                scaleScore = -root + cRoot * y * DistributionNumerics.Exprel(k * y);
                kScore = y * root - cRoot * y * y * SecondExponentialRelative(k * y);
            }
            else
            {
                scaleScore = -root + (cExponentialRoot - cRoot) / k;
                kScore = y * root + ((1 + k * y) * cRoot - cExponentialRoot) / k / k;
            }
            double hScore;
            if (Math.Abs(h * u) < .1 || h == 0)
                hScore = u * root - (1 - h) * u * u * root * SecondExponentialRelative(h * u);
            else hScore = (u * root - Math.Exp(logOneMinusH + logA + logroot)) / h;
            return new[] { cExponentialRoot, scaleScore, kScore, hScore };
        }

        /// <summary>Pairs the full lower and upper p tails in a common finite quadrature coordinate.</summary>
        private static double[] Integrand(double coordinate, double k, double h, int count, int power)
        {
            double logCoordinate = Math.Log(coordinate);
            double logtail = power * logCoordinate - Math.Log(2);
            double logother = DistributionNumerics.Log1mExp(logtail);
            double logroot = .5 * (Math.Log(power / 2d) + (power - 1) * logCoordinate);
            double[] lower = WeightedScores(logtail, logother, k, h, logroot);
            double[] upper = WeightedScores(logother, logtail, k, h, logroot);
            var values = new double[count + count * (count + 1) / 2];
            double root = Math.Exp(logroot);
            for (int i = 0; i < count; i++) values[i] = (lower[i] + upper[i]) * root;
            int index = count;
            for (int i = 0; i < count; i++)
            for (int j = i; j < count; j++, index++) values[index] = lower[i] * lower[j] + upper[i] * upper[j];
            foreach (double value in values)
                if (!DistributionNumerics.IsFinite(value))
                    throw new InvalidOperationException("A complete-tail information integrand could not be represented finitely.");
            return values;
        }

        /// <summary>Stores a subinterval's Kronrod estimate and conservative local error indicators.</summary>
        private sealed class Interval
        {
            internal double Lower, Upper;
            internal double[] Values = Array.Empty<double>();
            internal double[] Errors = Array.Empty<double>();
        }

        /// <summary>G10K21 rule with absolute-deviation rescaling and a floating-point roundoff floor.</summary>
        private static Interval Evaluate(double a, double b, double k, double h, int count, int power)
        {
            double center = .5 * (a + b), half = .5 * (b - a);
            var nodes = new double[21][];
            nodes[20] = Integrand(center, k, h, count, power);
            int dimension = nodes[20].Length;
            var kronrod = new double[dimension];
            var gauss = new double[dimension];
            var absolute = new double[dimension];
            for (int j = 0; j < dimension; j++)
            {
                kronrod[j] = KronrodWeights[10] * nodes[20][j];
                absolute[j] = KronrodWeights[10] * Math.Abs(nodes[20][j]);
            }
            for (int i = 0; i < 10; i++)
            {
                nodes[2 * i] = Integrand(center - half * Nodes[i], k, h, count, power);
                nodes[2 * i + 1] = Integrand(center + half * Nodes[i], k, h, count, power);
                for (int j = 0; j < dimension; j++)
                {
                    double sum = nodes[2 * i][j] + nodes[2 * i + 1][j];
                    kronrod[j] += KronrodWeights[i] * sum;
                    absolute[j] += KronrodWeights[i] * (Math.Abs(nodes[2 * i][j]) + Math.Abs(nodes[2 * i + 1][j]));
                    if (i % 2 == 1) gauss[j] += GaussWeights[i / 2] * sum;
                }
            }
            var errors = new double[dimension];
            for (int j = 0; j < dimension; j++)
            {
                double mean = kronrod[j] / 2;
                double deviation = KronrodWeights[10] * Math.Abs(nodes[20][j] - mean);
                for (int i = 0; i < 10; i++)
                    deviation += KronrodWeights[i] * (Math.Abs(nodes[2 * i][j] - mean) + Math.Abs(nodes[2 * i + 1][j] - mean));
                deviation *= half;
                double error = Math.Abs(kronrod[j] - gauss[j]) * half;
                if (deviation != 0 && error != 0) error = deviation * Math.Min(1, Math.Pow(200 * error / deviation, 1.5));
                errors[j] = Math.Max(error, 50 * Roundoff * half * absolute[j]);
                kronrod[j] *= half;
            }
            return new Interval { Lower = a, Upper = b, Values = kronrod, Errors = errors };
        }

        /// <summary>Refines the interval with the largest normalized error until every component converges.</summary>
        private static Interval Integrate(double k, double h, int count, int power)
        {
            var initial = Evaluate(0, 1, k, h, count, power);
            var intervals = new List<Interval> { initial };
            var total = new Interval { Values = (double[])initial.Values.Clone(), Errors = (double[])initial.Errors.Clone() };
            while (true)
            {
                bool success = true;
                for (int j = 0; j < total.Values.Length; j++)
                    success &= total.Errors[j] <= AbsoluteTolerance + RelativeTolerance * Math.Abs(total.Values[j]);
                if (success) return total;
                if (intervals.Count >= MaximumIntervals)
                    throw new InvalidOperationException("Complete-tail expected-information quadrature did not meet its error tolerances.");
                double worst = -1;
                int selected = 0;
                for (int i = 0; i < intervals.Count; i++)
                for (int j = 0; j < total.Values.Length; j++)
                {
                    double ratio = intervals[i].Errors[j] / (AbsoluteTolerance + RelativeTolerance * Math.Abs(total.Values[j]));
                    if (ratio > worst) { worst = ratio; selected = i; }
                }
                Interval old = intervals[selected];
                double center = .5 * (old.Lower + old.Upper);
                if (center == old.Lower || center == old.Upper)
                    throw new InvalidOperationException("Expected-information quadrature exhausted floating-point subdivision resolution.");
                Interval left = Evaluate(old.Lower, center, k, h, count, power);
                Interval right = Evaluate(center, old.Upper, k, h, count, power);
                intervals[selected] = left;
                intervals.Add(right);
                // Re-sum error estimates to avoid a negative error from subtracting parent estimates.
                Array.Clear(total.Values, 0, total.Values.Length);
                Array.Clear(total.Errors, 0, total.Errors.Length);
                foreach (Interval interval in intervals)
                for (int j = 0; j < total.Values.Length; j++)
                {
                    total.Values[j] += interval.Values[j];
                    total.Errors[j] += interval.Errors[j];
                }
            }
        }

        /// <summary>Inverts an SPD matrix with diagonal scaling and checked Cholesky solves.</summary>
        private static double[,] InvertInformation(double[,] information)
        {
            int count = information.GetLength(0);
            var scale = new double[count];
            var normalized = new double[count, count];
            var lower = new double[count, count];
            for (int i = 0; i < count; i++)
            {
                if (!(information[i, i] > 0) || !DistributionNumerics.IsFinite(information[i, i]))
                    throw new InvalidOperationException("Expected information has a nonpositive or nonfinite diagonal.");
                scale[i] = Math.Sqrt(information[i, i]);
            }
            for (int i = 0; i < count; i++)
            for (int j = 0; j < count; j++) normalized[i, j] = information[i, j] / scale[i] / scale[j];
            for (int i = 0; i < count; i++)
            for (int j = 0; j <= i; j++)
            {
                double value = normalized[i, j];
                for (int m = 0; m < j; m++) value -= lower[i, m] * lower[j, m];
                if (i == j)
                {
                    if (!(value > 0) || !DistributionNumerics.IsFinite(value))
                        throw new InvalidOperationException("Expected information is not numerically positive definite.");
                    lower[i, j] = Math.Sqrt(value);
                }
                else lower[i, j] = value / lower[j, j];
            }
            var inverse = new double[count, count];
            for (int column = 0; column < count; column++)
            {
                var solution = new double[count];
                for (int i = 0; i < count; i++)
                {
                    double value = i == column ? 1 : 0;
                    for (int j = 0; j < i; j++) value -= lower[i, j] * solution[j];
                    solution[i] = value / lower[i, i];
                }
                for (int i = count - 1; i >= 0; i--)
                {
                    double value = solution[i];
                    for (int j = i + 1; j < count; j++) value -= lower[j, i] * solution[j];
                    solution[i] = value / lower[i, i];
                    inverse[i, column] = solution[i];
                }
                for (int i = 0; i < count; i++)
                {
                    double value = 0;
                    for (int j = 0; j < count; j++) value += normalized[i, j] * solution[j];
                    if (!DistributionNumerics.IsFinite(value) || Math.Abs(value - (i == column ? 1 : 0)) > 1E-9)
                        throw new InvalidOperationException("The expected-information inverse failed its residual check.");
                }
            }
            for (int i = 0; i < count; i++)
            for (int j = i; j < count; j++)
            {
                double value = .5 * (inverse[i, j] + inverse[j, i]) / scale[i] / scale[j];
                if (!DistributionNumerics.IsFinite(value)) throw new InvalidOperationException("Expected-information inversion was nonfinite.");
                inverse[i, j] = inverse[j, i] = value;
            }
            return inverse;
        }
    }
}
