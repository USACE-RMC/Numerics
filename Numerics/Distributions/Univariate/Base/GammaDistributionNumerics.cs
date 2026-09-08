using System;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{
    internal static partial class DistributionNumerics
    {
        /// <summary>Trigamma with an exact recurrence and a sufficiently large asymptotic argument for covariance work.</summary>
        internal static double AccurateTrigamma(double shape)
        {
            ValidateGammaShape(shape);
            double sum = 0;
            while (shape < 32) { double inverse = 1 / shape; sum += inverse * inverse; shape++; }
            double r = 1 / shape, s = r * r;
            return sum + r + .5 * s + r * s * (1.0 / 6 + s * (-1.0 / 30 + s * (1.0 / 42 + s * (-1.0 / 30 + s * (5.0 / 66 - s * 691.0 / 2730)))));
        }

        /// <summary>Log regularized lower gamma integral P(a,x).</summary>
        internal static double GammaLogCDF(double shape, double x) => GammaLogTail(shape, x, false, out _);

        /// <summary>Log regularized upper gamma integral Q(a,x).</summary>
        internal static double GammaLogSurvival(double shape, double x) => GammaLogTail(shape, x, true, out _);

        /// <summary>Log unit-scale gamma density, including one-sided endpoint limits.</summary>
        internal static double GammaLogDensity(double shape, double x)
        {
            ValidateGammaShape(shape);
            if (double.IsNaN(x)) return double.NaN;
            if (x < 0 || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
            if (x == 0) return shape == 1 ? 0 : shape < 1 ? double.PositiveInfinity : double.NegativeInfinity;
            return shape == 1 ? -x : GammaLogKernel(shape, x) - Math.Log(x);
        }

        /// <summary>Implicit shape derivative of the actual unit-scale gamma quantile at its value.</summary>
        /// <remarks>Differentiates the convergent lower series or upper continued fraction together
        /// with the probability. The large-shape expansion is differentiated analytically. No
        /// frequency-factor approximation or perturbation across the shape boundary is used.</remarks>
        internal static double GammaQuantileShapeDerivative(double shape, double unitQuantile)
        {
            ValidateGammaShape(shape);
            if (unitQuantile == 0) return 0;
            if (!(unitQuantile > 0) || !IsFinite(unitQuantile))
                throw new ArgumentOutOfRangeException(nameof(unitQuantile));
            bool upper = unitQuantile >= shape;
            double log = GammaLogTail(shape, unitQuantile, upper, out double derivative);
            if (derivative == 0) return 0;
            return (upper ? Math.Sign(derivative) : -Math.Sign(derivative))
                * Math.Exp(Math.Log(Math.Abs(derivative)) + log - GammaLogDensity(shape, unitQuantile));
        }

        /// <summary>Inverts a gamma tail directly, retaining tiny upper and lower probabilities.</summary>
        /// <remarks>Uses a tail-aware bracketed Newton solve; a normal approximation supplies only
        /// the initial point. Tiny quantiles are solved in logarithmic coordinates.</remarks>
        internal static double GammaInverseCDF(double shape, double probability, bool upperTail = false)
        {
            ValidateGammaShape(shape);
            if (!(probability >= 0 && probability <= 1)) throw new ArgumentOutOfRangeException(nameof(probability));
            if (probability == 0) return upperTail ? double.PositiveInfinity : 0;
            if (probability == 1) return upperTail ? 0 : double.PositiveInfinity;
            if (probability > 0.5) { probability = 1 - probability; upperTail = !upperTail; }
            double target = Math.Log(probability);
            if (shape == 1) return upperTail ? -target : -Tools.Log1p(-probability);
            double smallLog = (target + LogGammaOnePlus(shape)) / shape;
            if (!upperTail && smallLog < -36) return Math.Exp(smallLog);
            double root = Math.Sqrt(shape);
            double z = Normal.StandardZ(probability) * (upperTail ? -1 : 1);
            double w = 1 - 1 / (9 * shape) + z / (3 * root);
            double guess = w > 0 ? shape * w * w * w : Math.Exp(smallLog);
            if (!(guess > 0) || !IsFinite(guess)) guess = Math.Max(shape, 1);
            double lower = 0, upper = Math.Max(Math.Max(shape, 1), guess);
            bool Below(double value)
            {
                double log = GammaLogTail(shape, value, upperTail, out _);
                return upperTail ? log > target : log < target;
            }
            while (Below(upper))
            {
                lower = upper;
                upper = upper < double.MaxValue / 2 ? upper * 2 : double.MaxValue;
                if (lower == upper) return double.PositiveInfinity;
            }
            double x = guess > lower && guess < upper ? guess : lower + (upper - lower) / 2;
            for (int iteration = 0; iteration < 100; iteration++)
            {
                double log = GammaLogTail(shape, x, upperTail, out _);
                double residual = log - target;
                if (Math.Abs(residual) <= 8E-15 * Math.Max(1, Math.Abs(target))) return x;
                if (upperTail ? residual > 0 : residual < 0) lower = x; else upper = x;
                double slope = Math.Exp(GammaLogDensity(shape, x) - log) * (upperTail ? -1 : 1);
                double next = x - residual / slope;
                if (!(next > lower && next < upper) || !IsFinite(next)) next = lower + (upper - lower) / 2;
                if (next == x || next == lower || next == upper) return next;
                x = next;
            }
            throw new InvalidOperationException("The gamma quantile solve did not converge.");
        }

        /// <summary>Requires a positive finite gamma shape.</summary>
        private static void ValidateGammaShape(double shape)
        {
            if (!(shape > 0) || !IsFinite(shape)) throw new ArgumentOutOfRangeException(nameof(shape), "Gamma shape must be positive and finite.");
        }

        /// <summary>Evaluates a gamma log tail and its fixed-observation shape derivative.</summary>
        private static double GammaLogTail(double a, double x, bool upper, out double derivative)
        {
            ValidateGammaShape(a);
            derivative = 0;
            if (double.IsNaN(x)) { derivative = double.NaN; return double.NaN; }
            if (x <= 0) return upper ? 0 : double.NegativeInfinity;
            if (double.IsPositiveInfinity(x)) return upper ? double.NegativeInfinity : 0;
            double delta = (x - a) / a;
            if (a >= 10000 && Math.Abs(delta) < 0.1) return GammaTemme(a, delta, upper, out derivative);

            if (a < 1 && x <= 1 && (upper || a * Math.Log(x) - LogGammaOnePlus(a) > -0.6931471805599453))
            {
                double logQ = GammaSmallUpper(a, x, out double dlogQ);
                if (upper) { derivative = dlogQ; return logQ; }
                double logP = Log1mExp(logQ);
                derivative = -dlogQ * Math.Exp(logQ - logP);
                return logP;
            }

            if (x < a + 1)
            {
                double sum = 1, term = 1, dsum = 0, dterm = 0;
                for (int n = 1; n <= 100000; n++)
                {
                    double ratio = x / (a + n);
                    dterm = dterm * ratio - term * ratio / (a + n);
                    term *= ratio;
                    sum += term;
                    dsum += dterm;
                    if (term <= sum * 2E-16 && Math.Abs(dterm) <= Math.Max(1, Math.Abs(dsum)) * 2E-16)
                    {
                        double log = a < 16 ? a * Math.Log(x) - x - LogGammaOnePlus(a) + Math.Log(sum)
                            : GammaLogKernel(a, x) - Math.Log(a) + Math.Log(sum);
                        double dlog = a < 16 ? Math.Log(x) - Gamma.Digamma(a + 1) + dsum / sum
                            : GammaKernelShapeDerivative(a, x) - 1 / a + dsum / sum;
                        if (!upper) { derivative = dlog; return log; }
                        double complement = Log1mExp(log);
                        derivative = -dlog * Math.Exp(log - complement);
                        return complement;
                    }
                }
                throw new InvalidOperationException("The lower gamma series did not converge.");
            }
            else
            {
                // Modified Lentz recurrence for DLMF 8.9.2, carrying its derivative.
                const double tiny = 1E-300;
                double b = x + 1 - a, c = 1 / tiny, dc = 0;
                double d = 1 / b, dd = d * d, h = d, dh = dd;
                for (int n = 1; n <= 100000; n++)
                {
                    double an = n * (a - n);
                    b += 2;
                    double denom = an * d + b;
                    double ddenom = n * d + an * dd - 1;
                    double nextC = b + an / c;
                    double nextDc = -1 + n / c - an * dc / c / c;
                    if (Math.Abs(denom) < tiny) denom = denom < 0 ? -tiny : tiny;
                    if (Math.Abs(nextC) < tiny) nextC = nextC < 0 ? -tiny : tiny;
                    d = 1 / denom;
                    dd = -ddenom * d * d;
                    c = nextC; dc = nextDc;
                    double factor = d * c, dfactor = dd * c + d * dc;
                    dh = dh * factor + h * dfactor;
                    h *= factor;
                    if (Math.Abs(factor - 1) <= 4E-16 && Math.Abs(dfactor) <= 4E-16)
                    {
                        double log = GammaLogKernel(a, x) + Math.Log(h);
                        double dlog = GammaKernelShapeDerivative(a, x) + dh / h;
                        if (upper) { derivative = dlog; return log; }
                        double complement = Log1mExp(log);
                        derivative = -dlog * Math.Exp(log - complement);
                        return complement;
                    }
                }
                throw new InvalidOperationException("The upper gamma continued fraction did not converge.");
            }
        }

        /// <summary>Direct Q series at small shape/argument, avoiding a near-one lower-tail subtraction.</summary>
        private static double GammaSmallUpper(double a, double x, out double derivative)
        {
            double sum = 0, dsum = 0, power = 1;
            for (int n = 1; n <= 1000; n++)
            {
                power *= -x / n;
                double term = power / (a + n);
                sum += term;
                dsum -= term / (a + n);
                if (Math.Abs(term) < Math.Abs(sum) * 2E-16) break;
            }
            double logx = Math.Log(x), u = a * logx - LogGammaOnePlus(a);
            double leading = Math.Exp(u), factor = a * leading;
            double q = -Tools.Expm1(u) - factor * sum;
            double du = logx - Gamma.Digamma(a + 1);
            derivative = (-leading * du - factor * (du * sum + dsum) - leading * sum) / q;
            return Math.Log(q);
        }

        /// <summary>Log(x^a exp(-x)/Gamma(a)) without subtracting large near-equal terms.</summary>
        private static double GammaLogKernel(double a, double x)
        {
            if (a < 16) return a * Math.Log(x) - x - Gamma.LogGamma(a);
            double delta = (x - a) / a;
            double deviance = Math.Abs(delta) < 0.25 ? a * Log1pMinusX(delta)
                : a * (Math.Log(x) - Math.Log(a)) - (x - a);
            return deviance + 0.5 * Math.Log(a) - Tools.LogSqrt2PI - StirlingRemainder(a);
        }

        /// <summary>Fixed-x derivative of the log kernel, retaining its small residual at large shape.</summary>
        private static double GammaKernelShapeDerivative(double a, double x)
        {
            if (a < 16) return Math.Log(x) - Gamma.Digamma(a);
            double inverse = 1 / a, square = inverse * inverse;
            double residual = 0.5 * inverse + square * (1.0 / 12 - square * (1.0 / 120 - square * (1.0 / 252 - square / 240)));
            double delta = (x - a) / a;
            return (Math.Abs(delta) < 0.25 ? Tools.Log1p(delta) : Math.Log(x) - Math.Log(a)) + residual;
        }

        /// <summary>Stirling log-gamma remainder for arguments at least sixteen.</summary>
        private static double StirlingRemainder(double a)
        {
            double r = 1 / a, r2 = r * r;
            return r * (1.0 / 12 + r2 * (-1.0 / 360 + r2 * (1.0 / 1260 + r2 * (-1.0 / 1680 + r2 * (1.0 / 1188 - r2 * 691.0 / 360360)))));
        }

        /// <summary>Cancellation-free log(1+x)-x.</summary>
        internal static double Log1pMinusX(double x)
        {
            if (Math.Abs(x) >= 0.25) return Tools.Log1p(x) - x;
            double power = -x * x, sum = power / 2;
            for (int n = 3; n < 100; n++)
            {
                power *= -x;
                double term = power / n;
                sum += term;
                if (Math.Abs(term) <= Math.Abs(sum) * 1E-17) break;
            }
            return sum;
        }

        /// <summary>Log Gamma(1+a) with the zeta Taylor series at small a.</summary>
        internal static double LogGammaOnePlus(double a)
        {
            if (a > 0.5) return Gamma.LogGamma(a + 1);
            double result = -0.57721566490153286061 * a, power = a;
            for (int n = 2; n < 60; n++)
            {
                power *= -a;
                double term = -power * ZetaInteger(n) / n;
                result += term;
                if (Math.Abs(term) <= Math.Abs(result) * 1E-17) break;
            }
            return result;
        }

        /// <summary>Integer zeta constants for the log Gamma(1+a) series.</summary>
        internal static double ZetaInteger(int n)
        {
            double[] values = { 1.6449340668482264365, 1.2020569031595942854, 1.0823232337111381915,
                1.0369277551433699263, 1.0173430619844491397, 1.0083492773819228268,
                1.0040773561979443394, 1.0020083928260822144, 1.0009945751278180853,
                1.0004941886041194646, 1.0002460865533080483, 1.0001227133475784891,
                1.0000612481350587048, 1.0000305882363070205, 1.0000152822594086519 };
            if (n <= 16) return values[n - 2];
            double sum = 1;
            for (int k = 2; k <= 32; k++) sum += Math.Pow(k, -n);
            return sum;
        }

        // Exact rational coefficients generated by docs/distributions/oracles/generate-gamma-temme-coefficients.py.
        // DLMF 8.12.9-11. Four inverse-shape terms suffice in the a>=10000, |x/a-1|<.1 region.
        private static readonly double[][] TemmeCoefficients =
        {
            new double[] { -.3333333333333333333,.0833333333333333333,-.01481481481481481481,.001157407407407407407,.0003527336860670194,-.0001787551440329218,3.919263178522438E-5,-2.185448510679992E-6,-1.85406221071516E-6,8.296711340953087E-7,-1.766595273682608E-7,6.707853543401498E-9,1.026180978424031E-8,-4.382036018453353E-9,9.14769958223679E-10,-2.551419399494625E-11,-5.830772132550426E-11,2.436194802066742E-11 },
            new double[] { -.001851851851851851852,-.003472222222222222222,.002645502645502645503,-.0009902263374485596,.000205761316872428,-4.018775720164609E-7,-1.809855033448998E-5,7.64916091608111E-6,-1.612090089456345E-6,4.647127802807434E-9,1.378633446915721E-7,-5.752545603517705E-8,1.195162859977815E-8,-1.754324171974765E-11,-1.009154371060041E-9,4.162792991842583E-10 },
            new double[] { .004133597883597883598,-.00268132716049382716,.0007716049382716049,2.009387860082305E-6,-.0001073665322636516,5.292344882912013E-5,-1.276063518861873E-5,3.423578734096138E-8,1.372195730906293E-6,-6.298992138380055E-7,1.428061420606424E-7,-2.047709842199087E-10,-1.409252991086752E-8,6.228974084922022E-9 },
            new double[] { .0006494341563786008,.0002294720936213992,-.0004691894943952557,.0002677206320628389,-7.561801671883977E-5,-2.396505113867297E-7,1.10826541153473E-5,-5.674952826991597E-6,1.423090073243588E-6,-2.786108029152814E-11,-1.695840409193028E-7,8.099464905388083E-8 }
        };

        /// <summary>Uniform gamma expansion with analytical fixed-x shape differentiation.</summary>
        private static double GammaTemme(double a, double delta, bool upper, out double derivative)
        {
            double eta = delta == 0 ? 0 : Math.Sign(delta) * Math.Sqrt(-2 * Log1pMinusX(delta));
            double root = Math.Sqrt(a), z = eta * root;
            double etaDerivative = delta == 0 ? -1 / a : -delta / eta / a;
            double zDerivative = eta / (2 * root) + root * etaDerivative;
            double series = 0, seriesEta = 0, seriesA = 0, inversePower = 1;
            for (int k = 0; k < TemmeCoefficients.Length; k++)
            {
                var coefficients = TemmeCoefficients[k];
                double value = coefficients[coefficients.Length - 1], dvalue = 0;
                for (int n = coefficients.Length - 2; n >= 0; n--) { dvalue = dvalue * eta + value; value = value * eta + coefficients[n]; }
                series += inversePower * value;
                seriesEta += inversePower * dvalue;
                seriesA -= k * inversePower / a * value;
                inversePower /= a;
            }
            double normal = upper ? NormalLogSurvival(z) : NormalLogCDF(z);
            double logPhi = -(0.5 * z) * z - Tools.LogSqrt2PI;
            double signedCorrection = (upper ? 1 : -1) * series * Math.Exp(logPhi - 0.5 * Math.Log(a) - normal);
            double log = normal + Tools.Log1p(signedCorrection);
            double correctionDerivative = (seriesEta * etaDerivative + seriesA - (z * zDerivative + 0.5 / a) * series) / root;
            derivative = Math.Exp(logPhi - log) * (upper ? -zDerivative + correctionDerivative : zDerivative - correctionDerivative);
            return log;
        }
    }
}
