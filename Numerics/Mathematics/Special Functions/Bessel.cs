using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.SpecialFunctions
{

    /// <summary>
    /// Contains methods for evaluating Bessel functions of the first kind (J), second kind (Y),
    /// modified Bessel functions of the first kind (I), and modified Bessel functions of the second kind (K).
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
    /// Bessel functions are solutions to Bessel's differential equation:
    /// <code>
    ///     x² y'' + x y' + (x² - n²) y = 0
    /// </code>
    /// They arise naturally in problems with cylindrical or spherical symmetry, including heat conduction,
    /// wave propagation, electrostatics, and fluid dynamics. In statistics, the modified Bessel functions
    /// I₀ and I₁ appear in the von Mises distribution for circular data.
    /// </para>
    /// <para>
    /// This class provides the following functions:
    /// <list type="bullet">
    /// <item><description><b>J₀(x), J₁(x), Jₙ(x)</b>: Bessel functions of the first kind. Solutions regular at the origin.</description></item>
    /// <item><description><b>Y₀(x), Y₁(x), Yₙ(x)</b>: Bessel functions of the second kind (Neumann functions). Singular at the origin.</description></item>
    /// <item><description><b>I₀(x), I₁(x), Iₙ(x)</b>: Modified Bessel functions of the first kind. Exponentially growing solutions of the modified equation.</description></item>
    /// <item><description><b>I₀e(x), I₁e(x)</b>: Exponentially scaled exp(-|x|) I₀(x) and exp(-|x|) I₁(x), accurate to double precision and finite for every argument (Cephes Chebyshev expansions).</description></item>
    /// <item><description><b>K₀(x), K₁(x), Kₙ(x)</b>: Modified Bessel functions of the second kind. Exponentially decaying solutions of the modified equation.</description></item>
    /// </list>
    /// </para>
    /// <para>
    /// For orders 0 and 1, rational polynomial and polynomial approximations are used directly. For arbitrary
    /// integer orders, stable recurrence relations are employed: forward recurrence for Y and K (which are
    /// stable in the forward direction), and Miller's downward recurrence for J and I (which are stable in the
    /// backward direction), with normalization against the known order-0 values.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Abramowitz, M. and Stegun, I.A. (1972). "Handbook of Mathematical Functions."
    /// National Bureau of Standards, Applied Mathematics Series 55. Sections 9.1-9.8.
    /// </description></item>
    /// <item><description>
    /// Press, W.H., Teukolsky, S.A., Vetterling, W.T. and Flannery, B.P. (2007).
    /// "Numerical Recipes: The Art of Scientific Computing," 3rd ed. Cambridge University Press.
    /// Sections 6.5-6.7.
    /// </description></item>
    /// <item><description>
    /// Olver, F.W.J. et al. (2010). "NIST Handbook of Mathematical Functions."
    /// Cambridge University Press. Chapter 10.
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public static class Bessel
    {

        /// <summary>
        /// Accuracy parameter for Miller's downward recurrence. Larger values give more accurate results
        /// at the cost of more iterations. A value of 40 gives approximately double precision accuracy.
        /// </summary>
        private const int IACC = 40;

        /// <summary>
        /// Rescaling threshold for preventing overflow during Miller's downward recurrence.
        /// </summary>
        private const double BIGNO = 1.0e10;

        /// <summary>
        /// Inverse of BIGNO for rescaling during Miller's downward recurrence.
        /// </summary>
        private const double BIGNI = 1.0e-10;

        #region Modified Bessel Functions of the First Kind (I)

        /// <summary>
        /// Computes the modified Bessel function of the first kind, order 0: I₀(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate I₀. Can be any real number.</param>
        /// <returns>The value of I₀(x). Always positive. I₀(x) = I₀(-x).</returns>
        /// <remarks>
        /// <para>
        /// I₀(x) = exp(|x|) I₀e(x), where I₀e is the double-precision Chebyshev expansion of <see cref="I0e(double)"/>
        /// (Cephes i0, Moshier 1989). The exponential is split for |x| &gt; 709, so the result overflows only
        /// where I₀ itself exceeds the double range (|x| above about 713.9). Measured against mpmath at 50 digits
        /// over |x| in [1E-8, 700] the worst relative error is about 1E-15 (see
        /// docs/mathematics/oracles/special-functions).
        /// </para>
        /// <para>
        /// I₀(x) is the zeroth-order modified Bessel function of the first kind, satisfying:
        /// <code>
        ///     x² y'' + x y' - x² y = 0
        /// </code>
        /// It is related to J₀ by I₀(x) = J₀(ix), and appears as the normalization constant
        /// in the von Mises distribution: f(x|mu,kappa) = exp(kappa*cos(x-mu)) / (2*pi*I₀(kappa)).
        /// </para>
        /// </remarks>
        public static double I0(double x)
        {
            double ax = Math.Abs(x);
            if (ax <= 8.0)
                return Math.Exp(ax) * Evaluate.Chebyshev(I0eSmall, ax / 2.0 - 2.0);
            return ScaledExp(ax, Evaluate.Chebyshev(I0eLarge, 32.0 / ax - 2.0));
        }

        /// <summary>
        /// Evaluates exp(ax) / sqrt(ax) * p without overflowing exp(ax) before the product does.
        /// </summary>
        /// <param name="ax">The non-negative argument.</param>
        /// <param name="p">The polynomial factor of the asymptotic expansion.</param>
        /// <returns>exp(ax) / sqrt(ax) * p, or +∞ only when the product itself overflows.</returns>
        /// <remarks>
        /// Below 709 the direct product is used unchanged; above it exp(ax) is split into two halves.
        /// </remarks>
        private static double ScaledExp(double ax, double p)
        {
            if (double.IsPositiveInfinity(ax))
                return double.PositiveInfinity;
            if (ax <= 709.0)
                return (Math.Exp(ax) / Math.Sqrt(ax)) * p;
            double h = Math.Exp(0.5 * ax);
            return h * (h / Math.Sqrt(ax) * p);
        }

        /// <summary>
        /// Computes the modified Bessel function of the first kind, order 1: I₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate I₁. Can be any real number.</param>
        /// <returns>The value of I₁(x). I₁(-x) = -I₁(x) (odd function).</returns>
        /// <remarks>
        /// <para>
        /// I₁(x) = exp(|x|) I₁e(x), where I₁e is the double-precision Chebyshev expansion of <see cref="I1e(double)"/>
        /// (Cephes i1, Moshier 1989). The exponential is split for |x| &gt; 709, so the result overflows only
        /// where I₁ itself exceeds the double range. Measured against mpmath at 50 digits over |x| in [1E-8, 700]
        /// the worst relative error is about 2E-15.
        /// </para>
        /// <para>
        /// I₁(x) appears in circular statistics: the mean resultant length of the von Mises distribution
        /// is A(kappa) = I₁(kappa)/I₀(kappa), and the circular variance is V = 1 - A(kappa).
        /// </para>
        /// </remarks>
        public static double I1(double x)
        {
            double ax = Math.Abs(x);
            double result = ax <= 8.0
                ? Math.Exp(ax) * (Evaluate.Chebyshev(I1eSmall, ax / 2.0 - 2.0) * ax)
                : ScaledExp(ax, Evaluate.Chebyshev(I1eLarge, 32.0 / ax - 2.0));
            return x < 0 ? -result : result;
        }

        /// <summary>
        /// Chebyshev coefficients for exp(-x) I₀(x) on [0, 8], order-zero term first (Cephes i0.c, array A).
        /// </summary>
        private static readonly double[] I0eSmall =
        {
            6.76795274409476084995E-1, -3.04682672343198398683E-1, 1.71620901522208775349E-1,
            -9.49010970480476444210E-2, 4.93052842396707084878E-2, -2.37374148058994688156E-2,
            1.05464603945949983183E-2, -4.32430999505057594430E-3, 1.63947561694133579842E-3,
            -5.76375574538582365885E-4, 1.88502885095841655729E-4, -5.75419501008210370398E-5,
            1.64484480707288970893E-5, -4.41673835845875056359E-6, 1.11738753912010371815E-6,
            -2.67079385394061173391E-7, 6.04699502254191894932E-8, -1.30002500998624804212E-8,
            2.65982372468238665035E-9, -5.18979560163526290666E-10, 9.67580903537323691224E-11,
            -1.72682629144155570723E-11, 2.95505266312963983461E-12, -4.85644678311192946090E-13,
            7.67618549860493561688E-14, -1.16853328779934516808E-14, 1.71539128555513303061E-15,
            -2.43127984654795469359E-16, 3.33079451882223809783E-17, -4.41534164647933937950E-18
        };

        /// <summary>
        /// Chebyshev coefficients for exp(-x) sqrt(x) I₀(x) on the inverted interval [8, ∞), order-zero term
        /// first (Cephes i0.c, array B).
        /// </summary>
        private static readonly double[] I0eLarge =
        {
            8.04490411014108831608E-1, 3.36911647825569408990E-3, 6.88975834691682398426E-5,
            2.89137052083475648297E-6, 2.04891858946906374183E-7, 2.26666899049817806459E-8,
            3.39623202570838634515E-9, 4.94060238822496958910E-10, 1.18891471078464383424E-11,
            -3.14991652796324136454E-11, -1.32158118404477131188E-11, -1.79417853150680611778E-12,
            7.18012445138366623367E-13, 3.85277838274214270114E-13, 1.54008621752140982691E-14,
            -4.15056934728722208663E-14, -9.55484669882830764870E-15, 3.81168066935262242075E-15,
            1.77256013305652638360E-15, -3.42548561967721913462E-16, -2.82762398051658348494E-16,
            3.46122286769746109310E-17, 4.46562142029675999901E-17, -4.83050448594418207126E-18,
            -7.23318048787475395456E-18
        };

        /// <summary>
        /// Chebyshev coefficients for exp(-x) I₁(x) / x on [0, 8], order-zero term first (Cephes i1.c, array A).
        /// </summary>
        private static readonly double[] I1eSmall =
        {
            2.52587186443633654823E-1, -1.76416518357834055153E-1, 1.02643658689847095384E-1,
            -5.29459812080949914269E-2, 2.47264490306265168283E-2, -1.05640848946261981558E-2,
            4.15642294431288815669E-3, -1.51357245063125314899E-3, 5.12285956168575772895E-4,
            -1.61760815825896745588E-4, 4.78156510755005422638E-5, -1.32731636560394358279E-5,
            3.47025130813767847674E-6, -8.56872026469545474066E-7, 2.00329475355213526229E-7,
            -4.44505912879632808065E-8, 9.38153738649577178388E-9, -1.88724975172282928790E-9,
            3.62559028155211703701E-10, -6.66348972350202774223E-11, 1.17361862988909016308E-11,
            -1.98397439776494371520E-12, 3.22379336594557470981E-13, -5.04218550472791168711E-14,
            7.60068429473540693410E-15, -1.10559694773538630805E-15, 1.55363195773620046921E-16,
            -2.11142121435816608115E-17, 2.77791411276104639959E-18
        };

        /// <summary>
        /// Chebyshev coefficients for exp(-x) sqrt(x) I₁(x) on the inverted interval [8, ∞), order-zero term
        /// first (Cephes i1.c, array B).
        /// </summary>
        private static readonly double[] I1eLarge =
        {
            7.78576235018280120474E-1, -9.76109749136146840777E-3, -1.10588938762623716291E-4,
            -3.88256480887769039346E-6, -2.51223623787020892529E-7, -2.63146884688951950684E-8,
            -3.83538038596423702205E-9, -5.58974346219658380687E-10, -1.89749581235054123450E-11,
            3.25260358301548823856E-11, 1.41258074366137813316E-11, 2.03562854414708950722E-12,
            -7.19855177624590851209E-13, -4.08355111109219731823E-13, -2.10154184277266431302E-14,
            4.27244001671195135429E-14, 1.04202769841288027642E-14, -3.81440307243700780478E-15,
            -1.88035477551078244854E-15, 3.30820231092092828324E-16, 2.96262899764595013876E-16,
            -3.20952592199342395980E-17, -4.65030536848935832153E-17, 4.41434832307170791151E-18,
            7.51729631084210481353E-18
        };

        /// <summary>
        /// Computes the exponentially scaled modified Bessel function of the first kind, order 0:
        /// I₀e(x) = exp(-|x|) I₀(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate I₀e. Can be any real number.</param>
        /// <returns>
        /// The value of exp(-|x|) I₀(x), in (0, 1]. I₀e(0) = 1, I₀e(-x) = I₀e(x), I₀e(±∞) = 0, and NaN for NaN input.
        /// </returns>
        /// <remarks>
        /// <para>
        /// The scaling removes the exponential growth of I₀, so the function stays finite for every argument and
        /// ratios such as I₁(x)/I₀(x) = I₁e(x)/I₀e(x) can be formed without overflow. <see cref="I0(double)"/> is this
        /// expansion multiplied by exp(|x|).
        /// </para>
        /// <para>
        /// The range is split at |x| = 8. On [0, 8] a 30-term Chebyshev expansion of exp(-x) I₀(x) in x/2 - 2 is
        /// summed; above 8 a 25-term expansion of exp(-x) sqrt(x) I₀(x) in 32/x - 2 is summed and divided by sqrt(x).
        /// The coefficients and the method are those of i0e in the Cephes Math Library (Moshier, 1989; Cephes 2.8,
        /// 2000), whose stated peak relative error on [0, 30] is 5.4E-16. The coefficients are identical to netlib
        /// cephes/bessel.tgz (i0.c), and the results agree with scipy.special.i0e and with mpmath
        /// besseli(0, x) exp(-x) to within 1E-15 relative over x in [0, 1E4].
        /// </para>
        /// <para>
        /// References: Moshier, S. L. (1989). Methods and Programs for Mathematical Functions. Ellis Horwood,
        /// Chichester. Cephes Math Library, <see href="https://www.netlib.org/cephes/"/>.
        /// </para>
        /// </remarks>
        public static double I0e(double x)
        {
            double ax = Math.Abs(x);
            if (ax <= 8.0)
                return Evaluate.Chebyshev(I0eSmall, ax / 2.0 - 2.0);
            return Evaluate.Chebyshev(I0eLarge, 32.0 / ax - 2.0) / Math.Sqrt(ax);
        }

        /// <summary>
        /// Computes the exponentially scaled modified Bessel function of the first kind, order 1:
        /// I₁e(x) = exp(-|x|) I₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate I₁e. Can be any real number.</param>
        /// <returns>
        /// The value of exp(-|x|) I₁(x). I₁e(0) = 0, I₁e(-x) = -I₁e(x), I₁e(±∞) = ±0, and NaN for NaN input.
        /// </returns>
        /// <remarks>
        /// <para>
        /// The scaling removes the exponential growth of I₁, so the function stays finite for every argument.
        /// The ratio I₁e(x)/I₀e(x) equals I₁(x)/I₀(x), the mean resultant length A(x) of the von Mises distribution.
        /// <see cref="I1(double)"/> is this expansion multiplied by exp(|x|).
        /// </para>
        /// <para>
        /// The range is split at |x| = 8. On [0, 8] a 29-term Chebyshev expansion of exp(-x) I₁(x) / x in x/2 - 2
        /// is summed and multiplied by |x|; above 8 a 25-term expansion of exp(-x) sqrt(x) I₁(x) in 32/x - 2 is
        /// summed and divided by sqrt(x). The coefficients and the method are those of i1e in the Cephes Math
        /// Library (Moshier, 1989; Cephes 2.8, 2000), whose stated peak relative error on [0, 30] is 2.0E-15. The
        /// coefficients are identical to netlib cephes/bessel.tgz (i1.c), and the results agree with
        /// scipy.special.i1e to within 1E-15 relative and with mpmath besseli(1, x) exp(-x) to within 2E-15
        /// relative over x in [0, 1E4].
        /// </para>
        /// <para>
        /// References: Moshier, S. L. (1989). Methods and Programs for Mathematical Functions. Ellis Horwood,
        /// Chichester. Cephes Math Library, <see href="https://www.netlib.org/cephes/"/>.
        /// </para>
        /// </remarks>
        public static double I1e(double x)
        {
            double z = Math.Abs(x);
            if (z <= 8.0)
                z = Evaluate.Chebyshev(I1eSmall, z / 2.0 - 2.0) * z;
            else
                z = Evaluate.Chebyshev(I1eLarge, 32.0 / z - 2.0) / Math.Sqrt(z);
            return x < 0.0 ? -z : z;
        }

        /// <summary>
        /// Computes the modified Bessel function of the first kind for integer order n: Iₙ(x).
        /// </summary>
        /// <param name="n">The order of the function. Must be non-negative (n &gt;= 0).</param>
        /// <param name="x">The argument at which to evaluate Iₙ. Can be any real number.</param>
        /// <returns>
        /// The value of Iₙ(x). For odd n, Iₙ(-x) = -Iₙ(x). For even n, Iₙ(-x) = Iₙ(x).
        /// </returns>
        /// <remarks>
        /// <para>
        /// For n = 0 and n = 1, delegates to <see cref="I0(double)"/> and <see cref="I1(double)"/> respectively.
        /// For n &gt;= 2, uses Miller's downward recurrence algorithm, which is numerically stable:
        /// <code>
        ///     I_{n-1}(x) = (2n/x) I_n(x) + I_{n+1}(x)
        /// </code>
        /// The result is normalized against the known value of I₀(x). The starting order is
        /// 2(N + sqrt(40 N)) with N = max(n, ⌈|x|⌉), so the recurrence also converges when |x| &gt; n.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if n &lt; 0 or n &gt; 1,000,000.</exception>
        public static double In(int n, double x)
        {
            if (n < 0 || n > 1_000_000)
                throw new ArgumentOutOfRangeException(nameof(n), "The order n must be between 0 and 1,000,000.");
            if (n == 0) return I0(x);
            if (n == 1) return I1(x);
            if (x == 0.0) return 0.0;

            double ax = Math.Abs(x);
            double tox = 2.0 / ax;
            double bip = 0.0, bi = 1.0;
            double ans = 0.0;

            // Starting index for downward recurrence. The recurrence must start well above both n
            // and |x|, so the base order grows with |x| when |x| > n. Beyond |x| = 1000 the I0
            // normalization overflows anyway, so the base is capped there to bound the work.
            int nb = n;
            if (ax > n)
                nb = (int)Math.Ceiling(Math.Min(ax, Math.Max(1000.0, n)));
            int m = 2 * (nb + (int)Math.Sqrt(IACC * nb));

            for (int j = m; j > 0; j--)
            {
                double bim = bip + j * tox * bi;
                bip = bi;
                bi = bim;

                // Rescale to prevent overflow
                if (Math.Abs(bi) > BIGNO)
                {
                    ans *= BIGNI;
                    bi *= BIGNI;
                    bip *= BIGNI;
                }

                if (j == n) ans = bip;
            }

            // Normalize using I0
            ans *= I0(ax) / bi;
            return (x < 0.0 && (n & 1) != 0) ? -ans : ans;
        }

        #endregion

        #region Modified Bessel Functions of the Second Kind (K)

        /// <summary>
        /// Computes the modified Bessel function of the second kind, order 0: K₀(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate K₀. Must be positive (x &gt; 0).</param>
        /// <returns>The value of K₀(x). Always positive, monotonically decreasing.</returns>
        /// <remarks>
        /// <para>
        /// The method of Cephes k0 (Moshier, 1989) is used. For 0 &lt; x &lt;= 2, K₀(x) = A(x²) - ln(x/2) I₀(x),
        /// where A is a 10-term Chebyshev expansion in x² - 2. For x &gt; 2, K₀(x) = exp(-x) B(x) / sqrt(x),
        /// where B is a 26-term Chebyshev expansion of exp(x) sqrt(x) K₀(x) in 8/x - 2. The exponential is
        /// split above x = 700, so the result underflows gradually rather than through a subnormal exp(-x).
        /// </para>
        /// <para>
        /// The Chebyshev coefficients were computed with mpmath at 50 digits by
        /// docs/mathematics/oracles/special-functions/bessel_chebyshev_coefficients.py. Measured against
        /// mpmath over x in [1E-8, 700] the worst relative error is about 1E-15.
        /// </para>
        /// <para>
        /// K₀(x) decays exponentially as x → ∞ and diverges logarithmically as x → 0⁺.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if x &lt;= 0.</exception>
        public static double K0(double x)
        {
            if (x <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(x), "The argument x must be positive for K0.");

            if (x <= 2.0)
                return Evaluate.Chebyshev(K0Small, x * x - 2.0) - Math.Log(0.5 * x) * I0(x);
            return DecayingExp(x, Evaluate.Chebyshev(K0Large, 8.0 / x - 2.0));
        }

        /// <summary>
        /// Computes the modified Bessel function of the second kind, order 1: K₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate K₁. Must be positive (x &gt; 0).</param>
        /// <returns>The value of K₁(x). Always positive, monotonically decreasing.</returns>
        /// <remarks>
        /// <para>
        /// The method of Cephes k1 (Moshier, 1989) is used. For 0 &lt; x &lt;= 2, K₁(x) = ln(x/2) I₁(x) + A(x²) / x,
        /// where A is an 11-term Chebyshev expansion in x² - 2. For x &gt; 2, K₁(x) = exp(-x) B(x) / sqrt(x),
        /// where B is a 26-term Chebyshev expansion of exp(x) sqrt(x) K₁(x) in 8/x - 2. The coefficients were
        /// computed with mpmath at 50 digits. Measured against mpmath over x in [1E-8, 700] the worst relative
        /// error is about 1E-15.
        /// </para>
        /// <para>
        /// K₁(x) decays exponentially as x → ∞ and diverges as 1/x as x → 0⁺.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if x &lt;= 0.</exception>
        public static double K1(double x)
        {
            if (x <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(x), "The argument x must be positive for K1.");

            if (x <= 2.0)
                return Math.Log(0.5 * x) * I1(x) + Evaluate.Chebyshev(K1Small, x * x - 2.0) / x;
            return DecayingExp(x, Evaluate.Chebyshev(K1Large, 8.0 / x - 2.0));
        }

        /// <summary>
        /// Evaluates exp(-x) / sqrt(x) * b without a subnormal intermediate exp(-x).
        /// </summary>
        /// <param name="x">The argument, greater than 2.</param>
        /// <param name="b">The value of the scaled expansion exp(x) sqrt(x) K(x).</param>
        /// <returns>exp(-x) / sqrt(x) * b, rounded into the subnormal range only when the product itself is.</returns>
        private static double DecayingExp(double x, double b)
        {
            if (x <= 700.0)
                return Math.Exp(-x) * b / Math.Sqrt(x);
            double h = Math.Exp(-0.5 * x);
            return h * (h * b / Math.Sqrt(x));
        }

        /// <summary>
        /// Computes the modified Bessel function of the second kind for integer order n: Kₙ(x).
        /// </summary>
        /// <param name="n">The order of the function. Must be non-negative (n &gt;= 0).</param>
        /// <param name="x">The argument at which to evaluate Kₙ. Must be positive (x &gt; 0).</param>
        /// <returns>The value of Kₙ(x). Always positive for x &gt; 0.</returns>
        /// <remarks>
        /// <para>
        /// For n = 0 and n = 1, delegates to <see cref="K0(double)"/> and <see cref="K1(double)"/> respectively.
        /// For n &gt;= 2, uses the forward recurrence relation, which is numerically stable for K:
        /// <code>
        ///     K_{n+1}(x) = K_{n-1}(x) + (2n/x) K_n(x)
        /// </code>
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if n &lt; 0 or x &lt;= 0.</exception>
        public static double Kn(int n, double x)
        {
            if (n < 0)
                throw new ArgumentOutOfRangeException(nameof(n), "The order n must be non-negative.");
            if (x <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(x), "The argument x must be positive for Kn.");
            if (n == 0) return K0(x);
            if (n == 1) return K1(x);

            double tox = 2.0 / x;
            double bkm = K0(x);
            double bk = K1(x);

            for (int j = 1; j < n; j++)
            {
                double bkp = bkm + j * tox * bk;
                bkm = bk;
                bk = bkp;
            }

            return bk;
        }

        #endregion

        #region Bessel Functions of the First Kind (J)

        /// <summary>
        /// The argument above which Hankel's asymptotic series is summed for J₀, J₁, Y₀ and Y₁.
        /// </summary>
        private const double HankelThreshold = 25.0;

        /// <summary>
        /// Evaluates Hankel's auxiliary functions P and Q of order 0 or 1 for x &gt;= 8.
        /// </summary>
        /// <param name="order">The order, 0 or 1.</param>
        /// <param name="x">The argument, at least 8.</param>
        /// <param name="p">On return, P(x).</param>
        /// <param name="q">On return, Q(x).</param>
        /// <remarks>
        /// <para>
        /// P and Q are defined by J(x) = sqrt(2/(πx)) (P cos χ - Q sin χ) and Y(x) = sqrt(2/(πx)) (P sin χ + Q cos χ)
        /// with χ = x - (2ν + 1)π/4 (Abramowitz and Stegun 9.2.5 and 9.2.6). For 8 &lt;= x &lt;= 25 they are
        /// Chebyshev expansions in u = 8/x, computed with mpmath at 50 digits. Above 25 Hankel's asymptotic series
        /// (Abramowitz and Stegun 9.2.9 and 9.2.10; NIST DLMF 10.17.1 to 10.17.4) is summed term by term with
        /// a_k(ν) / x^k = a_(k-1)(ν) / x^(k-1) (4ν² - (2k - 1)²) / (8kx), until a term falls below 1E-18. The
        /// smallest term at x = 25 is about exp(-50), so the truncation is far below double precision.
        /// </para>
        /// </remarks>
        private static void HankelPQ(int order, double x, out double p, out double q)
        {
            if (x <= HankelThreshold)
            {
                double u = 8.0 / x;
                double t = (2.0 * u - (1.0 + 8.0 / HankelThreshold)) / (0.5 * (1.0 - 8.0 / HankelThreshold));
                p = Evaluate.Chebyshev(order == 0 ? P0Mid : P1Mid, t);
                q = u * Evaluate.Chebyshev(order == 0 ? Q0Mid : Q1Mid, t);
                return;
            }

            double mu = 4.0 * order * order;
            double term = 1.0;
            p = 1.0;
            q = 0.0;
            for (int k = 1; k <= 60; k++)
            {
                double odd = 2.0 * k - 1.0;
                term *= (mu - odd * odd) / (8.0 * k * x);
                // The k-th term enters P (even k) or Q (odd k) with the sign (-1)^floor(k/2).
                double signed = ((k >> 1) & 1) == 0 ? term : -term;
                if ((k & 1) == 0) p += signed;
                else q += signed;
                if (Math.Abs(term) < 1E-18) break;
            }
        }

        /// <summary>
        /// Computes the Bessel function of the first kind, order 0: J₀(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate J₀. Can be any real number.</param>
        /// <returns>The value of J₀(x). J₀(0) = 1. J₀(x) = J₀(-x) (even function). Range: [-0.4028, 1].</returns>
        /// <remarks>
        /// <para>
        /// For |x| &lt;= 8, J₀(x) = (x - j₁)(x + j₁)(x - j₂)(x + j₂) g(x²), where j₁ and j₂ are the first two
        /// zeros, each stored as the sum of two doubles, and g is a 15-term Chebyshev expansion in x² on [0, 64].
        /// Factoring the zeros out, as Cephes j0 does (Moshier, 1989), keeps the relative error small next to them.
        /// For |x| &gt; 8, J₀(x) = (P (cos x + sin x) - Q (sin x - cos x)) / sqrt(πx) with Hankel's P and Q, from
        /// Chebyshev expansions in 8/x up to x = 25 and from Hankel's asymptotic series above. The coefficients
        /// were computed with mpmath at 50 digits by
        /// docs/mathematics/oracles/special-functions/bessel_chebyshev_coefficients.py.
        /// </para>
        /// <para>
        /// Measured against mpmath, the relative error is about 1E-15 for |x| &lt;= 8, including within 1E-6 of
        /// the zeros, and the error above 8 is about 1E-15 of the envelope sqrt(2/(πx)), as for scipy.special.j0.
        /// Above 8 the relative error grows only next to a zero.
        /// </para>
        /// <para>
        /// J₀(x) is the unique solution of Bessel's equation of order 0 that is finite at x = 0.
        /// It oscillates with decreasing amplitude: J₀(x) ~ sqrt(2/(pi*x)) cos(x - pi/4) for large x.
        /// </para>
        /// </remarks>
        public static double J0(double x)
        {
            double ax = Math.Abs(x);
            if (ax <= 8.0)
            {
                double factor = (ax - J0Zero1Hi - J0Zero1Lo) * (ax + J0Zero1Hi)
                    * ((ax - J0Zero2Hi - J0Zero2Lo) * (ax + J0Zero2Hi));
                return factor * Evaluate.Chebyshev(J0Small, ax * ax / 16.0 - 2.0);
            }
            HankelPQ(0, ax, out double p, out double q);
            double s = Math.Sin(ax), c = Math.Cos(ax);
            return (p * (c + s) - q * (s - c)) / Math.Sqrt(Math.PI * ax);
        }

        /// <summary>
        /// Computes the Bessel function of the first kind, order 1: J₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate J₁. Can be any real number.</param>
        /// <returns>The value of J₁(x). J₁(0) = 0. J₁(-x) = -J₁(x) (odd function). Range: [-0.5819, 0.5819].</returns>
        /// <remarks>
        /// <para>
        /// For |x| &lt;= 8, J₁(x) = x (x - j₁)(x + j₁)(x - j₂)(x + j₂) g(x²), where j₁ and j₂ are the first two
        /// positive zeros, each stored as the sum of two doubles, and g is a 15-term Chebyshev expansion in x² on
        /// [0, 64]. For |x| &gt; 8, J₁(x) = (P (sin x - cos x) + Q (sin x + cos x)) / sqrt(πx) with Hankel's P and
        /// Q of order 1. The coefficients were computed with mpmath at 50 digits.
        /// </para>
        /// <para>
        /// Measured against mpmath, the relative error is about 1E-15 for |x| &lt;= 8, including next to the
        /// zeros, and the error above 8 is about 1E-15 of the envelope sqrt(2/(πx)), as for scipy.special.j1.
        /// </para>
        /// <para>
        /// J₁(x) oscillates with decreasing amplitude: J₁(x) ~ sqrt(2/(pi*x)) cos(x - 3*pi/4) for large x.
        /// </para>
        /// </remarks>
        public static double J1(double x)
        {
            double ax = Math.Abs(x);
            if (ax <= 8.0)
            {
                double factor = (ax - J1Zero1Hi - J1Zero1Lo) * (ax + J1Zero1Hi)
                    * ((ax - J1Zero2Hi - J1Zero2Lo) * (ax + J1Zero2Hi));
                // The factor x (not |x|) carries the odd symmetry, including the sign of zero.
                return x * factor * Evaluate.Chebyshev(J1Small, ax * ax / 16.0 - 2.0);
            }
            HankelPQ(1, ax, out double p, out double q);
            double s = Math.Sin(ax), c = Math.Cos(ax);
            double ans = (p * (s - c) + q * (s + c)) / Math.Sqrt(Math.PI * ax);
            return x < 0.0 ? -ans : ans;
        }

        /// <summary>
        /// Computes the Bessel function of the first kind for integer order n: Jₙ(x).
        /// </summary>
        /// <param name="n">The order of the function. Must be non-negative (n &gt;= 0).</param>
        /// <param name="x">The argument at which to evaluate Jₙ. Can be any real number.</param>
        /// <returns>
        /// The value of Jₙ(x). For odd n, Jₙ(-x) = -Jₙ(x). For even n, Jₙ(-x) = Jₙ(x).
        /// Jₙ(0) = 0 for n &gt; 0, and J₀(0) = 1.
        /// </returns>
        /// <remarks>
        /// <para>
        /// For n = 0 and n = 1, delegates to <see cref="J0(double)"/> and <see cref="J1(double)"/> respectively.
        /// For n &gt;= 2, uses two strategies depending on the relative magnitudes of n and |x|:
        /// </para>
        /// <para>
        /// When |x| &gt; n, forward recurrence from J₀ and J₁ is numerically stable:
        /// <code>
        ///     J_{n+1}(x) = (2n/x) J_n(x) - J_{n-1}(x)
        /// </code>
        /// </para>
        /// <para>
        /// When |x| &lt;= n, the ratio Jₙ(x)/Jₙ₋₁(x) is first found from its continued fraction
        /// (Abramowitz and Stegun 9.1.73), summed by the modified Lentz method, and the recurrence is run
        /// backward from orders n and n - 1 to orders 1 and 0, where it is normalized by <see cref="J0(double)"/>, or by
        /// <see cref="J1(double)"/> when |J₁(x)| &gt; |J₀(x)|, so a zero of J₀ never sets the scale. This is the scheme of
        /// Cephes jn (Moshier, 1989) and of Press et al. (2007), section 6.6. Measured against mpmath for
        /// n up to 50, the error is about 1E-15 relative where |x| &lt;= n and about 1E-15 of the envelope above.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if n &lt; 0 or n &gt; 1,000,000.</exception>
        public static double Jn(int n, double x)
        {
            if (n < 0 || n > 1_000_000)
                throw new ArgumentOutOfRangeException(nameof(n), "The order n must be between 0 and 1,000,000.");
            if (n == 0) return J0(x);
            if (n == 1) return J1(x);

            double ax = Math.Abs(x);
            if (ax == 0.0) return 0.0;

            double ans;
            double tox = 2.0 / ax;
            if (ax > (double)n)
            {
                // Forward recurrence from J0, J1 (stable for x > n)
                double bjm = J0(ax);
                double bj = J1(ax);
                for (int j = 1; j < n; j++)
                {
                    double bjp = j * tox * bj - bjm;
                    bjm = bj;
                    bj = bjp;
                }
                ans = bj;
            }
            else
            {
                // h = J_{n-1}/J_n = b_n - 1/(b_{n+1} - 1/(b_{n+2} - ...)), b_k = 2k/x, by the modified Lentz method.
                const double Tiny = 1E-300;
                double h = n * tox;
                double c = h;
                double d = 0.0;
                for (int k = n + 1; k < n + 100_000; k++)
                {
                    double b = k * tox;
                    d = b - d;
                    if (d == 0.0) d = Tiny;
                    c = b - 1.0 / c;
                    if (c == 0.0) c = Tiny;
                    d = 1.0 / d;
                    double delta = c * d;
                    h *= delta;
                    if (Math.Abs(delta - 1.0) < 1E-16) break;
                }

                // Backward recurrence J_{k-1} = (2k/x) J_k - J_{k+1} from J_n = 1, J_{n-1} = h.
                double scale = 1.0;
                double upper = 1.0;
                double lower = h;
                for (int k = n - 1; k >= 1; k--)
                {
                    double next = k * tox * lower - upper;
                    upper = lower;
                    lower = next;
                    if (Math.Abs(lower) > BIGNO)
                    {
                        lower *= BIGNI;
                        upper *= BIGNI;
                        scale *= BIGNI;
                    }
                }

                // lower and upper now hold J_0 and J_1 on the scale where J_n = scale.
                double j0 = J0(ax);
                double j1 = J1(ax);
                ans = Math.Abs(j0) >= Math.Abs(j1) ? scale * (j0 / lower) : scale * (j1 / upper);
            }

            return (x < 0.0 && (n & 1) != 0) ? -ans : ans;
        }

        #endregion

        #region Bessel Functions of the Second Kind (Y)

        /// <summary>
        /// Computes the Bessel function of the second kind, order 0: Y₀(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate Y₀. Must be positive (x &gt; 0).</param>
        /// <returns>The value of Y₀(x). Y₀(x) → -∞ as x → 0⁺. Oscillates for large x.</returns>
        /// <remarks>
        /// <para>
        /// For x &lt;= 0.5, Y₀(x) = R(x²) + (2/π) ln(x) J₀(x), with R a 7-term Chebyshev expansion in x², the
        /// form of Cephes y0 (Moshier, 1989). On (0.5, 8] the range is split at 1.5, 3 and 5.5, and on each piece
        /// that holds a zero y_k of Y₀ the expansion is of Y₀(x) / (x - y_k), with y_k stored as the sum of two
        /// doubles, so the relative error stays small next to the zeros 0.8936, 3.9577 and 7.0861. For x &gt; 8,
        /// Y₀(x) = (P (sin x - cos x) + Q (sin x + cos x)) / sqrt(πx), with Hankel's P and Q from Chebyshev
        /// expansions in 8/x up to x = 25 and from Hankel's asymptotic series (Abramowitz and Stegun 9.2.9 and
        /// 9.2.10) above. The coefficients were computed with mpmath at 50 digits by
        /// docs/mathematics/oracles/special-functions/bessel_chebyshev_coefficients.py.
        /// </para>
        /// <para>
        /// Measured against mpmath, the relative error is about 1E-15 for x &lt;= 8, including within 1E-6 of
        /// the zeros, and the error above 8 is about 1E-15 of the envelope sqrt(2/(πx)), as for scipy.special.y0.
        /// </para>
        /// <para>
        /// Y₀(x) is the second linearly independent solution of Bessel's equation of order 0.
        /// Unlike J₀, it is singular at x = 0. For large x: Y₀(x) ~ sqrt(2/(pi*x)) sin(x - pi/4).
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if x &lt;= 0.</exception>
        public static double Y0(double x)
        {
            if (x <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(x), "The argument x must be positive for Y0.");

            if (x <= 0.5)
                return Evaluate.Chebyshev(Y0Log, 16.0 * x * x - 2.0) + TwoOverPi * Math.Log(x) * J0(x);
            if (x <= 1.5)
                return (x - Y0Zero1Hi - Y0Zero1Lo) * Evaluate.Chebyshev(Y0Zero1, 4.0 * x - 4.0);
            if (x <= 3.0)
                return Evaluate.Chebyshev(Y0Mid, (4.0 * x - 9.0) / 1.5);
            if (x <= 5.5)
                return (x - Y0Zero2Hi - Y0Zero2Lo) * Evaluate.Chebyshev(Y0Zero2, (4.0 * x - 17.0) / 2.5);
            if (x <= 8.0)
                return (x - Y0Zero3Hi - Y0Zero3Lo) * Evaluate.Chebyshev(Y0Zero3, (4.0 * x - 27.0) / 2.5);

            HankelPQ(0, x, out double p, out double q);
            double s = Math.Sin(x), c = Math.Cos(x);
            return (p * (s - c) + q * (s + c)) / Math.Sqrt(Math.PI * x);
        }

        /// <summary>
        /// Computes the Bessel function of the second kind, order 1: Y₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate Y₁. Must be positive (x &gt; 0).</param>
        /// <returns>The value of Y₁(x). Y₁(x) → -∞ as x → 0⁺. Oscillates for large x.</returns>
        /// <remarks>
        /// <para>
        /// For x &lt;= 1.5, Y₁(x) = x R(x²) + (2/π)(ln(x) J₁(x) - 1/x), with R a 9-term Chebyshev expansion in x²,
        /// the form of Cephes y1 (Moshier, 1989). On (1.5, 3.5] and (3.5, 8] the expansions are of Y₁(x) / (x - y_k)
        /// for the zeros y₁ = 2.1971 and y₂ = 5.4297, each stored as the sum of two doubles, so the relative error
        /// stays small next to them. For x &gt; 8, Y₁(x) = (Q (sin x - cos x) - P (sin x + cos x)) / sqrt(πx) with
        /// Hankel's P and Q of order 1. The coefficients were computed with mpmath at 50 digits.
        /// </para>
        /// <para>
        /// Measured against mpmath, the relative error is about 1E-15 for x &lt;= 8, including next to the zeros,
        /// and the error above 8 is about 1E-15 of the envelope sqrt(2/(πx)), as for scipy.special.y1.
        /// </para>
        /// <para>
        /// For large x: Y₁(x) ~ sqrt(2/(pi*x)) sin(x - 3*pi/4).
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if x &lt;= 0.</exception>
        public static double Y1(double x)
        {
            if (x <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(x), "The argument x must be positive for Y1.");

            if (x <= 1.5)
                return x * Evaluate.Chebyshev(Y1Log, x * x / 0.5625 - 2.0) + TwoOverPi * (J1(x) * Math.Log(x) - 1.0 / x);
            if (x <= 3.5)
                return (x - Y1Zero1Hi - Y1Zero1Lo) * Evaluate.Chebyshev(Y1Zero1, 2.0 * x - 5.0);
            if (x <= 8.0)
                return (x - Y1Zero2Hi - Y1Zero2Lo) * Evaluate.Chebyshev(Y1Zero2, (4.0 * x - 23.0) / 4.5);

            HankelPQ(1, x, out double p, out double q);
            double s = Math.Sin(x), c = Math.Cos(x);
            return (q * (s - c) - p * (s + c)) / Math.Sqrt(Math.PI * x);
        }

        /// <summary>
        /// Computes the Bessel function of the second kind for integer order n: Yₙ(x).
        /// </summary>
        /// <param name="n">The order of the function. Must be non-negative (n &gt;= 0).</param>
        /// <param name="x">The argument at which to evaluate Yₙ. Must be positive (x &gt; 0).</param>
        /// <returns>The value of Yₙ(x). Singular at x = 0.</returns>
        /// <remarks>
        /// <para>
        /// For n = 0 and n = 1, delegates to <see cref="Y0(double)"/> and <see cref="Y1(double)"/> respectively.
        /// For n &gt;= 2, uses the forward recurrence relation, which is numerically stable for Y:
        /// <code>
        ///     Y_{n+1}(x) = (2n/x) Y_n(x) - Y_{n-1}(x)
        /// </code>
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if n &lt; 0 or x &lt;= 0.</exception>
        public static double Yn(int n, double x)
        {
            if (n < 0)
                throw new ArgumentOutOfRangeException(nameof(n), "The order n must be non-negative.");
            if (x <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(x), "The argument x must be positive for Yn.");
            if (n == 0) return Y0(x);
            if (n == 1) return Y1(x);

            double tox = 2.0 / x;
            double bym = Y0(x);
            double by = Y1(x);

            for (int j = 1; j < n; j++)
            {
                double byp = j * tox * by - bym;
                bym = by;
                by = byp;
            }

            return by;
        }

        #endregion

        #region Chebyshev Expansions

        // The expansions below are summed by Evaluate.Chebyshev, order-zero term first. They
        // were computed with mpmath at 50 digits by
        // docs/mathematics/oracles/special-functions/bessel_chebyshev_coefficients.py, which reproduces them.

        /// <summary>The constant 2/π.</summary>
        private const double TwoOverPi = 0.63661977236758134;

        /// <summary>Leading double of zero 1 of J₀.</summary>
        private const double J0Zero1Hi = 2.404825557695773;

        /// <summary>Trailing double of zero 1 of J₀: the zero is J0Zero1Hi + J0Zero1Lo to about 32 digits.</summary>
        private const double J0Zero1Lo = -1.176691651530894E-16;

        /// <summary>Leading double of zero 2 of J₀.</summary>
        private const double J0Zero2Hi = 5.520078110286311;

        /// <summary>Trailing double of zero 2 of J₀: the zero is J0Zero2Hi + J0Zero2Lo to about 32 digits.</summary>
        private const double J0Zero2Lo = 8.088597146146722E-17;

        /// <summary>Leading double of zero 1 of J₁.</summary>
        private const double J1Zero1Hi = 3.8317059702075125;

        /// <summary>Trailing double of zero 1 of J₁: the zero is J1Zero1Hi + J1Zero1Lo to about 32 digits.</summary>
        private const double J1Zero1Lo = -1.5269184090088067E-16;

        /// <summary>Leading double of zero 2 of J₁.</summary>
        private const double J1Zero2Hi = 7.015586669815619;

        /// <summary>Trailing double of zero 2 of J₁: the zero is J1Zero2Hi + J1Zero2Lo to about 32 digits.</summary>
        private const double J1Zero2Lo = -9.414165653410389E-17;

        /// <summary>Leading double of zero 1 of Y₀.</summary>
        private const double Y0Zero1Hi = 0.8935769662791675;

        /// <summary>Trailing double of zero 1 of Y₀: the zero is Y0Zero1Hi + Y0Zero1Lo to about 32 digits.</summary>
        private const double Y0Zero1Lo = 2.6596231539720385E-17;

        /// <summary>Leading double of zero 2 of Y₀.</summary>
        private const double Y0Zero2Hi = 3.957678419314858;

        /// <summary>Trailing double of zero 2 of Y₀: the zero is Y0Zero2Hi + Y0Zero2Lo to about 32 digits.</summary>
        private const double Y0Zero2Lo = -1.0764340697562706E-16;

        /// <summary>Leading double of zero 3 of Y₀.</summary>
        private const double Y0Zero3Hi = 7.086051060301773;

        /// <summary>Trailing double of zero 3 of Y₀: the zero is Y0Zero3Hi + Y0Zero3Lo to about 32 digits.</summary>
        private const double Y0Zero3Lo = -8.835285723085408E-17;

        /// <summary>Leading double of zero 1 of Y₁.</summary>
        private const double Y1Zero1Hi = 2.197141326031017;

        /// <summary>Trailing double of zero 1 of Y₁: the zero is Y1Zero1Hi + Y1Zero1Lo to about 32 digits.</summary>
        private const double Y1Zero1Lo = -4.8259835876454966E-17;

        /// <summary>Leading double of zero 2 of Y₁.</summary>
        private const double Y1Zero2Hi = 5.429681040794135;

        /// <summary>Trailing double of zero 2 of Y₁: the zero is Y1Zero2Hi + Y1Zero2Lo to about 32 digits.</summary>
        private const double Y1Zero2Lo = 4.162514026670377E-16;

        /// <summary>J₀(x) / ((x² - j₀,₁²)(x² - j₀,₂²)) in z = x² on [0, 64].</summary>
        private static readonly double[] J0Small =
        {
            0.003980711000320297, -0.002620939138497025, 0.0008687721472062577,
            -0.00017040584951557357, 2.204908639454817E-05, -2.0248894980012834E-06,
            1.3898025173428905E-07, -7.406701254344167E-09, 3.1559927583611334E-10,
            -1.1003860827939697E-11, 3.1989232970233477E-13, -7.875041254709154E-15,
            1.6632588015211278E-16, -3.0475993756698376E-18, 4.8912290055811176E-20
        };

        /// <summary>J₁(x) / (x (x² - j₁,₁²)(x² - j₁,₂²)) in z = x² on [0, 64].</summary>
        private static readonly double[] J1Small =
        {
            0.0005579141977522621, -0.0003115684036643562, 8.553646588540538E-05,
            -1.4149356425623485E-05, 1.576470099225889E-06, -1.269075502369425E-07,
            7.747614077073057E-09, -3.716577146675772E-10, 1.4395417669312554E-11,
            -4.60007822413452E-13, 1.2341465589619718E-14, -2.8205466239524856E-16,
            5.558821065951572E-18, -9.546930865755219E-20, 1.4418291781315051E-21
        };

        /// <summary>Y₀(x) - (2/π) ln(x) J₀(x) in z = x² on [0, 0.25].</summary>
        private static readonly double[] Y0Log =
        {
            -0.1039553136934708, 0.02170236926617634, -0.00012400802744650475,
            2.6068033752962025E-07, -2.877948821859841E-10, 1.9647455203887995E-13,
            -9.127739625832616E-17
        };

        /// <summary>Y₀(x) / (x - y₀,₁) on [0.5, 1.5].</summary>
        private static readonly double[] Y0Zero1 =
        {
            1.7064582474256436, -0.24289162163307515, 0.02529773907170474,
            -0.0061511838930721095, 0.001442379524001996, -0.0003231933673867555,
            7.450008557158439E-05, -1.7531902451562595E-05, 4.185154006441054E-06,
            -1.010751529886327E-06, 2.46458968028983E-07, -6.057902497774891E-08,
            1.499139309290079E-08, -3.731479995658509E-09, 9.334707144988312E-10,
            -2.345444132620023E-10, 5.915988768647463E-11, -1.497333187899378E-11,
            3.801363955830564E-12, -9.67732674294546E-13, 2.4697376011228824E-13,
            -6.317221253288305E-14, 1.6191786126889323E-14, -4.157969341537637E-15,
            1.0695947958980984E-15, -2.755815382404065E-16, 7.11085386178469E-17,
            -1.8373325763192018E-17, 4.753424873541765E-18, -1.2312350823963516E-18
        };

        /// <summary>Y₀(x) on [1.5, 3].</summary>
        private static readonly double[] Y0Mid =
        {
            0.8991279610481316, -0.0072251978939073735, -0.07020068643706694,
            0.004410992111554405, 0.0002932436587864384, 1.3924633463273487E-05,
            -6.958599639988895E-06, 8.063758963881566E-07, -1.0943386372417973E-07,
            1.744087303008207E-08, -2.73757633923206E-09, 4.2860069489262145E-10,
            -6.766764218521532E-11, 1.0752413845846688E-11, -1.7173941476264612E-12,
            2.7556392294753683E-13, -4.4395740282392144E-14, 7.178473295041094E-15,
            -1.1644728006474257E-15, 1.8944870635738637E-16, -3.0902751336250493E-17,
            5.0528956206889594E-18, -8.279979584278538E-19
        };

        /// <summary>Y₀(x) / (x - y₀,₂) on [3, 5.5].</summary>
        private static readonly double[] Y0Zero2 =
        {
            -0.6886043790410038, 0.09089304683356025, 0.038022895048569844,
            -0.004245463228750419, -0.0005323017083957325, 4.962228357728892E-05,
            4.691432667256341E-06, -4.131964775294702E-07, -1.1484541798716266E-08,
            3.0928640440108364E-10, 2.2118553460788837E-10, -2.633701414837581E-11,
            2.98087307134756E-12, -4.3097295090885396E-13, 6.231203852035348E-14,
            -8.812518466205934E-15, 1.2507999699440095E-15, -1.7835333897144863E-16,
            2.5493866690739737E-17, -3.652086093952828E-18, 5.242853413923317E-19
        };

        /// <summary>Y₀(x) / (x - y₀,₃) on [5.5, 8].</summary>
        private static readonly double[] Y0Zero3 =
        {
            0.5296440696721009, 0.01517732810674114, -0.03618847560068274,
            9.288480506124671E-05, 0.0006762432410947991, -7.631168798553423E-06,
            -5.773053117312353E-06, 7.47985919873368E-08, 2.900442850206552E-08,
            -3.9291184116249246E-10, -9.446769451172849E-11, 1.1680346125267796E-12,
            2.277031977771271E-13, -3.2058214322535506E-15, -3.4690235689913347E-16,
            8.307611094390935E-19, 8.548377883989016E-19
        };

        /// <summary>(Y₁(x) - (2/π)(ln(x) J₁(x) - 1/x)) / x in z = x² on [0, 2.25].</summary>
        private static readonly double[] Y1Log =
        {
            -0.28055418641072555, 0.05403326126873243, -0.0017225901354540915,
            2.3956137315801513E-05, -1.8893203495591734E-07, 9.638814729253337E-10,
            -3.448274919187277E-12, 9.138761360821853E-15, -1.8673559354499098E-17
        };

        /// <summary>Y₁(x) / (x - y₁,₁) on [1.5, 3.5].</summary>
        private static readonly double[] Y1Zero1 =
        {
            0.9335979448957102, -0.13825807742548832, -0.014353805147792643,
            7.657866298969309E-05, 0.0006677202885425471, -0.00010938293859302674,
            1.860782644365838E-05, -4.113077401751908E-06, 8.83828033983229E-07,
            -1.8505988241531371E-07, 3.87643296680694E-08, -8.121444873930584E-09,
            1.699970944961634E-09, -3.556098188198514E-10, 7.43562864479522E-11,
            -1.5542381724315902E-11, 3.2479206876871926E-12, -6.785866442626836E-13,
            1.4175397281607625E-13, -2.9607984767732935E-14, 6.1835321903082136E-15,
            -1.2912987820633445E-15, 2.696408148566485E-16, -5.630129935672028E-17,
            1.1755178288076227E-17, -2.4542645695332483E-18, 5.123864317192076E-19
        };

        /// <summary>Y₁(x) / (x - y₁,₂) on [3.5, 8].</summary>
        private static readonly double[] Y1Zero2 =
        {
            -0.4520525274382522, 0.0903467035205593, 0.09395085473448606,
            -0.015493318422092322, -0.005092916759693723, 0.0007001087578538037,
            0.0001390897220923952, -1.6947654678076103E-05, -2.0323356852794283E-06,
            2.047987384297512E-07, 2.7304818500930765E-08, -3.159090049366525E-09,
            4.29713547384621E-11, -2.868404405873767E-11, 9.80574970899084E-12,
            -1.89522068070456E-12, 3.69220900540239E-13, -7.615029127595644E-14,
            1.56596638527455E-14, -3.2009908285096767E-15, 6.541475128490997E-16,
            -1.3367935946684902E-16, 2.7308293843618197E-17, -5.576827410389544E-18,
            1.1386041563580327E-18
        };

        /// <summary>Hankel's P0(x) in u = 8/x on [8, 25].</summary>
        private static readonly double[] P0Mid =
        {
            1.9989337250663841, -0.00048106486454967687, -5.957179199746459E-05,
            5.999912872020045E-07, 2.783056959969125E-08, -1.363509943981793E-09,
            -1.8653701506557085E-12, 3.3835667650797142E-12, -1.7970762324282873E-13,
            -1.4235537018948742E-15, 9.392644877101111E-16, -7.208382592343428E-17,
            1.0358205209104134E-18
        };

        /// <summary>Hankel's Q0(x) / u in u = 8/x on [8, 25].</summary>
        private static readonly double[] Q0Mid =
        {
            -0.03111323728974103, 6.125509272453135E-05, 7.315586031435822E-06,
            -1.4062095683918164E-07, -5.40790575850404E-09, 4.2963779389537513E-10,
            -4.832309037010757E-12, -1.1397631791404546E-12, 9.22602912859179E-14,
            -1.3741509207543389E-15, -3.990247056333184E-16, 4.616928146114218E-17,
            -1.922037696497648E-18, -1.4756824662818317E-19, 3.4272832568568806E-20
        };

        /// <summary>Hankel's P1(x) in u = 8/x on [8, 25].</summary>
        private static readonly double[] P1Mid =
        {
            2.001783820884021, 0.0008062650839909221, 0.00010074542735820923,
            -7.81826306394218E-07, -3.738330500932916E-08, 1.6509849421820302E-09,
            5.2102839973974545E-12, -4.005006474342676E-12, 2.0055166459337564E-13,
            2.1579104893942667E-15, -1.0671418689437694E-15, 7.852630015395434E-17,
            -9.178278825017994E-19
        };

        /// <summary>Hankel's Q1(x) / u in u = 8/x on [8, 25].</summary>
        private static readonly double[] Q1Mid =
        {
            0.09355777472904885, -8.626193692862056E-05, -1.0402551855536108E-05,
            1.7419468395008437E-07, 6.946605321526864E-09, -5.079820942147065E-10,
            4.782588984313255E-12, 1.3358061697190562E-12, -1.0292787750706666E-13,
            1.3023011334169809E-15, 4.543597568903915E-16, -5.049392984781845E-17,
            1.9892163883922136E-18, 1.7051784936992025E-19
        };

        /// <summary>K₀(x) + ln(x/2) I₀(x) in z = x² on [0, 4].</summary>
        private static readonly double[] K0Small =
        {
            -0.5353273932339028, 0.3442898999246285, 0.0359799365153615,
            0.001264615411446926, 2.286212103119452E-05, 2.5347910790261494E-07,
            1.904516377220209E-09, 1.0349695257633625E-11, 4.2598161427910826E-14,
            1.3744654358807508E-16
        };

        /// <summary>exp(x) sqrt(x) K₀(x) in t = 4/x - 1 on [2, ∞).</summary>
        private static readonly double[] K0Large =
        {
            2.4403030820659555, -0.0314481013119645, 0.0015698838857300533,
            -0.00012849549581627802, 1.39498137188765E-05, -1.8317555227191195E-06,
            2.766813639445015E-07, -4.660489897687948E-08, 8.574034017414225E-09,
            -1.6975345093890614E-09, 3.5773972814003283E-10, -7.957489244477396E-11,
            1.8559491149549264E-11, -4.514597883374519E-12, 1.1403405882073441E-12,
            -2.9800969231481784E-13, 8.032890775068375E-14, -2.2275133267462965E-14,
            6.340076476276646E-15, -1.848593377920907E-15, 5.5120559994043335E-16,
            -1.6782311257549006E-16, 5.2103917776435543E-17, -1.6475805939842632E-17,
            5.3004337711773354E-18, -1.7331712005821001E-18
        };

        /// <summary>x (K₁(x) - ln(x/2) I₁(x)) in z = x² on [0, 4].</summary>
        private static readonly double[] K1Small =
        {
            1.5253002273389478, -0.3531559607765449, -0.12261118082265715,
            -0.006975723859639864, -0.0001730288957513052, -2.4334061415659684E-06,
            -2.213387630734726E-08, -1.4114883926335278E-10, -6.666901694199329E-13,
            -2.427449850519366E-15, -7.023863479386288E-18
        };

        /// <summary>exp(x) sqrt(x) K₁(x) in t = 4/x - 1 on [2, ∞).</summary>
        private static readonly double[] K1Large =
        {
            2.7206261904844427, 0.10392373657681724, -0.002857816859622779,
            0.00019521551847135162, -1.936197974166083E-05, 2.406484947837217E-06,
            -3.5019606030878126E-07, 5.7410841254500495E-08, -1.0345762465678097E-08,
            2.0150497551970347E-09, -4.1903547593419254E-10, 9.218315187605315E-11,
            -2.129967838427791E-11, 5.139639673482343E-12, -1.2891739609498229E-12,
            3.348419666052243E-13, -8.976705182010146E-14, 2.4771544242195988E-14,
            -7.0198370892147685E-15, 2.038703166239861E-15, -6.057047270643018E-16,
            1.8380935752430455E-16, -5.689462849193648E-17, 1.7940510478863572E-17,
            -5.7567444820733025E-18, 1.8778651901623268E-18
        };

        #endregion

        #region List Overloads

        /// <summary>
        /// Computes the modified Bessel function of the first kind of order 0 for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="I0(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="I0(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] I0(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = I0(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the modified Bessel function of the first kind of order 1 for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="I1(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="I1(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] I1(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = I1(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the exponentially scaled modified Bessel function exp(-|x|) I0(x) for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="I0e(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="I0e(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] I0e(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = I0e(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the exponentially scaled modified Bessel function exp(-|x|) I1(x) for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="I1e(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="I1e(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] I1e(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = I1e(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the modified Bessel function of the first kind of integer order n for each value of a list.
        /// </summary>
        /// <param name="n">The order, as in <see cref="In(int, double)"/>, used for every value.</param>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="In(int, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="In(int, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] In(int n, IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = In(n, x[i]);
            return result;
        }

        /// <summary>
        /// Computes the modified Bessel function of the second kind of order 0 for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="K0(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="K0(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] K0(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = K0(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the modified Bessel function of the second kind of order 1 for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="K1(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="K1(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] K1(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = K1(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the modified Bessel function of the second kind of integer order n for each value of a list.
        /// </summary>
        /// <param name="n">The order, as in <see cref="Kn(int, double)"/>, used for every value.</param>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="Kn(int, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Kn(int, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Kn(int n, IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Kn(n, x[i]);
            return result;
        }

        /// <summary>
        /// Computes the Bessel function of the first kind of order 0 for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="J0(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="J0(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] J0(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = J0(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the Bessel function of the first kind of order 1 for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="J1(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="J1(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] J1(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = J1(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the Bessel function of the first kind of integer order n for each value of a list.
        /// </summary>
        /// <param name="n">The order, as in <see cref="Jn(int, double)"/>, used for every value.</param>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="Jn(int, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Jn(int, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Jn(int n, IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Jn(n, x[i]);
            return result;
        }

        /// <summary>
        /// Computes the Bessel function of the second kind of order 0 for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="Y0(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Y0(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Y0(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Y0(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the Bessel function of the second kind of order 1 for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="Y1(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Y1(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Y1(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Y1(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the Bessel function of the second kind of integer order n for each value of a list.
        /// </summary>
        /// <param name="n">The order, as in <see cref="Yn(int, double)"/>, used for every value.</param>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="Yn(int, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Yn(int, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Yn(int n, IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Yn(n, x[i]);
            return result;
        }

        #endregion

    }
}
