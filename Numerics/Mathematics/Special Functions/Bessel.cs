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
        /// Uses polynomial approximations from Abramowitz and Stegun (1972), sections 9.8.1 and 9.8.2.
        /// For |x| &lt;= 3.75, uses a polynomial in (x/3.75)² with maximum error |e| &lt; 1.6x10⁻⁷.
        /// For |x| &gt; 3.75, uses an asymptotic expansion with maximum error |e| &lt; 1.9x10⁻⁷.
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
            if (ax < 3.75)
            {
                double t = x / 3.75;
                t *= t;
                return 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492
                    + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
            }
            else
            {
                if (ax > 709) return double.PositiveInfinity;
                double t = 3.75 / ax;
                return (Math.Exp(ax) / Math.Sqrt(ax)) * (0.39894228 + t * (0.01328592
                    + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281 + t * (-0.02057706
                    + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377))))))));
            }
        }

        /// <summary>
        /// Computes the modified Bessel function of the first kind, order 1: I₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate I₁. Can be any real number.</param>
        /// <returns>The value of I₁(x). I₁(-x) = -I₁(x) (odd function).</returns>
        /// <remarks>
        /// <para>
        /// Uses polynomial approximations from Abramowitz and Stegun (1972), sections 9.8.3 and 9.8.4.
        /// For |x| &lt;= 3.75, uses a polynomial in (x/3.75)² with maximum error |e| &lt; 8x10⁻⁹.
        /// For |x| &gt; 3.75, uses an asymptotic expansion with maximum error |e| &lt; 2.2x10⁻⁷.
        /// </para>
        /// <para>
        /// I₁(x) appears in circular statistics: the mean resultant length of the von Mises distribution
        /// is A(kappa) = I₁(kappa)/I₀(kappa), and the circular variance is V = 1 - A(kappa).
        /// </para>
        /// </remarks>
        public static double I1(double x)
        {
            double ax = Math.Abs(x);
            double result;
            if (ax < 3.75)
            {
                double t = x / 3.75;
                t *= t;
                result = ax * (0.5 + t * (0.87890594 + t * (0.51498869 + t * (0.15084934
                    + t * (0.02658733 + t * (0.00301532 + t * 0.00032411))))));
            }
            else
            {
                double t = 3.75 / ax;
                result = (Math.Exp(ax) / Math.Sqrt(ax)) * (0.39894228 + t * (-0.03988024
                    + t * (-0.00362018 + t * (0.00163801 + t * (-0.01031555 + t * (0.02282967
                    + t * (-0.02895312 + t * (0.01787654 + t * (-0.00420059)))))))));
            }
            return x < 0 ? -result : result;
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
        /// For n = 0 and n = 1, delegates to <see cref="I0"/> and <see cref="I1"/> respectively.
        /// For n &gt;= 2, uses Miller's downward recurrence algorithm, which is numerically stable:
        /// <code>
        ///     I_{n-1}(x) = (2n/x) I_n(x) + I_{n+1}(x)
        /// </code>
        /// The result is normalized against the known value of I₀(x).
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if n &lt; 0.</exception>
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

            // Starting index for downward recurrence
            int m = 2 * (n + (int)Math.Sqrt(IACC * n));

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
        /// Uses polynomial approximations from Abramowitz and Stegun (1972), sections 9.8.5 and 9.8.6.
        /// For 0 &lt; x &lt;= 2, uses: K₀(x) = -ln(x/2)I₀(x) + polynomial in (x/2)², with |e| &lt; 1x10⁻⁸.
        /// For x &gt; 2, uses: K₀(x) = (1/sqrt(x))exp(-x) * polynomial in (2/x), with |e| &lt; 1.9x10⁻⁷.
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
            {
                double t = x * x / 4.0; // (x/2)^2
                return -Math.Log(x / 2.0) * I0(x) +
                    (-0.57721566 + t * (0.42278420 + t * (0.23069756
                    + t * (0.03488590 + t * (0.00262698 + t * (0.00010750
                    + t * 0.00000740))))));
            }
            else
            {
                double t = 2.0 / x;
                return (Math.Exp(-x) / Math.Sqrt(x)) *
                    (1.25331414 + t * (-0.07832358 + t * (0.02189568
                    + t * (-0.01062446 + t * (0.00587872 + t * (-0.00251540
                    + t * 0.00053208))))));
            }
        }

        /// <summary>
        /// Computes the modified Bessel function of the second kind, order 1: K₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate K₁. Must be positive (x &gt; 0).</param>
        /// <returns>The value of K₁(x). Always positive, monotonically decreasing.</returns>
        /// <remarks>
        /// <para>
        /// Uses polynomial approximations from Abramowitz and Stegun (1972), sections 9.8.7 and 9.8.8.
        /// For 0 &lt; x &lt;= 2, uses: x*K₁(x) = x*ln(x/2)*I₁(x) + polynomial in (x/2)², with |e| &lt; 8x10⁻⁹.
        /// For x &gt; 2, uses: K₁(x) = (1/sqrt(x))exp(-x) * polynomial in (2/x), with |e| &lt; 2.2x10⁻⁷.
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
            {
                double t = x * x / 4.0; // (x/2)^2
                return Math.Log(x / 2.0) * I1(x) + (1.0 / x) *
                    (1.0 + t * (0.15443144 + t * (-0.67278579
                    + t * (-0.18156897 + t * (-0.01919402 + t * (-0.00110404
                    + t * (-0.00004686)))))));
            }
            else
            {
                double t = 2.0 / x;
                return (Math.Exp(-x) / Math.Sqrt(x)) *
                    (1.25331414 + t * (0.23498619 + t * (-0.03655620
                    + t * (0.01504268 + t * (-0.00780353 + t * (0.00325614
                    + t * (-0.00068245)))))));
            }
        }

        /// <summary>
        /// Computes the modified Bessel function of the second kind for integer order n: Kₙ(x).
        /// </summary>
        /// <param name="n">The order of the function. Must be non-negative (n &gt;= 0).</param>
        /// <param name="x">The argument at which to evaluate Kₙ. Must be positive (x &gt; 0).</param>
        /// <returns>The value of Kₙ(x). Always positive for x &gt; 0.</returns>
        /// <remarks>
        /// <para>
        /// For n = 0 and n = 1, delegates to <see cref="K0"/> and <see cref="K1"/> respectively.
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
        /// Computes the Bessel function of the first kind, order 0: J₀(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate J₀. Can be any real number.</param>
        /// <returns>The value of J₀(x). J₀(0) = 1. J₀(x) = J₀(-x) (even function). Range: [-0.4028, 1].</returns>
        /// <remarks>
        /// <para>
        /// Uses rational polynomial approximations for |x| &lt; 8 and Hankel's asymptotic expansion
        /// for |x| &gt;= 8, following Press et al. (2007). Accuracy is approximately double precision
        /// for |x| &lt; 8 and approximately 10⁻⁸ for the asymptotic region.
        /// </para>
        /// <para>
        /// J₀(x) is the unique solution of Bessel's equation of order 0 that is finite at x = 0.
        /// It oscillates with decreasing amplitude: J₀(x) ~ sqrt(2/(pi*x)) cos(x - pi/4) for large x.
        /// </para>
        /// </remarks>
        public static double J0(double x)
        {
            double ax = Math.Abs(x);
            if (ax < 8.0)
            {
                double y = x * x;
                double ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                    + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
                double ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                    + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
                return ans1 / ans2;
            }
            else
            {
                double z = 8.0 / ax;
                double y = z * z;
                double xx = ax - 0.785398164; // ax - pi/4
                double p = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                    + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                double q = -0.1562499995e-1 + y * (0.1430488765e-3
                    + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                    + y * (-0.934945152e-7))));
                return Math.Sqrt(0.636619772 / ax) * (Math.Cos(xx) * p - z * Math.Sin(xx) * q);
            }
        }

        /// <summary>
        /// Computes the Bessel function of the first kind, order 1: J₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate J₁. Can be any real number.</param>
        /// <returns>The value of J₁(x). J₁(0) = 0. J₁(-x) = -J₁(x) (odd function). Range: [-0.5819, 0.5819].</returns>
        /// <remarks>
        /// <para>
        /// Uses rational polynomial approximations for |x| &lt; 8 and Hankel's asymptotic expansion
        /// for |x| &gt;= 8, following Press et al. (2007). Accuracy is approximately double precision
        /// for |x| &lt; 8 and approximately 10⁻⁸ for the asymptotic region.
        /// </para>
        /// <para>
        /// J₁(x) oscillates with decreasing amplitude: J₁(x) ~ sqrt(2/(pi*x)) cos(x - 3*pi/4) for large x.
        /// </para>
        /// </remarks>
        public static double J1(double x)
        {
            double ax = Math.Abs(x);
            double ans;
            if (ax < 8.0)
            {
                double y = x * x;
                double ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
                    + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
                double ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
                    + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
                ans = ans1 / ans2;
            }
            else
            {
                double z = 8.0 / ax;
                double y = z * z;
                double xx = ax - 2.356194491; // ax - 3*pi/4
                double p = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                    + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                double q = 0.04687499995 + y * (-0.2002690873e-3
                    + y * (0.8449199096e-5 + y * (-0.88228987e-6
                    + y * 0.105787412e-6)));
                ans = Math.Sqrt(0.636619772 / ax) * (Math.Cos(xx) * p - z * Math.Sin(xx) * q);
                if (x < 0.0) ans = -ans;
            }
            return ans;
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
        /// For n = 0 and n = 1, delegates to <see cref="J0"/> and <see cref="J1"/> respectively.
        /// For n &gt;= 2, uses two strategies depending on the relative magnitudes of n and |x|:
        /// </para>
        /// <para>
        /// When |x| &gt; n, forward recurrence from J₀ and J₁ is numerically stable:
        /// <code>
        ///     J_{n+1}(x) = (2n/x) J_n(x) - J_{n-1}(x)
        /// </code>
        /// </para>
        /// <para>
        /// When |x| &lt;= n, Miller's downward recurrence is used, starting from a large order m &gt;&gt; n
        /// and recursing down. The result is normalized using the identity:
        /// <code>
        ///     J₀(x) + 2 J₂(x) + 2 J₄(x) + ... = 1
        /// </code>
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if n &lt; 0.</exception>
        public static double Jn(int n, double x)
        {
            if (n < 0 || n > 1_000_000)
                throw new ArgumentOutOfRangeException(nameof(n), "The order n must be between 0 and 1,000,000.");
            if (n == 0) return J0(x);
            if (n == 1) return J1(x);

            double ax = Math.Abs(x);
            if (ax == 0.0) return 0.0;

            double ans;
            if (ax > (double)n)
            {
                // Forward recurrence from J0, J1 (stable for x > n)
                double tox = 2.0 / ax;
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
                // Miller's downward recurrence (stable for x <= n)
                double tox = 2.0 / ax;

                // Starting index - must be even and much larger than n
                int m = 2 * ((n + (int)Math.Sqrt(IACC * n)) / 2);
                bool jsum = false;
                double bjp = 0.0;
                ans = 0.0;
                double sum = 0.0;
                double bj = 1.0;

                for (int j = m; j > 0; j--)
                {
                    double bjm = j * tox * bj - bjp;
                    bjp = bj;
                    bj = bjm;

                    // Rescale to prevent overflow
                    if (Math.Abs(bj) > BIGNO)
                    {
                        bj *= BIGNI;
                        bjp *= BIGNI;
                        ans *= BIGNI;
                        sum *= BIGNI;
                    }

                    if (jsum) sum += bj;
                    jsum = !jsum;
                    if (j == n) ans = bjp;
                }

                // Normalize: J0 + 2*J2 + 2*J4 + ... = 1
                sum = 2.0 * sum - bj;
                ans /= sum;
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
        /// Uses rational polynomial approximations for x &lt; 8 and Hankel's asymptotic expansion
        /// for x &gt;= 8, following Press et al. (2007). For x &lt; 8, the approximation involves
        /// J₀(x)·ln(x), reflecting the logarithmic singularity at the origin.
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

            if (x < 8.0)
            {
                double y = x * x;
                double ans1 = -2957821389.0 + y * (7062834065.0 + y * (-512359803.6
                    + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
                double ans2 = 40076544269.0 + y * (745249964.8 + y * (7189466.438
                    + y * (47447.26470 + y * (226.1030244 + y * 1.0))));
                return (ans1 / ans2) + 0.636619772 * J0(x) * Math.Log(x);
            }
            else
            {
                double z = 8.0 / x;
                double y = z * z;
                double xx = x - 0.785398164; // x - pi/4
                double p = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                    + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                double q = -0.1562499995e-1 + y * (0.1430488765e-3
                    + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                    + y * (-0.934945152e-7))));
                return Math.Sqrt(0.636619772 / x) * (Math.Sin(xx) * p + z * Math.Cos(xx) * q);
            }
        }

        /// <summary>
        /// Computes the Bessel function of the second kind, order 1: Y₁(x).
        /// </summary>
        /// <param name="x">The argument at which to evaluate Y₁. Must be positive (x &gt; 0).</param>
        /// <returns>The value of Y₁(x). Y₁(x) → -∞ as x → 0⁺. Oscillates for large x.</returns>
        /// <remarks>
        /// <para>
        /// Uses rational polynomial approximations for x &lt; 8 and Hankel's asymptotic expansion
        /// for x &gt;= 8, following Press et al. (2007). For x &lt; 8, the approximation involves
        /// J₁(x)·ln(x) - 1/x, reflecting the singular behavior at the origin.
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

            if (x < 8.0)
            {
                double y = x * x;
                double ans1 = x * (-4900604943000.0 + y * (1275274390000.0
                    + y * (-51534381390.0 + y * (734926455.1
                    + y * (-4237922.726 + y * 8511.937935)))));
                double ans2 = 24995805700000.0 + y * (424441966400.0
                    + y * (3733650367.0 + y * (22459040.02
                    + y * (102042.6050 + y * (354.9632885 + y * 1.0)))));
                return (ans1 / ans2) + 0.636619772 * (J1(x) * Math.Log(x) - 1.0 / x);
            }
            else
            {
                double z = 8.0 / x;
                double y = z * z;
                double xx = x - 2.356194491; // x - 3*pi/4
                double p = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                    + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                double q = 0.04687499995 + y * (-0.2002690873e-3
                    + y * (0.8449199096e-5 + y * (-0.88228987e-6
                    + y * 0.105787412e-6)));
                return Math.Sqrt(0.636619772 / x) * (Math.Sin(xx) * p + z * Math.Cos(xx) * q);
            }
        }

        /// <summary>
        /// Computes the Bessel function of the second kind for integer order n: Yₙ(x).
        /// </summary>
        /// <param name="n">The order of the function. Must be non-negative (n &gt;= 0).</param>
        /// <param name="x">The argument at which to evaluate Yₙ. Must be positive (x &gt; 0).</param>
        /// <returns>The value of Yₙ(x). Singular at x = 0.</returns>
        /// <remarks>
        /// <para>
        /// For n = 0 and n = 1, delegates to <see cref="Y0"/> and <see cref="Y1"/> respectively.
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

    }
}
