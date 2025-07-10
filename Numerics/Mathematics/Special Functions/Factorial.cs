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

namespace Numerics.Mathematics.SpecialFunctions
{

    /// <summary>
    /// Factorial functions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// </para>
    /// <b> Description: </b>
    /// A factorial is a non negative integer n, denoted by n!, that is the product of all positive integers 
    /// less than or equal to n. For example, 5! equals (1 * 2 * 3 * 4 * 5) which is 120. It is also important to note
    /// that 0! equals 1.
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Factorial" />
    /// </description></item>
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Binomial_coefficient" />
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class Factorial
    {
        private const int FactorialMaxArgument = 170;
        private static readonly double[] _factorialCache = new double[171] { 1d, 1d, 2d, 6d, 24d, 120d, 720d, 5040d, 40320d, 362880d, 3628800d, 39916800d, 479001600d, 6227020800d, 87178291200d, 1307674368000d, 20922789888000d, 355687428096000d, 6402373705728000d, 1.21645100408832E+17d, 2.43290200817664E+18d, 5.109094217170944E+19d, 1.1240007277776077E+21d, 2.5852016738884978E+22d, 6.2044840173323941E+23d, 1.5511210043330986E+25d, 4.0329146112660565E+26d, 1.0888869450418352E+28d, 3.0488834461171384E+29d, 8.8417619937397008E+30d, 2.6525285981219103E+32d, 8.2228386541779224E+33d, 2.6313083693369352E+35d, 8.6833176188118859E+36d, 2.9523279903960412E+38d, 1.0333147966386144E+40d, 3.7199332678990118E+41d, 1.3763753091226343E+43d, 5.23022617466601E+44d, 2.0397882081197442E+46d, 8.1591528324789768E+47d, 3.3452526613163803E+49d, 1.4050061177528798E+51d, 6.0415263063373834E+52d, 2.6582715747884485E+54d, 1.1962222086548019E+56d, 5.5026221598120885E+57d, 2.5862324151116818E+59d, 1.2413915592536073E+61d, 6.0828186403426752E+62d, 3.0414093201713376E+64d, 1.5511187532873822E+66d, 8.0658175170943877E+67d, 4.2748832840600255E+69d, 2.3084369733924138E+71d, 1.2696403353658276E+73d, 7.1099858780486348E+74d, 4.0526919504877221E+76d, 2.3505613312828789E+78d, 1.3868311854568987E+80d, 8.3209871127413916E+81d, 5.0758021387722484E+83d, 3.1469973260387939E+85d, 1.98260831540444E+87d, 1.2688693218588417E+89d, 8.2476505920824715E+90d, 5.4434493907744307E+92d, 3.6471110918188683E+94d, 2.4800355424368305E+96d, 1.711224524281413E+98d, 1.197857166996989E+100d, 8.5047858856786218E+101d, 6.1234458376886077E+103d, 4.4701154615126834E+105d, 3.3078854415193856E+107d, 2.4809140811395391E+109d, 1.8854947016660499E+111d, 1.4518309202828584E+113d, 1.1324281178206295E+115d, 8.9461821307829729E+116d, 7.1569457046263779E+118d, 5.7971260207473655E+120d, 4.75364333701284E+122d, 3.9455239697206569E+124d, 3.314240134565352E+126d, 2.8171041143805494E+128d, 2.4227095383672724E+130d, 2.1077572983795269E+132d, 1.8548264225739836E+134d, 1.6507955160908453E+136d, 1.4857159644817607E+138d, 1.3520015276784023E+140d, 1.24384140546413E+142d, 1.1567725070816409E+144d, 1.0873661566567424E+146d, 1.0329978488239052E+148d, 9.916779348709491E+149d, 9.6192759682482062E+151d, 9.426890448883242E+153d, 9.33262154439441E+155d, 9.33262154439441E+157d, 9.4259477598383536E+159d, 9.6144667150351211E+161d, 9.9029007164861754E+163d, 1.0299016745145622E+166d, 1.0813967582402904E+168d, 1.1462805637347078E+170d, 1.2265202031961373E+172d, 1.3246418194518284E+174d, 1.4438595832024928E+176d, 1.5882455415227421E+178d, 1.7629525510902437E+180d, 1.9745068572210728E+182d, 2.2311927486598123E+184d, 2.5435597334721862E+186d, 2.9250936934930141E+188d, 3.3931086844518965E+190d, 3.969937160808719E+192d, 4.6845258497542883E+194d, 5.5745857612076033E+196d, 6.6895029134491239E+198d, 8.09429852527344E+200d, 9.8750442008335976E+202d, 1.2146304367025325E+205d, 1.5061417415111404E+207d, 1.8826771768889254E+209d, 2.3721732428800459E+211d, 3.0126600184576582E+213d, 3.8562048236258025E+215d, 4.9745042224772855E+217d, 6.4668554892204716E+219d, 8.4715806908788174E+221d, 1.1182486511960039E+224d, 1.4872707060906852E+226d, 1.9929427461615181E+228d, 2.6904727073180495E+230d, 3.6590428819525472E+232d, 5.01288874827499E+234d, 6.9177864726194859E+236d, 9.6157231969410859E+238d, 1.346201247571752E+241d, 1.89814375907617E+243d, 2.6953641378881614E+245d, 3.8543707171800706E+247d, 5.5502938327393013E+249d, 8.0479260574719866E+251d, 1.17499720439091E+254d, 1.7272458904546376E+256d, 2.5563239178728637E+258d, 3.8089226376305671E+260d, 5.7133839564458505E+262d, 8.6272097742332346E+264d, 1.3113358856834518E+267d, 2.0063439050956811E+269d, 3.0897696138473489E+271d, 4.7891429014633912E+273d, 7.47106292628289E+275d, 1.1729568794264138E+278d, 1.8532718694937338E+280d, 2.9467022724950369E+282d, 4.714723635992059E+284d, 7.5907050539472148E+286d, 1.2296942187394488E+289d, 2.0044015765453015E+291d, 3.2872185855342945E+293d, 5.423910666131586E+295d, 9.0036917057784329E+297d, 1.5036165148649983E+300d, 2.5260757449731969E+302d, 4.2690680090047027E+304d, 7.257415615307994E+306d };


        /// <summary>
        /// Computes the factorial of an integer.
        /// </summary>
        /// <param name="x">An integer greater than 0 </param>
        /// <returns>
        /// The product of all integers less than or equal to the given x
        /// </returns>
        public static double Function(int x)
        {
            if (x < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(x), "X must be non-negative.");
            }

            if (x < _factorialCache.Length)
                return _factorialCache[x];
            return double.PositiveInfinity;
        }

        /// <summary>
        /// Computes the factorial of an integer.
        /// </summary>
        /// <param name="x">An integer greater than 0 </param>
        /// <returns>
        /// The product of all integers less than or equal to the given x
        /// </returns>
        public static double Function(long x) 
        {
            if (x < 0L)
            {
                throw new ArgumentOutOfRangeException(nameof(x), "X must be non-negative.");
            }

            if (x == 0L)
                return 1.0d;
            long r = x;
            while (System.Threading.Interlocked.Decrement(ref x) > 1L)
                r *= x;
            return r;
        }

        /// <summary>
        /// Computes the logarithmic factorial function x = ln(x!) of an integer x > 0.
        /// </summary>
        /// <param name="x">An integer greater than 0 </param>
        /// <returns>
        /// The log of the product of all integers less than or equal to the given x
        /// </returns>
        public static double LogFactorial(int x)
        {
            if (x < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(x), "X must be non-negative.");
            }

            if (x <= 0)
                return 0.0d;
            if (x < _factorialCache.Length)
                return Math.Log(_factorialCache[x]);
            return Gamma.LogGamma(x + 1.0d);
        }

        /// <summary>
        /// Computes the binomial coefficient.
        /// </summary>
        /// <remarks>
        /// The equation for the binomial coefficient is :
        /// <code>
        ///     n chose k = n! / k!(n-k)!
        /// </code>
        /// </remarks>
        /// <param name="n">An integer greater than or equal to 0 </param>
        /// <param name="k">An integer greater than or equal to 0 </param>
        /// <returns>
        /// The binomial coefficient evaluated at n and k (n chose k)
        /// </returns>
        public static double BinomialCoefficient(int n, int k)
        {
            if (k < 0 || n < 0 || k > n)
                return 0.0d;
            return Math.Floor(0.5d + Math.Exp(LogFactorial(n) - LogFactorial(k) - LogFactorial(n - k)));
        }

       
        /// <summary>
        /// Helper function for FindCombinations() 
        /// </summary>
        /// <param name="buffer"> An array to store the current combination being constructed </param>
        /// <param name="done"> The index of where the next element of the combination should be placed </param>
        /// <param name="begin"> Index of where to start adding elements. Ensures each element is only included once</param>
        /// <param name="end"> The last index (exclusive) to be inlcuded in the current combination </param>
        private static IEnumerable<int[]> FindCombosRecursive(int[] buffer, int done, int begin, int end)
        {
            for (int i = begin; i < end; i++)
            {
                buffer[done] = i;

                if (done == buffer.Length - 1)
                    yield return buffer;
                else
                    foreach (int[] child in FindCombosRecursive(buffer, done + 1, i + 1, end))
                        yield return child;
            }
        }

        /// <summary>
        /// Finds each m combinations within n.
        /// </summary>
        /// <param name="m">The combination size.</param>
        /// <param name="n">The overall count.</param>
        public static IEnumerable<int[]> FindCombinations(int m, int n)
        {
            return FindCombosRecursive(new int[m], 0, 0, n);
        }

        /// <summary>
        /// Finds all combinations of m within n without replacement.           
        /// </summary>
        /// <param name="n">The overall count.</param>
        /// <returns>
        /// An array of all possible combinations, represented with 1s and 0s, where 1 is the occurrence
        /// of an event.
        /// </returns>
        public static int[,] AllCombinations(int n)
        {
            int f = (int)Math.Pow(2, n) - 1;
            var output = new int[f, n];
            int t = 0;
            for (int i = 1; i <= n; i++)
            {
                foreach (int[] c in FindCombinations(i, n))
                {
                    for (int j = 0; j < c.Length; j++)
                    {
                        output[t, c[j]] = 1;
                    }
                    t++;
                }
            }
            return output;
        }

    }
}