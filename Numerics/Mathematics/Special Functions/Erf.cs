using Numerics.Distributions;

namespace Numerics.Mathematics.SpecialFunctions
{
    /// <summary>
    /// The error function.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The error function, also known as the Gauss error function, is defined 
    /// as:
    /// </para>
    /// <code>
    ///                        X
    ///     erf(X) = 2/sqrt(𝜋) ∫  e^(-t^2) dt
    ///                        0           
    /// </code>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Error_function#cite_note-2"/>
    ///  </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public partial class Erf
    {
        /// <summary>
        /// Computes the error function
        /// </summary>
        /// <param name="X"> The upper bound </param>>
        /// <returns>
        /// The error function evaluated with the given upper bound
        /// </returns>
        public static double Function(double X)
        {
            if (X < 0d)
            {
                return -Gamma.LowerIncomplete(0.5d, X * X);
            }
            else
            {
                return Gamma.LowerIncomplete(0.5d, X * X);
            }
        }

        /// <summary>
        /// Computes the complement of the error function
        /// </summary>
        /// <param name="X"> The upper bound </param>>
        /// <returns>
        /// The complement of the error function evaluated with the given upper bound
        /// </returns>
        public static double Erfc(double X)
        {
            return 1d - Function(X);
        }

        /// <summary>
        /// Computes the inverse error function
        /// </summary>
        ///<param name="y"> The value to be evaluated (such that y = erf(erf^-1(y)) ) </param>
        /// <returns>
        /// The inverse error function evaluated at the given y
        /// </returns>
        public static double InverseErf(double y)
        {
            double s = Normal.StandardZ(0.5d * y + 0.5d);
            double r = s * Tools.Sqrt2 / 2.0d;
            return r;
        }

        /// <summary>
        /// Computes the inverse of the complement of the error function
        /// </summary>
        ///<param name="y"> The value to be evaluated (such that y = erf(erf^-1(y)) ) </param>
        /// <returns>
        /// The inverse of the complement of the error function evaluated at the given y
        /// </returns>
        public static double InverseErfc(double y)
        {
            double s = Normal.StandardZ(-0.5d * y + 1d);
            double r = s * Tools.Sqrt2 / 2.0d;
            return r;
        }
    }
}