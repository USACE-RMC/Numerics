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

using Numerics.Mathematics.Optimization;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Class for performing Yeo-Johnson transformation.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// <para>
    /// <b> Description: </b>
    /// This method transforms non-normal dependent variables into a normal shape.
    /// </para>
    /// </remarks>
    public class YeoJohnson
    {
        /// <summary>
        /// Fit the transform parameter using maximum likelihood estimation.
        /// </summary>
        /// <param name="values">The list of values to transform.</param>
        /// <param name="lambda">Output. The transformation exponent. Range -5 to +5.</param>
        public static void FitLambda(IList<double> values, out double lambda)
        {
            // Solve with Brent 
            var brent = new BrentSearch((x) => { return LogLikelihood(values, x); }, -5d, 5d);
            brent.Maximize();
            if (brent.Status == OptimizationStatus.Success)
            {
                lambda = brent.BestParameterSet.Values[0];
            }
            else
            {
                lambda = double.NaN;
            }
        }

        /// <summary>
        /// The log-likelihood function. The transformed observations are assumed to come from a
        /// normal distribution. The change of variable formula is used to write the log-likelihood function.
        /// </summary>
        /// <param name="values">The list of values to transform.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        /// <returns>
        /// The value of log-likelihood function evaluated at the given values and lambdas.
        /// </returns>
        public static double LogLikelihood(IList<double> values, double lambda)
        {
            int n = values.Count;
            var transformed = new double[n];
            double sum = 0d;
            double logJacobianSum = 0d;

            // Transform values and compute sum and log-Jacobian
            for (int i = 0; i < n; i++)
            {
                double xi = values[i];
                double yi = Transform(xi, lambda);
                transformed[i] = yi;
                sum += yi;

                // Compute derivative dT/dy for log-Jacobian
                double dTdy;
                if (xi >= 0)
                {
                    dTdy = Math.Pow(xi + 1, lambda - 1);
                }
                else
                {
                    dTdy = Math.Pow(-xi + 1, 1 - lambda);
                }

                // Avoid log of zero or negative values
                if (dTdy > 0)
                    logJacobianSum += Math.Log(dTdy);
                else
                    return double.NegativeInfinity; // log-likelihood undefined
            }

            // Compute mean and SSE
            double mu = sum / n;
            double sse = 0d;
            for (int i = 0; i < n; i++)
            {
                double resid = transformed[i] - mu;
                sse += resid * resid;
            }

            double sigmaSq = sse / n;
            if (sigmaSq <= 0 || double.IsNaN(sigmaSq))
                return double.NegativeInfinity;

            double logLikelihood = -n / 2.0 * Tools.LogSqrt2PI - n / 2.0 * Math.Log(sigmaSq) - sse / (2.0 * sigmaSq) + logJacobianSum;

            return logLikelihood;
        }

        /// <summary>
        /// Computes the Log-Jacobian used to adjust the log-likelihood function. 
        /// </summary>
        /// <param name="values">The list of values to transform.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        /// <returns>Returns the Log Jacobian.</returns>
        public static double LogJacobian(IList<double> values, double lambda)
        {
            double logJacobianSum = 0d;
            int n = values.Count;

            // Transform values and compute sum and log-Jacobian
            for (int i = 0; i < n; i++)
            {
                double xi = values[i];

                // Compute derivative dT/dy for log-Jacobian
                double dTdy;
                if (xi >= 0)
                {
                    dTdy = Math.Pow(xi + 1, lambda - 1);
                }
                else
                {
                    dTdy = Math.Pow(-xi + 1, 1 - lambda);
                }

                // Avoid log of zero or negative values
                if (dTdy > 0)
                    logJacobianSum += Math.Log(dTdy);
                else
                    return double.NegativeInfinity; // log-likelihood undefined
            }
            return logJacobianSum;
        }

        /// <summary>
        /// Returns the Yeo-Johnson transformation of the value.
        /// </summary>
        /// <param name="value">The value to transform.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        public static double Transform(double value, double lambda)
        {
            if (Math.Abs(lambda) > 5d) return double.NaN;
            if (value >= 0 && Math.Abs(lambda) >= 1E-8)
            {
                return (Math.Pow(value + 1, lambda) - 1) / lambda;
            }
            else if (value >= 0 && Math.Abs(lambda) < 1E-8)
            {
                return Math.Log(value + 1);
            }
            else if (value < 0 && lambda != 2)
            {
                return -(Math.Pow(-value + 1, 2 - lambda) - 1) / (2 - lambda);
            }
            else if (value < 0 && lambda == 2)
            {
                return -Math.Log(-value + 1);
            }
            else
            {
                return double.NaN;
            }
        }

        /// <summary>
        /// Returns the Yeo-Johnson transformation of each value in the list.
        /// </summary>
        /// <param name="values">The list of values to transform.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        public static List<double> Transform(IList<double> values, double lambda)
        {
            var newValues = new List<double>();
            for (int i = 0; i < values.Count; i++)
                newValues.Add(Transform(values[i], lambda));
            return newValues;
        }

        /// <summary>
        /// Returns the inverse of the Yeo-Johnson transformed value.
        /// </summary>
        /// <param name="value">The value to reverse transform.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        public static double InverseTransform(double value, double lambda)
        {
            if (Math.Abs(lambda) > 5d) return double.NaN;
            if (value >= 0 && Math.Abs(lambda) >= 1E-8)
            {
                return Math.Pow(lambda * value + 1, 1 / lambda) - 1;
            }
            else if (value >= 0 && Math.Abs(lambda) < 1E-8)
            {
                return Math.Exp(value) - 1;
            }
            else if (value < 0 && lambda != 2)
            {
                return 1 - (Math.Pow(1 - (2 - lambda) * value, 1 / (2 - lambda)) - 1);
            }
            else if (value < 0 && lambda == 2)
            {
                return 1 - Math.Exp(-value);
            }
            else
            {
                return double.NaN;
            }
        }

        /// <summary>
        /// Returns the inverse of each Yeo-Johnson transformed value in the list.
        /// </summary>
        /// <param name="values">The list of values to reverse transform.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        public static List<double> InverseTransform(IList<double> values, double lambda)
        {
            var newValues = new List<double>();
            for (int i = 0; i < values.Count; i++)
                newValues.Add(InverseTransform(values[i], lambda));
            return newValues;
        }
    }
}
