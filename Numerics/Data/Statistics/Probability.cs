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

using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace Numerics.Data.Statistics
{

    /// <summary>
    /// A class for performing probability calculations for risk analysis. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class Probability
    {

        /// <summary>
        /// Enumeration of dependency types. 
        /// </summary>
        public enum DependencyType
        {
            /// <summary>
            /// Statistically independent.
            /// </summary>
            Independent,
            /// <summary>
            /// Perfectly positively dependent.
            /// </summary>
            PerfectlyPositive,
            /// <summary>
            /// Perfectly negatively dependent.
            /// </summary>
            PerfectlyNegative,
            /// <summary>
            /// User-defined correlation matrix.
            /// </summary>
            CorrelationMatrix
        }

        // A large but finite z to hand to CDFs that might not love +∞ in all paths.
        private const double Z_MAX = 8.0; // Φ(8) ~ 0.9999999999999993

        #region Basic Probability Rules for Two Random Variables

        /// <summary>
        /// Returns the probability of intersection (or joint probability) of A and B, P(A and B). 
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double AAndB(double A, double B, double rho = 0d)
        {
            if (A == 0d || B == 0d) return 0d;
            if (A == 1d) return B;
            if (B == 1d) return A;
            if (rho <= -0.999) return Math.Max(0d, A + B - 1);
            if (rho >= 0.999) return Math.Min(A, B);
            if (Math.Abs(rho) <= 1E-3) return A * B;
            return Tools.Clamp(MultivariateNormal.BivariateCDF(Normal.StandardZ(1 - A), Normal.StandardZ(1 - B), rho), 0, 1);
        }

        /// <summary>
        /// Returns the probability of union of A and B, P(A or B).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double AOrB(double A, double B, double rho = 0d)
        {
            return Tools.Clamp(A + B - AAndB(A, B, rho), 0, 1);
        }

        /// <summary>
        /// Returns the probability of A and not B, P(A and not B).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double ANotB(double A, double B, double rho = 0d)
        {
            return Tools.Clamp(A - AAndB(A, B, rho), 0, 1);
        }

        /// <summary>
        /// Returns the probability of B and not A, P(B and not A).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double BNotA(double A, double B, double rho = 0d)
        {
            return Tools.Clamp(B - AAndB(A, B, rho), 0, 1);
        }

        /// <summary>
        /// Returns the probability of A given B, P(A|B).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double AGivenB(double A, double B, double rho = 0d)
        {
            return Tools.Clamp(AAndB(A, B, rho) / B, 0, 1);
        }

        /// <summary>
        /// Returns the probability of B given A, P(B|A).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double BGivenA(double A, double B, double rho = 0d)
        {
            return Tools.Clamp(AAndB(A, B, rho) / A, 0, 1);
        }

        #endregion

        #region Joint Probability

        /// <summary>
        /// Returns the joint probability.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="dependency">The dependency type. Default = Independent.</param>
        public static double JointProbability(IList<double> probabilities, DependencyType dependency = DependencyType.Independent)
        {
            if (dependency == DependencyType.Independent)
            {
                return IndependentJointProbability(probabilities);
            }
            else if (dependency == DependencyType.PerfectlyPositive)
            {
                return PositiveJointProbability(probabilities);
            }
            else if (dependency == DependencyType.PerfectlyNegative)
            {
                return NegativeJointProbability(probabilities);
            }
            return double.NaN;
        }

        /// <summary>
        /// Computes the joint probability of multiple events with dependency, using the Product of Conditional Marginals (PCM) method.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency. Default = null.</param>
        /// <param name="dependency">The dependency type. Default = Correlation matrix.</param>
        /// <returns>The joint probability.</returns>
        public static double JointProbability(IList<double> probabilities, int[] indicators, double[,]? correlationMatrix = null, DependencyType dependency = DependencyType.CorrelationMatrix)
        {
            if (dependency == DependencyType.CorrelationMatrix && correlationMatrix != null)
            {
                return JointProbabilityHPCM(probabilities, indicators, correlationMatrix);
            }
            else if (dependency == DependencyType.Independent)
            {
                return IndependentJointProbability(probabilities, indicators);
            }
            else if (dependency == DependencyType.PerfectlyPositive)
            {
                return PositiveJointProbability(probabilities, indicators);
            }
            else if (dependency == DependencyType.PerfectlyNegative)
            {
                return NegativeJointProbability(probabilities, indicators);
            }
            return double.NaN;
        }

        /// <summary>
        /// Returns the joint probability assuming perfect independence. 
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        public static double IndependentJointProbability(IList<double> probabilities)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));

            double p = 1;
            for (int i = 0; i < probabilities.Count; i++)
            {
                p *= probabilities[i];
                if (p == 0d) return 0d;
            }
            return Tools.Clamp(p, 0, 1);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect independence. 
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double IndependentJointProbability(IList<double> probabilities, int[] indicators)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.Length == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.Length)
                throw new ArgumentException("The probabilities and indicators arrays must have the same length.", nameof(probabilities));
            return Tools.Clamp(Tools.Product(probabilities, indicators), 0, 1);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect positive dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double PositiveJointProbability(IList<double> probabilities)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            return Tools.Clamp(Tools.Min(probabilities), 0, 1);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect positive dependence. 
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double PositiveJointProbability(IList<double> probabilities, int[] indicators)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.Length == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.Length)
                throw new ArgumentException("The probabilities and indicators arrays must have the same length.", nameof(probabilities));
            return Tools.Clamp(Tools.Min(probabilities, indicators), 0, 1);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect negative dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double NegativeJointProbability(IList<double> probabilities)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            return Math.Max(0, Math.Min(1, Tools.Sum(probabilities)) - 1);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect negative dependence. 
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double NegativeJointProbability(IList<double> probabilities, int[] indicators)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.Length == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.Length)
                throw new ArgumentException("The probabilities and indicators arrays must have the same length.", nameof(probabilities));
            return Math.Max(0, Math.Min(1, Tools.Sum(probabilities, indicators)) - 1);
        }

        /// <summary>
        /// Computes the joint probability of multiple events with dependency, using Haden Smith's modification of Pandey's method for the Product of Conditional Marginals (PCM).
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event. These values should lie between 0 and 1, representing the marginal probabilities of individual events.</param>
        /// <param name="indicators">An array of indicators, where 0 means the event did not occur, and 1 means the event did occur. This is used to adjust the conditional probabilities during the calculation of the joint probability.</param>
        /// <param name="correlationMatrix">A 2D array representing the correlation matrix that defines the dependencies between the events. It is assumed that the matrix is symmetric and valid, meaning the correlation coefficient between each pair of events lies between -1 and 1.</param>
        /// <param name="conditionalProbabilities">An optional output array that, if provided, will store the conditional probabilities for each event after the calculation. This array must have the same length as the number of events.</param>
        /// <returns>
        /// The joint probability of the events, adjusted for dependencies as defined by the correlation matrix. The return value is between 0 and 1.
        /// </returns>
        /// <remarks>
        /// This method utilizes a modified version of Pandey's PCM method.
        /// </remarks>
        public static double JointProbabilityHPCM(IList<double> probabilities, int[] indicators, double[,] correlationMatrix, double[]? conditionalProbabilities = null)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.Length == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.Length)
                throw new ArgumentException("The probabilities and indicators arrays must have the same length.", nameof(probabilities));
            if (correlationMatrix == null)
                throw new ArgumentNullException(nameof(correlationMatrix), "The correlation matrix cannot be null.");

            int n = probabilities.Count;
            int rows = correlationMatrix.GetLength(0);
            int cols = correlationMatrix.GetLength(1);

            if (rows != n || cols != n)
            {
                throw new ArgumentException("The correlation matrix must be a square matrix with dimensions equal to the length of the probabilities array.", nameof(correlationMatrix));
            }

            // Get z-values
            double zMin = -9, zMax = 9;
            var R = new double[n, n];
            Array.Copy(correlationMatrix, R, correlationMatrix.Length);
            int i, j, k, ir, ic;
            double pdf, cdf, A, B, z1, z2, r12, z21, p21, jp;
            for (i = 0; i < n; i++)
            {
                if (indicators[i] == 0)
                {
                    R[i, i] = Tools.Clamp(Normal.StandardZ(1), zMin, zMax);
                }
                else
                {
                    R[i, i] = Tools.Clamp(Normal.StandardZ(probabilities[i]), zMin, zMax);
                }
            }
            // Update the conditional correlation matrix
            // First cycle
            z1 = R[0, 0];
            pdf = Normal.StandardPDF(z1);
            cdf = Normal.StandardCDF(z1);
            //if (cdf < 1e-300) cdf = 1e-300;
            A = pdf / cdf;
            B = A * (z1 + A);
            for (k = 1; k < n; k++)
            {
                z2 = R[k, k];
                r12 = R[0, k];
                r12 = Math.Abs(r12) < 1E-3 ? 0: r12;
                p21 = MultivariateNormal.BivariateCDF(-z1, -z2, r12) / cdf;
                p21 = Math.Max(0, Math.Min(1, p21));
                z21 = Tools.Clamp(Normal.StandardZ(p21), zMin, zMax);
                R[k, 0] = z21;
            }
            // update condition correlation r[ir|ic] and store them in R[ir,ic]
            for (ir = 1; ir < n - 1; ir++)
            {
                for (ic = ir + 1; ic < n; ic++)
                {
                   R[ir, ic] = (R[ir, ic] - R[0, ir] * R[0, ic] * B) / Math.Sqrt((1d - R[0, ir] * R[0, ir] * B) * (1d - R[0, ic] * R[0, ic] * B));
                }
            }
            // Remaining cycles
            for (j = 1; j < n - 1; j++)
            {
                z1 = R[j, j - 1];
                pdf = Normal.StandardPDF(z1);
                cdf = Normal.StandardCDF(z1);
                if (cdf < 1e-300) cdf = 1e-300;
                A = pdf / cdf;
                B = A * (z1 + A);
                for (k = j + 1; k < n; k++)
                {
                    z2 = R[k, j - 1];
                    r12 = R[j, k];
                    r12 = Math.Abs(r12) < 1E-3 ? 0 : r12;
                    p21 = MultivariateNormal.BivariateCDF(-z1, -z2, r12) / cdf;
                    p21 = Math.Max(0, Math.Min(1, p21));
                    z21 = Tools.Clamp(Normal.StandardZ(p21), zMin, zMax);
                    R[k, j] = z21;

                }
                for (ir = j + 1; ir < n - 1; ir++)
                {
                    for (ic = ir + 1; ic < n; ic++)
                    {
                        R[ir, ic] = (R[ir, ic] - R[j, ir] * R[j, ic] * B) / Math.Sqrt((1d - R[j, ir] * R[j, ir] * B) * (1d - R[j, ic] * R[j, ic] * B));
                    }
                }
            }

            // Calculate the product of conditional marginals (PCM)
            jp = Math.Log(Normal.StandardCDF(R[0, 0]));
            if (conditionalProbabilities != null && conditionalProbabilities.Length == n)
                conditionalProbabilities[0] = Normal.StandardCDF(R[0, 0]);
            for (i = 1; i < n; i++)
            {
                jp += Math.Log(Normal.StandardCDF(R[i, i - 1]));
                if (conditionalProbabilities != null && conditionalProbabilities.Length == n)
                    conditionalProbabilities[i] = Normal.StandardCDF(R[i, i - 1]);
            }
            jp = Math.Exp(jp);
            jp = Math.Min(1, Math.Max(0, jp));
            if (double.IsNaN(jp)) jp = 0;
            return jp;
        }

        /// <summary>
        /// Computes the joint probability of multiple events with dependency, using Pandey's original method for the Product of Conditional Marginals (PCM).
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event. These values should lie between 0 and 1, representing the marginal probabilities of individual events.</param>
        /// <param name="indicators">An array of indicators, where 0 means the event did not occur, and 1 means the event did occur. This is used to adjust the conditional probabilities during the calculation of the joint probability.</param>
        /// <param name="correlationMatrix">A 2D array representing the correlation matrix that defines the dependencies between the events. It is assumed that the matrix is symmetric and valid, meaning the correlation coefficient between each pair of events lies between -1 and 1.</param>
        /// <param name="conditionalProbabilities">An optional output array that, if provided, will store the conditional probabilities for each event after the calculation. This array must have the same length as the number of events.</param>
        /// <returns>
        /// The joint probability of the events, adjusted for dependencies as defined by the correlation matrix. The return value is between 0 and 1.
        /// </returns>
        public static double JointProbabilityPCM(IList<double> probabilities, int[] indicators, double[,] correlationMatrix, double[]? conditionalProbabilities = null)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.Length == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.Length)
                throw new ArgumentException("The probabilities and indicators arrays must have the same length.", nameof(probabilities));
            if (correlationMatrix == null)
                throw new ArgumentNullException(nameof(correlationMatrix), "The correlation matrix cannot be null.");

            int n = probabilities.Count;
            int rows = correlationMatrix.GetLength(0);
            int cols = correlationMatrix.GetLength(1);

            if (rows != n || cols != n)
            {
                throw new ArgumentException("The correlation matrix must be a square matrix with dimensions equal to the length of the probabilities array.", nameof(correlationMatrix));
            }

            // Get z-values
            double zMin = -9, zMax = 9;
            var R = new double[n, n];
            Array.Copy(correlationMatrix, R, correlationMatrix.Length);
            int i, j, k, ir, ic;
            double A, B, z1, z2, z21, r12;
            for (i = 0; i < n; i++)
            {
                if (indicators[i] == 0)
                {
                    R[i, i] = Tools.Clamp(Normal.StandardZ(1), zMin, zMax);
                }
                else
                {
                    R[i, i] = Tools.Clamp(Normal.StandardZ(probabilities[i]), zMin, zMax);
                }
            }
            // Update the conditional correlation matrix
            // First cycle
            z1 = R[0, 0];
            var pdf = Normal.StandardPDF(z1);
            var cdf = Normal.StandardCDF(z1);
            if (cdf < 1e-300) cdf = 1e-300;
            A = pdf / cdf;
            B = A * (z1 + A);
            // calculate z[k|0] and store them in R[k,0], k = 1,...,n
            for (k = 1; k < n; k++)
            {
                z2 = R[k, k];
                r12 = R[0, k];
                z21 = (z2 + r12 * A) / Math.Sqrt(1d - r12 * r12 * B);
                R[k, 0] = z21;
            }
            // update r[ir|ic] and store them in R[ir,ic]
            for (ir = 1; ir < n - 1; ir++)
            {
                for (ic = ir + 1; ic < n; ic++)
                {
                    R[ir, ic] = (R[ir, ic] - R[0, ir] * R[0, ic] * B) / Math.Sqrt((1d - R[0, ir] * R[0, ir] * B) * (1d - R[0, ic] * R[0, ic] * B));
                }
            }
            // Remaining cycles
            for (j = 1; j < n - 1; j++)
            {
                z1 = R[j, j - 1];
                pdf = Normal.StandardPDF(z1);
                cdf = Normal.StandardCDF(z1);
                if (cdf < 1e-300) cdf = 1e-300;
                A = pdf / cdf;
                B = A * (z1 + A);
                for (k = j + 1; k < n; k++)
                {
                    z2 = R[k, j - 1];
                    r12 = R[j, k];
                    z21 = (z2 + r12 * A) / Math.Sqrt(1d - r12 * r12 * B);
                    R[k, j] = z21;
                }
                for (ir = j + 1; ir < n - 1; ir++)
                {
                    for (ic = ir + 1; ic < n; ic++)
                    {
                        R[ir, ic] = (R[ir, ic] - R[j, ir] * R[j, ic] * B) / Math.Sqrt((1d - R[j, ir] * R[j, ir] * B) * (1d - R[j, ic] * R[j, ic] * B));
                    }
                }
            }
            // Calculate the product of conditional marginals (PCM)
            double jp = Math.Log(Normal.StandardCDF(R[0, 0]));
            if (conditionalProbabilities != null && conditionalProbabilities.Length == n)
                conditionalProbabilities[0] = Normal.StandardCDF(R[0, 0]);
            for (i = 1; i < n; i++)
            {
                jp += Math.Log(Normal.StandardCDF(R[i, i - 1]));
                if (conditionalProbabilities != null && conditionalProbabilities.Length == n)
                    conditionalProbabilities[i] = Normal.StandardCDF(R[i, i - 1]);
            }
            jp = Math.Exp(jp);
            jp = Math.Min(1, Math.Max(0, jp));
            if (double.IsNaN(jp)) jp = 0;
            return jp;
        }

        /// <summary>
        /// Returns an array of joint probabilities of multiple events with dependency, using the Product of Conditional Marginals (PCM) method.
        /// </summary>
        /// <param name="probabilities">And array of probabilities for each event.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        public static double[] JointProbabilitiesPCM(IList<double> probabilities, int[,] indicators, double[,] correlationMatrix)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must be non-null and contain at least one row.");
            if (correlationMatrix == null)
                throw new ArgumentNullException(nameof(correlationMatrix), "The correlation matrix cannot be null.");
            
            var result = new double[indicators.GetLength(0)];

            Parallel.For(0, indicators.GetLength(0), idx =>
            {
                if (idx < probabilities.Count)
                {
                    result[idx] = probabilities[idx];
                }
                else
                {
                    result[idx] = JointProbabilityPCM(probabilities, indicators.GetRow(idx), correlationMatrix);
                }
            });
            return result;
        }

        /// <summary>
        /// Returns the joint probability of multiple events with dependency. 
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution for computing the joint probability.</param>
        public static double JointProbabilityMVN(IList<double> probabilities, int[] indicators, MultivariateNormal multivariateNormal)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must be non-null and contain at least one row.");
            if (multivariateNormal == null)
                throw new ArgumentException("The multivariateNormal distribution must be non-null.");

            var zVals = new double[indicators.Length];
            for (int i = 0; i < indicators.Length; i++)
            {
                if (indicators[i] == 0)
                {
                    zVals[i] = double.PositiveInfinity;
                }
                else
                {
                    zVals[i] = Normal.StandardZ(probabilities[i]);
                }
            }
            var p = multivariateNormal.CDF(zVals);
            p = Math.Max(0, Math.Min(1, p));
            return p;
        }

        /// <summary>
        /// Returns an array of joint probabilities of multiple events with dependency.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution for computing the joint probability.</param>
        public static double[] JointProbabilitiesMVN(IList<double> probabilities, int[,] indicators, MultivariateNormal multivariateNormal)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must be non-null and contain at least one row.");
            if (multivariateNormal == null)
                throw new ArgumentException("The multivariateNormal distribution must be non-null.");

            var result = new double[indicators.GetLength(0)];

            Parallel.For(0, indicators.GetLength(0), idx =>
            {
                if (idx < probabilities.Count)
                {
                    result[idx] = probabilities[idx];
                }
                else
                {
                    result[idx] = JointProbabilityMVN(probabilities, indicators.GetRow(idx), (MultivariateNormal)multivariateNormal.Clone());
                }           
            });
            return result;
        }

        #endregion

        #region Probability of Union

        /// <summary>
        /// Compute the probability of union.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="dependency">The dependency type. Default = Independent.</param>
        public static double Union(IList<double> probabilities, DependencyType dependency = DependencyType.Independent)
        {
            if (dependency == DependencyType.Independent)
            {
                return IndependentUnion(probabilities);
            }
            else if (dependency == DependencyType.PerfectlyPositive)
            {
                return PositivelyDependentUnion(probabilities);
            }
            else if (dependency == DependencyType.PerfectlyNegative)
            {
                return NegativelyDependentUnion(probabilities);
            }
            return double.NaN;
        }

        /// <summary>
        /// Returns the probability of union assuming independence (De Morgan's Rule).
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double IndependentUnion(IList<double> probabilities)
        {
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (probabilities.Count == 1) return probabilities[0];

            double numerator = 1d;
            for (int i = 0; i < probabilities.Count; i++)
            {
                var q = 1d - probabilities[i];
                if (q == 0d) return 1d; // any event certain -> union = 1
                numerator *= q;
            }
            return 1d - numerator;
        }

        /// <summary>
        /// Returns the unimodal bound for the probability of union assuming perfect positive dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double PositivelyDependentUnion(IList<double> probabilities)
        {
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (probabilities.Count == 1) return probabilities[0];
            return Tools.Clamp(Tools.Max(probabilities), 0, 1);
        }

        /// <summary>
        /// Returns the unimodal bound for the probability of union assuming perfect negative dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double NegativelyDependentUnion(IList<double> probabilities)
        {
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (probabilities.Count == 1) return probabilities[0];
            return Tools.Clamp(Tools.Sum(probabilities), 0, 1);
        }

        /// <summary>
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        public static double UnionPCM(IList<double> probabilities, double[,] correlationMatrix, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            // Ensure the input is valid
            if (probabilities == null || probabilities.Count == 0 || correlationMatrix == null)
            {
                throw new ArgumentException("Input arrays must be non-empty and correlation matrix must not be null.");
            }

            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var binomialCombinations = new int[N];
            for (int i = 1; i <= N; i++)
            {
                binomialCombinations[i - 1] = (int)Factorial.BinomialCoefficient(N, i);
            }

            // Get combination indicators
            var indicators = Factorial.AllCombinations(N);

            // Return Union
            return UnionPCM(probabilities, binomialCombinations, indicators, correlationMatrix, absoluteTolerance, relativeTolerance);
         
        }

        /// <summary>
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the PCM method.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <returns>The probability of the union of the events using the inclusion-exclusion method.</returns>
        public static double UnionPCM(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, double[,] correlationMatrix, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            // Ensure the input is valid
            if (probabilities == null || probabilities.Count == 0 || binomialCombinations.Length == 0 || indicators.GetLength(0) == 0 || correlationMatrix == null)
            {
                throw new ArgumentException("Input arrays must be non-empty and correlation matrix must not be null.");
            }

            double result = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];
            double inc = double.NaN;
            double exc = double.NaN;
            int numIndicators = indicators.GetLength(0);  // Reduce redundant lookups for array lengths

            for (int i = 0; i < numIndicators; i++)
            {
                // Adjust result for inclusion-exclusion based on binomial combinations
                if (i == c)
                {
                    if (j > 0)
                    {
                        // Set inc and exc when transitioning
                        if (s == 1) inc = result;
                        else if (s == -1) exc = result;
                    }

                    // Check for convergence based on tolerance
                    double diff = Math.Abs(inc - exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        return result + 0.5 * diff; // Converged, return the result with half of the difference
                    }

                    s *= -1; // Alternate sign for inclusion-exclusion
                    j++;

                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j]; // Update the next binomial combination index
                    }

                }

                // Add the contribution of the current event, either from the probability or joint probability calculation
                if (i < probabilities.Count)
                {
                    result += s * probabilities[i];
                }
                else
                {
                    result += s * JointProbability(probabilities, indicators.GetRow(i), correlationMatrix);
                }

            }

            return result;
        }

        /// <summary>
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the PCM method.
        /// </summary>
        /// <param name="probabilities">List of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">A 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency between events.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities.</param>
        /// <param name="eventIndicators">Output. A list of exclusive event indicators that were evaluated.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <returns>The probability of the union of the events using the inclusion-exclusion method.</returns>
        public static double UnionPCM(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, double[,] correlationMatrix, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            // Ensure the input is valid
            if (probabilities == null || probabilities.Count == 0 || binomialCombinations.Length == 0 || indicators.GetLength(0) == 0 || correlationMatrix == null)
            {
                throw new ArgumentException("Input arrays must be non-empty and correlation matrix must not be null.");
            }

            // Initialize output lists
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();

            double union = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];
            double inc = double.NaN;
            double exc = double.NaN;

            int numIndicators = indicators.GetLength(0);  // Cache the number of rows for efficiency

            // Loop through all combinations of events
            for (int i = 0; i < numIndicators; i++)
            {
                if (i == c)
                {
                    // Set inc and exc when transitioning between inclusion-exclusion steps
                    if (j > 0)
                    {
                        if (s == 1) inc = union;
                        else if (s == -1) exc = union;
                    }

                    // Check for convergence based on tolerance
                    double diff = Math.Abs(inc - exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        eventIndicators.Add(indicators.GetRow(indicators.GetLength(0) - 1));  // Add the last row for event indicators
                        eventProbabilities.Add(0.5 * diff);  // Add the averaged difference
                        return union + 0.5 * diff; // Converged, return the result with half of the difference
                    }

                    s *= -1; // Alternate the sign for inclusion-exclusion
                    j++;

                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j]; // Update the next binomial combination index
                    }

                }

                // Record indicators for the current event
                eventIndicators.Add(indicators.GetRow(i));

                // Add the contribution of the current event, either from the probability or joint probability calculation
                if (i < probabilities.Count)
                {
                    union += s * probabilities[i];  // If the event is within the range of probabilities, add directly
                    eventProbabilities.Add(probabilities[i]);  // Store the probability
                }
                else
                {
                    var jp = JointProbability(probabilities, indicators.GetRow(i), correlationMatrix);  // Otherwise, calculate the joint probability
                    union += s * jp;  // Add the joint probability contribution
                    eventProbabilities.Add(jp);  // Store the joint probability
                }

            }

            return union;
        }

        /// <summary>
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        public static double UnionMVN(IList<double> probabilities, MultivariateNormal multivariateNormal)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (multivariateNormal == null)
                throw new ArgumentException("The multivariateNormal distribution must be non-null.");

            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var binomialCombinations = new int[N];
            for (int i = 1; i <= N; i++)
            {
                binomialCombinations[i - 1] = (int)Factorial.BinomialCoefficient(N, i);
            }

            // Get combination indicators
            var indicators = Factorial.AllCombinations(N);

            // Return result
            return UnionMVN(probabilities, binomialCombinations, indicators, multivariateNormal);

        }

        /// <summary>
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        public static double UnionMVN(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, MultivariateNormal multivariateNormal)
        {
            
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (binomialCombinations == null || binomialCombinations.Length == 0)
                throw new ArgumentException("The binomialCombinations array must be non-null and contain at least one element.");
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must be non-null and contain at least one row.");
            if (multivariateNormal == null)
                throw new ArgumentException("The multivariateNormal distribution must be non-null.");

            // Ensure that the length of binomialCombinations matches the number of rows in the indicators array
            if (binomialCombinations.Length != indicators.GetLength(1))
                throw new ArgumentException("The length of binomialCombinations must match the number of columns in the indicators array.");


            double result = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];

            // Loop through all possible combinations of events
            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                // Adjust the inclusion-exclusion signs based on the binomial combinations
                if (i == c)
                {
                    s *= -1; // Alternate the sign for inclusion-exclusion
                    j++;
                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j]; // Update the combination index for the next subset
                    }
                }

                // If we are within the range of probabilities, simply add the probability with the current sign
                if (i < probabilities.Count)
                {
                    result += s * probabilities[i];
                }
                else
                {
                    // If beyond the range of probabilities, calculate the joint probability using the MVN distribution
                    result += s * JointProbabilityMVN(probabilities, indicators.GetRow(i), multivariateNormal);
                }
            }

            return result;
        }

        #endregion

        #region Exclusive Probability of all Combinations of Events

        #region Independent

        /// <summary>
        /// Returns the exclusive probability of multiple events occurring assuming independence.
        /// This method calculates the probability of the intersection of events occurring assuming independence. 
        /// It uses the indicator array to determine which events are considered in the calculation.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event. Each element represents the probability of an individual event occurring.</param>
        /// <param name="indicators">An array of indicators, where 0 means the event did not occur and 1 means the event did occur.</param>
        /// <returns>The exclusive probability of the events occurring, calculated as the product of the individual event probabilities, considering the indicator array.</returns>
        /// <exception cref="ArgumentException">Thrown if the probabilities array is null, empty, or if the lengths of the probabilities and indicators arrays do not match.</exception>
        public static double IndependentExclusive(IList<double> probabilities, int[] indicators)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.Length == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.Length)
                throw new ArgumentException("The probabilities and indicators arrays must have the same length.", nameof(probabilities));

            double result = 1;
            for (int i = 0; i < probabilities.Count; i++)
            {
                if (double.IsNaN(probabilities[i])) return double.NaN;
                if (indicators[i] == 1)
                {
                    result *= probabilities[i];
                }
                else
                {
                    result *= (1 - probabilities[i]);
                }
            }
            return result;
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events occurring assuming independence.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">A 2D array of indicators, where each row represents a combination of events, and 0 means the event did not occur and 1 means the event did occur.</param>
        /// <returns>An array of exclusive probabilities for each combination of events, calculated as the product of the event probabilities considering the corresponding indicator array.</returns>
        /// <exception cref="ArgumentException">Thrown if the probabilities array is null or empty, or if the indicators array is null or has invalid dimensions.</exception>
        public static double[] IndependentExclusive(IList<double> probabilities, int[,] indicators)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must be non-null and contain at least one row.");

            var result = new double[indicators.GetLength(0)];
            for (int i = 0; i < indicators.GetLength(0); i++)
            {             
                result[i] = IndependentExclusive(probabilities, indicators.GetRow(i));
            }
            return result;         
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events occurring assuming independence.
        /// This method calculates the exclusive probabilities for all possible combinations of events using the probabilities array.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <returns>An array of exclusive probabilities for all possible combinations of the events, assuming independence.</returns>
        /// <exception cref="ArgumentException">Thrown if the probabilities array is null or empty.</exception>
        public static double[] IndependentExclusive(IList<double> probabilities)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));

            int n = probabilities.Count;
            int f = (int)Math.Pow(2, n) - 1; // Number of non-empty subsets
            var result = new double[f];
            int t = 0;

            // Loop through all possible combinations of events
            for (int i = 1; i <= n; i++)
            {
                foreach (int[] combos in Factorial.FindCombinations(i, n))
                {
                    var indicators = new int[n];
                    for (int j = 0; j < combos.Length; j++)
                    {
                        indicators[combos[j]] = 1; // Mark the event as occurring in the combination
                    }
                    result[t] = IndependentExclusive(probabilities, indicators); // Calculate the exclusive probability for the current combination
                    t++;
                }
            }
            return result;
        }

        /// <summary>
        /// Returns a list of exclusive probabilities of multiple events occurring assuming independence.
        /// This method calculates the exclusive probabilities for all possible event combinations using the provided probabilities and indicators arrays.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event. Each element represents the probability of an individual event occurring.</param>
        /// <param name="binomialCombinations">An array of binomial combinations that define the number of events to consider for each calculation.</param>
        /// <param name="indicators">A 2D array of indicators, where each row represents a combination of events, and 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities for each combination, calculated based on the event indicators.</param>
        /// <param name="eventIndicators">Output. A list of event indicators that correspond to the probabilities in the eventProbabilities list.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <exception cref="ArgumentException">Thrown if the probabilities array is null, empty, or if the lengths of the probabilities and indicators arrays do not match.</exception>
        /// <remarks>
        /// This method uses the inclusion-exclusion principle to compute the exclusive probability of each event combination.
        /// The result is added to a list, and convergence is monitored using the specified tolerances.
        /// </remarks>
        public static void IndependentExclusive(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.GetLength(1))
                throw new ArgumentException("The probabilities array and the indicator array must have the same length.", nameof(probabilities));

            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();

            double union = 0;
            double s = 1; // Sign for inclusion-exclusion
            int j = 0; // Binomial combination index
            int c = binomialCombinations[j]; // Current combination limit
            double inc = double.NaN; // Temporary variable for inclusion value
            double exc = double.NaN; // Temporary variable for exclusion value
            int numIndicators = indicators.GetLength(0);  // Cache the number of rows for efficiency

            // Loop through each row of indicators (each event combination)
            for (int i = 0; i < numIndicators; i++)
            {
                if (i == c)
                {

                    // Set inc and exc when transitioning between inclusion-exclusion steps
                    if (j > 0)
                    {
                        if (s == 1) inc = union;
                        else if (s == -1) exc = union;
                    }

                    // Check for convergence based on tolerance
                    double diff = Math.Abs(inc - exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        eventIndicators.Add(indicators.GetRow(indicators.GetLength(0) - 1)); // Add last indicator row
                        eventProbabilities.Add(0.5 * diff); // Add the average of the difference to the event probabilities
                        return; // Exit early when convergence is reached
                    }

                    // Flip the sign for the next inclusion-exclusion term
                    s *= -1;
                    j++;
                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j]; // Update the current combination limit
                    }

                }

                // Record the current indicators
                eventIndicators.Add(indicators.GetRow(i));

                // Compute the exclusive event probability and add to the list
                eventProbabilities.Add(IndependentExclusive(probabilities, eventIndicators.Last()));

                // Calculate the union of probabilities (inclusion-exclusion)
                if (i < probabilities.Count)
                {
                    union += s * probabilities[i];
                }
                else
                {
                    union += s * IndependentJointProbability(probabilities, eventIndicators.Last());
                }

            }

        }

        #endregion

        #region Positively Dependent

        /// <summary>
        /// Returns the exclusive probability of multiple events occurring assuming perfect positive dependence.
        /// This method computes the exclusive probability by finding the minimum probability for events that occur and the maximum for events that do not occur.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event. Each element represents the probability of an individual event occurring.</param>
        /// <param name="indicators">An array of indicators where 0 means the event did not occur and 1 means the event did occur.</param>
        /// <returns>The exclusive probability of multiple events occurring with perfect positive dependence.</returns>
        /// <exception cref="ArgumentException">Thrown if the probabilities array is null, empty, or if the lengths of the probabilities and indicators arrays do not match.</exception>
        public static double PositivelyDependentExclusive(IList<double> probabilities, int[] indicators)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.Length == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.Length)
                throw new ArgumentException("The probabilities and indicators arrays must have the same length.", nameof(probabilities));

            double min = 1.0;
            double max = 0.0;
            for (int i = 0; i < probabilities.Count; i++)
            {
                if (double.IsNaN(probabilities[i])) return double.NaN;
                if (indicators[i] == 1)
                {
                    if (probabilities[i] < min) min = probabilities[i];
                }
                else
                {
                    if (probabilities[i] > max) max = probabilities[i];
                }
            }
            return Math.Max(min - max, 0);
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events occurring assuming perfect positive dependence.
        /// This method computes the exclusive probabilities for each event combination based on the given indicator array.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">A 2D array of indicators, where 0 means the event did not occur and 1 means the event did occur.</param>
        /// <returns>An array of exclusive probabilities for each event combination, assuming perfect positive dependence.</returns>
        /// <exception cref="ArgumentException">Thrown if the probabilities array is null, empty, or if the lengths of the probabilities and indicators arrays do not match.</exception>
        public static double[] PositivelyDependentExclusive(IList<double> probabilities, int[,] indicators)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.GetLength(1))
                throw new ArgumentException("The probabilities array and the indicator array must have the same length.", nameof(probabilities));
            
            var result = new double[indicators.GetLength(0)];
            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                result[i] = PositivelyDependentExclusive(probabilities, indicators.GetRow(i));
            }
            return result;
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events occurring assuming perfect positive dependence.
        /// This method calculates the exclusive probabilities for all combinations of events.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <returns>An array of exclusive probabilities for each event combination, assuming perfect positive dependence.</returns>
        /// <exception cref="ArgumentException">Thrown if the probabilities array is null, empty, or if any event combination is not valid.</exception>
        public static double[] PositivelyDependentExclusive(IList<double> probabilities)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));

            int n = probabilities.Count;
            int f = (int)Math.Pow(2, n) - 1;
            var result = new double[f];
            int t = 0;
            for (int i = 1; i <= n; i++)
            {
                foreach (int[] combos in Factorial.FindCombinations(i, n))
                {
                    var indicators = new int[n];
                    for (int j = 0; j < combos.Length; j++)
                    {
                        indicators[combos[j]] = 1;
                    }
                    result[t] = PositivelyDependentExclusive(probabilities, indicators);
                    t++;
                }
            }
            return result;
        }

        /// <summary>
        /// Returns a list of exclusive probabilities of multiple events occurring assuming perfect positive dependence.
        /// The method calculates the exclusive probabilities for each event combination using the inclusion-exclusion principle and convergence checks.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations used to define the event groupings.</param>
        /// <param name="indicators">A 2D array of indicators where 0 means the event did not occur and 1 means the event did occur.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities for each combination.</param>
        /// <param name="eventIndicators">Output. A list of event indicators corresponding to each combination.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for convergence of the inclusion-exclusion algorithm. Default is 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for convergence of the inclusion-exclusion algorithm. Default is 1E-4.</param>
        /// <exception cref="ArgumentException">Thrown if any array is null or contains invalid values.</exception>
        /// <remarks>
        /// This method uses the inclusion-exclusion principle to evaluate the union of event combinations, with convergence checks to avoid unnecessary calculations.
        /// </remarks>
        public static void PositivelyDependentExclusive(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.GetLength(1))
                throw new ArgumentException("The probabilities array and the indicator array must have the same length.", nameof(probabilities));

            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();

            double union = 0;
            double s = 1; // Sign for inclusion-exclusion
            int j = 0; // Binomial combination index
            int c = binomialCombinations[j]; // Current combination limit
            double inc = double.NaN; // Temporary variable for inclusion value
            double exc = double.NaN; // Temporary variable for exclusion value
            int numIndicators = indicators.GetLength(0);  // Cache the number of rows for efficiency

            // Loop through each row of indicators (each event combination)
            for (int i = 0; i < numIndicators; i++)
            {
                // Check if the current index matches the binomial combination threshold
                if (i == c)
                {
                    // Set inc and exc when transitioning between inclusion-exclusion steps
                    if (j > 0)
                    {
                        if (s == 1) inc = union;
                        else if (s == -1) exc = union;
                    }

                    // Check for convergence based on tolerance
                    double diff = Math.Abs(inc - exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        eventIndicators.Add(indicators.GetRow(indicators.GetLength(0) - 1)); // Add last indicator row
                        eventProbabilities.Add(0.5 * diff); // Add the average of the difference to the event probabilities
                        return; // Exit early when convergence is reached
                    }

                    // Flip the sign for the next inclusion-exclusion term
                    s *= -1;
                    j++;
                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j]; // Update the current combination limit
                    }

                }

                // Record the current indicators
                eventIndicators.Add(indicators.GetRow(i));

                // Compute the exclusive event probability and add to the list
                eventProbabilities.Add(PositivelyDependentExclusive(probabilities, eventIndicators.Last()));

                // Calculate the union of probabilities (inclusion-exclusion)
                if (i < probabilities.Count)
                {
                    union += s * probabilities[i];
                }
                else
                {
                    union += s * PositiveJointProbability(probabilities, eventIndicators.Last());
                }

            }

        }

        #endregion

        #region Any Dependency

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method.
        /// Dependence between events is captured with the PCM method.
        /// </summary>
        /// <param name="probabilities">A list of probabilities for each event. Each probability represents the likelihood of the corresponding event occurring.</param>
        /// <param name="correlationMatrix">The correlation matrix that defines the dependency between the events.</param>
        /// <returns>An array of exclusive probabilities for each combination of events, assuming the dependence described by the correlation matrix.</returns>
        /// <exception cref="ArgumentException">Thrown if the probabilities or correlation matrix is null, or if their lengths do not match.</exception>
        public static double[] ExclusivePCM(IList<double> probabilities, double[,] correlationMatrix)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (correlationMatrix == null || correlationMatrix.GetLength(0) != probabilities.Count)
                throw new ArgumentException("The correlation matrix must be square and match the length of the probabilities array.", nameof(correlationMatrix));

            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var binomialCombinations = new int[N];
            for (int i = 1; i <= N; i++)
            {
                binomialCombinations[i - 1] = (int)Factorial.BinomialCoefficient(N, i);
            }

            // Get combination indicators
            var indicators = Factorial.AllCombinations(N);

            // Call the second method to calculate the exclusive probabilities
            return ExclusivePCM(probabilities, binomialCombinations, indicators, correlationMatrix);
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method.
        /// Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">A list of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations representing the number of possible event combinations for each subset.</param>
        /// <param name="indicators">A 2D array of indicators where 0 means the event did not occur and 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency between the events.</param>
        /// <returns>An array of exclusive probabilities of the events based on the inclusion-exclusion method and dependency captured by the correlation matrix.</returns>
        /// <exception cref="ArgumentException">Thrown if the probabilities, binomialCombinations, indicators, or correlationMatrix are invalid.</exception>
        public static double[] ExclusivePCM(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, double[,] correlationMatrix)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (binomialCombinations == null || binomialCombinations.Length == 0)
                throw new ArgumentException("The binomial combinations array must have a length greater than 0.", nameof(binomialCombinations));
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (correlationMatrix == null || correlationMatrix.GetLength(0) != probabilities.Count)
                throw new ArgumentException("The correlation matrix must be square and match the length of the probabilities array.", nameof(correlationMatrix));

            // Calculate cumulative combinations
            int N = probabilities.Count;
            var cumCombos = new int[N - 1];
            cumCombos[0] = binomialCombinations[0];
            for (int i = 1; i < N - 1; i++)
            {
                cumCombos[i] = cumCombos[i - 1] + binomialCombinations[i];
            }

            // Get joint probabilities
            var pVals = JointProbabilitiesPCM(probabilities, indicators, correlationMatrix);

            var result = new double[indicators.GetLength(0)];
            int j = 0;
            int c = binomialCombinations[j];

            // Inclusion-exclusion loop
            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                if (i == c)
                {
                    j++;
                    if (j < binomialCombinations.Length) c += binomialCombinations[j];
                }

                result[i] = pVals[i];
                double s = 1;
                for (int k = j; k < cumCombos.Length; k++)
                {
                    s *= -1;
                    int c1 = cumCombos[k];
                    int c2 = k == cumCombos.Length - 1 ? cumCombos[k] + 1 : cumCombos[k + 1];
                    var sP = SumSearch(pVals, indicators.GetRow(i), indicators, c1, c2);
                    result[i] += s * sP;
                }

                // Correct for floating point issues
                if (result[i] < 0d) result[i] = 0d;
            }

            return result;
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method.
        /// Dependence between events is captured with the PCM method.
        /// This method includes tolerance checks for early termination of calculations if convergence is reached.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">A 2D array of indicators, where 0 means the event did not occur and 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency between events.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities for each event combination.</param>
        /// <param name="eventIndicators">Output. A list of event indicators corresponding to each event combination.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for convergence of the inclusion-exclusion algorithm. Default is 1E-8.</param>
        /// <param name="relativeTolerance">The relative tolerance for convergence of the inclusion-exclusion algorithm. Default is 1E-4.</param>
        /// <returns>A list of exclusive probabilities of the events based on the inclusion-exclusion method with early convergence checks.</returns>
        /// <exception cref="ArgumentException">Thrown if any array is invalid or if their lengths do not match.</exception>
        public static void ExclusivePCM(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, double[,] correlationMatrix,
                                        out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            var jointProbabilities = new List<double>();
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();
            var union = UnionPCM(probabilities, binomialCombinations, indicators, correlationMatrix, out jointProbabilities, out eventIndicators, absoluteTolerance, relativeTolerance);

            // Validate input arrays
            int N = probabilities.Count;
            var cumCombos = new List<int>();
            cumCombos.Add(binomialCombinations[0]);
            for (int i = 1; i < N - 1; i++)
            {
                cumCombos.Add(cumCombos[i - 1] + binomialCombinations[i]);
                if (cumCombos[i] > eventIndicators.Count)
                {
                    cumCombos.RemoveAt(i);
                    break;
                }
            }

            // Get joint probabilities
            var pVals = jointProbabilities.ToArray();

            var result = new double[eventIndicators.Count];
            int j = 0;
            int c = binomialCombinations[j];

            // Loop through event indicators
            for (int i = 0; i < eventIndicators.Count; i++)
            {
                if (i == c)
                {
                    j++;
                    if (j < binomialCombinations.Length) c += binomialCombinations[j];
                }

                result[i] = pVals[i];
                double s = 1;
                for (int k = j; k < cumCombos.Count; k++)
                {
                    s *= -1;
                    int c1 = cumCombos[k];
                    int c2 = k == cumCombos.Count - 1 ? cumCombos[k] + 1 : cumCombos[k + 1];
                    var sP = SumSearch(jointProbabilities, eventIndicators[i], eventIndicators, c1, c2);
                    result[i] += s * sP;
                }

                // Correct for floating point issues
                if (result[i] < 0d) result[i] = 0d;
            }
            eventProbabilities = result.ToList();
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. 
        /// Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities for each event.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        /// <returns>An array of exclusive probabilities for multiple events assuming inclusion-exclusion with MVN dependence.</returns>
        /// <exception cref="ArgumentException">Thrown if any input parameter is invalid.</exception>
        public static double[] ExclusiveMVN(IList<double> probabilities, MultivariateNormal multivariateNormal)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must not be null or empty.", nameof(probabilities));
            if (multivariateNormal == null)
                throw new ArgumentException("The multivariate normal distribution must not be null.", nameof(multivariateNormal));

            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var binomialCombinations = new int[N];
            var cumCombos = new int[N - 1];
            binomialCombinations[0] = (int)Factorial.BinomialCoefficient(N, 1);
            cumCombos[0] = binomialCombinations[0];

            // Compute binomial combinations and cumulative combinations
            for (int i = 1; i < N - 1; i++)
            {
                binomialCombinations[i] = (int)Factorial.BinomialCoefficient(N, i + 1);
                cumCombos[i] = cumCombos[i - 1] + binomialCombinations[i];
            }

            // Get combination indicators (all possible event combinations)
            var indicators = Factorial.AllCombinations(N);

            // Get joint probabilities for each combination
            var pVals = JointProbabilitiesMVN(probabilities, indicators, multivariateNormal);

            var result = new double[indicators.GetLength(0)];
            int j = 0;
            int c = binomialCombinations[j];

            // Iterate through the indicator combinations and apply inclusion-exclusion method
            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                if (i == c)
                {
                    j++;
                    if (j < binomialCombinations.Length) c += binomialCombinations[j];
                }

                result[i] = pVals[i];
                double s = 1;
                for (int k = j; k < cumCombos.Length; k++)
                {
                    s *= -1;
                    int c1 = cumCombos[k];
                    int c2 = k == cumCombos.Length - 1 ? cumCombos[k] + 1 : cumCombos[k + 1];
                    var sP = SumSearch(pVals, indicators.GetRow(i), indicators, c1, c2);
                    result[i] += s * sP;
                }

                // Correct small negative values due to floating point precision issues
                if (result[i] < 0d) result[i] = 0d;
            }

            return result;
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. 
        /// Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations (number of possible combinations per subset).</param>
        /// <param name="indicators">A 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        /// <returns>An array of exclusive probabilities for each combination of events, applying the inclusion-exclusion method with MVN dependence.</returns>
        /// <exception cref="ArgumentException">Thrown if any input parameter is invalid.</exception>
        public static double[] ExclusiveMVN(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, MultivariateNormal multivariateNormal)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must not be null or empty.", nameof(probabilities));
            if (binomialCombinations == null || binomialCombinations.Length == 0)
                throw new ArgumentException("The binomial combinations array must not be null or empty.", nameof(binomialCombinations));
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must not be null or empty.", nameof(indicators));
            if (multivariateNormal == null)
                throw new ArgumentException("The multivariate normal distribution must not be null.", nameof(multivariateNormal));

            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var cumCombos = new int[N - 1];
            cumCombos[0] = binomialCombinations[0];
            for (int i = 1; i < N - 1; i++)
            {
                cumCombos[i] = cumCombos[i - 1] + binomialCombinations[i];
            }

            // Get joint probabilities for each combination
            var pVals = JointProbabilitiesMVN(probabilities, indicators, multivariateNormal);

            var result = new double[indicators.GetLength(0)];
            int j = 0;
            int c = binomialCombinations[j];

            // Iterate through the indicator combinations and apply inclusion-exclusion method
            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                if (i == c)
                {
                    j++;
                    if (j < binomialCombinations.Length) c += binomialCombinations[j];
                }

                result[i] = pVals[i];
                double s = 1;
                for (int k = j; k < cumCombos.Length; k++)
                {
                    s *= -1;
                    int c1 = cumCombos[k];
                    int c2 = k == cumCombos.Length - 1 ? cumCombos[k] + 1 : cumCombos[k + 1];
                    var sP = SumSearch(pVals, indicators.GetRow(i), indicators, c1, c2);
                    result[i] += s * sP;
                }

                // Correct small negative values due to floating point precision issues
                if (result[i] < 0d) result[i] = 0d;
            }

            return result;
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. 
        /// Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">A list of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations for event subsets.</param>
        /// <param name="indicators">A 2D array of indicator values (0 or 1) for each event combination.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute joint probabilities.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities.</param>
        /// <param name="eventIndicators">Output. A list of exclusive event indicators that were evaluated.</param>
        /// <param name="absoluteTol">The absolute tolerance for convergence evaluation. Default is 1E-8.</param>
        /// <param name="relativeTol">The relative tolerance for convergence evaluation. Default is 1E-4.</param>
        /// <exception cref="ArgumentException">Thrown if any input parameter is invalid.</exception>
        public static void ExclusiveMVN(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, MultivariateNormal multivariateNormal, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTol = 1E-4, double relativeTol = 1E-4)
        {
            // Validate input parameters
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must not be null or empty.", nameof(probabilities));
            if (binomialCombinations == null || binomialCombinations.Length == 0)
                throw new ArgumentException("The binomial combinations array must not be null or empty.", nameof(binomialCombinations));
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must not be null or empty.", nameof(indicators));
            if (multivariateNormal == null)
                throw new ArgumentException("The multivariate normal distribution must not be null.", nameof(multivariateNormal));

            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var cumCombos = new int[N - 1];
            cumCombos[0] = binomialCombinations[0];
            for (int i = 1; i < N - 1; i++)
            {
                cumCombos[i] = cumCombos[i - 1] + binomialCombinations[i];
            }

            // Initialize output lists
            var jointProbabilities = new List<double>();
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();

            double union = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];
            double inc = double.NaN;
            double exc = double.NaN;

            // Iterate through each indicator and calculate the union probabilities
            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                if (i == c)
                {
                    if (j > 0 && s == 1)
                    {
                        inc = union;
                    }
                    else if (j > 0 && s == -1)
                    {
                        exc = union;
                    }

                    // Check for convergence
                    double diff = Math.Abs(inc - exc);
                    double tol = absoluteTol + relativeTol * Math.Min(inc, exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= tol)
                    {
                        eventIndicators.Add(indicators.GetRow(indicators.GetLength(0) - 1));
                        jointProbabilities.Add(0.5 * diff);
                        goto Exclusive;
                    }

                    s *= -1;
                    j++;
                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j];
                    }
                }

                // Record event indicators
                eventIndicators.Add(indicators.GetRow(i));

                // Compute union probability
                if (i < probabilities.Count)
                {
                    jointProbabilities.Add(probabilities[i]);
                    union += s * jointProbabilities.Last();
                }
                else
                {
                    jointProbabilities.Add(JointProbabilityMVN(probabilities, eventIndicators.Last(), multivariateNormal));
                    union += s * jointProbabilities.Last();
                }
            }

        Exclusive:

            // Recalculate exclusive event probabilities
            j = 0;
            c = binomialCombinations[j];

            for (int i = 0; i < eventIndicators.Count; i++)
            {
                if (i == c)
                {
                    j++;
                    if (j < binomialCombinations.Length) c += binomialCombinations[j];
                }

                double prob = jointProbabilities[i];
                s = 1;
                for (int k = j; k < cumCombos.Length; k++)
                {
                    s *= -1;
                    int c1 = cumCombos[k];
                    int c2 = k == cumCombos.Length - 1 ? cumCombos[k] + 1 : cumCombos[k + 1];
                    if (c2 >= eventIndicators.Count - 1)
                    {
                        break;
                    }
                    else
                    {
                        var sP = SumSearch(jointProbabilities, eventIndicators[i], eventIndicators, c1, c2);
                        prob += s * sP;
                    }
                }

                // Correct small negative values due to floating point precision issues
                if (prob < 0d) prob = 0d;
                eventProbabilities.Add(prob);
            }
        }

        /// <summary>
        /// Computes the sum of the joint probabilities based on a subset of indicators and their indices.
        /// </summary>
        /// <param name="probabilityValues">An array of probability values to sum over.</param>
        /// <param name="indicatorValues">An array of indicator values representing the current event combination.</param>
        /// <param name="indicators">A 2D array of all indicators representing different event combinations.</param>
        /// <param name="startIndex">The start index for summing the joint probabilities.</param>
        /// <param name="endIndex">The end index for summing the joint probabilities.</param>
        /// <returns>The summed joint probabilities for the specified range of event combinations.</returns>
        private static double SumSearch(double[] probabilityValues, int[] indicatorValues, int[,] indicators, int startIndex, int endIndex)
        {
            double result = 0;
            var indices = new List<int>(indicatorValues.Length);
            // Collect indices of events that occurred
            for (int i = 0; i < indicatorValues.Length; i++)
            {
                if (indicatorValues[i] == 1)
                {
                    indices.Add(i);
                }
            }
            // Iterate through the specified range and sum the joint probabilities
            for (int i = startIndex; i < endIndex; i++)
            {
                bool inclusive = true;
                for (int j = 0; j < indices.Count; j++)
                {
                    if (indicators[i, indices[j]] == 0)
                    {
                        inclusive = false;
                        break;
                    }
                }
                if (inclusive) result += probabilityValues[i];
            }
            return result;
        }

        /// <summary>
        /// Computes the sum of joint probabilities based on a subset of indicators and their indices for event combinations.
        /// </summary>
        /// <param name="probabilityValues">A list of probability values for each event combination.</param>
        /// <param name="indicatorValues">An array of indicator values representing the current event combination (1 for event occurrence, 0 for non-occurrence).</param>
        /// <param name="indicators">A list of arrays representing all possible indicator combinations for the events.</param>
        /// <param name="startIndex">The start index for summing the joint probabilities.</param>
        /// <param name="endIndex">The end index for summing the joint probabilities.</param>
        /// <returns>The summed joint probabilities for the specified range of event combinations.</returns>
        private static double SumSearch(List<double> probabilityValues, int[] indicatorValues, List<int[]> indicators, int startIndex, int endIndex)
        {
            double result = 0;

            // Create a list of indices for the events that are marked as 1 (indicating occurrence)
            var indices = new List<int>(indicatorValues.Length);
            for (int i = 0; i < indicatorValues.Length; i++)
            {
                if (indicatorValues[i] == 1)
                {
                    indices.Add(i);
                }
            }

            // Iterate over the range specified by startIndex and endIndex
            for (int i = startIndex; i < endIndex; i++)
            {
                bool inclusive = true;

                // Check if the current combination includes all the required indicators
                for (int j = 0; j < indices.Count; j++)
                {
                    if (indicators[i][indices[j]] == 0)
                    {
                        inclusive = false;
                        break;
                    }
                }

                // If all required indicators are present, add the corresponding probability value to the result
                if (inclusive) result += probabilityValues[i];
                
            }
            return result;
        }

        #endregion

        #region Common Cause Adjustment

        /// <summary>
        /// Computes the common cause adjustment factor.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <returns>The common cause adjustment factor.</returns>
        public static double CommonCauseAdjustment(IList<double> probabilities)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (probabilities.Count == 1) return 1d;

            double numerator = 1d;
            double denominator = 0d;
            for (int i = 0; i < probabilities.Count; i++)
            {
                numerator *= 1d - probabilities[i];
                denominator += probabilities[i];
            }
            if (denominator == 0) return 1d;
            return (1d - numerator) / denominator;
        }

        /// <summary>
        /// Computes the common cause adjustment factor.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="dependency">The dependency type. Default = Correlation matrix.</param>
        /// <returns>The common cause adjustment factor.</returns>
        public static double CommonCauseAdjustment(IList<double> probabilities, double[,]? correlationMatrix = null, DependencyType dependency = DependencyType.CorrelationMatrix)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (correlationMatrix == null)
                throw new ArgumentNullException(nameof(correlationMatrix), "The correlation matrix cannot be null.");
            if (probabilities.Count == 1) return 1d;

            var indicators = new int[probabilities.Count];
            var complement = new double[probabilities.Count];
            double denominator = 0;
            for (int i = 0;i < probabilities.Count; i++)
            {
                indicators[i] = 1;
                complement[i] = 1 - probabilities[i];
                denominator += probabilities[i];
            }
            if (denominator == 0) return 1d;
            double numerator = JointProbability(complement, indicators, correlationMatrix, dependency);
            return (1d - numerator) / denominator;
        }

        /// <summary>
        /// Computes the mutually exclusive adjustment factor.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <returns>The mutually exclusive adjustment factor.</returns> 
        public static double MutuallyExclusiveAdjustment(IList<double> probabilities)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (probabilities.Count == 1) return 1d;

            double numerator = 1d;
            double denominator = 0d;
            for (int i = 0; i < probabilities.Count; i++)
                denominator += probabilities[i];
            if (denominator <= 1) return 1d;
            return numerator / denominator;
        }

        #endregion

        #endregion

    }
}
