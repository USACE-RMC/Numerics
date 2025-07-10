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

        #region Basic Probability Rules for Two Random Variables

        /// <summary>
        /// Returns the probability of intersection (or joint probability) of A and B, P(A and B). 
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double AAndB(double A, double B, double rho = 0d)
        {
            if (Math.Abs(rho) <= 1E-3) return A * B;
            if (rho >= 0.999) return Math.Min(A, B);
            return MultivariateNormal.BivariateCDF(Normal.StandardZ(1 - A), Normal.StandardZ(1 - B), rho);
        }

        /// <summary>
        /// Returns the probability of union of A and B, P(A or B).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double AOrB(double A, double B, double rho = 0d)
        {
            return A + B - AAndB(A, B, rho);
        }

        /// <summary>
        /// Returns the probability of A and not B, P(A and not B).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double ANotB(double A, double B, double rho = 0d)
        {
            return A - AAndB(A, B, rho);
        }

        /// <summary>
        /// Returns the probability of B and not A, P(B and not A).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double BNotA(double A, double B, double rho = 0d)
        {
            return B - AAndB(A, B, rho);
        }

        /// <summary>
        /// Returns the probability of A given B, P(A|B).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double AGivenB(double A, double B, double rho = 0d)
        {
            return AAndB(A, B, rho) / B;
        }

        /// <summary>
        /// Returns the probability of B given A, P(B|A).
        /// </summary>
        /// <param name="A">Marginal probability of A.</param>
        /// <param name="B">Marginal probability of B.</param>
        /// <param name="rho">Pearson's correlation coefficient. Default = 0.</param>
        public static double BGivenA(double A, double B, double rho = 0d)
        {
            return AAndB(A, B, rho) / A;
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
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="dependency">The dependency type. Default = Correlation matrix.</param>
        /// <returns>The joint probability.</returns>
        public static double JointProbability(IList<double> probabilities, int[] indicators, double[,] correlationMatrix = null, DependencyType dependency = DependencyType.CorrelationMatrix)
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
            double p = 1;
            for (int i = 0; i < probabilities.Count; i++)
            {
                p *= probabilities[i];
            }
            return p;
        }

        /// <summary>
        /// Returns the joint probability assuming perfect independence. 
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double IndependentJointProbability(IList<double> probabilities, int[] indicators)
        {
            return Tools.Product(probabilities, indicators);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect positive dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double PositiveJointProbability(IList<double> probabilities)
        {
            return Tools.Min(probabilities);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect positive dependence. 
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double PositiveJointProbability(IList<double> probabilities, int[] indicators)
        {
            return Tools.Min(probabilities, indicators);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect negative dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double NegativeJointProbability(IList<double> probabilities)
        {
            return Math.Max(0, Math.Min(1, Tools.Sum(probabilities)) - 1);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect negative dependence. 
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double NegativeJointProbability(IList<double> probabilities, int[] indicators)
        {
            return Math.Max(0, Math.Min(1, Tools.Sum(probabilities, indicators)) - 1);
        }

        /// <summary>
        /// Returns the joint probability of multiple events with dependency, using the Product of Conditional Marginals (PCM) method.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="conditionalProbabilities">Returns the array of conditional probabilities for each event.</param>
        public static double JointProbabilityHPCM(IList<double> probabilities, int[] indicators, double[,] correlationMatrix, double[] conditionalProbabilities = null)
        {
            // Get z-values
            int n = probabilities.Count;
            var R = new double[n, n];
            Array.Copy(correlationMatrix, R, correlationMatrix.Length);
            int i, j, k, ir, ic;
            double pdf, cdf, A, B, z1, z2, r12, z21, p21, jp;
            for (i = 0; i < n; i++)
            {
                if (indicators[i] == 0)
                {
                    R[i, i] = Normal.StandardZ(1);
                }
                else
                {
                    R[i, i] = Normal.StandardZ(probabilities[i]);
                }
            }
            // Update the conditional correlation matrix
            // First cycle
            z1 = R[0, 0];
            pdf = Normal.StandardPDF(z1);
            cdf = Normal.StandardCDF(z1);
            A = pdf / cdf;
            B = A * (z1 + A);
            for (k = 1; k < n; k++)
            {
                z2 = R[k, k];
                r12 = R[0, k];
                r12 = Math.Abs(r12) < 1E-3 ? 0: r12;
                p21 = MultivariateNormal.BivariateCDF(-z1, -z2, r12) / cdf;
                p21 = Math.Max(0, Math.Min(1, p21));
                z21 = Normal.StandardZ(p21);
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
                A = pdf / cdf;
                B = A * (z1 + A);
                for (k = j + 1; k < n; k++)
                {
                    z2 = R[k, j - 1];
                    r12 = R[j, k];
                    r12 = Math.Abs(r12) < 1E-3 ? 0 : r12;
                    p21 = MultivariateNormal.BivariateCDF(-z1, -z2, r12) / cdf;
                    p21 = Math.Max(0, Math.Min(1, p21));
                    z21 = Normal.StandardZ(p21);
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
        /// Returns the joint probability of multiple events with dependency, using the Product of Conditional Marginals (PCM) method.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="conditionalProbabilities">Returns the array of conditional probabilities for each event.</param>
        public static double JointProbabilityPCM(IList<double> probabilities, int[] indicators, double[,] correlationMatrix, double[] conditionalProbabilities = null)
        {
            // Get z-values
            int n = probabilities.Count;
            var R = new double[n, n];
            Array.Copy(correlationMatrix, R, correlationMatrix.Length);
            int i, j, k, ir, ic;
            double A, B, z1, z2, z21, r12;
            for (i = 0; i < n; i++)
            {
                if (indicators[i] == 0)
                {
                    R[i, i] = Normal.StandardZ(1);
                }
                else
                {
                    R[i, i] = Normal.StandardZ(probabilities[i]);
                }
            }
            // Update the conditional correlation matrix
            // First cycle
            z1 = R[0, 0];
            A = Normal.StandardPDF(z1) / Normal.StandardCDF(z1);
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
                A = Normal.StandardPDF(z1) /  Normal.StandardCDF(z1);
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
        public static double JointProbability(IList<double> probabilities, int[] indicators, MultivariateNormal multivariateNormal)
        {
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
        public static double[] JointProbabilities(IList<double> probabilities, int[,] indicators, MultivariateNormal multivariateNormal)
        {
            var result = new double[indicators.GetLength(0)];

            Parallel.For(0, indicators.GetLength(0), idx =>
            {
                if (idx < probabilities.Count)
                {
                    result[idx] = probabilities[idx];
                }
                else
                {
                    result[idx] = JointProbability(probabilities, indicators.GetRow(idx), (MultivariateNormal)multivariateNormal.Clone());
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
            if (probabilities.Count == 1) return probabilities[0];
            double numerator = 1d;
            for (int i = 0; i < probabilities.Count; i++)
                numerator *= 1d - probabilities[i];
            return 1d - numerator;
        }

        /// <summary>
        /// Returns the unimodal bound for the probability of union assuming perfect positive dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double PositivelyDependentUnion(IList<double> probabilities)
        {
            return Tools.Max(probabilities);
        }

        /// <summary>
        /// Returns the unimodal bound for the probability of union assuming perfect negative dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double NegativelyDependentUnion(IList<double> probabilities)
        {
            return Math.Min(1, Tools.Sum(probabilities));
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
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        public static double UnionPCM(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, double[,] correlationMatrix, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {

            double result = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];
            double inc = double.NaN;
            double exc = double.NaN;

            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                if (i == c)
                {
                    if (j > 0 && s == 1)
                    {
                        inc = result;
                    }
                    else if (j > 0 && s == -1)
                    {
                        exc = result;
                    }

                    // Check for convergence
                    // Add the tolerance to the last joint event
                    double diff = Math.Abs(inc - exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        return result + 0.5 * diff;
                    }

                    s *= -1;
                    j++;
                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j];
                    }

                }
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
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities.</param>
        /// <param name="eventIndicators">Output. A list of exclusive event indicators that were evaluated.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        public static double UnionPCM(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, double[,] correlationMatrix, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();

            double union = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];
            double inc = double.NaN;
            double exc = double.NaN;

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
                    // Add the tolerance to the last joint event
                    double diff = Math.Abs(inc - exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        eventIndicators.Add(indicators.GetRow(indicators.GetLength(0) - 1));
                        eventProbabilities.Add(0.5 * diff);
                        return union + 0.5 * diff;
                    }

                    s *= -1;
                    j++;
                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j];
                    }

                }

                // Record indicators
                eventIndicators.Add(indicators.GetRow(i));
                if (i < probabilities.Count)
                {
                    union += s * probabilities[i];
                    eventProbabilities.Add(probabilities[i]);
                }
                else
                {
                    var jp = JointProbability(probabilities, indicators.GetRow(i), correlationMatrix);
                    union += s * jp;
                    eventProbabilities.Add(jp);
                }

            }

            return union;
        }

        /// <summary>
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        public static double Union(IList<double> probabilities, MultivariateNormal multivariateNormal)
        {

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
            return Union(probabilities, binomialCombinations, indicators, multivariateNormal);

        }

        /// <summary>
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        public static double Union(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, MultivariateNormal multivariateNormal)
        {

            double result = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];

            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                if (i == c)
                {
                    s *= -1;
                    j++;
                    if (j < binomialCombinations.Length) c += binomialCombinations[j];
                }
                if (i < probabilities.Count)
                {
                    result += s * probabilities[i];
                }
                else
                {
                    result += s * JointProbability(probabilities, indicators.GetRow(i), multivariateNormal);
                }
            }

            return result;
        }

        #endregion

        #region Exclusive Probability of all Combinations of Events

        #region Independent
        
        /// <summary>
        /// Returns the exclusive probability of multiple events occurring assuming independence.
        /// </summary>
        /// <param name="probabilities">And array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double IndependentExclusive(IList<double> probabilities, int[] indicators)
        {
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
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double[] IndependentExclusive(IList<double> probabilities, int[,] indicators)
        {
            var result = new double[indicators.GetLength(0)];
            for (int i = 0; i < indicators.GetLength(0); i++)
            {             
                result[i] = IndependentExclusive(probabilities, indicators.GetRow(i));
            }
            return result;         
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events occurring assuming independence.
        /// </summary>
        /// <param name="probabilities">And array of probabilities for each event.</param>
        public static double[] IndependentExclusive(IList<double> probabilities)
        {
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
                    result[t] = IndependentExclusive(probabilities, indicators);
                    t++;
                }
            }
            return result;
        }

        /// <summary>
        /// Returns a list of exclusive probabilities of multiple events occurring assuming independence.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities.</param>
        /// <param name="eventIndicators">Output. A list of exclusive event indicators that were evaluated.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        public static void IndependentExclusive(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();

            double union = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];
            double inc = double.NaN;
            double exc = double.NaN;

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
                    // Add the tolerance to the last joint event
                    double diff = Math.Abs(inc - exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        eventIndicators.Add(indicators.GetRow(indicators.GetLength(0) - 1));
                        eventProbabilities.Add(0.5 * diff);
                        return;
                    }

                    s *= -1;
                    j++;
                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j];
                    }

                }

                // Record indicators
                eventIndicators.Add(indicators.GetRow(i));
                // Compute the exclusive event probability
                eventProbabilities.Add(IndependentExclusive(probabilities, eventIndicators.Last()));
                // Compute the union
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
        /// </summary>
        /// <param name="probabilities">And array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double PositivelyDependentExclusive(IList<double> probabilities, int[] indicators)
        {
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
        /// </summary>
        /// <param name="probabilities">And array of probabilities for each event.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        public static double[] PositivelyDependentExclusive(IList<double> probabilities, int[,] indicators)
        {
            var result = new double[indicators.GetLength(0)];
            for (int i = 0; i < indicators.GetLength(0); i++)
            {
                result[i] = PositivelyDependentExclusive(probabilities, indicators.GetRow(i));
            }
            return result;
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events occurring assuming perfect positive dependence.
        /// </summary>
        /// <param name="probabilities">And array of probabilities for each event.</param>
        public static double[] PositivelyDependentExclusive(IList<double> probabilities)
        {
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
        /// Returns a list of exclusive probabilities of multiple events occurring assuming independence.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities.</param>
        /// <param name="eventIndicators">Output. A list of exclusive event indicators that were evaluated.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        public static void PositivelyDependentExclusive(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();

            double union = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];
            double inc = double.NaN;
            double exc = double.NaN;

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
                    // Add the tolerance to the last joint event
                    double diff = Math.Abs(inc - exc);
                    if (j > 0 && j < binomialCombinations.Length && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        eventIndicators.Add(indicators.GetRow(indicators.GetLength(0) - 1));
                        eventProbabilities.Add(0.5 * diff);
                        return;
                    }

                    s *= -1;
                    j++;
                    if (j < binomialCombinations.Length)
                    {
                        c += binomialCombinations[j];
                    }

                }

                // Record indicators
                eventIndicators.Add(indicators.GetRow(i));
                // Compute the exclusive event probability
                eventProbabilities.Add(PositivelyDependentExclusive(probabilities, eventIndicators.Last()));
                // Compute the union
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
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        public static double[] ExclusivePCM(IList<double> probabilities, double[,] correlationMatrix)
        {
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
            return ExclusivePCM(probabilities, binomialCombinations, indicators, correlationMatrix);

        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        public static double[] ExclusivePCM(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, double[,] correlationMatrix)
        {

            // Get number of unique combinations by subset
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
                // Due to floating point issues, this can sometimes be very small, but negative. 
                if (result[i] < 0d) result[i] = 0d;
            }

            return result;
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities.</param>
        /// <param name="eventIndicators">Output. A list of exclusive event indicators that were evaluated.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-8.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        public static void ExclusivePCM(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, double[,] correlationMatrix, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            var jointProbabilities = new List<double>();
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();
            var union = UnionPCM(probabilities, binomialCombinations, indicators, correlationMatrix, out jointProbabilities,  out eventIndicators, absoluteTolerance, relativeTolerance);

            // Get number of unique combinations by subset
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
                // Due to floating point issues, this can sometimes be very small, but negative. 
                if (result[i] < 0d) result[i] = 0d;
            }
            eventProbabilities = result.ToList();
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        public static double[] Exclusive(IList<double> probabilities, MultivariateNormal multivariateNormal)
        {
            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var binomialCombinations = new int[N];
            var cumCombos = new int[N - 1];
            binomialCombinations[0] = (int)Factorial.BinomialCoefficient(N, 1);
            cumCombos[0] = binomialCombinations[0];
            for (int i = 1; i < N - 1; i++)
            {
                binomialCombinations[i] = (int)Factorial.BinomialCoefficient(N, i + 1);
                cumCombos[i] = cumCombos[i - 1] + binomialCombinations[i];
            }

            // Get combination indicators
            var indicators = Factorial.AllCombinations(N);

            // Get joint probabilities
            var pVals = JointProbabilities(probabilities, indicators, multivariateNormal);

            var result = new double[indicators.GetLength(0)];

            int j = 0;
            int c = binomialCombinations[j];

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
                // Due to floating point issues, this can sometimes be very small, but negative. 
                if (result[i] < 0d) result[i] = 0d;
            }

            return result;

        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        public static double[] Exclusive(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, MultivariateNormal multivariateNormal)
        {

            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var cumCombos = new int[N - 1];
            cumCombos[0] = binomialCombinations[0];
            for (int i = 1; i < N - 1; i++)
            {
                cumCombos[i] = cumCombos[i - 1] + binomialCombinations[i];
            }

            // Get joint probabilities
            var pVals = JointProbabilities(probabilities, indicators, multivariateNormal);

            var result = new double[indicators.GetLength(0)];

            int j = 0;
            int c = binomialCombinations[j];

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
                // Due to floating point issues, this can sometimes be very small, but negative. 
                if (result[i] < 0d) result[i] = 0d;
            }

            return result;
        }

        /// <summary>
        /// Returns an array of exclusive probabilities of multiple events using the inclusion-exclusion method. Dependence between events is captured with the multivariate normal distribution.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="binomialCombinations">An array of binomial combinations.</param>
        /// <param name="indicators">An 2D array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <param name="multivariateNormal">The multivariate normal distribution used to compute the joint probabilities.</param>
        /// <param name="eventProbabilities">Output. A list of exclusive event probabilities.</param>
        /// <param name="eventIndicators">Output. A list of exclusive event indicators that were evaluated.</param>
        /// <param name="absoluteTol">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-8.</param>
        /// <param name="relativeTol">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        public static void Exclusive(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, MultivariateNormal multivariateNormal, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTol = 1E-4, double relativeTol = 1E-4)
        {


            // Get number of unique combinations by subset
            int N = probabilities.Count;
            var cumCombos = new int[N - 1];
            cumCombos[0] = binomialCombinations[0];
            for (int i = 1; i < N - 1; i++)
            {
                cumCombos[i] = cumCombos[i - 1] + binomialCombinations[i];
            }

            var jointProbabilities = new List<double>();
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();

            double union = 0;
            double s = 1;
            int j = 0;
            int c = binomialCombinations[j];
            double inc = double.NaN;
            double exc = double.NaN;

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
                    // Add the tolerance to the last joint event
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

                // Record indicators
                eventIndicators.Add(indicators.GetRow(i));

                // Compute the union
                if (i < probabilities.Count)
                {
                    jointProbabilities.Add(probabilities[i]);
                    union += s * jointProbabilities.Last();
                }
                else
                {
                    jointProbabilities.Add(JointProbability(probabilities, eventIndicators.Last(), multivariateNormal));
                    union += s * jointProbabilities.Last();
                }

            }

        Exclusive:

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
                // Due to floating point issues, this can sometimes be very small, but negative. 
                if (prob < 0d) prob = 0d;
                eventProbabilities.Add(prob);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="probabilityValues"></param>
        /// <param name="indicatorValues"></param>
        /// <param name="indicators"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        /// <returns></returns>
        private static double SumSearch(double[] probabilityValues, int[] indicatorValues, int[,] indicators, int startIndex, int endIndex)
        {
            double result = 0;
            var indeces = new List<int>();
            for (int i = 0; i < indicatorValues.Length; i++)
            {
                if (indicatorValues[i] == 1)
                {
                    indeces.Add(i);
                }
            }
            for (int i = startIndex; i < endIndex; i++)
            {
                bool inclusive = true;
                for (int j = 0; j < indeces.Count; j++)
                {
                    if (indicators[i, indeces[j]] == 0)
                    {
                        inclusive = false;
                        break;
                    }
                }
                if (inclusive == true) result += probabilityValues[i];
            }
            return result;
        }

        private static double SumSearch(List<double> probabilityValues, int[] indicatorValues, List<int[]> indicators, int startIndex, int endIndex)
        {
            double result = 0;
            var indeces = new List<int>();
            for (int i = 0; i < indicatorValues.Length; i++)
            {
                if (indicatorValues[i] == 1)
                {
                    indeces.Add(i);
                }
            }
            for (int i = startIndex; i < endIndex; i++)
            {
                bool inclusive = true;
                for (int j = 0; j < indeces.Count; j++)
                {
                    if (indicators[i][indeces[j]] == 0)
                    {
                        inclusive = false;
                        break;
                    }
                }
                if (inclusive == true) result += probabilityValues[i];
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
        public static double CommonCauseAdjustment(IList<double> probabilities, double[,] correlationMatrix = null, DependencyType dependency = DependencyType.CorrelationMatrix)
        {
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
