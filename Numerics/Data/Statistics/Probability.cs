using Numerics.Distributions;
using Numerics.Mathematics.Integration;
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
            if (A == 1d) return Tools.Clamp(B, 0d, 1d);
            if (B == 1d) return Tools.Clamp(A, 0d, 1d);
            if (rho <= -0.999) return Tools.Clamp(A + B - 1d, 0d, 1d);
            if (rho >= 0.999) return Tools.Clamp(Math.Min(A, B), 0d, 1d);
            if (Math.Abs(rho) <= 1E-3) return Tools.Clamp(A * B, 0d, 1d);
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
        /// <returns>The Fréchet–Hoeffding lower bound max(0, Σpᵢ − (n − 1)), where n is the number of events.</returns>
        /// <remarks>
        /// For two events with probabilities 0.8 and 0.9 the joint probability is 0.7: under perfect
        /// negative dependence the events overlap only by the amount their total probability exceeds one.
        /// When the probabilities sum to one or less, perfectly negatively dependent events are disjoint
        /// and the joint probability is zero.
        /// </remarks>
        public static double NegativeJointProbability(IList<double> probabilities)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            return Tools.Clamp(Tools.Sum(probabilities) - (probabilities.Count - 1d), 0d, 1d);
        }

        /// <summary>
        /// Returns the joint probability assuming perfect negative dependence.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event.</param>
        /// <param name="indicators">An array of indicators, 0 means the event did not occur, 1 means the event did occur.</param>
        /// <returns>The Fréchet–Hoeffding lower bound max(0, Σpᵢ − (k − 1)) over the indicated events, where k is the number of indicated events.</returns>
        /// <remarks>
        /// Only the events whose indicator is 1 participate, matching the other joint-probability
        /// overloads. With no indicated events the joint probability of the empty intersection is one.
        /// </remarks>
        public static double NegativeJointProbability(IList<double> probabilities, int[] indicators)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.Length == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.Length)
                throw new ArgumentException("The probabilities and indicators arrays must have the same length.", nameof(probabilities));
            int indicated = 0;
            for (int i = 0; i < indicators.Length; i++)
                if (indicators[i] == 1) indicated++;
            return Tools.Clamp(Tools.Sum(probabilities, indicators) - (indicated - 1d), 0d, 1d);
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
            const double minimumCdf = 1E-300;
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
            if (cdf < minimumCdf) cdf = minimumCdf;
            A = pdf / cdf;
            B = A * (z1 + A);
            for (k = 1; k < n; k++)
            {
                z2 = R[k, k];
                r12 = R[0, k];
                r12 = Math.Abs(r12) < 1E-3 ? 0: r12;
                p21 = MultivariateNormal.BivariateCDF(-z1, -z2, r12) / cdf;
                p21 = Tools.Clamp(p21, 0d, 1d);
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
                if (cdf < minimumCdf) cdf = minimumCdf;
                A = pdf / cdf;
                B = A * (z1 + A);
                for (k = j + 1; k < n; k++)
                {
                    z2 = R[k, j - 1];
                    r12 = R[j, k];
                    r12 = Math.Abs(r12) < 1E-3 ? 0 : r12;
                    p21 = MultivariateNormal.BivariateCDF(-z1, -z2, r12) / cdf;
                    p21 = Tools.Clamp(p21, 0d, 1d);
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
                conditionalProbabilities[0] = Tools.Clamp(Normal.StandardCDF(R[0, 0]), 0d, 1d);
            for (i = 1; i < n; i++)
            {
                jp += Math.Log(Normal.StandardCDF(R[i, i - 1]));
                if (conditionalProbabilities != null && conditionalProbabilities.Length == n)
                    conditionalProbabilities[i] = Tools.Clamp(Normal.StandardCDF(R[i, i - 1]), 0d, 1d);
            }
            jp = Math.Exp(jp);
            jp = Tools.Clamp(jp, 0d, 1d);
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
                conditionalProbabilities[0] = Tools.Clamp(Normal.StandardCDF(R[0, 0]), 0d, 1d);
            for (i = 1; i < n; i++)
            {
                jp += Math.Log(Normal.StandardCDF(R[i, i - 1]));
                if (conditionalProbabilities != null && conditionalProbabilities.Length == n)
                    conditionalProbabilities[i] = Tools.Clamp(Normal.StandardCDF(R[i, i - 1]), 0d, 1d);
            }
            jp = Math.Exp(jp);
            jp = Tools.Clamp(jp, 0d, 1d);
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
                    result[idx] = Tools.Clamp(probabilities[idx], 0d, 1d);
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
            p = Tools.Clamp(p, 0d, 1d);
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
                    result[idx] = Tools.Clamp(probabilities[idx], 0d, 1d);
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
            if (probabilities.Count == 1) return Tools.Clamp(probabilities[0], 0d, 1d);

            double numerator = 1d;
            for (int i = 0; i < probabilities.Count; i++)
            {
                var q = 1d - probabilities[i];
                if (q == 0d) return 1d; // any event certain -> union = 1
                numerator *= q;
            }
            return Tools.Clamp(1d - numerator, 0d, 1d);
        }

        /// <summary>
        /// Returns the unimodal bound for the probability of union assuming perfect positive dependence. 
        /// </summary>
        /// <param name="probabilities">List of probabilities.</param>
        public static double PositivelyDependentUnion(IList<double> probabilities)
        {
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            if (probabilities.Count == 1) return Tools.Clamp(probabilities[0], 0d, 1d);
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
            if (probabilities.Count == 1) return Tools.Clamp(probabilities[0], 0d, 1d);
            return Tools.Clamp(Tools.Sum(probabilities), 0, 1);
        }

        /// <summary>
        /// Returns the probability of union under an equicorrelated single-factor Gaussian
        /// dependence structure, where the events are conditionally independent given one shared
        /// standard normal factor.
        /// </summary>
        /// <param name="probabilities">List of marginal event probabilities.</param>
        /// <param name="rho">The common correlation of the underlying Gaussian variates, within [0, 1].</param>
        /// <param name="relativeTolerance">Optional. The relative tolerance of the quadrature. Default = 1E-8.</param>
        /// <remarks>
        /// <para>
        /// With thresholds b(i) = StandardZ(p(i)) and loading sqrt(rho) on the shared factor z, the
        /// union is P = Integral of phi(z) * [1 - Product(1 - Phi((b(i) - sqrt(rho)*z)/sqrt(1-rho)))] dz.
        /// Conditional independence makes each node O(n), rather than the 2^n inclusion-exclusion of
        /// <see cref="UnionPCM(IList{double}, double[,], double, double)"/> and
        /// <see cref="UnionMVN(IList{double}, MultivariateNormal)"/>, so the method scales to very wide event
        /// sets. The integral is evaluated in probability space u = Phi(z) with adaptive
        /// Gauss-Kronrod quadrature over [1e-16, 1 - 1e-16] (the excluded endpoint slivers bound the
        /// truncation error by 2e-16), and the conditional survival product is accumulated in log
        /// space through <see cref="Tools.Log1p"/> and recovered through <see cref="Tools.Expm1"/>,
        /// so rare-event unions retain relative accuracy. The result is deterministic.
        /// </para>
        /// <para>
        /// The equicorrelated Gaussian correlation matrix is positive semi-definite down to
        /// -1/(n-1), but a negative correlation has no real single-factor loading, so this method
        /// requires rho within [0, 1]; use <see cref="UnionPCM(IList{double}, double[,], double, double)"/> or
        /// <see cref="UnionMVN(IList{double}, MultivariateNormal)"/> for negative dependence. At rho = 1 the comonotone limit, the maximum marginal probability,
        /// is returned analytically; at rho = 0 the events are independent. Zero-probability events
        /// carry no mass and any certain event returns one.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentException">Thrown when the probabilities list is null or empty.</exception>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when a probability is outside [0, 1], the correlation is outside [0, 1], or the
        /// relative tolerance is outside the quadrature's accepted range of [1E-15, 1].
        /// </exception>
        public static double UnionSingleFactor(IList<double> probabilities, double rho, double relativeTolerance = 1E-8)
        {
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities list must be non-null and contain at least one element.");
            for (int i = 0; i < probabilities.Count; i++)
            {
                if (!Tools.IsFinite(probabilities[i]) || probabilities[i] < 0d || probabilities[i] > 1d)
                    throw new ArgumentOutOfRangeException(nameof(probabilities), "Probabilities must be finite and within [0, 1].");
            }
            if (!Tools.IsFinite(rho) || rho < 0d || rho > 1d)
                throw new ArgumentOutOfRangeException(nameof(rho), "The common correlation must be within [0, 1]. A negative equicorrelation (positive semi-definite down to -1/(n-1)) has no real single-factor loading; use UnionPCM or UnionMVN for negative dependence.");

            // Zero-probability events never occur; certainty and the comonotone limit are analytic.
            double maxP = 0d;
            int m = 0;
            for (int i = 0; i < probabilities.Count; i++)
            {
                if (probabilities[i] >= 1d) return 1d;
                if (probabilities[i] > maxP) maxP = probabilities[i];
                if (probabilities[i] > 0d) m++;
            }
            if (m == 0) return 0d;
            if (rho == 1d) return maxP;

            var thresholds = new double[m];
            int j = 0;
            for (int i = 0; i < probabilities.Count; i++)
            {
                if (probabilities[i] > 0d) thresholds[j++] = Normal.StandardZ(probabilities[i]);
            }

            double sqrtRho = Math.Sqrt(rho);
            double sqrtComplement = Math.Sqrt(1d - rho);
            var conditional = new double[m];
            Func<double, double> integrand = u =>
            {
                double z = Normal.StandardZ(u);
                SingleFactorConditionalCore(thresholds, sqrtRho, sqrtComplement, z, conditional);
                double logSurvival = 0d;
                for (int i = 0; i < m; i++)
                    logSurvival += Tools.Log1p(-conditional[i]);
                return -Tools.Expm1(logSurvival);
            };
            var quadrature = new AdaptiveGaussKronrod(integrand, 1E-16, 1d - 1E-16)
            {
                RelativeTolerance = relativeTolerance,
                // The union can be arbitrarily small, and the quadrature accepts an interval on the
                // absolute OR the relative criterion - pin the absolute tolerance at the framework
                // floor so the relative criterion always governs.
                AbsoluteTolerance = 1E-15,
                MinDepth = 2,
                ReportFailure = true
            };
            quadrature.Integrate();
            return Tools.Clamp(quadrature.Result, 0d, 1d);
        }

        /// <summary>
        /// Fills a caller-owned buffer with the conditional event probabilities of the
        /// equicorrelated single-factor Gaussian structure at a given factor value:
        /// Phi((b(i) - sqrt(rho)*z)/sqrt(1-rho)) for each threshold b(i).
        /// </summary>
        /// <param name="normalThresholds">The standard normal thresholds b(i) = StandardZ(p(i)), precomputed by the caller.</param>
        /// <param name="rho">The common correlation of the underlying Gaussian variates, within [0, 1).</param>
        /// <param name="z">The shared standard normal factor value.</param>
        /// <param name="conditional">The caller-owned output buffer, at least as long as the thresholds. Entries beyond the threshold count are left untouched.</param>
        /// <remarks>
        /// Conditional on the factor the events are independent, so downstream combination kernels
        /// built for independent probabilities can run per factor node with an outer quadrature
        /// around them. The method allocates nothing. At rho = 1 the conditionals degenerate to
        /// indicators of z against each threshold; handle that comonotone limit analytically rather
        /// than through this method. Non-finite thresholds propagate through the normal CDF without
        /// validation - the buffer fill is a hot path.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either array is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the buffer is shorter than the thresholds.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the correlation is outside [0, 1) or the factor value is not finite.</exception>
        public static void SingleFactorConditionalProbabilities(double[] normalThresholds, double rho, double z, double[] conditional)
        {
            if (normalThresholds == null) throw new ArgumentNullException(nameof(normalThresholds));
            if (conditional == null) throw new ArgumentNullException(nameof(conditional));
            if (conditional.Length < normalThresholds.Length)
                throw new ArgumentException("The conditional buffer must be at least as long as the thresholds.", nameof(conditional));
            if (!Tools.IsFinite(rho) || rho < 0d || rho >= 1d)
                throw new ArgumentOutOfRangeException(nameof(rho), "The common correlation must be within [0, 1). At rho = 1 the conditional probabilities degenerate to indicators; handle the comonotone limit analytically.");
            if (!Tools.IsFinite(z))
                throw new ArgumentOutOfRangeException(nameof(z), "The factor value must be finite.");
            SingleFactorConditionalCore(normalThresholds, Math.Sqrt(rho), Math.Sqrt(1d - rho), z, conditional);
        }

        /// <summary>
        /// The unguarded conditional-probability fill shared by the union quadrature and the public
        /// buffer helper.
        /// </summary>
        /// <param name="normalThresholds">The standard normal thresholds.</param>
        /// <param name="sqrtRho">The square root of the common correlation.</param>
        /// <param name="sqrtComplement">The square root of one minus the common correlation.</param>
        /// <param name="z">The shared standard normal factor value.</param>
        /// <param name="conditional">The output buffer.</param>
        private static void SingleFactorConditionalCore(double[] normalThresholds, double sqrtRho, double sqrtComplement, double z, double[] conditional)
        {
            for (int i = 0; i < normalThresholds.Length; i++)
                conditional[i] = Normal.StandardCDF((normalThresholds[i] - sqrtRho * z) / sqrtComplement);
        }

        /// <summary>
        /// Returns the probability of union using the inclusion-exclusion method. Dependence between events is captured with the PCM method.
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

            return UnionPCMLazy(probabilities, correlationMatrix, out _, absoluteTolerance, relativeTolerance);
        }

        /// <summary>
        /// Lazily returns the probability of union using inclusion-exclusion with Product of Conditional Marginals (PCM) dependence.
        /// </summary>
        /// <param name="probabilities">List of marginal event probabilities.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="status">The enumeration completion status.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <returns>The probability of the union of the events.</returns>
        /// <exception cref="ArgumentException">Thrown when the probability collection is null or empty, or the correlation matrix is null.</exception>
        /// <remarks>
        /// Combinations are generated in the subset-size and lexicographic order used by
        /// <see cref="Factorial.AllCombinations(int)"/>. The PCM calculation, alternating
        /// inclusion-exclusion signs, dual convergence predicate, and half-gap closure are
        /// identical to the dense overload.
        /// </remarks>
        public static double UnionPCMLazy(IList<double> probabilities, double[,] correlationMatrix,
            out ExclusiveEnumerationStatus status, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            ValidateLazyPCMInputs(probabilities, correlationMatrix, absoluteTolerance, relativeTolerance);

            int n = probabilities.Count;
            var row = new int[n];
            double union = 0d;
            double sign = 1d;
            double inclusion = double.NaN;
            double exclusion = double.NaN;
            int previousSubsetSize = 0;

            foreach (int[] combination in Factorial.AllCombinationsLazy(n))
            {
                int subsetSize = combination.Length;
                if (subsetSize != previousSubsetSize)
                {
                    previousSubsetSize = subsetSize;
                    if (subsetSize >= 2)
                    {
                        int block = subsetSize - 2;
                        if (block > 0)
                        {
                            if (sign == 1d) inclusion = union;
                            else if (sign == -1d) exclusion = union;
                        }

                        double difference = Math.Abs(inclusion - exclusion);
                        if (block > 0 && block < n &&
                            difference <= absoluteTolerance &&
                            difference <= relativeTolerance * Math.Min(inclusion, exclusion))
                        {
                            status = ExclusiveEnumerationStatus.Converged;
                            return Tools.Clamp(union + 0.5d * difference, 0d, 1d);
                        }

                        sign *= -1d;
                    }
                }

                Array.Clear(row, 0, row.Length);
                for (int i = 0; i < combination.Length; i++) row[combination[i]] = 1;
                double jointProbability = subsetSize == 1
                    ? probabilities[combination[0]]
                    : JointProbability(probabilities, row, correlationMatrix);
                union += sign * jointProbability;
            }

            status = ExclusiveEnumerationStatus.Complete;
            return Tools.Clamp(union, 0d, 1d);
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
                        return Tools.Clamp(result + 0.5d * diff, 0d, 1d); // Converged, return the result with half of the difference
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

            return Tools.Clamp(result, 0d, 1d);
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
                        eventProbabilities.Add(Tools.Clamp(0.5d * diff, 0d, 1d));  // Add the averaged difference
                        return Tools.Clamp(union + 0.5d * diff, 0d, 1d); // Converged, return the result with half of the difference
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
                    eventProbabilities.Add(Tools.Clamp(probabilities[i], 0d, 1d));  // Store the probability
                }
                else
                {
                    var jp = JointProbability(probabilities, indicators.GetRow(i), correlationMatrix);  // Otherwise, calculate the joint probability
                    union += s * jp;  // Add the joint probability contribution
                    eventProbabilities.Add(Tools.Clamp(jp, 0d, 1d));  // Store the joint probability
                }

            }

            return Tools.Clamp(union, 0d, 1d);
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

            return Tools.Clamp(result, 0d, 1d);
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
            return Tools.Clamp(result, 0d, 1d);
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
        /// <exception cref="ArgumentNullException">Thrown when the binomial-combination metadata is null.</exception>
        /// <exception cref="ArgumentException">Thrown when probabilities are missing or the combination and indicator metadata is structurally inconsistent.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a probability or tolerance is outside its valid range.</exception>
        /// <exception cref="OverflowException">Thrown when a required combination count exceeds <see cref="int.MaxValue"/>.</exception>
        /// <remarks>
        /// This method uses the inclusion-exclusion principle to compute the exclusive probability of each event combination.
        /// The result is added to a list, and convergence is monitored using the specified tolerances.
        /// </remarks>
        public static void IndependentExclusive(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, out List<double> eventProbabilities, out List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            eventProbabilities = new List<double>();
            eventIndicators = new List<int[]>();
            IndependentExclusive(probabilities, binomialCombinations, indicators, eventProbabilities, eventIndicators, absoluteTolerance, relativeTolerance);
        }

        /// <summary>
        /// Computes independent exclusive probabilities into caller-owned output collections.
        /// Existing indicator arrays with the required length are refilled and reused.
        /// </summary>
        /// <param name="probabilities">The probability of each event.</param>
        /// <param name="binomialCombinations">The number of combinations for each subset size.</param>
        /// <param name="indicators">The event-indicator rows in subset-size order.</param>
        /// <param name="eventProbabilities">The output probabilities; cleared and refilled.</param>
        /// <param name="eventIndicators">The output indicator rows; matching arrays are reused.</param>
        /// <param name="absoluteTolerance">The non-negative absolute convergence tolerance.</param>
        /// <param name="relativeTolerance">The non-negative relative convergence tolerance.</param>
        /// <returns>
        /// <see langword="true"/> when convergence stops the inclusion-exclusion expansion
        /// before every combination is enumerated; otherwise, <see langword="false"/>.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when an output collection or metadata array is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the probability or indicator metadata is structurally inconsistent.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a probability or tolerance is outside its valid range.</exception>
        /// <remarks>
        /// After the output lists reach their required capacity, repeated calls with the same
        /// dimensions reuse the indicator arrays and do not allocate output rows.
        /// </remarks>
        public static bool IndependentExclusive(IList<double> probabilities, int[] binomialCombinations, int[,] indicators, List<double> eventProbabilities, List<int[]> eventIndicators, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            // Validation Checks
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (indicators == null || indicators.GetLength(0) == 0)
                throw new ArgumentException("The indicators array must have at least one row.", nameof(indicators));
            if (probabilities.Count != indicators.GetLength(1))
                throw new ArgumentException("The probabilities array and the indicator array must have the same length.", nameof(probabilities));
            if (eventProbabilities == null) throw new ArgumentNullException(nameof(eventProbabilities));
            if (eventIndicators == null) throw new ArgumentNullException(nameof(eventIndicators));
            ValidatePooledExclusiveMetadata(probabilities, binomialCombinations, indicators, absoluteTolerance, relativeTolerance);

            int n = probabilities.Count;
            int used = 0; // Number of output rows used by this call.
            eventProbabilities.Clear();

            // Copies an indicator row into the pooled output slot, reusing an existing row
            // array of matching length in place (the steady-state zero-allocation path).
            int[] PlaceRow(int rowIndex)
            {
                int[] row;
                if (used < eventIndicators.Count && eventIndicators[used] != null && eventIndicators[used].Length == n)
                {
                    row = eventIndicators[used];
                }
                else
                {
                    row = new int[n];
                    if (used < eventIndicators.Count) eventIndicators[used] = row;
                    else eventIndicators.Add(row);
                }
                for (int column = 0; column < n; column++) row[column] = indicators[rowIndex, column];
                used++;
                return row;
            }

            // Trims stale rows from a previous, larger call.
            void TrimToUsed()
            {
                while (eventIndicators.Count > used) eventIndicators.RemoveAt(eventIndicators.Count - 1);
            }

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
                        PlaceRow(indicators.GetLength(0) - 1); // Add last indicator row
                        eventProbabilities.Add(Tools.Clamp(0.5d * diff, 0d, 1d)); // Add the average of the difference to the event probabilities
                        TrimToUsed();
                        return true; // Report that convergence ended the expansion early.
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
                var currentRow = PlaceRow(i);

                // Compute the exclusive event probability and add to the list
                eventProbabilities.Add(Tools.Clamp(IndependentExclusive(probabilities, currentRow), 0d, 1d));

                // Calculate the union of probabilities (inclusion-exclusion)
                if (i < probabilities.Count)
                {
                    union += s * probabilities[i];
                }
                else
                {
                    union += s * IndependentJointProbability(probabilities, currentRow);
                }

            }

            TrimToUsed();
            return false;
        }

        /// <summary>
        /// Validates the common inputs for lazy PCM enumeration without changing the legacy tolerance or marginal-probability handling.
        /// </summary>
        /// <param name="probabilities">The marginal event probabilities.</param>
        /// <param name="correlationMatrix">The PCM correlation matrix.</param>
        /// <param name="absoluteTolerance">The absolute convergence tolerance.</param>
        /// <param name="relativeTolerance">The relative convergence tolerance.</param>
        /// <exception cref="ArgumentException">Thrown when the probability collection is null or empty, or the correlation matrix is null.</exception>
        private static void ValidateLazyPCMInputs(IList<double> probabilities, double[,] correlationMatrix,
            double absoluteTolerance, double relativeTolerance)
        {
            if (probabilities == null || probabilities.Count == 0 || correlationMatrix == null)
                throw new ArgumentException("Input arrays must be non-empty and correlation matrix must not be null.");

            // The established dense PCM overloads do not reject negative or non-finite tolerances.
            // Preserve that validation behavior; these parameters are accepted here only to keep
            // the common lazy call contract explicit.
            _ = absoluteTolerance;
            _ = relativeTolerance;
        }
        /// <summary>
        /// Validates the structural metadata used by the pooled exclusive-probability overload.
        /// </summary>
        /// <param name="probabilities">The event probabilities.</param>
        /// <param name="binomialCombinations">The expected combination count for each subset size.</param>
        /// <param name="indicators">The materialized indicator rows.</param>
        /// <param name="absoluteTolerance">The absolute convergence tolerance.</param>
        /// <param name="relativeTolerance">The relative convergence tolerance.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="binomialCombinations"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the combination counts or indicator dimensions are inconsistent.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a probability or tolerance is outside its valid range.</exception>
        /// <exception cref="OverflowException">Thrown when a required combination count exceeds <see cref="int.MaxValue"/>.</exception>
        private static void ValidatePooledExclusiveMetadata(IList<double> probabilities, int[] binomialCombinations,
            int[,] indicators, double absoluteTolerance, double relativeTolerance)
        {
            if (binomialCombinations == null)
                throw new ArgumentNullException(nameof(binomialCombinations));
            if (binomialCombinations.Length != probabilities.Count)
                throw new ArgumentException("The binomial metadata must contain one count for each subset size.", nameof(binomialCombinations));
            if (!Tools.IsFinite(absoluteTolerance) || absoluteTolerance < 0d)
                throw new ArgumentOutOfRangeException(nameof(absoluteTolerance), "The absolute tolerance must be finite and non-negative.");
            if (!Tools.IsFinite(relativeTolerance) || relativeTolerance < 0d)
                throw new ArgumentOutOfRangeException(nameof(relativeTolerance), "The relative tolerance must be finite and non-negative.");

            long rowCount = 0;
            for (int i = 0; i < probabilities.Count; i++)
            {
                if (!Tools.IsFinite(probabilities[i]) || probabilities[i] < 0d || probabilities[i] > 1d)
                    throw new ArgumentOutOfRangeException(nameof(probabilities), "Probabilities must be finite and between zero and one.");

                int expected = checked((int)Factorial.BinomialCoefficient(probabilities.Count, i + 1));
                if (binomialCombinations[i] != expected)
                    throw new ArgumentException("The binomial metadata does not match the probability count.", nameof(binomialCombinations));
                rowCount += expected;
            }
            if (rowCount != indicators.GetLength(0))
                throw new ArgumentException("The indicator row count does not match the binomial metadata.", nameof(indicators));
        }
        /// <summary>
        /// The outcome of a lazily enumerated exclusive-probability expansion.
        /// </summary>
        public enum ExclusiveEnumerationStatus
        {
            /// <summary>Every combination was enumerated.</summary>
            Complete,

            /// <summary>The inclusion-exclusion expansion converged; the deepest combinations were not enumerated.</summary>
            Converged,

            /// <summary>The emitted-combination cap was reached; the remaining combinations were not enumerated.</summary>
            Capped,
        }

        /// <summary>
        /// The lazily enumerated form of
        /// <see cref="IndependentExclusive(IList{double}, int[], int[,], List{double}, List{int[]}, double, double)"/>:
        /// the combinations are generated in <see cref="Factorial.AllCombinations(int)"/> order
        /// rather than read from a materialized <c>n·(2^n − 1)</c> matrix, so the caller never
        /// allocates one.
        /// </summary>
        /// <param name="probabilities">An array of probabilities for each event; n = Count.</param>
        /// <param name="eventProbabilities">The caller-owned output list of exclusive event probabilities; cleared and refilled.</param>
        /// <param name="eventIndicators">The caller-owned output list of indicator rows; rows of matching length are refilled in place, and the list is trimmed to the produced count.</param>
        /// <param name="includeNoEventRow">True to emit the all-zero (no-event) combination first, carrying <c>Π(1 − pᵢ)</c>. Excluded from the inclusion-exclusion bracket.</param>
        /// <param name="maxEmittedCombinations">The cap on emitted rows, excluding any closing row; non-positive means no cap.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence of the inclusion-exclusion algorithm. Default = 1E-4.</param>
        /// <returns>Whether the expansion completed, converged early, or hit the cap.</returns>
        /// <exception cref="ArgumentNullException">Thrown when either output list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the probabilities list is null or empty.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a probability or tolerance is outside its valid range.</exception>
        /// <remarks>
        /// Emits the same rows, in the same order, with the same probabilities as the dense
        /// overload. On convergence it closes with the same half-gap pseudo-row; at the cap it
        /// closes with the exact residual <c>1 − Σ(emitted)</c>, which for independent events is
        /// the mass of everything not enumerated. Both closing rows are attributed to the all-ones
        /// combination.
        /// </remarks>
        public static ExclusiveEnumerationStatus IndependentExclusiveLazy(IList<double> probabilities,
            List<double> eventProbabilities, List<int[]> eventIndicators, bool includeNoEventRow = false,
            long maxEmittedCombinations = 0, double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (eventProbabilities == null) throw new ArgumentNullException(nameof(eventProbabilities));
            if (eventIndicators == null) throw new ArgumentNullException(nameof(eventIndicators));
            for (int i = 0; i < probabilities.Count; i++)
            {
                if (!Tools.IsFinite(probabilities[i]) || probabilities[i] < 0d || probabilities[i] > 1d)
                    throw new ArgumentOutOfRangeException(nameof(probabilities), "Probabilities must be finite and between zero and one.");
            }

            int n = probabilities.Count;
            int used = 0;
            eventProbabilities.Clear();

            int[] Row()
            {
                int[] row;
                if (used < eventIndicators.Count && eventIndicators[used] != null && eventIndicators[used].Length == n)
                {
                    row = eventIndicators[used];
                    Array.Clear(row, 0, n);
                }
                else
                {
                    row = new int[n];
                    if (used < eventIndicators.Count) eventIndicators[used] = row;
                    else eventIndicators.Add(row);
                }
                used++;
                return row;
            }

            void TrimToUsed()
            {
                while (eventIndicators.Count > used) eventIndicators.RemoveAt(eventIndicators.Count - 1);
            }

            // Closes the expansion on an all-ones row carrying the supplied mass.
            ExclusiveEnumerationStatus Close(double mass, ExclusiveEnumerationStatus status)
            {
                var row = Row();
                for (int column = 0; column < n; column++) row[column] = 1;
                eventProbabilities.Add(Tools.Clamp(mass, 0d, 1d));
                TrimToUsed();
                return status;
            }

            double noEventMass = 1d;
            for (int i = 0; i < n; i++) noEventMass *= 1d - probabilities[i];
            double totalOutputMass = includeNoEventRow ? 1d : 1d - noEventMass;
            double emittedMass = 0d;
            if (includeNoEventRow)
            {
                Row();
                eventProbabilities.Add(Tools.Clamp(noEventMass, 0d, 1d));
                emittedMass += noEventMass;
            }

            double union = 0;
            double s = 1;
            double inc = double.NaN;
            double exc = double.NaN;
            long emitted = 0;

            for (int k = 1; k <= n; k++)
            {
                // The size-block transition, mirroring the dense form: the bracket is read at the
                // first row of each size from two upward, and the sign alternates per block.
                if (k >= 2)
                {
                    int block = k - 2;
                    if (block > 0)
                    {
                        if (s == 1) inc = union;
                        else if (s == -1) exc = union;
                    }
                    double diff = Math.Abs(inc - exc);
                    if (block > 0 && block < n && diff <= absoluteTolerance && diff <= relativeTolerance * Math.Min(inc, exc))
                    {
                        return Close(0.5 * diff, ExclusiveEnumerationStatus.Converged);
                    }
                    s *= -1;
                }

                var combination = new int[k];
                for (int t = 0; t < k; t++) combination[t] = t;
                do
                {
                    if (maxEmittedCombinations > 0 && emitted >= maxEmittedCombinations)
                    {
                        return Close(Tools.Clamp(totalOutputMass - emittedMass, 0d, 1d), ExclusiveEnumerationStatus.Capped);
                    }

                    var row = Row();
                    for (int t = 0; t < k; t++) row[combination[t]] = 1;

                    double exclusive = IndependentExclusive(probabilities, row);
                    eventProbabilities.Add(Tools.Clamp(exclusive, 0d, 1d));
                    emittedMass += exclusive;
                    emitted++;

                    union += s * (k == 1 ? probabilities[combination[0]] : IndependentJointProbability(probabilities, row));
                }
                while (Factorial.NextCombinationUnchecked(combination, n));
            }

            TrimToUsed();
            return ExclusiveEnumerationStatus.Complete;
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
            return Tools.Clamp(min - max, 0d, 1d);
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
                        eventProbabilities.Add(Tools.Clamp(0.5d * diff, 0d, 1d)); // Add the average of the difference to the event probabilities
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
                eventProbabilities.Add(Tools.Clamp(PositivelyDependentExclusive(probabilities, eventIndicators.Last()), 0d, 1d));

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

        /// <summary>
        /// Lazily enumerates mutually exclusive event probabilities under perfect positive dependence.
        /// </summary>
        /// <param name="probabilities">The marginal event probabilities.</param>
        /// <param name="eventProbabilities">The caller-owned output list of exclusive event probabilities; cleared and refilled.</param>
        /// <param name="eventIndicators">The caller-owned output list of indicator rows; matching rows are reused.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence. Default = 1E-4.</param>
        /// <returns>Whether every row was enumerated or the established convergence rule closed the expansion early.</returns>
        /// <exception cref="ArgumentException">Thrown when the probability collection is null or empty.</exception>
        /// <exception cref="ArgumentNullException">Thrown when an output list is null.</exception>
        /// <remarks>
        /// Row ordering, joint-probability association, sign transitions, the dual convergence
        /// predicate, and the all-ones half-gap closing row match the dense overload exactly.
        /// </remarks>
        public static ExclusiveEnumerationStatus PositivelyDependentExclusiveLazy(IList<double> probabilities,
            List<double> eventProbabilities, List<int[]> eventIndicators,
            double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            if (probabilities == null || probabilities.Count == 0)
                throw new ArgumentException("The probabilities array must have a length greater than 0.", nameof(probabilities));
            if (eventProbabilities == null) throw new ArgumentNullException(nameof(eventProbabilities));
            if (eventIndicators == null) throw new ArgumentNullException(nameof(eventIndicators));

            int n = probabilities.Count;
            int used = 0;
            eventProbabilities.Clear();

            int[] Row()
            {
                int[] row;
                if (used < eventIndicators.Count && eventIndicators[used] != null && eventIndicators[used].Length == n)
                {
                    row = eventIndicators[used];
                    Array.Clear(row, 0, n);
                }
                else
                {
                    row = new int[n];
                    if (used < eventIndicators.Count) eventIndicators[used] = row;
                    else eventIndicators.Add(row);
                }

                used++;
                return row;
            }

            void TrimToUsed()
            {
                while (eventIndicators.Count > used) eventIndicators.RemoveAt(eventIndicators.Count - 1);
            }

            ExclusiveEnumerationStatus Close(double mass)
            {
                int[] row = Row();
                for (int i = 0; i < n; i++) row[i] = 1;
                eventProbabilities.Add(Tools.Clamp(mass, 0d, 1d));
                TrimToUsed();
                return ExclusiveEnumerationStatus.Converged;
            }

            double union = 0d;
            double sign = 1d;
            double inclusion = double.NaN;
            double exclusion = double.NaN;

            for (int subsetSize = 1; subsetSize <= n; subsetSize++)
            {
                if (subsetSize >= 2)
                {
                    int block = subsetSize - 2;
                    if (block > 0)
                    {
                        if (sign == 1d) inclusion = union;
                        else if (sign == -1d) exclusion = union;
                    }

                    double difference = Math.Abs(inclusion - exclusion);
                    if (block > 0 && block < n &&
                        difference <= absoluteTolerance &&
                        difference <= relativeTolerance * Math.Min(inclusion, exclusion))
                    {
                        return Close(0.5d * difference);
                    }

                    sign *= -1d;
                }

                var combination = new int[subsetSize];
                for (int i = 0; i < subsetSize; i++) combination[i] = i;
                do
                {
                    int[] row = Row();
                    for (int i = 0; i < combination.Length; i++) row[combination[i]] = 1;
                    eventProbabilities.Add(Tools.Clamp(PositivelyDependentExclusive(probabilities, row), 0d, 1d));
                    union += sign * (subsetSize == 1
                        ? probabilities[combination[0]]
                        : PositiveJointProbability(probabilities, row));
                }
                while (Factorial.NextCombinationUnchecked(combination, n));
            }

            TrimToUsed();
            return ExclusiveEnumerationStatus.Complete;
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
        /// Dependence between events is captured with the PCM method.
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
                result[i] = Tools.Clamp(result[i], 0d, 1d);
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
        /// <param name="absoluteTolerance">The absolute tolerance for convergence of the inclusion-exclusion algorithm. Default is 1E-4.</param>
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
                result[i] = Tools.Clamp(result[i], 0d, 1d);
            }
            eventProbabilities = result.ToList();
        }

        /// <summary>
        /// Lazily enumerates mutually exclusive event probabilities using PCM joint probabilities and inclusion-exclusion.
        /// </summary>
        /// <param name="probabilities">The marginal event probabilities.</param>
        /// <param name="correlationMatrix">The correlation matrix defining the dependency.</param>
        /// <param name="eventProbabilities">The caller-owned output list of exclusive event probabilities; cleared and refilled.</param>
        /// <param name="eventIndicators">The caller-owned output list of indicator rows; matching rows are reused.</param>
        /// <param name="absoluteTolerance">The absolute tolerance for evaluation convergence. Default = 1E-4.</param>
        /// <param name="relativeTolerance">The relative tolerance for evaluation convergence. Default = 1E-4.</param>
        /// <returns>Whether every row was enumerated or the established convergence rule closed the expansion early.</returns>
        /// <exception cref="ArgumentException">Thrown when the probability collection is null or empty, or the correlation matrix is null.</exception>
        /// <exception cref="ArgumentNullException">Thrown when an output list is null.</exception>
        /// <remarks>
        /// This method streams the same joint rows as the dense overload, then applies its
        /// exclusive inclusion-exclusion transform over the emitted rows. It does not renormalize
        /// the resulting probabilities.
        /// </remarks>
        public static ExclusiveEnumerationStatus ExclusivePCMLazy(IList<double> probabilities,
            double[,] correlationMatrix, List<double> eventProbabilities, List<int[]> eventIndicators,
            double absoluteTolerance = 1E-4, double relativeTolerance = 1E-4)
        {
            ValidateLazyPCMInputs(probabilities, correlationMatrix, absoluteTolerance, relativeTolerance);
            if (eventProbabilities == null) throw new ArgumentNullException(nameof(eventProbabilities));
            if (eventIndicators == null) throw new ArgumentNullException(nameof(eventIndicators));

            int n = probabilities.Count;
            int used = 0;
            eventProbabilities.Clear();
            var jointProbabilities = new List<double>();
            var cumulativeCombinations = new List<int>();

            int[] Row()
            {
                int[] row;
                if (used < eventIndicators.Count && eventIndicators[used] != null && eventIndicators[used].Length == n)
                {
                    row = eventIndicators[used];
                    Array.Clear(row, 0, n);
                }
                else
                {
                    row = new int[n];
                    if (used < eventIndicators.Count) eventIndicators[used] = row;
                    else eventIndicators.Add(row);
                }

                used++;
                return row;
            }

            void TrimToUsed()
            {
                while (eventIndicators.Count > used) eventIndicators.RemoveAt(eventIndicators.Count - 1);
            }

            ExclusiveEnumerationStatus status = ExclusiveEnumerationStatus.Complete;
            double union = 0d;
            double sign = 1d;
            double inclusion = double.NaN;
            double exclusion = double.NaN;

            for (int subsetSize = 1; subsetSize <= n; subsetSize++)
            {
                if (subsetSize >= 2)
                {
                    int block = subsetSize - 2;
                    if (block > 0)
                    {
                        if (sign == 1d) inclusion = union;
                        else if (sign == -1d) exclusion = union;
                    }

                    double difference = Math.Abs(inclusion - exclusion);
                    if (block > 0 && block < n &&
                        difference <= absoluteTolerance &&
                        difference <= relativeTolerance * Math.Min(inclusion, exclusion))
                    {
                        int[] closingRow = Row();
                        for (int i = 0; i < n; i++) closingRow[i] = 1;
                        jointProbabilities.Add(Tools.Clamp(0.5d * difference, 0d, 1d));
                        status = ExclusiveEnumerationStatus.Converged;
                        break;
                    }

                    sign *= -1d;
                }

                var combination = new int[subsetSize];
                for (int i = 0; i < subsetSize; i++) combination[i] = i;
                do
                {
                    int[] row = Row();
                    for (int i = 0; i < combination.Length; i++) row[combination[i]] = 1;
                    double jointProbability = subsetSize == 1
                        ? probabilities[combination[0]]
                        : JointProbability(probabilities, row, correlationMatrix);
                    jointProbabilities.Add(Tools.Clamp(jointProbability, 0d, 1d));
                    union += sign * jointProbability;
                }
                while (Factorial.NextCombinationUnchecked(combination, n));

                if (subsetSize < n) cumulativeCombinations.Add(jointProbabilities.Count);
            }

            int combinationBlock = 0;
            int nextBlockStart = cumulativeCombinations.Count == 0
                ? int.MaxValue
                : cumulativeCombinations[0];

            for (int i = 0; i < used; i++)
            {
                if (i == nextBlockStart)
                {
                    combinationBlock++;
                    nextBlockStart = combinationBlock < cumulativeCombinations.Count
                        ? cumulativeCombinations[combinationBlock]
                        : int.MaxValue;
                }

                double exclusiveProbability = jointProbabilities[i];
                double association = 1d;
                for (int block = combinationBlock; block < cumulativeCombinations.Count; block++)
                {
                    association *= -1d;
                    int start = cumulativeCombinations[block];
                    int end = block == cumulativeCombinations.Count - 1
                        ? cumulativeCombinations[block] + 1
                        : cumulativeCombinations[block + 1];
                    exclusiveProbability += association *
                        SumSearch(jointProbabilities, eventIndicators[i], eventIndicators, start, end);
                }

                eventProbabilities.Add(Tools.Clamp(exclusiveProbability, 0d, 1d));
            }

            TrimToUsed();
            return status;
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
                result[i] = Tools.Clamp(result[i], 0d, 1d);
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
                result[i] = Tools.Clamp(result[i], 0d, 1d);
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
        /// <param name="absoluteTol">The absolute tolerance for convergence evaluation. Default is 1E-4.</param>
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
                        jointProbabilities.Add(0.5d * diff);
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
                prob = Tools.Clamp(prob, 0d, 1d);
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
