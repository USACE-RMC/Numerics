using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;

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
            lambda = double.NaN;
            if (!CanFitLambda(values))
                return;

            // Keep BrentSearch arithmetic finite even when the profile likelihood is undefined.
            var brent = new BrentSearch((x) =>
            {
                double logLikelihood = LogLikelihood(values, x);
                return Tools.IsFinite(logLikelihood) ? logLikelihood : -double.MaxValue;
            }, -5d, 5d)
            {
                ReportFailure = false,
                ComputeHessian = false,
                RecordTraces = false
            };

            brent.Maximize();
            if (brent.Status != OptimizationStatus.Success ||
                brent.BestParameterSet.Values == null ||
                brent.BestParameterSet.Values.Length == 0)
                return;

            double candidate = brent.BestParameterSet.Values[0];
            if (!Tools.IsFinite(candidate) ||
                Math.Abs(candidate) > 5d ||
                !Tools.IsFinite(LogLikelihood(values, candidate)))
                return;

            lambda = candidate;
        }

        /// <summary>
        /// Fit the transform parameter using maximum likelihood estimation.
        /// </summary>
        /// <param name="values">The list of values to transform.</param>
        /// <returns>The fitted transformation exponent.</returns>
        public static double FitLambda(IList<double> values)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Count < 2) throw new ArgumentException("At least 2 values are required to fit lambda.", nameof(values));

            FitLambda(values, out double lambda);
            return lambda;
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
            if (!CanFitLambda(values) || !Tools.IsFinite(lambda) || Math.Abs(lambda) > 5d)
                return double.NegativeInfinity;

            int n = values.Count;
            var transformed = new double[n];
            double sum = 0d;
            double logJacobianSum = 0d;

            // Transform values and compute sum and log-Jacobian
            for (int i = 0; i < n; i++)
            {
                double xi = values[i];
                double yi = Transform(xi, lambda);
                if (!Tools.IsFinite(yi))
                    return double.NegativeInfinity;

                transformed[i] = yi;
                sum += yi;
                if (!Tools.IsFinite(sum))
                    return double.NegativeInfinity;

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
                if (!Tools.IsFinite(logJacobianSum))
                    return double.NegativeInfinity;
            }

            // Compute mean and SSE
            double mu = sum / n;
            double sse = 0d;
            for (int i = 0; i < n; i++)
            {
                double resid = transformed[i] - mu;
                sse += resid * resid;
            }
            if (!Tools.IsFinite(sse) || sse <= 0d)
                return double.NegativeInfinity;

            double sigmaSq = sse / n;
            if (!Tools.IsFinite(sigmaSq) || sigmaSq <= 0)
                return double.NegativeInfinity;

            double logLikelihood = -n / 2.0 * Tools.LogSqrt2PI - n / 2.0 * Math.Log(sigmaSq) - sse / (2.0 * sigmaSq) + logJacobianSum;

            return Tools.IsFinite(logLikelihood) ? logLikelihood : double.NegativeInfinity;
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
                logJacobianSum += Math.Log(Math.Abs(dTdy));
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
            else if (value < 0 && Math.Abs(lambda - 2d) >= 1E-8)
            {
                return -(Math.Pow(-value + 1, 2 - lambda) - 1) / (2 - lambda);
            }
            else if (value < 0 && Math.Abs(lambda - 2d) < 1E-8)
            {
                return -Math.Log(-value + 1);
            }
            else
            {
                return double.NaN;
            }
        }

        /// <summary>
        /// Returns the derivative of the Yeo-Johnson transformation with respect to the original value.
        /// </summary>
        /// <param name="value">The value at which to evaluate the derivative.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        public static double Derivative(double value, double lambda)
        {
            if (Math.Abs(lambda) > 5d) return double.NaN;
            return value >= 0d
                ? Math.Pow(value + 1d, lambda - 1d)
                : Math.Pow(-value + 1d, 1d - lambda);
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
            else if (value < 0 && Math.Abs(lambda - 2d) >= 1E-8)
            {
                return 1 - (Math.Pow(1 - (2 - lambda) * value, 1 / (2 - lambda)));
            }
            else if (value < 0 && Math.Abs(lambda - 2d) < 1E-8)
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

        /// <summary>
        /// Determines whether a sample can support Yeo-Johnson lambda fitting.
        /// </summary>
        /// <param name="values">The sample values to inspect.</param>
        /// <returns><see langword="true"/> when the sample is finite and non-degenerate.</returns>
        private static bool CanFitLambda(IList<double> values)
        {
            if (values == null || values.Count < 2)
                return false;

            double first = values[0];
            if (!Tools.IsFinite(first))
                return false;

            bool hasDifferentValue = false;
            for (int i = 1; i < values.Count; i++)
            {
                double value = values[i];
                if (!Tools.IsFinite(value))
                    return false;
                if (value != first)
                    hasDifferentValue = true;
            }

            return hasDifferentValue;
        }
    }
}
