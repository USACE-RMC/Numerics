using System;
using System.Collections.Generic;
using Numerics.Mathematics.Optimization;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Class for performing Box-Cox transformation.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// <para>
    /// <b> Description: </b>
    /// This method transforms non-normal dependent variables into a normal shape.
    /// </para>
    /// </remarks>
    public class BoxCox
    {

        /// <summary>
        /// Fit the transformation parameters using maximum likelihood estimation.
        /// </summary>
        /// <param name="values">The list of values to transform.</param>
        /// <param name="lambda">Output. The transformation exponent. Range -5 to +5.</param>
        /// <remarks>
        /// https://www.rdocumentation.org/packages/EnvStats/versions/2.4.0/topics/boxcox
        /// </remarks>
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
        /// The log-likelihood function. The transformed observations are assumed to come from a
        /// normal distribution. The change of variable formula is used to write the log-likelihood function.
        /// </summary>
        /// <param name="values">The list of values to transform.</param>
        /// <param name="lambda1">The transformation exponent. Range -5 to +5.</param>
        /// <returns>
        /// The value of log-likelihood function evaluated at the given values and lambdas.
        /// </returns>
        public static double LogLikelihood(IList<double> values, double lambda1)
        {
            if (!CanFitLambda(values) || !Tools.IsFinite(lambda1) || Math.Abs(lambda1) > 5d)
                return double.NegativeInfinity;

            int n = values.Count;
            var y = new double[n];
            double mu = 0d;
            var sumX = 0d;
            for (int i = 0; i < n; i++)
            {
                y[i] = Transform(values[i], lambda1);
                if (!Tools.IsFinite(y[i]))
                    return double.NegativeInfinity;

                mu += y[i];
                sumX += Math.Log(values[i]);
            }
            if (!Tools.IsFinite(mu) || !Tools.IsFinite(sumX))
                return double.NegativeInfinity;

            mu = mu / n;
            double sse = 0d;
            for (int i = 0; i < n; i++)
                sse += Math.Pow(y[i] - mu, 2d);
            if (!Tools.IsFinite(sse) || sse <= 0d)
                return double.NegativeInfinity;

            double sigma = Math.Sqrt(sse / n);
            if (!Tools.IsFinite(sigma) || sigma <= 0d)
                return double.NegativeInfinity;

            double ll = -n / 2.0d * Tools.LogSqrt2PI - n / 2.0d * Math.Log(sigma * sigma) - 1.0d / (2d * sigma * sigma) * sse + (lambda1 - 1d) * sumX;
            return Tools.IsFinite(ll) ? ll : double.NegativeInfinity;
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

            for (int i = 0; i < n; i++)
            {
                double xi = values[i];

                if (xi <= 0)
                    return double.NegativeInfinity; // Box-Cox undefined for non-positive values

                // For Box-Cox: log derivative is (lambda - 1) * log(x)
                logJacobianSum += (lambda - 1d) * Math.Log(xi);
            }

            return logJacobianSum;
        }

        /// <summary>
        /// Returns the Box-Cox transformation of the value.
        /// </summary>
        /// <param name="value">The value to transform.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        public static double Transform(double value, double lambda)
        {
            if (value <= 0) return double.NaN;
            if (Math.Abs(lambda) > 5d) return double.NaN;
            if (Math.Abs(lambda) < 1e-8) return Math.Log(value);
            return (Math.Pow(value, lambda) - 1.0d) / lambda;
        }

        /// <summary>
        /// Returns the Box-Cox transformation of each value in the list.
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
        /// Returns the reverse of the Box-Cox transformed value.
        /// </summary>
        /// <param name="value">The value to reverse transform.</param>
        /// <param name="lambda">The transformation exponent. Range -5 to +5.</param>
        public static double InverseTransform(double value, double lambda)
        {
            if (Math.Abs(lambda) > 5d) return double.NaN;
            if (Math.Abs(lambda) < 1e-8) return Math.Exp(value);
            return Math.Pow(value * lambda + 1.0d, 1.0d / lambda);
        }

        /// <summary>
        /// Returns the inverse of each Box-Cox transformed value in the list.
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
        /// Determines whether a sample can support Box-Cox lambda fitting.
        /// </summary>
        /// <param name="values">The sample values to inspect.</param>
        /// <returns><see langword="true"/> when the sample is finite, positive, and non-degenerate.</returns>
        private static bool CanFitLambda(IList<double> values)
        {
            if (values == null || values.Count < 2)
                return false;

            double first = values[0];
            if (!Tools.IsFinite(first) || first <= 0d)
                return false;

            bool hasDifferentValue = false;
            for (int i = 1; i < values.Count; i++)
            {
                double value = values[i];
                if (!Tools.IsFinite(value) || value <= 0d)
                    return false;
                if (value != first)
                    hasDifferentValue = true;
            }

            return hasDifferentValue;
        }
    }
}