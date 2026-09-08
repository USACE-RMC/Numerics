using System;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;

namespace Numerics.Distributions
{
    /// <summary>Checked full-support central-moment integration for composite distributions.</summary>
    internal static class DistributionMomentIntegration
    {
        /// <summary>Integrates a normalized log density about a local reference without raw-moment subtraction.</summary>
        /// <remarks>Maps each side of the reference through x=center±scale*t/(1-t), including infinite
        /// endpoints. Integrates the mean offset first, then directly integrates centered powers.
        /// No probability tails are truncated and no estimated mass is silently normalized.</remarks>
        internal static double[] Compute(Func<double, double> logDensity, double minimum, double maximum, double center, double scale)
        {
            if (!DistributionNumerics.IsFinite(center) || !(scale > 0) || !DistributionNumerics.IsFinite(scale))
                throw new InvalidOperationException("A finite reference and positive local scale are required for composite moment integration.");
            double lower = Limit(center - minimum, scale), upper = Limit(maximum - center, scale);
            double logScale = Math.Log(scale);
            double Moment(int order, double offset)
            {
                double Side(bool positive, double limit)
                {
                    if (limit == 0) return 0;
                    double Function(double t)
                    {
                        double magnitude = t / (1 - t);
                        double coordinate = positive ? magnitude : -magnitude;
                        double x = center + scale * coordinate;
                        double log = logDensity(x);
                        if (double.IsNegativeInfinity(log)) return 0;
                        double centered = coordinate - offset;
                        if (order != 0)
                        {
                            if (centered == 0) return 0;
                            log += order * Math.Log(Math.Abs(centered));
                        }
                        double value = Math.Exp(log + logScale - 2 * Tools.Log1p(-t));
                        if (!DistributionNumerics.IsFinite(value))
                            throw new InvalidOperationException("A composite moment is divergent or cannot be resolved numerically.");
                        return (order % 2 != 0 && centered < 0) ? -value : value;
                    }
                    var integrator = new AdaptiveGaussKronrod(Function, 0, limit)
                    {
                        RelativeTolerance = 1E-8, AbsoluteTolerance = 1E-10,
                        MaxFunctionEvaluations = 200000, ReportFailure = true
                    };
                    integrator.Integrate();
                    if (integrator.Status != IntegrationStatus.Success || !DistributionNumerics.IsFinite(integrator.Result)
                        || !DistributionNumerics.IsFinite(integrator.StandardError)
                        || integrator.StandardError > Math.Max(1E-10, Math.Abs(integrator.Result) * 1E-8))
                        throw new InvalidOperationException("Composite moment integration did not meet its error tolerance.");
                    return integrator.Result;
                }
                return Side(false, lower) + Side(true, upper);
            }
            double mass = Moment(0, 0);
            if (Math.Abs(mass - 1) > 5E-8) throw new InvalidOperationException("Composite moment integration did not recover unit probability mass.");
            double offset = Moment(1, 0);
            double mean = center + scale * offset;
            double variance = Moment(2, offset);
            if (!(variance > 0)) throw new InvalidOperationException("Composite moment integration did not produce a positive variance.");
            double third = Moment(3, offset), fourth = Moment(4, offset);
            return new[] { mean, scale * Math.Sqrt(variance), third / variance / Math.Sqrt(variance), fourth / variance / variance };
        }

        /// <summary>Maps a finite or infinite one-sided support width to [0,1].</summary>
        private static double Limit(double width, double scale)
        {
            if (width <= 0) return 0;
            if (double.IsPositiveInfinity(width)) return 1;
            double ratio = width / scale;
            return double.IsPositiveInfinity(ratio) ? 1 : ratio / (1 + ratio);
        }
    }
}
