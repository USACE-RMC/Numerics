using System;

namespace Numerics.Distributions
{
    /// <summary>Provides overflow-resistant transforms used by Hosking-form distribution tails.</summary>
    internal static partial class DistributionNumerics
    {
        /// <summary>Evaluates the Hosking shape transform without losing a finite logarithm to affine or product overflow.</summary>
        /// <param name="x">The observation in physical coordinates.</param>
        /// <param name="location">The finite location.</param>
        /// <param name="scale">The finite positive scale.</param>
        /// <param name="shape">The finite Hosking shape.</param>
        /// <returns>Minus log(1-shape*(x-location)/scale) divided by shape, with its continuous zero-shape limit.</returns>
        /// <remarks>The caller validates parameters and support. Ordinary finite products retain log1p
        /// arithmetic. If standardization or multiplication overflows, the same support expression is
        /// assembled from signed physical differences and logarithms.</remarks>
        internal static double HoskingShapeTransform(double x, double location, double scale, double shape)
        {
            double standardized = Standardize(x, location, scale);
            if (shape == 0 || double.IsNaN(standardized)) return standardized;
            double product = shape * standardized;
            if (Tools.IsFinite(product)) return product == 0 ? standardized : -Tools.Log1p(-product) / shape;
            if (!Tools.IsFinite(x) || !Tools.IsFinite(location)) return -Tools.Log1p(-product) / shape;

            double difference = x - location;
            double logDifference;
            if (Tools.IsFinite(difference)) logDifference = Math.Log(Math.Abs(difference));
            else
            {
                double magnitude = Math.Max(Math.Abs(x), Math.Abs(location));
                logDifference = Math.Log(magnitude) + Math.Log(Math.Abs(x / magnitude - location / magnitude));
            }
            double logProduct = Math.Log(Math.Abs(shape)) + logDifference - Math.Log(scale);
            bool negativeProduct = shape > 0 ? x < location : x > location;
            double logSupport = negativeProduct ? LogSum(0, logProduct) : Log1mExp(logProduct);
            return -logSupport / shape;
        }
    }
}
