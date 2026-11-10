using System;
using System.Collections.Generic;

namespace Numerics.Distributions.Regional
{
    /// <summary>
    /// The regional growth curve of the index-flood procedure: a distribution fitted to the regional average L-moment
    /// ratios with mean 1, whose quantiles scale by each site's index flood.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Hosking and Wallis (1997), section 6.2. The regional L-moments {1, t^R, t3^R, t4^R} are those of the region
    /// scaled to l1 = 1, so the fitted distribution is the growth curve q(F), whose mean is 1, and the quantile of
    /// site i at non-exceedance probability F is Q_i(F) = l1^(i) q(F), with the site mean l1^(i) as index flood.
    /// </para>
    /// <para>
    /// The distribution is fitted by its own ParametersFromLinearMoments: the generalized logistic, generalized
    /// extreme value, generalized normal, Pearson type III and generalized Pareto read (1, t^R, t3^R); the
    /// four-parameter kappa also reads t4^R.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// Hosking, J. R. M. and Wallis, J. R. (1997). Regional Frequency Analysis: An Approach Based on L-Moments.
    /// Cambridge University Press. https://doi.org/10.1017/CBO9780511529443
    /// </para>
    /// </remarks>
    public class RegionalGrowthCurve
    {
        private readonly UnivariateDistributionBase _distribution;

        /// <summary>
        /// Constructs the growth curve of a distribution fitted to regional L-moment ratios.
        /// </summary>
        /// <param name="distributionType">The distribution: <see cref="UnivariateDistributionType.GeneralizedLogistic"/>,
        /// <see cref="UnivariateDistributionType.GeneralizedExtremeValue"/>, <see cref="UnivariateDistributionType.GeneralizedNormal"/>,
        /// <see cref="UnivariateDistributionType.PearsonTypeIII"/>, <see cref="UnivariateDistributionType.GeneralizedPareto"/>
        /// or <see cref="UnivariateDistributionType.KappaFour"/>.</param>
        /// <param name="lCV">The regional L-CV t^R, positive.</param>
        /// <param name="tau3">The regional L-skewness t3^R.</param>
        /// <param name="tau4">Optional. The regional L-kurtosis t4^R; required only by the kappa distribution.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the distribution is not one of the listed
        /// families, the L-CV is not positive and finite, t3 is not finite, t4 is not finite for the kappa
        /// distribution, or the fitted parameters are invalid, including a Pearson type III skew beyond ±6.</exception>
        /// <exception cref="ArgumentException">Thrown when the family's L-moment fit rejects the ratios.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the kappa fit does not converge.</exception>
        public RegionalGrowthCurve(UnivariateDistributionType distributionType, double lCV, double tau3, double tau4 = double.NaN)
        {
            if (!IsSupported(distributionType))
                throw new ArgumentOutOfRangeException(nameof(distributionType), "The growth curve supports GLO, GEV, GNO, PE3, GPA and the four-parameter kappa.");
            if (!Tools.IsFinite(lCV) || lCV <= 0d) throw new ArgumentOutOfRangeException(nameof(lCV), "The regional L-CV must be positive and finite.");
            if (!Tools.IsFinite(tau3)) throw new ArgumentOutOfRangeException(nameof(tau3), "The regional L-skewness must be finite.");
            if (distributionType == UnivariateDistributionType.KappaFour && !Tools.IsFinite(tau4))
                throw new ArgumentOutOfRangeException(nameof(tau4), "The kappa distribution needs a finite regional L-kurtosis.");
            DistributionType = distributionType;
            _distribution = Fit(distributionType, [1d, lCV, tau3, Tools.IsFinite(tau4) ? tau4 : 0d]);
        }

        /// <summary>
        /// Constructs the growth curve of a distribution fitted to the record-length weighted regional average
        /// L-moments of a region.
        /// </summary>
        /// <param name="distributionType">The distribution, as for
        /// <see cref="RegionalGrowthCurve(UnivariateDistributionType, double, double, double)"/>.</param>
        /// <param name="data">The regional data; the ratios are <see cref="RegionalData.RegionalAverage()"/>.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown as for the ratio constructor.</exception>
        /// <exception cref="ArgumentException">Thrown when the family's L-moment fit rejects the ratios.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the kappa fit does not converge.</exception>
        public RegionalGrowthCurve(UnivariateDistributionType distributionType, RegionalData data)
            : this(distributionType, RegionalAverageOf(data))
        {
            Data = data;
        }

        /// <summary>
        /// Constructs the growth curve from regional average L-moments {1, t^R, t3^R, t4^R, ...}.
        /// </summary>
        /// <param name="distributionType">The distribution family.</param>
        /// <param name="regional">The regional average L-moments.</param>
        private RegionalGrowthCurve(UnivariateDistributionType distributionType, double[] regional)
            : this(distributionType, regional[1], regional[2], regional[3])
        {
        }

        /// <summary>
        /// Gets the fitted distribution family.
        /// </summary>
        public UnivariateDistributionType DistributionType { get; }

        /// <summary>
        /// Gets the regional data the curve was fitted to, or null when it was fitted to ratios.
        /// </summary>
        public RegionalData? Data { get; }

        /// <summary>
        /// Gets the parameters of the growth curve distribution, in the family's own order.
        /// </summary>
        public double[] Parameters => _distribution.GetParameters;

        /// <summary>
        /// Returns a copy of the growth curve distribution, whose mean is 1.
        /// </summary>
        /// <returns>A clone of the fitted distribution.</returns>
        public UnivariateDistributionBase GetDistribution() => _distribution.Clone();

        /// <summary>
        /// Returns the growth factor q(F) at a non-exceedance probability.
        /// </summary>
        /// <param name="probability">The non-exceedance probability F, in [0, 1].</param>
        /// <returns>The regional growth factor.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the probability is outside [0, 1].</exception>
        public double GrowthFactor(double probability) => _distribution.InverseCDF(probability);

        /// <summary>
        /// Returns the growth factors at several non-exceedance probabilities.
        /// </summary>
        /// <param name="probabilities">The non-exceedance probabilities, in [0, 1].</param>
        /// <returns>The growth factors, one per probability.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="probabilities"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a probability is outside [0, 1].</exception>
        public double[] GrowthFactors(IList<double> probabilities)
        {
            if (probabilities == null) throw new ArgumentNullException(nameof(probabilities));
            return _distribution.InverseCDF(probabilities);
        }

        /// <summary>
        /// Returns site quantiles: the index flood times the growth factors.
        /// </summary>
        /// <param name="indexFlood">The site's index flood, its mean l1.</param>
        /// <param name="probabilities">The non-exceedance probabilities, in [0, 1].</param>
        /// <returns>Q(F) = indexFlood q(F), one per probability.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="probabilities"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the index flood is not positive and finite, or
        /// a probability is outside [0, 1].</exception>
        public double[] SiteQuantiles(double indexFlood, IList<double> probabilities)
        {
            if (!Tools.IsFinite(indexFlood) || indexFlood <= 0d)
                throw new ArgumentOutOfRangeException(nameof(indexFlood), "The index flood must be positive and finite.");
            var q = GrowthFactors(probabilities);
            for (int j = 0; j < q.Length; j++) q[j] *= indexFlood;
            return q;
        }

        /// <summary>
        /// Returns the quantiles of every site of a region.
        /// </summary>
        /// <param name="data">The regional data; each site's index flood is its <see cref="RegionalData.Mean"/>.</param>
        /// <param name="probabilities">The non-exceedance probabilities, in [0, 1].</param>
        /// <returns>The quantiles, indexed [site, probability].</returns>
        /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a site mean is not positive, or a probability is
        /// outside [0, 1].</exception>
        public double[,] SiteQuantiles(RegionalData data, IList<double> probabilities)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            var q = GrowthFactors(probabilities);
            var result = new double[data.SiteCount, q.Length];
            for (int i = 0; i < data.SiteCount; i++)
            {
                double index = data.Mean[i];
                if (!(index > 0d)) throw new ArgumentOutOfRangeException(nameof(data), $"The mean of site {i} must be positive.");
                for (int j = 0; j < q.Length; j++) result[i, j] = index * q[j];
            }
            return result;
        }

        /// <summary>
        /// Returns whether a family can serve as a growth curve.
        /// </summary>
        /// <param name="type">The distribution family.</param>
        /// <returns>True for GLO, GEV, GNO, PE3, GPA and the four-parameter kappa.</returns>
        internal static bool IsSupported(UnivariateDistributionType type)
        {
            return type == UnivariateDistributionType.GeneralizedLogistic ||
                   type == UnivariateDistributionType.GeneralizedExtremeValue ||
                   type == UnivariateDistributionType.GeneralizedNormal ||
                   type == UnivariateDistributionType.PearsonTypeIII ||
                   type == UnivariateDistributionType.GeneralizedPareto ||
                   type == UnivariateDistributionType.KappaFour;
        }

        /// <summary>
        /// Fits a family to the L-moments {l1, l2, t3, t4} by its ParametersFromLinearMoments.
        /// </summary>
        /// <param name="type">A family for which <see cref="IsSupported"/> is true.</param>
        /// <param name="moments">The L-moments {l1, l2, t3, t4}.</param>
        /// <returns>The fitted distribution.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the fitted parameters are invalid.</exception>
        internal static UnivariateDistributionBase Fit(UnivariateDistributionType type, double[] moments)
        {
            var distribution = UnivariateDistributionFactory.CreateDistribution(type);
            var parameters = ((ILinearMomentEstimation)distribution).ParametersFromLinearMoments(moments);
            distribution.ValidateParameters(parameters, true);
            distribution.SetParameters(parameters);
            return distribution;
        }

        /// <summary>
        /// Returns the regional average of the data, refusing null.
        /// </summary>
        /// <param name="data">The regional data.</param>
        /// <returns>{1, t^R, t3^R, t4^R, t5^R}.</returns>
        private static double[] RegionalAverageOf(RegionalData data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            return data.RegionalAverage();
        }
    }
}
