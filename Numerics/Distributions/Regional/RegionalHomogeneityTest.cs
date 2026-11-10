using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Threading.Tasks;
using Numerics.Data.Statistics;
using Numerics.Sampling;

namespace Numerics.Distributions.Regional
{
    /// <summary>
    /// The heterogeneity measures H1, H2 and H3 and the goodness-of-fit measures Z^DIST of a region, from Monte Carlo
    /// simulation of a homogeneous region with a four-parameter kappa parent.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The statistics of Hosking's REGTST routine (LMOMENTS Fortran package) and of lmomRFA's regtst, defined in
    /// Hosking and Wallis (1993) and Hosking and Wallis (1997), sections 3.2.3, 4.3.3 and 5.2.3.
    /// </para>
    /// <para>
    /// A four-parameter kappa distribution is fitted to the regional average L-moments {1, t^R, t3^R, t4^R}
    /// (Hosking, 1994). When that fit fails, which happens when (t3^R, t4^R) lies above the generalized logistic line,
    /// the parent is the generalized logistic fitted to {1, t^R, t3^R}, the kappa with h = −1, as REGTST does. Each
    /// simulated region draws, independently at every site, as many values from the parent as the site's record
    /// length, and computes the sites' sample L-moment ratios.
    /// </para>
    /// <para>
    /// The dispersion statistics are V1 = (Σ n_i (t^(i) − t^R)² / Σ n_i)^½,
    /// V2 = Σ n_i ((t^(i) − t^R)² + (t3^(i) − t3^R)²)^½ / Σ n_i and
    /// V3 = Σ n_i ((t3^(i) − t3^R)² + (t4^(i) − t4^R)²)^½ / Σ n_i, each simulated region measured about its own
    /// regional averages. H_j = (V_j,obs − μ_V) / σ_V with μ_V and σ_V the mean and sample standard deviation of the
    /// simulated V_j. Hosking and Wallis call a region acceptably homogeneous for H &lt; 1, possibly heterogeneous for
    /// 1 ≤ H &lt; 2 and definitely heterogeneous for H ≥ 2; H1 is the recommended measure.
    /// </para>
    /// <para>
    /// For each candidate distribution fitted to {1, t^R, t3^R}, Z^DIST = (τ4^DIST − t4^R + B4) / σ4, with τ4^DIST the
    /// fitted distribution's L-kurtosis, B4 the mean simulated regional t4^R less the observed t4^R, and σ4 the sample
    /// standard deviation of the simulated regional t4^R. A fit is acceptable at the 10% level when |Z| ≤ 1.64.
    /// The Pearson type III candidate is limited by <see cref="PearsonTypeIII"/> to |skew| ≤ 6 (t3 up to about 0.772);
    /// beyond that its τ4 and Z are NaN, where Hosking's REGTST still reports them.
    /// </para>
    /// <para>
    /// Simulation m uses its own <see cref="MersenneTwister"/>, seeded by the m-th of the integers drawn from a
    /// MersenneTwister seeded with <see cref="Seed"/>, so the results depend only on the seed and the number of
    /// simulations, not on the number of threads. The uniform variates are on the open interval (0, 1). The results
    /// match Hosking's within Monte Carlo error, not exactly, because the random number generators differ.
    /// </para>
    /// <para>
    /// The sites are simulated as independent. Intersite correlation inflates the variance of the regional averages
    /// and is not represented.
    /// </para>
    /// <para>
    /// A simulated site whose sample mean is zero or negative, possible when the parent has a large L-CV and a lower
    /// bound below zero, has an infinite or negative L-CV, which corrupts that region's V statistics and therefore H.
    /// REGTST behaves the same way, and the algorithm keeps that behavior; <see cref="NonpositiveMeanSimulations"/>
    /// counts the simulated regions affected so that a caller can detect it. When the simulated V_j or the simulated
    /// regional t4 have no spread (σ_V or σ4 equal to 0), H_j or Z is infinite or NaN; no error is raised.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Hosking, J. R. M. (1994). The four-parameter kappa distribution. IBM Journal of Research and Development,
    /// 38(3), 251-258. https://doi.org/10.1147/rd.383.0251
    /// </description></item>
    /// <item><description>
    /// Hosking, J. R. M. and Wallis, J. R. (1993). Some statistics useful in regional frequency analysis. Water
    /// Resources Research, 29(2), 271-281. https://doi.org/10.1029/92WR01980
    /// </description></item>
    /// <item><description>
    /// Hosking, J. R. M. and Wallis, J. R. (1997). Regional Frequency Analysis: An Approach Based on L-Moments.
    /// Cambridge University Press. https://doi.org/10.1017/CBO9780511529443
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class RegionalHomogeneityTest
    {
        /// <summary>
        /// The candidate distributions of the goodness-of-fit measure, in the order of REGTST.
        /// </summary>
        private static readonly UnivariateDistributionType[] Candidates =
        [
            UnivariateDistributionType.GeneralizedLogistic,
            UnivariateDistributionType.GeneralizedExtremeValue,
            UnivariateDistributionType.GeneralizedNormal,
            UnivariateDistributionType.PearsonTypeIII,
            UnivariateDistributionType.GeneralizedPareto,
        ];

        private readonly double[] _regionalLinearMoments, _discordancy, _kappaParameters, _observedDispersion;
        private readonly double[] _simulatedDispersionMean, _simulatedDispersionStandardDeviation, _heterogeneity;
        private readonly double[] _distributionTau4, _z;
        private readonly UnivariateDistributionType[] _distributions;

        /// <summary>
        /// Runs the simulations and computes the test statistics of a region.
        /// </summary>
        /// <param name="data">The regional data, with at least two sites, each with a record length of at least 4.</param>
        /// <param name="simulations">Optional. The number of simulated regions, at least 2. Default = 500, the number
        /// Hosking and Wallis use.</param>
        /// <param name="seed">Optional. The master seed of the simulations. Default = 12345.</param>
        /// <param name="maxDegreeOfParallelism">Optional. The largest number of simulations run at once; −1 (the
        /// default) leaves the choice to the scheduler. The results do not depend on it.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the region has fewer than two sites, or a site has a record
        /// length below 4.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="simulations"/> is below 2, or
        /// <paramref name="maxDegreeOfParallelism"/> is 0 or below −1.</exception>
        /// <exception cref="InvalidOperationException">Thrown when neither the kappa nor the generalized logistic can
        /// be fitted to the regional L-moments.</exception>
        public RegionalHomogeneityTest(RegionalData data, int simulations = 500, int seed = 12345, int maxDegreeOfParallelism = -1)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.SiteCount < 2) throw new ArgumentException("The heterogeneity measures need at least two sites.", nameof(data));
            if (simulations < 2) throw new ArgumentOutOfRangeException(nameof(simulations), "The number of simulations must be at least 2.");
            if (maxDegreeOfParallelism == 0 || maxDegreeOfParallelism < -1)
                throw new ArgumentOutOfRangeException(nameof(maxDegreeOfParallelism), "The degree of parallelism must be positive or -1.");
            int sites = data.SiteCount;
            var n = new int[sites];
            double sumLength = 0d;
            for (int i = 0; i < sites; i++)
            {
                n[i] = data.RecordLengths[i];
                if (n[i] < 4) throw new ArgumentException($"Site {i} has a record length below 4.", nameof(data));
                sumLength += n[i];
            }

            Data = data;
            Simulations = simulations;
            Seed = seed;
            _regionalLinearMoments = data.RegionalAverage();
            _discordancy = data.Discordancy();
            DiscordancyCriticalValue = RegionalData.DiscordancyCriticalValue(sites);
            var rmom = _regionalLinearMoments;
            double t4R = rmom[3];

            // The simulated parent.
            _kappaParameters = FitParent(rmom, out bool fallback);
            IsGeneralizedLogisticParent = fallback;
            var parent = new KappaFour(_kappaParameters[0], _kappaParameters[1], _kappaParameters[2], _kappaParameters[3]);

            // The observed dispersion about the regional averages.
            var lcv = new double[sites];
            var t3 = new double[sites];
            var t4 = new double[sites];
            for (int i = 0; i < sites; i++)
            {
                lcv[i] = data.LCV[i];
                t3[i] = data.Tau3[i];
                t4[i] = data.Tau4[i];
            }
            _observedDispersion = Dispersion(n, sumLength, lcv, t3, t4, rmom[1], rmom[2], rmom[3]);

            // The simulated regions.
            var seeds = new MersenneTwister(seed).NextIntegers(1, int.MaxValue, simulations);
            var v = new double[simulations, 4];
            var nonpositive = new bool[simulations];
            var options = new ParallelOptions { MaxDegreeOfParallelism = maxDegreeOfParallelism };
            Parallel.For(0, simulations, options, m =>
            {
                var random = new MersenneTwister(seeds[m]);
                var simCV = new double[sites];
                var sim3 = new double[sites];
                var sim4 = new double[sites];
                double sum2 = 0d, sum3 = 0d, sum4 = 0d;
                for (int i = 0; i < sites; i++)
                {
                    var x = new double[n[i]];
                    for (int j = 0; j < x.Length; j++) x[j] = parent.InverseCDF(random.GenRandReal3());
                    var moments = Statistics.LinearMoments(x, 4);
                    if (!(moments[0] > 0d)) nonpositive[m] = true;
                    simCV[i] = moments[1] / moments[0];
                    sim3[i] = moments[2];
                    sim4[i] = moments[3];
                    sum2 += n[i] * simCV[i];
                    sum3 += n[i] * sim3[i];
                    sum4 += n[i] * sim4[i];
                }
                sum2 /= sumLength;
                sum3 /= sumLength;
                sum4 /= sumLength;
                var dispersion = Dispersion(n, sumLength, simCV, sim3, sim4, sum2, sum3, sum4);
                v[m, 0] = dispersion[0];
                v[m, 1] = dispersion[1];
                v[m, 2] = dispersion[2];
                v[m, 3] = sum4;
            });

            int affected = 0;
            for (int m = 0; m < simulations; m++)
                if (nonpositive[m]) affected++;
            NonpositiveMeanSimulations = affected;

            // Means and sample standard deviations over the simulations, in simulation order.
            var mean = new double[4];
            var sd = new double[4];
            for (int c = 0; c < 4; c++)
            {
                double s = 0d;
                for (int m = 0; m < simulations; m++) s += v[m, c];
                mean[c] = s / simulations;
                double ss = 0d;
                for (int m = 0; m < simulations; m++)
                {
                    double dev = v[m, c] - mean[c];
                    ss += dev * dev;
                }
                sd[c] = Math.Sqrt(ss / (simulations - 1));
            }
            _simulatedDispersionMean = [mean[0], mean[1], mean[2]];
            _simulatedDispersionStandardDeviation = [sd[0], sd[1], sd[2]];
            _heterogeneity = new double[3];
            for (int j = 0; j < 3; j++)
                _heterogeneity[j] = (_observedDispersion[j] - mean[j]) / sd[j];

            // Goodness of fit.
            Tau4Bias = mean[3] - t4R;
            Tau4StandardDeviation = sd[3];
            _distributions = (UnivariateDistributionType[])Candidates.Clone();
            _distributionTau4 = new double[Candidates.Length];
            _z = new double[Candidates.Length];
            for (int d = 0; d < Candidates.Length; d++)
            {
                _distributionTau4[d] = CandidateTau4(Candidates[d], rmom);
                _z[d] = (_distributionTau4[d] - t4R + Tau4Bias) / Tau4StandardDeviation;
            }
        }

        /// <summary>
        /// Gets the regional data tested.
        /// </summary>
        public RegionalData Data { get; }

        /// <summary>
        /// Gets the number of simulated regions.
        /// </summary>
        public int Simulations { get; }

        /// <summary>
        /// Gets the master seed of the simulations.
        /// </summary>
        public int Seed { get; }

        /// <summary>
        /// Gets the record-length weighted regional average L-moments {1, t^R, t3^R, t4^R, t5^R}.
        /// </summary>
        public ReadOnlyCollection<double> RegionalLinearMoments => Array.AsReadOnly(_regionalLinearMoments);

        /// <summary>
        /// Gets the discordancy measure D_i of each site (<see cref="RegionalData.Discordancy()"/>).
        /// </summary>
        public ReadOnlyCollection<double> Discordancy => Array.AsReadOnly(_discordancy);

        /// <summary>
        /// Gets the 10% critical value of the discordancy measure for the region's number of sites.
        /// </summary>
        public double DiscordancyCriticalValue { get; }

        /// <summary>
        /// Gets the parameters {ξ, α, k, h} of the simulated kappa parent, in Hosking's parameterization; h = −1 when
        /// the generalized logistic stands in.
        /// </summary>
        public ReadOnlyCollection<double> KappaParameters => Array.AsReadOnly(_kappaParameters);

        /// <summary>
        /// Gets whether the kappa fit failed and the generalized logistic (the kappa with h = −1) is the parent.
        /// </summary>
        public bool IsGeneralizedLogisticParent { get; }

        /// <summary>
        /// Gets the observed dispersion statistics {V1, V2, V3}.
        /// </summary>
        public ReadOnlyCollection<double> ObservedDispersion => Array.AsReadOnly(_observedDispersion);

        /// <summary>
        /// Gets the mean of each dispersion statistic over the simulated regions.
        /// </summary>
        public ReadOnlyCollection<double> SimulatedDispersionMean => Array.AsReadOnly(_simulatedDispersionMean);

        /// <summary>
        /// Gets the sample standard deviation of each dispersion statistic over the simulated regions.
        /// </summary>
        public ReadOnlyCollection<double> SimulatedDispersionStandardDeviation => Array.AsReadOnly(_simulatedDispersionStandardDeviation);

        /// <summary>
        /// Gets the heterogeneity measures {H1, H2, H3}.
        /// </summary>
        public ReadOnlyCollection<double> Heterogeneity => Array.AsReadOnly(_heterogeneity);

        /// <summary>
        /// Gets the number of simulated regions in which at least one site had a zero or negative sample mean, and so
        /// an infinite or negative simulated L-CV. A nonzero count means the heterogeneity measures are unreliable.
        /// </summary>
        public int NonpositiveMeanSimulations { get; }

        /// <summary>
        /// Gets the bias B4 of the regional L-kurtosis: the mean simulated t4^R less the observed t4^R.
        /// </summary>
        public double Tau4Bias { get; }

        /// <summary>
        /// Gets σ4, the sample standard deviation of the simulated regional L-kurtosis.
        /// </summary>
        public double Tau4StandardDeviation { get; }

        /// <summary>
        /// Gets the candidate distributions of the goodness-of-fit measure: GLO, GEV, GNO, PE3 and GPA.
        /// </summary>
        public ReadOnlyCollection<UnivariateDistributionType> Distributions => Array.AsReadOnly(_distributions);

        /// <summary>
        /// Gets the L-kurtosis τ4 of each candidate fitted to {1, t^R, t3^R}; NaN where the fit fails.
        /// </summary>
        public ReadOnlyCollection<double> DistributionTau4 => Array.AsReadOnly(_distributionTau4);

        /// <summary>
        /// Gets the goodness-of-fit measure Z^DIST of each candidate, in the order of <see cref="Distributions"/>;
        /// NaN where the fit fails.
        /// </summary>
        public ReadOnlyCollection<double> Z => Array.AsReadOnly(_z);

        /// <summary>
        /// Returns the goodness-of-fit measure of one candidate distribution.
        /// </summary>
        /// <param name="distributionType">One of <see cref="Distributions"/>.</param>
        /// <returns>Z^DIST.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the distribution is not a candidate.</exception>
        public double GetZ(UnivariateDistributionType distributionType)
        {
            int index = Array.IndexOf(_distributions, distributionType);
            if (index < 0) throw new ArgumentOutOfRangeException(nameof(distributionType), "The goodness-of-fit measure covers GLO, GEV, GNO, PE3 and GPA.");
            return _z[index];
        }

        /// <summary>
        /// Returns the candidates whose goodness-of-fit measure is acceptable, ordered by |Z|.
        /// </summary>
        /// <param name="criticalValue">Optional. The largest acceptable |Z|. Default = 1.64, the 10% level.</param>
        /// <returns>The acceptable candidates, best fit first.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the critical value is not positive and finite.</exception>
        public UnivariateDistributionType[] AcceptableDistributions(double criticalValue = 1.64)
        {
            if (!Tools.IsFinite(criticalValue) || criticalValue <= 0d)
                throw new ArgumentOutOfRangeException(nameof(criticalValue), "The critical value must be positive and finite.");
            var accepted = new List<int>();
            for (int d = 0; d < _z.Length; d++)
                if (Math.Abs(_z[d]) <= criticalValue) accepted.Add(d);
            accepted.Sort((a, b) => Math.Abs(_z[a]).CompareTo(Math.Abs(_z[b])));
            var result = new UnivariateDistributionType[accepted.Count];
            for (int k = 0; k < result.Length; k++) result[k] = _distributions[accepted[k]];
            return result;
        }

        /// <summary>
        /// Fits the simulation parent: the kappa when its L-moment fit succeeds, else the generalized logistic as the
        /// kappa with h = −1.
        /// </summary>
        /// <param name="rmom">The regional average L-moments {1, t^R, t3^R, t4^R, ...}.</param>
        /// <param name="fallback">True when the generalized logistic stands in.</param>
        /// <returns>The kappa parameters {ξ, α, k, h}.</returns>
        /// <exception cref="InvalidOperationException">Thrown when neither fit succeeds.</exception>
        private static double[] FitParent(double[] rmom, out bool fallback)
        {
            double[] moments = [rmom[0], rmom[1], rmom[2], rmom[3]];
            try
            {
                var kappa = new KappaFour().ParametersFromLinearMoments(moments);
                fallback = false;
                return kappa;
            }
            catch (Exception ex) when (ex is ArgumentException || ex is InvalidOperationException)
            {
                // REGTST replaces any kappa failure by the generalized logistic, the kappa with h = -1.
                fallback = true;
                try
                {
                    var glo = new GeneralizedLogistic().ParametersFromLinearMoments(moments);
                    return [glo[0], glo[1], glo[2], -1d];
                }
                catch (ArgumentException inner)
                {
                    throw new InvalidOperationException("Neither the kappa nor the generalized logistic distribution fits the regional L-moments.",
                        new AggregateException(ex, inner));
                }
            }
        }

        /// <summary>
        /// Returns the L-kurtosis of a candidate fitted to {1, t^R, t3^R}, or NaN when the fit fails.
        /// </summary>
        /// <param name="type">The candidate family.</param>
        /// <param name="rmom">The regional average L-moments.</param>
        /// <returns>τ4 of the fitted candidate.</returns>
        private static double CandidateTau4(UnivariateDistributionType type, double[] rmom)
        {
            try
            {
                var fitted = RegionalGrowthCurve.Fit(type, [rmom[0], rmom[1], rmom[2], rmom[3]]);
                return ((ILinearMomentEstimation)fitted).LinearMomentsFromParameters(fitted.GetParameters)[3];
            }
            catch (Exception ex) when (ex is ArgumentException || ex is InvalidOperationException)
            {
                // A candidate that cannot match the regional L-skewness has no goodness-of-fit measure.
                return double.NaN;
            }
        }

        /// <summary>
        /// Returns the dispersion statistics {V1, V2, V3} of sites about given regional averages.
        /// </summary>
        /// <param name="n">The record lengths.</param>
        /// <param name="sumLength">The sum of the record lengths.</param>
        /// <param name="lcv">The sites' L-CV.</param>
        /// <param name="t3">The sites' L-skewness.</param>
        /// <param name="t4">The sites' L-kurtosis.</param>
        /// <param name="rcv">The regional L-CV.</param>
        /// <param name="r3">The regional L-skewness.</param>
        /// <param name="r4">The regional L-kurtosis.</param>
        /// <returns>{V1, V2, V3}.</returns>
        private static double[] Dispersion(int[] n, double sumLength, double[] lcv, double[] t3, double[] t4, double rcv, double r3, double r4)
        {
            double v1 = 0d, v2 = 0d, v3 = 0d;
            for (int i = 0; i < n.Length; i++)
            {
                double d2 = (lcv[i] - rcv) * (lcv[i] - rcv);
                double d3 = (t3[i] - r3) * (t3[i] - r3);
                double d4 = (t4[i] - r4) * (t4[i] - r4);
                v1 += n[i] * d2;
                v2 += n[i] * Math.Sqrt(d2 + d3);
                v3 += n[i] * Math.Sqrt(d3 + d4);
            }
            return [Math.Sqrt(v1 / sumLength), v2 / sumLength, v3 / sumLength];
        }
    }
}
