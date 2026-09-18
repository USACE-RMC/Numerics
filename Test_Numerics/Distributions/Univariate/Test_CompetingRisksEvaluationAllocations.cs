using Numerics.Data.Statistics;
using Numerics.Distributions;

namespace Distributions
{
    /// <summary>Protects live evaluation from temporary collection-wrapper allocations.</summary>
    [TestClass]
    public class Test_CompetingRisksEvaluationAllocations
    {
#if NET6_0_OR_GREATER
        /// <summary>Known out-of-support densities and CDFs allocate only the existing bounds enumeration.</summary>
        [TestMethod]
        public void OutsideSupport_AllocatesOnlyBoundsEnumeration()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Weibull(5d, 1.4d), new Weibull(12d, 2.3d)
            })
            {
                Dependency = Probability.DependencyType.CorrelationMatrix,
                CorrelationMatrix = new[,] { { 1d, .35d }, { .35d, 1d } }
            };
            for (int i = 0; i < 100; i++) { _ = distribution.LogPDF(-1d); _ = distribution.CDF(-1d); }
            long before = GC.GetAllocatedBytesForCurrentThread();
            for (int i = 0; i < 200; i++) _ = distribution.Minimum;
            long boundsAllocation = GC.GetAllocatedBytesForCurrentThread() - before;
            before = GC.GetAllocatedBytesForCurrentThread();
            for (int i = 0; i < 100; i++) { _ = distribution.LogPDF(-1d); _ = distribution.CDF(-1d); }
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
            Assert.AreEqual(boundsAllocation, allocated, "Live validation and collection access should add no allocations to the bounds enumeration.");
            Assert.AreEqual(double.NegativeInfinity, distribution.LogPDF(-1d));
            Assert.AreEqual(0d, distribution.CDF(-1d));

            // Allocation removal must retain validation before the support shortcut.
            ((Weibull)distribution.Distributions[1]).Kappa = double.NaN;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.LogPDF(-1d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.CDF(-1d));
        }
#endif

        /// <summary>Derived owners retain every live support read even when their components are built-in Weibulls.</summary>
        [TestMethod]
        public void DerivedOwner_RetainsLiveSupportReads()
        {
            var distribution = new ObservedSupport(new UnivariateDistributionBase[]
            {
                new Weibull(5d, 1.4d), new Weibull(12d, 2.3d)
            })
            {
                Dependency = Probability.DependencyType.CorrelationMatrix,
                CorrelationMatrix = new[,] { { 1d, .35d }, { .35d, 1d } }
            };
            _ = distribution.LogPDF(3d);
            _ = distribution.LogPDF(3d);
            distribution.MinimumReads = distribution.MaximumReads = 0;
            distribution.CdfCalls = 0;
            _ = distribution.LogPDF(3d);
            Assert.AreEqual(5, distribution.MinimumReads);
            Assert.AreEqual(5, distribution.MaximumReads);
            Assert.AreEqual(2, distribution.CdfCalls);
        }

        /// <summary>Observes the virtual support extension points of a derived composite.</summary>
        private sealed class ObservedSupport : CompetingRisks
        {
            internal int MinimumReads;
            internal int MaximumReads;
            internal int CdfCalls;
            internal ObservedSupport(UnivariateDistributionBase[] distributions) : base(distributions) { }

            /// <inheritdoc/>
            public override double Minimum { get { MinimumReads++; return base.Minimum; } }

            /// <inheritdoc/>
            public override double Maximum { get { MaximumReads++; return base.Maximum; } }

            /// <inheritdoc/>
            public override double CDF(double x) { CdfCalls++; return base.CDF(x); }
        }
    }
}
