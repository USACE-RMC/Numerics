using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Functions;
using System;

namespace Functions
{
    /// <summary>
    /// Unit tests for the <see cref="FisherZLink"/> class.
    /// </summary>
    [TestClass]
    public class Test_FisherZLink
    {
        private const double RoundTripTol = 1e-12;
        private const double DeltaH = 1e-7;
        private const double DerivativeTol = 1e-6;

        /// <summary>
        /// Verifies known link values.
        /// </summary>
        [TestMethod]
        public void Link_KnownValues()
        {
            var link = new FisherZLink();

            Assert.AreEqual(0d, link.Link(0d), 1e-12);
            Assert.AreEqual(0.5d * Math.Log(3d), link.Link(0.5d), 1e-12);
            Assert.AreEqual(-0.5d * Math.Log(3d), link.Link(-0.5d), 1e-12);
        }

        /// <summary>
        /// Verifies inverse-link values.
        /// </summary>
        [TestMethod]
        public void InverseLink_KnownValues()
        {
            var link = new FisherZLink();

            Assert.AreEqual(0d, link.InverseLink(0d), 1e-12);
            Assert.AreEqual(0.5d, link.InverseLink(0.5d * Math.Log(3d)), 1e-12);
            Assert.AreEqual(-0.5d, link.InverseLink(-0.5d * Math.Log(3d)), 1e-12);
        }

        /// <summary>
        /// Verifies round-trip behavior over the correlation domain.
        /// </summary>
        [TestMethod]
        public void RoundTrip_RecoversInput()
        {
            var link = new FisherZLink();
            foreach (double x in new[] { -0.99d, -0.75d, -0.25d, 0d, 0.25d, 0.75d, 0.99d })
            {
                Assert.AreEqual(x, link.InverseLink(link.Link(x)), RoundTripTol);
            }
        }

        /// <summary>
        /// Verifies derivative values and finite-difference consistency.
        /// </summary>
        [TestMethod]
        public void DLink_FiniteDifferenceConsistency()
        {
            var link = new FisherZLink();
            foreach (double x in new[] { -0.8d, -0.25d, 0d, 0.25d, 0.8d })
            {
                double expected = 1d / (1d - x * x);
                double finiteDifference = (link.Link(x + DeltaH) - link.Link(x - DeltaH)) / (2d * DeltaH);
                Assert.AreEqual(expected, link.DLink(x), 1e-12);
                Assert.AreEqual(finiteDifference, link.DLink(x), DerivativeTol);
            }
        }

        /// <summary>
        /// Verifies domain checks.
        /// </summary>
        [TestMethod]
        public void Link_OutsideOpenInterval_Throws()
        {
            var link = new FisherZLink();

            Assert.Throws<ArgumentOutOfRangeException>(() => link.Link(-1d));
            Assert.Throws<ArgumentOutOfRangeException>(() => link.Link(1d));
            Assert.Throws<ArgumentOutOfRangeException>(() => link.DLink(-1d));
            Assert.Throws<ArgumentOutOfRangeException>(() => link.DLink(1d));
        }

        /// <summary>
        /// Verifies factory and XML construction.
        /// </summary>
        [TestMethod]
        public void Factory_CreatesFisherZLink()
        {
            Assert.IsInstanceOfType(LinkFunctionFactory.Create(LinkFunctionType.FisherZ), typeof(FisherZLink));
            Assert.IsInstanceOfType(LinkFunctionFactory.CreateFromXElement(new FisherZLink().ToXElement()), typeof(FisherZLink));
        }
    }
}
