using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Functions;
using System;
using System.Xml.Linq;

namespace Functions
{
    /// <summary>
    /// Unit tests for the <see cref="YeoJohnsonLink"/> class.
    /// </summary>
    [TestClass]
    public class Test_YeoJohnsonLink
    {
        private const double RoundTripTol = 1e-9;
        private const double DeltaH = 1e-7;
        private const double DerivativeTol = 1e-4;

        /// <summary>
        /// Verifies the default constructor uses the identity-transform lambda.
        /// </summary>
        [TestMethod]
        public void Constructor_Default_LambdaIsOne()
        {
            var link = new YeoJohnsonLink();

            Assert.AreEqual(1d, link.Lambda, 1e-12);
        }

        /// <summary>
        /// Verifies the lambda constructor stores the requested value.
        /// </summary>
        [TestMethod]
        public void Constructor_Lambda_StoresValue()
        {
            var link = new YeoJohnsonLink(0.5d);

            Assert.AreEqual(0.5d, link.Lambda, 1e-12);
        }

        /// <summary>
        /// Verifies the values constructor rejects null.
        /// </summary>
        [TestMethod]
        public void Constructor_Values_Null_Throws()
        {
            Assert.Throws<ArgumentNullException>(() => new YeoJohnsonLink((double[])null!));
        }

        /// <summary>
        /// Verifies the values constructor rejects degenerate samples.
        /// </summary>
        [TestMethod]
        public void Constructor_Values_SingleElement_Throws()
        {
            Assert.Throws<ArgumentException>(() => new YeoJohnsonLink(new[] { 1d }));
        }

        /// <summary>
        /// Verifies lambda fitting produces a finite value for a non-degenerate sample.
        /// </summary>
        [TestMethod]
        public void Constructor_Values_ProducesFiniteLambda()
        {
            var link = new YeoJohnsonLink(new[] { -2d, -1d, -0.25d, 0d, 0.5d, 1d, 3d });

            Assert.IsTrue(IsFinite(link.Lambda));
        }

        /// <summary>
        /// Verifies the XML constructor rejects null.
        /// </summary>
        [TestMethod]
        public void Constructor_XElement_Null_Throws()
        {
            Assert.Throws<ArgumentNullException>(() => new YeoJohnsonLink((XElement)null!));
        }

        /// <summary>
        /// Verifies lambda equal to 1 is the identity transform.
        /// </summary>
        [TestMethod]
        public void Link_LambdaOne_IsIdentity()
        {
            var link = new YeoJohnsonLink(1d);
            foreach (double x in new[] { -5d, -1d, 0d, 0.5d, 1d, 5d })
            {
                Assert.AreEqual(x, link.Link(x), 1e-10);
                Assert.AreEqual(1d, link.DLink(x), 1e-10);
            }
        }

        /// <summary>
        /// Verifies round-trip behavior for positive and negative values.
        /// </summary>
        [TestMethod]
        public void RoundTrip_PositiveAndNegativeValues()
        {
            var link = new YeoJohnsonLink(0.5d);
            foreach (double x in new[] { -10d, -5d, -1d, -0.1d, 0d, 0.1d, 1d, 5d, 10d })
            {
                double recovered = link.InverseLink(link.Link(x));
                Assert.AreEqual(x, recovered, RoundTripTol);
            }
        }

        /// <summary>
        /// Verifies the derivative against finite differences.
        /// </summary>
        [TestMethod]
        public void DLink_FiniteDifferenceConsistency()
        {
            var link = new YeoJohnsonLink(0.5d);
            foreach (double x in new[] { -5d, -2d, -0.5d, 0.5d, 1d, 5d })
            {
                double finiteDifference = (link.Link(x + DeltaH) - link.Link(x - DeltaH)) / (2d * DeltaH);
                Assert.AreEqual(finiteDifference, link.DLink(x), DerivativeTol);
            }
        }

        /// <summary>
        /// Verifies XML serialization preserves lambda.
        /// </summary>
        [TestMethod]
        public void XmlRoundTrip_PreservesLambda()
        {
            var original = new YeoJohnsonLink(0.75d);
            var restored = new YeoJohnsonLink(original.ToXElement());

            Assert.AreEqual(original.Lambda, restored.Lambda, 1e-12);
        }

        /// <summary>
        /// Verifies factory construction for Yeo-Johnson links.
        /// </summary>
        [TestMethod]
        public void Factory_CreatesYeoJohnsonLink()
        {
            Assert.IsInstanceOfType(LinkFunctionFactory.Create(LinkFunctionType.YeoJohnson), typeof(YeoJohnsonLink));

            var restored = LinkFunctionFactory.CreateFromXElement(new YeoJohnsonLink(0.25d).ToXElement());
            Assert.IsInstanceOfType(restored, typeof(YeoJohnsonLink));
            Assert.AreEqual(0.25d, ((YeoJohnsonLink)restored).Lambda, 1e-12);
        }

        private static bool IsFinite(double value)
        {
            return !double.IsNaN(value) && !double.IsInfinity(value);
        }
    }
}
