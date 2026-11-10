using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Functions;

namespace Functions
{
    /// <summary>
    /// Unit tests for the list overloads in <see cref="LinkFunctionExtensions"/>.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// </list>
    /// </remarks>
    [TestClass]
    public class Test_LinkFunctionExtensions
    {
        /// <summary>
        /// Every link class in the library, each with real-space values inside its domain.
        /// </summary>
        private static IEnumerable<(ILinkFunction Link, double[] X)> Links()
        {
            var unit = new[] { 1e-12, 0.05, 0.3, 0.5, 0.975, 1 - 1e-12 };
            yield return (new IdentityLink(), new[] { -3.5, 0.0, 2.25, 1e300 });
            yield return (new LogLink(), new[] { 1e-300, 0.5, 1.0, Math.E, 1e10 });
            yield return (new LogitLink(), unit);
            yield return (new ProbitLink(), unit);
            yield return (new ComplementaryLogLogLink(), unit);
            yield return (new FisherZLink(), new[] { -0.999, -0.5, 0.0, 0.5, 0.999 });
            yield return (new YeoJohnsonLink(0.4), new[] { -4.0, -0.5, 0.0, 0.5, 8.0 });
        }

        /// <summary>
        /// The list overloads give the same bits as the scalar members, element by element.
        /// </summary>
        [TestMethod]
        public void Test_ListOverloads_MatchScalarBits()
        {
            foreach (var (link, x) in Links())
            {
                string name = link.GetType().Name;
                var eta = link.Link(x);
                var back = link.InverseLink(eta);
                var slope = link.DLink(x);
                for (int i = 0; i < x.Length; i++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(link.Link(x[i])), BitConverter.DoubleToInt64Bits(eta[i]), name);
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(link.InverseLink(eta[i])), BitConverter.DoubleToInt64Bits(back[i]), name);
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(link.DLink(x[i])), BitConverter.DoubleToInt64Bits(slope[i]), name);
                }
            }
        }

        /// <summary>
        /// Known values: scipy.stats.norm.ppf(0.975) = 1.959963984540054, log(-log(0.7)) = -1.0309304331587228
        /// and arctanh(0.5) = 0.5493061443340549.
        /// </summary>
        [TestMethod]
        public void Test_ListOverloads_KnownValues()
        {
            Assert.AreEqual(1.959963984540054, new ProbitLink().Link(new[] { 0.975 })[0], 1e-14);
            Assert.AreEqual(-1.0309304331587228, new ComplementaryLogLogLink().Link(new List<double> { 0.3 })[0], 1e-14);
            Assert.AreEqual(0.5493061443340549, new FisherZLink().Link(new[] { 0.5 })[0], 1e-15);
            CollectionAssert.AreEqual(new[] { 1.0, 0.5 }, new LogLink().DLink(new[] { 1.0, 2.0 }));
            Assert.IsEmpty(new LogLink().InverseLink(new double[0]));
        }

        /// <summary>
        /// Null inputs throw.
        /// </summary>
        [TestMethod]
        public void Test_ListOverloads_InvalidInputs()
        {
            ILinkFunction link = new LogLink();
            IList<double> none = null!;
            Assert.Throws<ArgumentNullException>(() => link.Link(none));
            Assert.Throws<ArgumentNullException>(() => link.InverseLink(none));
            Assert.Throws<ArgumentNullException>(() => link.DLink(none));
            Assert.Throws<ArgumentNullException>(() => LinkFunctionExtensions.Link(null!, new double[1]));
            Assert.Throws<ArgumentNullException>(() => LinkFunctionExtensions.InverseLink(null!, new double[1]));
            Assert.Throws<ArgumentNullException>(() => LinkFunctionExtensions.DLink(null!, new double[1]));
        }
    }
}