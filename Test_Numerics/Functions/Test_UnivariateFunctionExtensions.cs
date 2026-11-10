using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Distributions;
using Numerics.Functions;

namespace Functions
{
    /// <summary>
    /// Unit tests for the list overloads in <see cref="UnivariateFunctionExtensions"/>.
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
    public class Test_UnivariateFunctionExtensions
    {
        /// <summary>
        /// Returns one instance of every function class in the library.
        /// </summary>
        private static IEnumerable<IUnivariateFunction> Functions()
        {
            yield return new LinearFunction(2.0, 0.5);
            yield return new PowerFunction(3.0, 1.5, 1.0);
            yield return new SegmentedPowerFunction(1);
            var opd = new UncertainOrderedPairedData(new double[] { 50, 100, 150, 200, 250 },
                new UnivariateDistributionBase[] { new Deterministic(100), new Deterministic(200), new Deterministic(300), new Deterministic(400), new Deterministic(500) },
                true, SortOrder.Ascending, true, SortOrder.Ascending, UnivariateDistributionType.Deterministic);
            yield return new TabularFunction(opd);
            yield return new CompositeFunction(new IUnivariateFunction[] { new LinearFunction(2.0, 0.5), new LinearFunction(1.0, 1.0) });
        }

        /// <summary>
        /// The list overloads give the same bits as the scalar members, element by element.
        /// </summary>
        [TestMethod]
        public void Test_ListOverloads_MatchScalarBits()
        {
            var x = new double[] { 60, 75.5, 120, 199.99, 240 };
            foreach (var function in Functions())
            {
                string name = function.GetType().Name;
                var values = function.Function(x);
                Assert.HasCount(x.Length, values);
                for (int i = 0; i < x.Length; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(function.Function(x[i])), BitConverter.DoubleToInt64Bits(values[i]), name);

                var inverse = function.InverseFunction(values);
                for (int i = 0; i < x.Length; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(function.InverseFunction(values[i])), BitConverter.DoubleToInt64Bits(inverse[i]), name);
            }
        }

        /// <summary>
        /// Known values of the linear function y = 2 + 0.5 x and its inverse; an empty list gives an empty array.
        /// </summary>
        [TestMethod]
        public void Test_ListOverloads_KnownValues()
        {
            IUnivariateFunction linear = new LinearFunction(2.0, 0.5);
            CollectionAssert.AreEqual(new[] { 2.0, 2.5, 7.0 }, linear.Function(new[] { 0.0, 1.0, 10.0 }));
            CollectionAssert.AreEqual(new[] { 0.0, 1.0, 10.0 }, linear.InverseFunction(new List<double> { 2.0, 2.5, 7.0 }));
            Assert.IsEmpty(linear.Function(new double[0]));
        }

        /// <summary>
        /// Null inputs throw.
        /// </summary>
        [TestMethod]
        public void Test_ListOverloads_InvalidInputs()
        {
            IUnivariateFunction linear = new LinearFunction(2.0, 0.5);
            IList<double> none = null!;
            Assert.Throws<ArgumentNullException>(() => linear.Function(none));
            Assert.Throws<ArgumentNullException>(() => linear.InverseFunction(none));
            Assert.Throws<ArgumentNullException>(() => UnivariateFunctionExtensions.Function(null!, new double[1]));
            Assert.Throws<ArgumentNullException>(() => UnivariateFunctionExtensions.InverseFunction(null!, new double[1]));
        }
    }
}