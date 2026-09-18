using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions.Copulas;
using System;
using System.Collections.Generic;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Unit tests for behavior supplied by the bivariate and Archimedean copula base classes.
    /// </summary>
    [TestClass]
    public class Test_BivariateCopula
    {
        /// <summary>
        /// A subclass that implements only the required product-copula primitives receives the
        /// finite-difference conditional and array-delegating inverse fallbacks from the base class.
        /// </summary>
        [TestMethod]
        public void Test_BaseConditionalFallbacks()
        {
            var copula = new ProductFallbackCopula();

            Assert.AreEqual(0.63d, copula.ConditionalCDF(0.37d, 0.63d), 1E-10);
            Assert.AreEqual(0.63d, copula.ConditionalCDF(0d, 0.63d), 1E-10);
            Assert.AreEqual(0.63d, copula.ConditionalCDF(1d, 0.63d), 1E-10);
            Assert.AreEqual(0.42d, copula.InverseConditionalCDF(0.37d, 0.42d), 0d);
        }

        /// <summary>
        /// Every generic Archimedean h-function returns the exact dependent-variable boundaries
        /// for an interior conditioning probability, without evaluating singular generators.
        /// </summary>
        [TestMethod]
        public void Test_ArchimedeanConditionalCDF_DependentVariableBoundaries()
        {
            var copulas = new ArchimedeanCopula[]
            {
                new AMHCopula(0.5d),
                new ClaytonCopula(2d),
                new FrankCopula(5d),
                new GumbelCopula(2d),
                new JoeCopula(2d)
            };

            foreach (ArchimedeanCopula copula in copulas)
            {
                Assert.AreEqual(0d, copula.ConditionalCDF(0.4d, 0d), 0d, copula.DisplayName);
                Assert.AreEqual(1d, copula.ConditionalCDF(0.4d, 1d), 0d, copula.DisplayName);
            }
        }

        /// <summary>
        /// Minimal product copula used to exercise the virtual base implementations directly.
        /// </summary>
        private sealed class ProductFallbackCopula : BivariateCopula
        {
            /// <inheritdoc/>
            public override CopulaType Type => CopulaType.Independence;

            /// <inheritdoc/>
            public override double ThetaMinimum => 0d;

            /// <inheritdoc/>
            public override double ThetaMaximum => 0d;

            /// <inheritdoc/>
            public override string[,] ParameterToString => new string[0, 2];

            /// <inheritdoc/>
            public override string ParameterNameShortForm => string.Empty;

            /// <inheritdoc/>
            public override string DisplayName => "Product fallback";

            /// <inheritdoc/>
            public override string ShortDisplayName => "Product";

            /// <inheritdoc/>
            public override double UpperTailDependence => 0d;

            /// <inheritdoc/>
            public override double LowerTailDependence => 0d;

            /// <inheritdoc/>
            public override int NumberOfCopulaParameters => 0;

            /// <inheritdoc/>
            public override double[] GetCopulaParameters => Array.Empty<double>();

            /// <inheritdoc/>
            public override double PDF(double u, double v) => 1d;

            /// <inheritdoc/>
            public override double CDF(double u, double v) => u * v;

            /// <inheritdoc/>
            public override double[] InverseCDF(double u, double v) => new double[] { u, v };

            /// <inheritdoc/>
            public override void SetCopulaParameters(double[] parameters)
            {
            }

            /// <inheritdoc/>
            public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
            {
                return new double[0, 2];
            }

            /// <inheritdoc/>
            public override ArgumentOutOfRangeException ValidateParameter(double parameter, bool throwException)
            {
                return null;
            }

            /// <inheritdoc/>
            public override BivariateCopula Clone()
            {
                return new ProductFallbackCopula();
            }
        }
    }
}
