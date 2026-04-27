using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Mathematics.Integration
{
    /// <summary>
    /// Test integration methods.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass()]
    public class Test_Integration
    {
        /// <summary>
        /// Test Gauss-Legendre.
        /// </summary>
        [TestMethod()]
        public void Test_GaussLegendre()
        {
            double e = Numerics.Mathematics.Integration.Integration.GaussLegendre(Integrands.FX3, 0d, 1d);
            double val = 0.25d;
            Assert.AreEqual(e, val, 1E-3);
        }

        /// <summary>
        /// Test 20-point Gauss-Legendre with polynomial, trigonometric, and logarithmic integrands.
        /// </summary>
        [TestMethod()]
        public void Test_GaussLegendre20()
        {
            // x^3 over [0,1] = 0.25 (exact for degree <= 39)
            double e1 = Numerics.Mathematics.Integration.Integration.GaussLegendre20(Integrands.FX3, 0d, 1d);
            Assert.AreEqual(0.25d, e1, 1E-14);

            // Cosine over [0,1] = sin(1)
            double e2 = Numerics.Mathematics.Integration.Integration.GaussLegendre20(Math.Cos, 0d, 1d);
            Assert.AreEqual(Math.Sin(1d), e2, 1E-14);

            // log(x) over [1,2] = 2*ln(2) - 1
            double e3 = Numerics.Mathematics.Integration.Integration.GaussLegendre20(Math.Log, 1d, 2d);
            Assert.AreEqual(2d * Math.Log(2d) - 1d, e3, 1E-14);
        }

        /// <summary>
        /// Test Trapezoidal Rule Method.
        /// </summary>
        [TestMethod()]
        public void Test_TrapezoidalRule()
        {
            double e = Numerics.Mathematics.Integration.Integration.TrapezoidalRule(Integrands.FX3, 0d, 1d, 1000);
            double val = 0.25d;
            Assert.AreEqual(e, val, 1E-3);

        }


        /// <summary>
        /// Test Simpsons Rule Method. 
        /// </summary>
        [TestMethod()]
        public void Test_SimpsonsRule()
        {
            double e = Numerics.Mathematics.Integration.Integration.SimpsonsRule(Integrands.FX3, 0d, 1d, 1000);
            double val = 0.25d;
            Assert.AreEqual(e, val, 1E-3);

        }

        /// <summary>
        /// Midpoint Method. 
        /// </summary>
        [TestMethod()]
        public void Test_MidPoint()
        {
            double e = Numerics.Mathematics.Integration.Integration.Midpoint(Integrands.FX3, 0d, 1d, 1000);
            double val = 0.25d;
            Assert.AreEqual(e, val, 1E-3);

        }
    }
}