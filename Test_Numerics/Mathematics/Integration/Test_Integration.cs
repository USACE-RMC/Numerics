/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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