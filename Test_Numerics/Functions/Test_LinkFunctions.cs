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

using System;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Functions;

namespace Functions
{
    /// <summary>
    /// Unit tests for the link function classes: ILinkFunction implementations,
    /// LinkController, LinkFunctionFactory, and LinkFunctionType.
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
    public class Test_LinkFunctions
    {
        /// <summary>
        /// Finite-difference step size for derivative verification.
        /// </summary>
        private const double DeltaH = 1e-7;

        /// <summary>
        /// Tolerance for round-trip (Link then InverseLink) identity tests.
        /// </summary>
        private const double RoundTripTol = 1E-10;

        /// <summary>
        /// Tolerance for derivative (finite-difference vs analytic) tests.
        /// </summary>
        private const double DerivativeTol = 1E-5;

        // ──────────────────────────────────────────────
        //  IdentityLink
        // ──────────────────────────────────────────────

        /// <summary>
        /// Test that IdentityLink.Link returns x unchanged.
        /// </summary>
        [TestMethod]
        public void Test_IdentityLink_Link()
        {
            var link = new IdentityLink();
            double[] values = { -100.0, -1.5, 0.0, 1.5, 100.0, double.MaxValue, double.MinValue };
            foreach (double x in values)
            {
                Assert.AreEqual(x, link.Link(x), 0.0);
            }
        }

        /// <summary>
        /// Test that IdentityLink.InverseLink returns eta unchanged.
        /// </summary>
        [TestMethod]
        public void Test_IdentityLink_InverseLink()
        {
            var link = new IdentityLink();
            double[] values = { -100.0, -1.5, 0.0, 1.5, 100.0 };
            foreach (double eta in values)
            {
                Assert.AreEqual(eta, link.InverseLink(eta), 0.0);
            }
        }

        /// <summary>
        /// Test that IdentityLink.DLink always returns 1.
        /// </summary>
        [TestMethod]
        public void Test_IdentityLink_DLink()
        {
            var link = new IdentityLink();
            double[] values = { -100.0, 0.0, 100.0, 42.0 };
            foreach (double x in values)
            {
                Assert.AreEqual(1.0, link.DLink(x), 0.0);
            }
        }

        /// <summary>
        /// Test round-trip: InverseLink(Link(x)) == x for IdentityLink.
        /// </summary>
        [TestMethod]
        public void Test_IdentityLink_RoundTrip()
        {
            var link = new IdentityLink();
            double[] values = { -1000.0, -1.0, 0.0, 1.0, 1000.0 };
            foreach (double x in values)
            {
                double recovered = link.InverseLink(link.Link(x));
                Assert.AreEqual(x, recovered, RoundTripTol);
            }
        }

        // ──────────────────────────────────────────────
        //  LogLink
        // ──────────────────────────────────────────────

        /// <summary>
        /// Test LogLink.Link against known values.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_Link_KnownValues()
        {
            var link = new LogLink();
            Assert.AreEqual(0.0, link.Link(1.0), 1E-12);
            Assert.AreEqual(Math.Log(2.0), link.Link(2.0), 1E-12);
            Assert.AreEqual(Math.Log(10.0), link.Link(10.0), 1E-12);
            Assert.AreEqual(Math.Log(0.5), link.Link(0.5), 1E-12);
        }

        /// <summary>
        /// Test LogLink.InverseLink against known values.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_InverseLink_KnownValues()
        {
            var link = new LogLink();
            Assert.AreEqual(1.0, link.InverseLink(0.0), 1E-12);
            Assert.AreEqual(Math.E, link.InverseLink(1.0), 1E-12);
            Assert.AreEqual(Math.Exp(2.0), link.InverseLink(2.0), 1E-12);
            Assert.AreEqual(Math.Exp(-1.0), link.InverseLink(-1.0), 1E-12);
        }

        /// <summary>
        /// Test LogLink.DLink against known values: h'(x) = 1/x.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_DLink_KnownValues()
        {
            var link = new LogLink();
            Assert.AreEqual(1.0, link.DLink(1.0), 1E-12);
            Assert.AreEqual(0.5, link.DLink(2.0), 1E-12);
            Assert.AreEqual(0.1, link.DLink(10.0), 1E-12);
            Assert.AreEqual(100.0, link.DLink(0.01), 1E-6);
        }

        /// <summary>
        /// Test round-trip: InverseLink(Link(x)) == x for LogLink over positive domain.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_RoundTrip()
        {
            var link = new LogLink();
            double[] values = { 0.001, 0.1, 1.0, 10.0, 100.0, 1e6 };
            foreach (double x in values)
            {
                double recovered = link.InverseLink(link.Link(x));
                Assert.AreEqual(x, recovered, x * 1E-10);
            }
        }

        /// <summary>
        /// Test LogLink.Link throws for non-positive x.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_Link_ThrowsForZero()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new LogLink().Link(0.0));
        }

        /// <summary>
        /// Test LogLink.Link throws for negative x.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_Link_ThrowsForNegative()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new LogLink().Link(-1.0));
        }

        /// <summary>
        /// Test LogLink.DLink throws for non-positive x.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_DLink_ThrowsForZero()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new LogLink().DLink(0.0));
        }

        /// <summary>
        /// Test LogLink derivative via finite difference.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_DerivativeConsistency()
        {
            var link = new LogLink();
            double[] testPoints = { 0.01, 0.5, 1.0, 5.0, 100.0 };
            foreach (double x in testPoints)
            {
                double finiteDiff = (link.Link(x + DeltaH) - link.Link(x - DeltaH)) / (2 * DeltaH);
                Assert.AreEqual(finiteDiff, link.DLink(x), DerivativeTol);
            }
        }

        // ──────────────────────────────────────────────
        //  LogitLink
        // ──────────────────────────────────────────────

        /// <summary>
        /// Test LogitLink.Link against known values: logit(0.5) = 0.
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_Link_KnownValues()
        {
            var link = new LogitLink();
            Assert.AreEqual(0.0, link.Link(0.5), 1E-12);
            // logit(0.731) = log(0.731/0.269) ≈ 0.99894
            Assert.AreEqual(Math.Log(0.731 / 0.269), link.Link(0.731), 1E-6);
            // logit(0.1) = log(1/9) ≈ -2.19722
            Assert.AreEqual(Math.Log(1.0 / 9.0), link.Link(0.1), 1E-10);
        }

        /// <summary>
        /// Test LogitLink.InverseLink (sigmoid) against known values.
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_InverseLink_KnownValues()
        {
            var link = new LogitLink();
            Assert.AreEqual(0.5, link.InverseLink(0.0), 1E-12);
            // sigmoid(large positive) ≈ 1
            Assert.AreEqual(1.0, link.InverseLink(100.0), 1E-10);
            // sigmoid(large negative) ≈ 0
            Assert.AreEqual(0.0, link.InverseLink(-100.0), 1E-10);
            // sigmoid(1) = 1/(1+e^{-1}) ≈ 0.7310586
            Assert.AreEqual(1.0 / (1.0 + Math.Exp(-1.0)), link.InverseLink(1.0), 1E-10);
        }

        /// <summary>
        /// Test LogitLink.DLink against known values: h'(x) = 1/(x(1-x)).
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_DLink_KnownValues()
        {
            var link = new LogitLink();
            // At x=0.5: h'(0.5) = 1/(0.5*0.5) = 4.0
            Assert.AreEqual(4.0, link.DLink(0.5), 1E-12);
            // At x=0.1: h'(0.1) = 1/(0.1*0.9) = 100/9 ≈ 11.1111
            Assert.AreEqual(1.0 / (0.1 * 0.9), link.DLink(0.1), 1E-10);
        }

        /// <summary>
        /// Test round-trip: InverseLink(Link(x)) == x for LogitLink over (0,1).
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_RoundTrip()
        {
            var link = new LogitLink();
            double[] values = { 0.001, 0.1, 0.25, 0.5, 0.75, 0.9, 0.999 };
            foreach (double x in values)
            {
                double recovered = link.InverseLink(link.Link(x));
                Assert.AreEqual(x, recovered, RoundTripTol);
            }
        }

        /// <summary>
        /// Test LogitLink.Link throws for x outside (0,1).
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_Link_ThrowsForZero()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new LogitLink().Link(0.0));
        }

        /// <summary>
        /// Test LogitLink.Link throws for x >= 1.
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_Link_ThrowsForOne()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new LogitLink().Link(1.0));
        }

        /// <summary>
        /// Test LogitLink.Link throws for negative x.
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_Link_ThrowsForNegative()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new LogitLink().Link(-0.5));
        }

        /// <summary>
        /// Test LogitLink sigmoid numerical stability for extreme eta values.
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_InverseLink_ExtremeValues()
        {
            var link = new LogitLink();
            // Very large positive eta: sigmoid should be very close to 1 without overflow
            double high = link.InverseLink(710.0);
            Assert.IsTrue(high > 0.999 && high <= 1.0);
            // Very large negative eta: sigmoid should be very close to 0 without underflow
            double low = link.InverseLink(-710.0);
            Assert.IsTrue(low >= 0.0 && low < 0.001);
        }

        /// <summary>
        /// Test LogitLink derivative via finite difference.
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_DerivativeConsistency()
        {
            var link = new LogitLink();
            double[] testPoints = { 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99 };
            foreach (double x in testPoints)
            {
                double finiteDiff = (link.Link(x + DeltaH) - link.Link(x - DeltaH)) / (2 * DeltaH);
                Assert.AreEqual(finiteDiff, link.DLink(x), DerivativeTol);
            }
        }

        // ──────────────────────────────────────────────
        //  ProbitLink
        // ──────────────────────────────────────────────

        /// <summary>
        /// Test ProbitLink.Link against known values: Phi^{-1}(0.5) = 0.
        /// </summary>
        [TestMethod]
        public void Test_ProbitLink_Link_KnownValues()
        {
            var link = new ProbitLink();
            Assert.AreEqual(0.0, link.Link(0.5), 1E-10);
            // Phi^{-1}(0.975) ≈ 1.95996
            Assert.AreEqual(1.95996, link.Link(0.975), 1E-4);
            // Phi^{-1}(0.025) ≈ -1.95996
            Assert.AreEqual(-1.95996, link.Link(0.025), 1E-4);
            // Phi^{-1}(0.8413) ≈ 1.0 (CDF at z=1)
            Assert.AreEqual(1.0, link.Link(0.8413), 1E-3);
        }

        /// <summary>
        /// Test ProbitLink.InverseLink against known values: Phi(0) = 0.5.
        /// </summary>
        [TestMethod]
        public void Test_ProbitLink_InverseLink_KnownValues()
        {
            var link = new ProbitLink();
            Assert.AreEqual(0.5, link.InverseLink(0.0), 1E-10);
            // Phi(1.96) ≈ 0.975
            Assert.AreEqual(0.975, link.InverseLink(1.96), 1E-3);
            // Phi(-1.96) ≈ 0.025
            Assert.AreEqual(0.025, link.InverseLink(-1.96), 1E-3);
        }

        /// <summary>
        /// Test round-trip: InverseLink(Link(x)) == x for ProbitLink over (0,1).
        /// </summary>
        [TestMethod]
        public void Test_ProbitLink_RoundTrip()
        {
            var link = new ProbitLink();
            double[] values = { 0.001, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.999 };
            foreach (double x in values)
            {
                double recovered = link.InverseLink(link.Link(x));
                Assert.AreEqual(x, recovered, RoundTripTol);
            }
        }

        /// <summary>
        /// Test ProbitLink.Link throws for x outside (0,1).
        /// </summary>
        [TestMethod]
        public void Test_ProbitLink_Link_ThrowsForZero()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new ProbitLink().Link(0.0));
        }

        /// <summary>
        /// Test ProbitLink.Link throws for x >= 1.
        /// </summary>
        [TestMethod]
        public void Test_ProbitLink_Link_ThrowsForOne()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new ProbitLink().Link(1.0));
        }

        /// <summary>
        /// Test ProbitLink derivative via finite difference.
        /// </summary>
        [TestMethod]
        public void Test_ProbitLink_DerivativeConsistency()
        {
            var link = new ProbitLink();
            double[] testPoints = { 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95 };
            foreach (double x in testPoints)
            {
                double finiteDiff = (link.Link(x + DeltaH) - link.Link(x - DeltaH)) / (2 * DeltaH);
                Assert.AreEqual(finiteDiff, link.DLink(x), DerivativeTol);
            }
        }

        // ──────────────────────────────────────────────
        //  ComplementaryLogLogLink
        // ──────────────────────────────────────────────

        /// <summary>
        /// Test CLogLog.Link against known values: h(x) = log(-log(1-x)).
        /// </summary>
        [TestMethod]
        public void Test_CLogLog_Link_KnownValues()
        {
            var link = new ComplementaryLogLogLink();
            // h(1-1/e) = h(0.63212...) = log(-log(1/e)) = log(1) = 0
            double x0 = 1.0 - 1.0 / Math.E;
            Assert.AreEqual(0.0, link.Link(x0), 1E-10);
            // h(0.5) = log(-log(0.5)) = log(log(2)) ≈ -0.36651
            Assert.AreEqual(Math.Log(Math.Log(2.0)), link.Link(0.5), 1E-10);
        }

        /// <summary>
        /// Test CLogLog.InverseLink against known values: h^{-1}(eta) = 1-exp(-exp(eta)).
        /// </summary>
        [TestMethod]
        public void Test_CLogLog_InverseLink_KnownValues()
        {
            var link = new ComplementaryLogLogLink();
            // h^{-1}(0) = 1 - exp(-1) = 1 - 1/e ≈ 0.63212
            Assert.AreEqual(1.0 - 1.0 / Math.E, link.InverseLink(0.0), 1E-10);
            // For large positive eta: result approaches 1
            Assert.AreEqual(1.0, link.InverseLink(10.0), 1E-4);
            // For large negative eta: result approaches 0
            Assert.AreEqual(0.0, link.InverseLink(-10.0), 1E-4);
        }

        /// <summary>
        /// Test round-trip: InverseLink(Link(x)) == x for CLogLog over (0,1).
        /// </summary>
        [TestMethod]
        public void Test_CLogLog_RoundTrip()
        {
            var link = new ComplementaryLogLogLink();
            double[] values = { 0.001, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.999 };
            foreach (double x in values)
            {
                double recovered = link.InverseLink(link.Link(x));
                Assert.AreEqual(x, recovered, RoundTripTol);
            }
        }

        /// <summary>
        /// Test CLogLog.Link throws for x outside (0,1).
        /// </summary>
        [TestMethod]
        public void Test_CLogLog_Link_ThrowsForZero()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new ComplementaryLogLogLink().Link(0.0));
        }

        /// <summary>
        /// Test CLogLog.Link throws for x >= 1.
        /// </summary>
        [TestMethod]
        public void Test_CLogLog_Link_ThrowsForOne()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new ComplementaryLogLogLink().Link(1.0));
        }

        /// <summary>
        /// Test CLogLog derivative via finite difference.
        /// </summary>
        [TestMethod]
        public void Test_CLogLog_DerivativeConsistency()
        {
            var link = new ComplementaryLogLogLink();
            double[] testPoints = { 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95 };
            foreach (double x in testPoints)
            {
                double finiteDiff = (link.Link(x + DeltaH) - link.Link(x - DeltaH)) / (2 * DeltaH);
                Assert.AreEqual(finiteDiff, link.DLink(x), DerivativeTol);
            }
        }

        // ──────────────────────────────────────────────
        //  LinkFunctionFactory
        // ──────────────────────────────────────────────

        /// <summary>
        /// Test that factory creates the correct type for each enum value.
        /// </summary>
        [TestMethod]
        public void Test_Factory_CreatesCorrectTypes()
        {
            Assert.IsInstanceOfType(LinkFunctionFactory.Create(LinkFunctionType.Identity), typeof(IdentityLink));
            Assert.IsInstanceOfType(LinkFunctionFactory.Create(LinkFunctionType.Log), typeof(LogLink));
            Assert.IsInstanceOfType(LinkFunctionFactory.Create(LinkFunctionType.Logit), typeof(LogitLink));
            Assert.IsInstanceOfType(LinkFunctionFactory.Create(LinkFunctionType.Probit), typeof(ProbitLink));
            Assert.IsInstanceOfType(LinkFunctionFactory.Create(LinkFunctionType.ComplementaryLogLog), typeof(ComplementaryLogLogLink));
        }

        /// <summary>
        /// Test that factory throws for invalid enum value.
        /// </summary>
        [TestMethod]
        public void Test_Factory_ThrowsForInvalidType()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => LinkFunctionFactory.Create((LinkFunctionType)999));
        }

        /// <summary>
        /// Test that factory-created links produce correct round-trip results.
        /// </summary>
        [TestMethod]
        public void Test_Factory_RoundTripViaFactory()
        {
            // Identity: domain = all reals
            var identity = LinkFunctionFactory.Create(LinkFunctionType.Identity);
            Assert.AreEqual(5.0, identity.InverseLink(identity.Link(5.0)), RoundTripTol);

            // Log: domain = (0, inf)
            var log = LinkFunctionFactory.Create(LinkFunctionType.Log);
            Assert.AreEqual(2.5, log.InverseLink(log.Link(2.5)), RoundTripTol);

            // Logit: domain = (0, 1)
            var logit = LinkFunctionFactory.Create(LinkFunctionType.Logit);
            Assert.AreEqual(0.3, logit.InverseLink(logit.Link(0.3)), RoundTripTol);

            // Probit: domain = (0, 1)
            var probit = LinkFunctionFactory.Create(LinkFunctionType.Probit);
            Assert.AreEqual(0.7, probit.InverseLink(probit.Link(0.7)), RoundTripTol);

            // CLogLog: domain = (0, 1)
            var cloglog = LinkFunctionFactory.Create(LinkFunctionType.ComplementaryLogLog);
            Assert.AreEqual(0.4, cloglog.InverseLink(cloglog.Link(0.4)), RoundTripTol);
        }

        // ──────────────────────────────────────────────
        //  LinkController
        // ──────────────────────────────────────────────

        /// <summary>
        /// Test that empty LinkController acts as identity for all parameters.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_Empty_ActsAsIdentity()
        {
            var ctrl = new LinkController();
            Assert.AreEqual(0, ctrl.Count);

            double[] x = { 100.0, 0.5, -0.1 };
            double[] eta = ctrl.Link(x);
            Assert.AreEqual(100.0, eta[0], 1E-12);
            Assert.AreEqual(0.5, eta[1], 1E-12);
            Assert.AreEqual(-0.1, eta[2], 1E-12);

            double[] xBack = ctrl.InverseLink(eta);
            Assert.AreEqual(100.0, xBack[0], 1E-12);
            Assert.AreEqual(0.5, xBack[1], 1E-12);
            Assert.AreEqual(-0.1, xBack[2], 1E-12);
        }

        /// <summary>
        /// Test LinkController with single LogLink on first parameter.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_SingleLink()
        {
            var ctrl = new LinkController(new LogLink());
            Assert.AreEqual(1, ctrl.Count);

            double[] x = { 2.0, 5.0, -3.0 };
            double[] eta = ctrl.Link(x);
            // First element transformed by LogLink
            Assert.AreEqual(Math.Log(2.0), eta[0], 1E-12);
            // Remaining elements pass through (no link registered)
            Assert.AreEqual(5.0, eta[1], 1E-12);
            Assert.AreEqual(-3.0, eta[2], 1E-12);
        }

        /// <summary>
        /// Test LinkController with mixed links (null = identity).
        /// </summary>
        [TestMethod]
        public void Test_LinkController_MixedLinks()
        {
            var ctrl = new LinkController(null, new LogLink(), null);
            Assert.AreEqual(3, ctrl.Count);

            double[] x = { 42.0, 10.0, -7.0 };
            double[] eta = ctrl.Link(x);
            Assert.AreEqual(42.0, eta[0], 1E-12);      // null = identity
            Assert.AreEqual(Math.Log(10.0), eta[1], 1E-12); // LogLink
            Assert.AreEqual(-7.0, eta[2], 1E-12);       // null = identity

            double[] xBack = ctrl.InverseLink(eta);
            Assert.AreEqual(42.0, xBack[0], 1E-12);
            Assert.AreEqual(10.0, xBack[1], 1E-10);
            Assert.AreEqual(-7.0, xBack[2], 1E-12);
        }

        /// <summary>
        /// Test LinkController round-trip with location=null, scale=Log, shape=Logit.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_RoundTrip()
        {
            var ctrl = new LinkController(null, new LogLink(), new LogitLink());
            double[] x = { 100.0, 5.0, 0.3 };
            double[] eta = ctrl.Link(x);
            double[] xBack = ctrl.InverseLink(eta);
            Assert.AreEqual(x[0], xBack[0], RoundTripTol);
            Assert.AreEqual(x[1], xBack[1], RoundTripTol);
            Assert.AreEqual(x[2], xBack[2], RoundTripTol);
        }

        /// <summary>
        /// Test LinkController handles array longer than link count (extra elements pass through).
        /// </summary>
        [TestMethod]
        public void Test_LinkController_ArrayLongerThanLinkCount()
        {
            var ctrl = new LinkController(new LogLink());
            double[] x = { 3.0, 7.0, 11.0, 13.0 };
            double[] eta = ctrl.Link(x);
            Assert.AreEqual(Math.Log(3.0), eta[0], 1E-12);
            Assert.AreEqual(7.0, eta[1], 1E-12);
            Assert.AreEqual(11.0, eta[2], 1E-12);
            Assert.AreEqual(13.0, eta[3], 1E-12);
        }

        /// <summary>
        /// Test LinkController indexer returns correct links and null for out-of-range.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_Indexer()
        {
            var logLink = new LogLink();
            var ctrl = new LinkController(null, logLink);
            Assert.IsNull(ctrl[0]);
            Assert.AreSame(logLink, ctrl[1]);
            Assert.IsNull(ctrl[2]);   // out of range
            Assert.IsNull(ctrl[-1]);  // negative index
        }

        /// <summary>
        /// Test ForLocationScaleShape factory method.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_ForLocationScaleShape()
        {
            var logLink = new LogLink();
            var ctrl = LinkController.ForLocationScaleShape(scaleLink: logLink);
            Assert.AreEqual(3, ctrl.Count);
            Assert.IsNull(ctrl[0]);
            Assert.AreSame(logLink, ctrl[1]);
            Assert.IsNull(ctrl[2]);
        }

        /// <summary>
        /// Test LinkJacobian returns correct diagonal matrix.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_LinkJacobian()
        {
            var ctrl = new LinkController(null, new LogLink());
            double[] x = { 42.0, 5.0 };
            var J = ctrl.LinkJacobian(x);

            // First element: identity => diagonal = 1.0
            Assert.AreEqual(1.0, J[0, 0], 1E-12);
            // Second element: LogLink => h'(5) = 1/5 = 0.2
            Assert.AreEqual(0.2, J[1, 1], 1E-12);
            // Off-diagonal = 0
            Assert.AreEqual(0.0, J[0, 1], 1E-12);
            Assert.AreEqual(0.0, J[1, 0], 1E-12);
        }

        /// <summary>
        /// Test LogDetJacobian for a single LogLink.
        /// For LogLink: deta/dtheta = 1/theta, so log|det J^{-1}| = -sum(-log|1/theta|) = log(theta).
        /// </summary>
        [TestMethod]
        public void Test_LinkController_LogDetJacobian()
        {
            var ctrl = new LinkController(new LogLink());
            // phi = [log(5)] in link-space
            double[] phi = { Math.Log(5.0) };
            double logDetJ = ctrl.LogDetJacobian(phi);
            // InverseLink(log5) = 5; DLink(5) = 1/5 = 0.2
            // log|det J^{-1}| = -log|0.2| = log(5) ≈ 1.60944
            Assert.AreEqual(Math.Log(5.0), logDetJ, 1E-10);
        }

        /// <summary>
        /// Test LogDetJacobian with identity (empty controller) returns 0.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_LogDetJacobian_Identity()
        {
            var ctrl = new LinkController();
            double[] phi = { 1.0, 2.0, 3.0 };
            double logDetJ = ctrl.LogDetJacobian(phi);
            Assert.AreEqual(0.0, logDetJ, 1E-12);
        }

        /// <summary>
        /// Test LogDetJacobian with multiple links.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_LogDetJacobian_MultipleLinks()
        {
            var ctrl = new LinkController(null, new LogLink(), new LogLink());
            // phi = [100, log(3), log(7)]
            double[] phi = { 100.0, Math.Log(3.0), Math.Log(7.0) };
            double logDetJ = ctrl.LogDetJacobian(phi);
            // index 0: null => no contribution (0)
            // index 1: theta=3, deta/dtheta=1/3 => contribution = -log(1/3) = log(3)
            // index 2: theta=7, deta/dtheta=1/7 => contribution = -log(1/7) = log(7)
            double expected = Math.Log(3.0) + Math.Log(7.0);
            Assert.AreEqual(expected, logDetJ, 1E-10);
        }

        #region XElement Serialization

        /// <summary>
        /// Test IdentityLink round-trip through ToXElement and CreateFromXElement.
        /// </summary>
        [TestMethod]
        public void Test_IdentityLink_RoundTrip_ToXElement()
        {
            var original = new IdentityLink();
            var xml = original.ToXElement();
            Assert.AreEqual("IdentityLink", xml.Name.LocalName);

            var restored = LinkFunctionFactory.CreateFromXElement(xml);
            Assert.IsInstanceOfType(restored, typeof(IdentityLink));
            Assert.AreEqual(original.Link(5.0), restored.Link(5.0), 1E-15);
            Assert.AreEqual(original.InverseLink(-3.0), restored.InverseLink(-3.0), 1E-15);
        }

        /// <summary>
        /// Test LogLink round-trip through ToXElement and CreateFromXElement.
        /// </summary>
        [TestMethod]
        public void Test_LogLink_RoundTrip_ToXElement()
        {
            var original = new LogLink();
            var xml = original.ToXElement();
            Assert.AreEqual("LogLink", xml.Name.LocalName);

            var restored = LinkFunctionFactory.CreateFromXElement(xml);
            Assert.IsInstanceOfType(restored, typeof(LogLink));
            Assert.AreEqual(original.Link(2.5), restored.Link(2.5), 1E-15);
            Assert.AreEqual(original.InverseLink(1.0), restored.InverseLink(1.0), 1E-15);
        }

        /// <summary>
        /// Test LogitLink round-trip through ToXElement and CreateFromXElement.
        /// </summary>
        [TestMethod]
        public void Test_LogitLink_RoundTrip_ToXElement()
        {
            var original = new LogitLink();
            var xml = original.ToXElement();
            Assert.AreEqual("LogitLink", xml.Name.LocalName);

            var restored = LinkFunctionFactory.CreateFromXElement(xml);
            Assert.IsInstanceOfType(restored, typeof(LogitLink));
            Assert.AreEqual(original.Link(0.7), restored.Link(0.7), 1E-15);
            Assert.AreEqual(original.InverseLink(0.5), restored.InverseLink(0.5), 1E-15);
        }

        /// <summary>
        /// Test ProbitLink round-trip through ToXElement and CreateFromXElement.
        /// </summary>
        [TestMethod]
        public void Test_ProbitLink_RoundTrip_ToXElement()
        {
            var original = new ProbitLink();
            var xml = original.ToXElement();
            Assert.AreEqual("ProbitLink", xml.Name.LocalName);

            var restored = LinkFunctionFactory.CreateFromXElement(xml);
            Assert.IsInstanceOfType(restored, typeof(ProbitLink));
            Assert.AreEqual(original.Link(0.3), restored.Link(0.3), 1E-15);
            Assert.AreEqual(original.InverseLink(-1.0), restored.InverseLink(-1.0), 1E-15);
        }

        /// <summary>
        /// Test ComplementaryLogLogLink round-trip through ToXElement and CreateFromXElement.
        /// </summary>
        [TestMethod]
        public void Test_ComplementaryLogLogLink_RoundTrip_ToXElement()
        {
            var original = new ComplementaryLogLogLink();
            var xml = original.ToXElement();
            Assert.AreEqual("ComplementaryLogLogLink", xml.Name.LocalName);

            var restored = LinkFunctionFactory.CreateFromXElement(xml);
            Assert.IsInstanceOfType(restored, typeof(ComplementaryLogLogLink));
            Assert.AreEqual(original.Link(0.5), restored.Link(0.5), 1E-15);
            Assert.AreEqual(original.InverseLink(0.0), restored.InverseLink(0.0), 1E-15);
        }

        /// <summary>
        /// Test that CreateFromXElement throws NotSupportedException for unknown element names.
        /// </summary>
        [TestMethod]
        public void Test_LinkFunctionFactory_CreateFromXElement_UnknownType_Throws()
        {
            var xml = new XElement("UnknownLink");
            Assert.Throws<NotSupportedException>(() => LinkFunctionFactory.CreateFromXElement(xml));
        }

        /// <summary>
        /// Test LinkController round-trip serialization with an empty controller.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_RoundTrip_Empty()
        {
            var original = new LinkController();
            var xml = original.ToXElement();
            Assert.AreEqual("LinkController", xml.Name.LocalName);

            var restored = new LinkController(xml);
            Assert.AreEqual(0, restored.Count);
        }

        /// <summary>
        /// Test LinkController round-trip serialization with all three link slots populated.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_RoundTrip_AllLinks()
        {
            var original = new LinkController(new IdentityLink(), new LogLink(), new LogitLink());
            Assert.AreEqual(3, original.Count);

            var xml = original.ToXElement();
            var restored = new LinkController(xml);

            Assert.AreEqual(3, restored.Count);
            Assert.IsInstanceOfType(restored[0], typeof(IdentityLink));
            Assert.IsInstanceOfType(restored[1], typeof(LogLink));
            Assert.IsInstanceOfType(restored[2], typeof(LogitLink));

            // Verify functional equivalence
            double[] x = { 5.0, 2.5, 0.7 };
            double[] originalLinked = original.Link(x);
            double[] restoredLinked = restored.Link(x);
            for (int i = 0; i < x.Length; i++)
                Assert.AreEqual(originalLinked[i], restoredLinked[i], 1E-15);
        }

        /// <summary>
        /// Test LinkController round-trip with null slots preserved.
        /// </summary>
        [TestMethod]
        public void Test_LinkController_RoundTrip_NullSlots()
        {
            var original = new LinkController(null, new LogLink(), null);
            Assert.AreEqual(3, original.Count);

            var xml = original.ToXElement();
            var restored = new LinkController(xml);

            Assert.AreEqual(3, restored.Count);
            Assert.IsNull(restored[0]);
            Assert.IsInstanceOfType(restored[1], typeof(LogLink));
            Assert.IsNull(restored[2]);

            // Verify functional equivalence: index 0 and 2 pass through, index 1 uses log
            double[] x = { 5.0, 2.5, 0.7 };
            double[] originalLinked = original.Link(x);
            double[] restoredLinked = restored.Link(x);
            for (int i = 0; i < x.Length; i++)
                Assert.AreEqual(originalLinked[i], restoredLinked[i], 1E-15);
        }

        #endregion
    }
}
