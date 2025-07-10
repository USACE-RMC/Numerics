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
using Numerics;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using System;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the Probability class. 
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </remarks>
    [TestClass]
    public class Test_Probability
    {
        /// <summary>
        /// Test basic probability rules for two random variables using textbook solutions.
        /// </summary>
        [TestMethod]
        public void Test_BasicRules()
        {
            double A = 0.15;
            double B = 0.62;
            double r = 0.75;

            // Test intersection
            Assert.AreEqual(0.146934, Probability.AAndB(A, B, r), 1E-6);
            // Test union
            Assert.AreEqual(0.623066, Probability.AOrB(A, B, r), 1E-6);
            // Test A not B
            Assert.AreEqual(0.003066, Probability.ANotB(A, B, r), 1E-6);
            // Test B not A
            Assert.AreEqual(0.473066, Probability.BNotA(A, B, r), 1E-6);
            // Test A given B
            Assert.AreEqual(0.236991, Probability.AGivenB(A, B, r), 1E-6);
            // Test B given A
            Assert.AreEqual(0.979562, Probability.BGivenA(A, B, r), 1E-6);
        }

        #region Joint Probability

        /// <summary>
        /// Test joint probability of ABCD assuming independence. 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_Independent()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;
            var joint = Probability.IndependentJointProbability(new double[] { A, B, C, D });
            Assert.AreEqual(0.021875, joint, 1E-6);
        }

        /// <summary>
        /// Test joint probability of ABCD using assuming independence using the Multivariate Normal (MVN). 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_Independent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.0;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var joint = Probability.JointProbability(new double[] { A, B, C, D }, new int[] { 1, 1, 1, 1 }, mvn);
            Assert.AreEqual(0.021875, joint, 1E-6);

        }

        /// <summary>
        /// Test joint probability of ABCD using assuming independence using the Product of Conditional Marginals (PCM). 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_Independent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.0;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            var joint = Probability.JointProbabilityHPCM(new double[] { A, B, C, D }, new int[] { 1, 1, 1, 1 }, corr);
            Assert.AreEqual(0.021875, joint, 1E-6);

        }

        /// <summary>
        /// Test joint probability of ABCD assuming perfect positive dependence. 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_PositivelyDependent()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;
            var joint = Probability.PositiveJointProbability(new double[] { A, B, C, D });
            Assert.AreEqual(0.25, joint, 1E-6);
        }

        /// <summary>
        /// Test joint probability of ABCD using assuming perfect positive dependence using the Multivariate Normal (MVN). 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_PositivelyDependent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.999;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var joint = Probability.JointProbability(new double[] { A, B, C, D }, new int[] { 1, 1, 1, 1 }, mvn);
            Assert.AreEqual(0.25, joint, 1E-6);

        }

        /// <summary>
        /// Test joint probability of ABCD using assuming perfect positive dependence using the Product of Conditional Marginals (PCM). 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_PositivelyDependent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.999;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            var joint = Probability.JointProbabilityHPCM(new double[] { A, B, C, D }, new int[] { 1, 1, 1, 1 }, corr);
            Assert.AreEqual(0.25, joint, 1E-6);

        }

        /// <summary>
        /// Test joint probability of ABCD assuming perfect negative dependence. 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_NegativelyDependent()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;
            var joint = Probability.NegativeJointProbability(new double[] { A, B, C, D });
            Assert.AreEqual(0.0, joint, 1E-6);
        }

        /// <summary>
        /// Test joint probability of ABCD using assuming perfect negative dependence using the Multivariate Normal (MVN). 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_NegativelyDependent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = -0.33;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var joint = Probability.JointProbability(new double[] { A, B, C, D }, new int[] { 1, 1, 1, 1 }, mvn);
            // True result comes from Monte Carlo simulation
            Assert.AreEqual(6.27974351606123E-14, joint, 1E-6);

        }

        /// <summary>
        /// Test joint probability of ABCD using assuming perfect negative dependence using the Product of Conditional Marginals (PCM). 
        /// </summary>
        [TestMethod]
        public void Test_JointABCD_NegativelyDependent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = -0.33;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            var joint = Probability.JointProbabilityHPCM(new double[] { A, B, C, D }, new int[] { 1, 1, 1, 1 }, corr);
            // True result comes from Monte Carlo simulation
            // PCM is less precise than MVN
            Assert.AreEqual(6.27974351606123E-14, joint, 1E-5);

        }

        #endregion

        #region Probability of Union

        /// <summary>
        /// Test union of ABCD assuming independence. 
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_Independent()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;
            // Uses De Morgan's rule
            var union = Probability.IndependentUnion(new double[] { A, B, C, D });
            Assert.AreEqual(0.878125, union, 1E-6);

        }

        /// <summary>
        /// Test union of ABCD assuming independence using Inclusion-Exclusion with the Multivariate Normal (MVN). 
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_Independent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.0;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get union from MVN
            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var union = Probability.Union(new double[] { A, B, C, D }, mvn);
            Assert.AreEqual(0.878125, union, 1E-6);

        }

        /// <summary>
        /// Test union of ABCD assuming independence using Inclusion-Exclusion with the Product of Conditional Marginals (PCM). 
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_Independent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.0;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get union from PCM
            var union = Probability.UnionPCM(new double[] { A, B, C, D }, corr);
            Assert.AreEqual(0.878125, union, 1E-6);

        }

        /// <summary>
        /// Test union of ABCD assuming perfect positive dependence.
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_PositivelyDependent()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;
            // Estimates the unimodal bound
            var union = Probability.PositivelyDependentUnion(new double[] { A, B, C, D });
            Assert.AreEqual(0.5, union, 1E-6);
        }

        /// <summary>
        /// Test union of ABCD assuming perfect positive dependence using Inclusion-Exclusion with the Multivariate Normal (MVN).
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_PositivelyDependent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.999;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get union from MVN
            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var union = Probability.Union(new double[] { A, B, C, D }, mvn);
            Assert.AreEqual(0.5, union, 1E-2);
        }

        /// <summary>
        /// Test union of ABCD assuming perfect positive dependence using Inclusion-Exclusion with the Product of Conditional Marginals (PCM).
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_PositivelyDependent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.999;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get union from PCM
            var union = Probability.UnionPCM(new double[] { A, B, C, D }, corr);
            // PCM is less precise than MVN
            Assert.AreEqual(0.5, union, 1E-2);
        }

        /// <summary>
        /// Test union of ABCD assuming perfect negative dependence.
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_NegativelyDependent()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;
            // Estimates the unimodal bound
            var union = Probability.NegativelyDependentUnion(new double[] { A, B, C, D });
            Assert.AreEqual(1.0, union, 1E-6);
        }

        /// <summary>
        /// Test union of ABCD assuming perfect negative dependence using Inclusion-Exclusion with the Multivariate Normal (MVN). 
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_NegativelyDependent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = -0.33;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get union from MVN
            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var union = Probability.Union(new double[] { A, B, C, D }, mvn);

            // True result comes from Monte Carlo simulation
            Assert.AreEqual(0.985016583, union, 1E-5);
        }

        /// <summary>
        /// Test union of ABCD assuming perfect negative dependence using Inclusion-Exclusion with the Product of Conditional Marginals (PCM).
        /// </summary>
        [TestMethod]
        public void Test_UnionABCD_NegativelyDependent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = -0.33;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get union from PCM
            var union = Probability.UnionPCM(new double[] { A, B, C, D }, corr);

            // True result comes from Monte Carlo simulation
            // PCM is less precise than MVN
            Assert.AreEqual(0.985016583, union, 1E-3);
        }

        #endregion

        #region Exclusive Probability

        /// <summary>
        /// Test exclusive probabilities of ABCD assuming independence. 
        /// </summary>
        [TestMethod]
        public void Test_ExclusiveABCD_Independent()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;
            var exclusive = Probability.IndependentExclusive(new double[] { A, B, C, D });

            // Test exclusive
            var trueE = new double[] { 0.040625, 0.065625, 0.121875, 0.121875, 0.021875, 0.040625, 0.040625, 0.065625, 0.065625, 0.121875, 0.021875, 0.021875, 0.040625, 0.065625, 0.021875 };
            for (int i = 0; i < exclusive.Length; i++)
            {
                Assert.AreEqual(trueE[i], exclusive[i], 1E-6);
            }

            // Test that sum = union
            var union = Tools.Sum(exclusive);
            Assert.AreEqual(0.878125, union, 1E-6);
        }

        /// <summary>
        /// Test exclusive probabilities of ABCD assuming independence using Inclusion-Exclusion with the Multivariate Normal (MVN). 
        /// </summary>
        [TestMethod]
        public void Test_ExclusiveABCD_Independent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.0;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get exclusive probabilities from MVN
            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var exclusive = Probability.Exclusive(new double[] { A, B, C, D }, mvn);

            // Test exclusive
            var trueE = new double[] { 0.040625, 0.065625, 0.121875, 0.121875, 0.021875, 0.040625, 0.040625, 0.065625, 0.065625, 0.121875, 0.021875, 0.021875, 0.040625, 0.065625, 0.021875 };
            for (int i = 0; i < exclusive.Length; i++)
            {
                Assert.AreEqual(trueE[i], exclusive[i],  1E-6);
            }

            // Test that sum = union
            var union = Tools.Sum(exclusive);
            Assert.AreEqual(0.878125, union,  1E-6);
        }

        /// <summary>
        /// Test exclusive probabilities of ABCD assuming independence using Inclusion-Exclusion with the Product of Conditional Marginals (PCM).
        /// </summary>
        [TestMethod]
        public void Test_ExclusiveABCD_Independent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.0;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get exclusive probabilities from PCM
            var exclusive = Probability.ExclusivePCM(new double[] { A, B, C, D }, corr);

            // Test exclusive
            var trueE = new double[] { 0.040625, 0.065625, 0.121875, 0.121875, 0.021875, 0.040625, 0.040625, 0.065625, 0.065625, 0.121875, 0.021875, 0.021875, 0.040625, 0.065625, 0.021875 };
            for (int i = 0; i < exclusive.Length; i++)
            {
                Assert.AreEqual(trueE[i], exclusive[i], 1E-6);
            }

            // Test that sum = union
            var union = Tools.Sum(exclusive);
            Assert.AreEqual(0.878125, union, 1E-6);
        }

        /// <summary>
        /// Test exclusive probabilities of ABCD assuming perfect positive dependence. 
        /// </summary>
        [TestMethod]
        public void Test_ExclusiveABCD_PositivelyDependent()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;
            var exclusive = Probability.PositivelyDependentExclusive(new double[] { A, B, C, D });

            // Test exclusive
            var trueE = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15, 0, 0, 0, 0.1, 0.25 };
            for (int i = 0; i < exclusive.Length; i++)
            {
                Assert.AreEqual(trueE[i], exclusive[i], 1E-6);
            }

            // Test that sum = union
            var union = Tools.Sum(exclusive);
            Assert.AreEqual(0.5, union, 1E-6);
        }

        /// <summary>
        /// Test exclusive probabilities of ABCD assuming perfect positive dependence using Inclusion-Exclusion with the Multivariate Normal (MVN). 
        /// </summary>
        [TestMethod]
        public void Test_ExclusiveABCD_PositivelyDependent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.999;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get exclusive probabilities from MVN
            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var exclusive = Probability.Exclusive(new double[] { A, B, C, D }, mvn);

            // Both MVN and PCM have poor precision on this test
            // Test exclusive
            var trueE = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15, 0, 0, 0, 0.1, 0.25 };
            for (int i = 0; i < exclusive.Length; i++)
            {
                Assert.AreEqual(trueE[i], exclusive[i], 1E-2);
            }

            // Test that sum = union
            var union = Tools.Sum(exclusive);
            Assert.AreEqual(0.5, union, 1E-2);
        }

        /// <summary>
        /// Test exclusive probabilities of ABCD assuming perfect positive dependence using Inclusion-Exclusion with the Product of Conditional Marginals (PCM).
        /// </summary>
        [TestMethod]
        public void Test_ExclusiveABCD_PositivelyDependent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;

            // Set correlation matrix
            var r = 0.999;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get exclusive probabilities from PCM
            var exclusive = Probability.ExclusivePCM(new double[] { A, B, C, D }, corr);

            // Both MVN and PCM have poor precision on this test. PCM is less precise than MVN.
            // Test exclusive
            var trueE = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15, 0, 0, 0, 0.1, 0.25 };
            for (int i = 0; i < exclusive.Length; i++)
            {
                Assert.AreEqual(trueE[i], exclusive[i], 1E-1);
            }

            // Test that sum = union
            var union = Tools.Sum(exclusive);
            Assert.AreEqual(0.5, union, 1E-1);
        }

        /// <summary>
        /// Test exclusive probabilities of ABCD assuming perfect negative dependence using Inclusion-Exclusion with the Multivariate Normal (MVN). 
        /// </summary>
        [TestMethod]
        public void Test_ExclusiveABCD_NegativelyDependent_MVN()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;


            var r = -0.33;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get exclusive probabilities from MVN
            var mvn = new MultivariateNormal(new double[] { 0, 0, 0, 0 }, corr) { MVNUNI = new Random(12345) };
            var exclusive = Probability.Exclusive(new double[] { A, B, C, D }, mvn);

            // Test exclusive
            var trueE = new double[] { 0.0591144184207316,  0.0841585674381143, 0.137370396630316,  0.137376624193301,  0.0381933344506983, 0.0639822078713156, 0.063975145059157,  0.0962373338116287, 0.0962331282359124, 0.160388585792858,  0.00596032423933968,    0.00596536506423,   0.012809204894465,  0.0232519467600138, 6.27974351606123E-14 };
            for (int i = 0; i < exclusive.Length; i++)
            {
                Assert.AreEqual(exclusive[i], trueE[i], 1E-4);
            }

            // Test that sum = union
            var union = Tools.Sum(exclusive);
            Assert.AreEqual(0.985016582862144, union, 1E-3);
        }

        /// <summary>
        /// Test exclusive probabilities of ABCD assuming perfect negative dependence using Inclusion-Exclusion with the Product of Conditional Marginals (PCM).
        /// </summary>
        [TestMethod]
        public void Test_ExclusiveABCD_NegativelyDependent_PCM()
        {
            double A = 0.25;
            double B = 0.35;
            double C = 0.5;
            double D = 0.5;


            var r = -0.33;
            var corr = new double[,]
            { { 1, r, r, r },
              { r, 1, r, r },
              { r, r, 1, r },
              { r, r, r, 1 }};

            // Get exclusive probabilities from PCM
            var exclusive = Probability.ExclusivePCM(new double[] { A, B, C, D }, corr);

            // Test exclusive
            var trueE = new double[] { 0.0591144184207316, 0.0841585674381143, 0.137370396630316, 0.137376624193301, 0.0381933344506983, 0.0639822078713156, 0.063975145059157, 0.0962373338116287, 0.0962331282359124, 0.160388585792858, 0.00596032423933968, 0.00596536506423, 0.012809204894465, 0.0232519467600138, 6.27974351606123E-14 };
            for (int i = 0; i < exclusive.Length; i++)
            {
                Assert.AreEqual(exclusive[i], trueE[i], 1E-3);
            }

            // True result comes from Monte Carlo simulation
            // PCM is less precise than MVN
            // Test that sum = union
            var union = Tools.Sum(exclusive);
            Assert.AreEqual(0.985016582862144, union, 1E-2);
        }




        #endregion


    }
}
