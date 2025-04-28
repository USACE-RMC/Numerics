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
using System.Collections.Generic;
using System.Linq;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Contains functions for computing the Multiple Grubbs Beck low outlier test.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// Cohn, T. A., England, J. F., Berenbrock, C. E., Mason, R. R., Stedinger, J. R., and Lamontagne, J. R. (2013). 
    /// A generalized Grubbs-Beck test statistic for detecting multiple potentially influential low outliers in flood series. 
    /// Water Resources Research, 49(8), 5047-5058.
    /// </item>
    /// <item>
    /// This code converted from the FORTRAN source code for PeakfqSA, which can be downloaded at:
    /// <see href = "https://sites.google.com/a/alumni.colostate.edu/jengland/resources" />
    /// </item>
    /// </list>
    /// </para>
    /// </remarks>
    public sealed class MultipleGrubbsBeckTest
    {
        /// <summary>
        /// Generalized Grubbs-Beck Test.
        /// <remarks>
        /// <para>
        /// This test identifies the number of low outliers using the generalized Grubbs-Beck test.
        /// </para>
        /// </remarks>
        /// </summary>
        /// <param name="X">A vector of flood in real space.</param>
        /// <returns>The number of low outliers.</returns>
        public static int Function(double[] X)
        {
            double Alphaout = 0.005d;
            double Alphain = 0.0d;
            double Alphazeroin = 0.1d;
            double maxFracLO = 0.5d;
            int N = X.Count();
            var zt = new double[N];
            var pvaluew = new double[N];
            for (int i = 0; i < N; i++)
            {
                zt[i] = Math.Log10(Math.Max(1.0E-88d, X[i]));
                pvaluew[i] = -99.0d;
            }

            // sort log flows from smallest to largest
            Array.Sort(zt);

            // set starting point for MGBT search at approximate median position (1/2 N)
            int n2 = (int)(N * maxFracLO);
            double S1 = 0d;
            double S2 = 0d;
            for (int i = N; i >= n2 + 2; i -= 1)
            {
                S1 += zt[i - 1];
                S2 += Math.Pow(zt[i - 1], 2d);
            }

            int NC;
            double XV;
            double XM;
            var W = new double[N];
            for (int i = n2; i >= 1; i -= 1)
            {
                S1 += zt[i];
                S2 += Math.Pow(zt[i], 2d);
                NC = N - i;
                XM = S1 / NC;
                XV = (S2 - NC * Math.Pow(XM, 2d)) / (NC - 1);
                W[i - 1] = (zt[i - 1] - XM) / Math.Sqrt(XV);
                pvaluew[i - 1] = GGBCRITP(N, i, W[i - 1]);
            }

            // Determine Number of low outliers in 2 Or 3 steps.
            // Based on TAC original code And JRS recommendations.
            // 
            // Step 1.   Outward sweep from median (always done).
            // alpha level of test = Alphaout
            // number of outliers = J1
            // 
            // Step 2.   Inward sweep from largest low outlier identified in Step 1.
            // alpha level of test = Alphain
            // number of outliers = J2
            // 
            // Step 3.   Inward sweep from smallest observation
            // alpha level of test = Alphazeroin
            // number of outliers = J3

            // Initialize counters
            int J1 = 0;   // Outward sweep number of low outliers
            int J2 = 0;   // Inward sweep number of low outliers
            int J3 = 0;   // Oth Inward sweep number of low outliers

            // 1) Outward sweep check: Loop over low flows up to median

            for (int i = n2; i >= 1; i -= 1)
            {
                if (pvaluew[i - 1] < Alphaout)
                {
                    J1 = i;
                    break;
                }
            }

            // 2) Inward sweep check with Alphain
            J2 = J1;
            for (int i = J1 + 1; i <= n2; i++)
            {
                if (pvaluew[i - 1] >= Alphain)
                {
                    J2 = i - 1;
                    break;
                }
            }

            // 3) Inward sweep check with Alphazeroin
            for (int i = 1; i <= n2; i++)
            {
                if (pvaluew[i - 1] >= Alphazeroin)
                {
                    J3 = i - 1;
                    break;
                }
            }

            // Set number of low outliers as max of 3 sweeps
            int MGBTP = Math.Max(J1, Math.Max(J2, J3));
            return MGBTP;
        }

        private static int _nIn;
        private static int _rIn;
        private static double _etaIn;

        /// <summary>
        /// Auxiliary routine used to compute p-values (GGCRITP) for a Generalized Grubbs-Beck Test.
        /// </summary>
        /// <returns>P-Values.</returns>
        private static double GGBCRITP(int N, int R, double ETA)
        {
            if (N < 10 | R > N / 2d)
            {
                return 0.5d;
            }
            else
            {
                _nIn = N;
                _rIn = R;
                _etaIn = ETA;
            }

            // The original FORTRAN source code utilized a globally adaptive Gauss-Kronrod integration method. 
            // This implementation, however, employs the adaptive Simpson's rule for integration.
            // The number of low outliers computed by this method is consistent with the results from the FORTRAN code.
            // Further testing is necessary to identify any edge cases where the adaptive Simpson's rule might prove insufficient.
            var sr = new AdaptiveSimpsonsRule(FGGB, 1E-16, 1 - 1E-16);
            sr.MaxDepth = 25;
            sr.ReportFailure = false;
            sr.Integrate();
            return sr.Status != IntegrationStatus.Failure ? sr.Result : double.NaN;
        }

        /// <summary>
        /// Auxiliary routine used in GGBCRITP
        /// </summary>
        private static double FGGB(double PZR)
        {
            double df, MuM, MuS2, VarM, VarS2, CovMS2;
            double EX1, EX2, EX3, EX4;
            double CovMS, VarS, alpha, beta;
            double MuMP, EtaP, H, Lambda, MuS, ncp, q, VarMP, PR, ZR, N2;
            double ANS;
            int N = _nIn;
            int R = _rIn;
            double ETA = _etaIn;

            // Compute the value of the r-th smallest obs. based on its order statistic
            N2 = N - R;
            var betaDist = new BetaDistribution(R, N + 1 - R);
            PR = betaDist.InverseCDF(PZR);
            ZR = Normal.StandardZ(PR);

            // Calculate the expected values of M, S2, S and their variances/covariances
            H = Normal.StandardPDF(ZR) / Math.Max(0.0000000001d, 1.0d - PR);
            EX1 = H;
            EX2 = 1d + H * ZR;
            EX3 = 2d * EX1 + H * Math.Pow(ZR, 2d);
            EX4 = 3d * EX2 + H * Math.Pow(ZR, 3d);
            MuM = EX1;
            MuS2 = EX2 - Math.Pow(EX1, 2d);
            VarM = MuS2 / N2;
            VarS2 = (EX4 - 4d * EX3 * EX1 + 6d * EX2 * Math.Pow(EX1, 2d) - 3d * Math.Pow(EX1, 4d) - Math.Pow(MuS2, 2d)) / N2 + 2.0d / ((N2 - 1.0d) * N2) * Math.Pow(MuS2, 2d);
            alpha = MuS2 * MuS2 / VarS2;
            beta = MuS2 / alpha;
            CovMS2 = (EX3 - 3d * EX2 * EX1 + 2d * Math.Pow(EX1, 3d)) / Math.Sqrt(N2 * (N2 - 1.0d));
            MuS = Math.Sqrt(beta) * Math.Exp(Gamma.LogGamma(alpha + 0.5d) - Gamma.LogGamma(alpha));
            CovMS = CovMS2 / (2d * MuS);
            VarS = MuS2 - Math.Pow(MuS, 2d);
            Lambda = CovMS / VarS;
            EtaP = ETA + Lambda;
            MuMP = MuM - Lambda * MuS;
            VarMP = VarM - CovMS * CovMS / VarS;
            df = 2.0d * alpha;
            ncp = (MuMP - ZR) / Math.Sqrt(VarMP);
            q = -Math.Sqrt(MuS2 / VarMP) * EtaP;
            var NCTDist = new NoncentralT(df, ncp);
            ANS = 1.0d - NCTDist.CDF(q);
            return ANS;
        }

        /// <summary>
        /// The original Grubbs and Beck test for detection of outliers. Test is performed at the 10% significance level based on a normal distribution.
        /// </summary>
        /// <param name="sample">The data sample.</param>
        /// <param name="XHi">Output. Values greater are considered high outliers.</param>
        /// <param name="XLo">Output. Values lower are considered low outliers.</param>
        public static void GrubbsBeckTest(IList<double> sample, out double XHi, out double XLo)
        {
            // The following polynomial approximation proposed by Pilon et al. (1985)
            int n = sample.Count;
            var logSample = new double[n];
            for (int i = 0; i < n; i++) logSample[i] = Math.Log(sample[i]);
            double mean = Statistics.Mean(logSample);
            double sd = Statistics.StandardDeviation(logSample);
            double Kn = -3.62201 + 6.28446 * Math.Pow(n, 0.25) - 2.49835 * Math.Pow(n, 0.5) + 0.491436 * Math.Pow(n, 0.75) - 0.037911 * n;
            XHi = Math.Exp(mean + Kn * sd);
            XLo = Math.Exp(mean - Kn * sd);
        }
    }
}