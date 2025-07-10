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

using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Distributions.Copulas
{
    /// <summary>
    /// A class for estimating the parameters of a bivariate copula.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class BivariateCopulaEstimation
    {

        /// <summary>
        /// Estimate the bivariate copula.
        /// </summary>
        /// <param name="estimationMethod"></param>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        public static void Estimate(ref BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY, CopulaEstimationMethod estimationMethod)
        {
            switch (estimationMethod)
            {
                case CopulaEstimationMethod.PseudoLikelihood:
                    MPL(copula, sampleDataX, sampleDataY);
                    break;
                case CopulaEstimationMethod.InferenceFromMargins:
                    IFM(copula, sampleDataX, sampleDataY);
                    break;
                case CopulaEstimationMethod.FullLikelihood:
                    MLE(copula, sampleDataX, sampleDataY);
                    break;
            }
        }

        /// <summary>
        /// The maximum pseudo likelihood method. 
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.When estimating with pseudo likelihood, this should be the plotting positions of the data.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.When estimating with pseudo likelihood, this should be the plotting positions of the data.</param>
        private static void MPL(BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY)
        {
            // Get constraints
            var LU = copula.ParameterConstraints(sampleDataX, sampleDataY);

            // Solve using Brent method
            Func<double, double> func = (x) =>
            {
                var C = copula.Clone();
                C.Theta = x;
                return C.PseudoLogLikelihood(sampleDataX, sampleDataY);
            };
            var brent = new BrentSearch(func, LU[0], LU[1]);
            brent.Maximize();
            copula.Theta = brent.BestParameterSet.Values[0];
        }

        /// <summary>
        /// The inference from margins method. 
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        private static void IFM(BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY)
        {
            // Get constraints
            var LU = copula.ParameterConstraints(sampleDataX, sampleDataY);

            // Solve using Brent method
            Func<double, double> func = (x) =>
            {
                var C = copula.Clone();
                C.Theta = x;
                return C.IFMLogLikelihood(sampleDataX, sampleDataY);
            };
            var brent = new BrentSearch(func, LU[0], LU[1]);
            brent.Maximize();
            copula.Theta = brent.BestParameterSet.Values[0];
        }

        /// <summary>
        /// The maximum likelihood estimation method.
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        private static void MLE(BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY)
        {
            // See if marginals are estimable
            IMaximumLikelihoodEstimation margin1 = (IMaximumLikelihoodEstimation)copula.MarginalDistributionX;
            IMaximumLikelihoodEstimation margin2 = (IMaximumLikelihoodEstimation)copula.MarginalDistributionY;
            if (margin1 == null || margin2 == null) throw new ArgumentOutOfRangeException("marginal distributions", "There marginal distributions must implement the IMaximumLikelihoodEstimation interface to use this method.");

            int np1 = copula.MarginalDistributionX.NumberOfParameters;
            int np2 = copula.MarginalDistributionY.NumberOfParameters;

            // Get constraints     
            var initials = new double[1 + np1 + np2];
            var lowers = new double[1 + np1 + np2];
            var uppers = new double[1 + np1 + np2];

            // Theta
            // get ranks of data
            var rank1 = Statistics.RanksInPlace(sampleDataX.ToArray());
            var rank2 = Statistics.RanksInPlace(sampleDataY.ToArray());
            // get plotting positions
            for (int i = 0; i < rank1.Length; i++)
            {
                rank1[i] = rank1[i] / (rank1.Length + 1d);
                rank2[i] = rank2[i] / (rank2.Length + 1d);
            }

            // Get constraints
            var LU = copula.ParameterConstraints(sampleDataX, sampleDataY);
            lowers[0] = LU[0];
            uppers[0] = LU[1];
            // Estimate copula using MPL
            MPL(copula, rank1, rank2);
            initials[0] = copula.Theta;

            // Estimate marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(sampleDataX, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(sampleDataY, ParameterEstimationMethod.MaximumLikelihood);
            

            var con = margin1.GetParameterConstraints(sampleDataX);
            var parms = copula.MarginalDistributionX.GetParameters;
            for (int i = 0; i < np1; i++)
            {
                initials[i + 1] = parms[i];
                lowers[i + 1] = con.Item2[i];
                uppers[i + 1] = con.Item3[i];
            }
            con = margin2.GetParameterConstraints(sampleDataY);
            parms = copula.MarginalDistributionY.GetParameters;
            for (int i = 0; i < np2; i++)
            {
                initials[i + 1 + np1] = parms[i];
                lowers[i + 1 + np1] = con.Item2[i];
                uppers[i + 1 + np1] = con.Item3[i];
            }

            // Log-likelihood function
            Func<double[], double> logLH = (double[] x) => {
                // Set copula
                var C = copula.Clone();
                C.Theta = x[0];

                // marginal 1
                var m1 = ((UnivariateDistributionBase)copula.MarginalDistributionX).Clone();
                var p1 = new double[np1];
                for (int i = 0; i < np1; i++)
                    p1[i] = x[i + 1];
                m1.SetParameters(p1);

                // marginal 2
                var m2 = ((UnivariateDistributionBase)copula.MarginalDistributionY).Clone();
                var p2 = new double[np2];
                for (int i = 0; i < np2; i++)
                    p2[i] = x[i + 1 + np1];
                m2.SetParameters(p2);

                C.MarginalDistributionX = m1;
                C.MarginalDistributionY = m2;
                return C.LogLikelihood(sampleDataX, sampleDataY);
            };

            var solver = new NelderMead(logLH, lowers.Length, initials, lowers, uppers);
            solver.Maximize();

            // Set parameters for copula and marginals
            copula.Theta = solver.BestParameterSet.Values[0];

            var par = new double[np1];
            for (int i = 0; i < np1; i++)
                par[i] = solver.BestParameterSet.Values[i + 1];
            copula.MarginalDistributionX.SetParameters(par);

            par = new double[np1];
            for (int i = 0; i < np2; i++)
                par[i] = solver.BestParameterSet.Values[i + 1 + np1];
            copula.MarginalDistributionY.SetParameters(par);

        }

    }
}
