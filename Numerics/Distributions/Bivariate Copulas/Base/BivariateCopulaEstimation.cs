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
using Numerics.Sampling.MCMC;
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
    /// <para>
    /// Supports copulas with one or more parameters. For single-parameter copulas, BrentSearch (1D) is used.
    /// For multi-parameter copulas (e.g., Student's t with rho and nu), NelderMead is used.
    /// Bayesian estimation uses the Adaptive Random Walk Metropolis-Hastings (ARWMH) MCMC sampler
    /// with uniform priors on the parameter constraints by default.
    /// </para>
    /// </remarks>
    [Serializable]
    public class BivariateCopulaEstimation
    {

        /// <summary>
        /// Estimate the bivariate copula.
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        /// <param name="estimationMethod">The estimation method to use.</param>
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
                case CopulaEstimationMethod.Bayesian:
                    BayesianMPL(copula, sampleDataX, sampleDataY);
                    break;
            }
        }

        /// <summary>
        /// Estimate the copula using Bayesian MCMC with optional custom priors.
        /// Sets the MAP estimate on the copula and returns the MCMC sampler for posterior analysis.
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        /// <param name="priorDistributions">Optional. The prior distributions for each copula parameter. If null, uniform priors on the parameter constraints are used.</param>
        /// <returns>The MCMC sampler with posterior samples accessible via <see cref="MCMCSampler.Output"/>.</returns>
        public static MCMCSampler EstimateBayesian(BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY, List<IUnivariateDistribution> priorDistributions = null)
        {
            return BayesianMPL(copula, sampleDataX, sampleDataY, priorDistributions);
        }

        /// <summary>
        /// The maximum pseudo likelihood method. Automatically selects BrentSearch for single-parameter
        /// copulas or NelderMead for multi-parameter copulas.
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable. When estimating with pseudo likelihood, this should be the plotting positions of the data.</param>
        /// <param name="sampleDataY">The sample data for the Y variable. When estimating with pseudo likelihood, this should be the plotting positions of the data.</param>
        private static void MPL(BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var constraints = copula.ParameterConstraints(sampleDataX, sampleDataY);
            int nParams = copula.NumberOfCopulaParameters;

            if (nParams == 1)
            {
                // Use BrentSearch for 1D optimization
                Func<double, double> func = (x) =>
                {
                    var C = copula.Clone();
                    C.SetCopulaParameters(new double[] { x });
                    return C.PseudoLogLikelihood(sampleDataX, sampleDataY);
                };
                var brent = new BrentSearch(func, constraints[0, 0], constraints[0, 1]);
                brent.Maximize();
                copula.SetCopulaParameters(new double[] { brent.BestParameterSet.Values[0] });
            }
            else
            {
                // Use NelderMead for multi-dimensional optimization
                var initials = copula.GetCopulaParameters;
                var lowers = new double[nParams];
                var uppers = new double[nParams];
                for (int i = 0; i < nParams; i++)
                {
                    lowers[i] = constraints[i, 0];
                    uppers[i] = constraints[i, 1];
                    // Clamp initials to be within bounds
                    initials[i] = Math.Max(lowers[i], Math.Min(uppers[i], initials[i]));
                }

                Func<double[], double> func = (x) =>
                {
                    var C = copula.Clone();
                    C.SetCopulaParameters(x);
                    return C.PseudoLogLikelihood(sampleDataX, sampleDataY);
                };

                var solver = new NelderMead(func, nParams, initials, lowers, uppers);
                solver.Maximize();
                copula.SetCopulaParameters(solver.BestParameterSet.Values);
            }
        }

        /// <summary>
        /// The inference from margins method. Automatically selects BrentSearch for single-parameter
        /// copulas or NelderMead for multi-parameter copulas.
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        private static void IFM(BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var constraints = copula.ParameterConstraints(sampleDataX, sampleDataY);
            int nParams = copula.NumberOfCopulaParameters;

            if (nParams == 1)
            {
                // Use BrentSearch for 1D optimization
                Func<double, double> func = (x) =>
                {
                    var C = copula.Clone();
                    C.SetCopulaParameters(new double[] { x });
                    return C.IFMLogLikelihood(sampleDataX, sampleDataY);
                };
                var brent = new BrentSearch(func, constraints[0, 0], constraints[0, 1]);
                brent.Maximize();
                copula.SetCopulaParameters(new double[] { brent.BestParameterSet.Values[0] });
            }
            else
            {
                // Use NelderMead for multi-dimensional optimization
                var initials = copula.GetCopulaParameters;
                var lowers = new double[nParams];
                var uppers = new double[nParams];
                for (int i = 0; i < nParams; i++)
                {
                    lowers[i] = constraints[i, 0];
                    uppers[i] = constraints[i, 1];
                    initials[i] = Math.Max(lowers[i], Math.Min(uppers[i], initials[i]));
                }

                Func<double[], double> func = (x) =>
                {
                    var C = copula.Clone();
                    C.SetCopulaParameters(x);
                    return C.IFMLogLikelihood(sampleDataX, sampleDataY);
                };

                var solver = new NelderMead(func, nParams, initials, lowers, uppers);
                solver.Maximize();
                copula.SetCopulaParameters(solver.BestParameterSet.Values);
            }
        }

        /// <summary>
        /// The maximum likelihood estimation method. Jointly estimates copula parameters and marginal
        /// distribution parameters using NelderMead optimization.
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        private static void MLE(BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY)
        {
            // See if marginals are estimable
            IMaximumLikelihoodEstimation margin1 = (IMaximumLikelihoodEstimation)copula.MarginalDistributionX;
            IMaximumLikelihoodEstimation margin2 = (IMaximumLikelihoodEstimation)copula.MarginalDistributionY;
            if (margin1 == null || margin2 == null) throw new ArgumentOutOfRangeException("marginal distributions", "The marginal distributions must implement the IMaximumLikelihoodEstimation interface to use this method.");

            int nCopula = copula.NumberOfCopulaParameters;
            int np1 = copula.MarginalDistributionX.NumberOfParameters;
            int np2 = copula.MarginalDistributionY.NumberOfParameters;
            int totalParams = nCopula + np1 + np2;

            var initials = new double[totalParams];
            var lowers = new double[totalParams];
            var uppers = new double[totalParams];

            // Get ranks and plotting positions for initial MPL estimate
            var rank1 = Statistics.RanksInPlace(sampleDataX.ToArray());
            var rank2 = Statistics.RanksInPlace(sampleDataY.ToArray());
            for (int i = 0; i < rank1.Length; i++)
            {
                rank1[i] = rank1[i] / (rank1.Length + 1d);
                rank2[i] = rank2[i] / (rank2.Length + 1d);
            }

            // Get copula parameter constraints and initial estimates via MPL
            var copulaConstraints = copula.ParameterConstraints(sampleDataX, sampleDataY);
            MPL(copula, rank1, rank2);
            var copulaParams = copula.GetCopulaParameters;
            for (int i = 0; i < nCopula; i++)
            {
                initials[i] = copulaParams[i];
                lowers[i] = copulaConstraints[i, 0];
                uppers[i] = copulaConstraints[i, 1];
            }

            // Estimate marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(sampleDataX, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(sampleDataY, ParameterEstimationMethod.MaximumLikelihood);

            var con = margin1.GetParameterConstraints(sampleDataX);
            var parms = copula.MarginalDistributionX.GetParameters;
            for (int i = 0; i < np1; i++)
            {
                initials[nCopula + i] = parms[i];
                lowers[nCopula + i] = con.Item2[i];
                uppers[nCopula + i] = con.Item3[i];
            }
            con = margin2.GetParameterConstraints(sampleDataY);
            parms = copula.MarginalDistributionY.GetParameters;
            for (int i = 0; i < np2; i++)
            {
                initials[nCopula + np1 + i] = parms[i];
                lowers[nCopula + np1 + i] = con.Item2[i];
                uppers[nCopula + np1 + i] = con.Item3[i];
            }

            // Log-likelihood function
            Func<double[], double> logLH = (double[] x) =>
            {
                // Set copula parameters
                var C = copula.Clone();
                var copulaVals = new double[nCopula];
                Array.Copy(x, 0, copulaVals, 0, nCopula);
                C.SetCopulaParameters(copulaVals);

                // Marginal 1
                var m1 = ((UnivariateDistributionBase)copula.MarginalDistributionX).Clone();
                var p1 = new double[np1];
                Array.Copy(x, nCopula, p1, 0, np1);
                m1.SetParameters(p1);

                // Marginal 2
                var m2 = ((UnivariateDistributionBase)copula.MarginalDistributionY).Clone();
                var p2 = new double[np2];
                Array.Copy(x, nCopula + np1, p2, 0, np2);
                m2.SetParameters(p2);

                C.MarginalDistributionX = m1;
                C.MarginalDistributionY = m2;
                return C.LogLikelihood(sampleDataX, sampleDataY);
            };

            var solver = new NelderMead(logLH, totalParams, initials, lowers, uppers);
            solver.Maximize();

            // Set copula parameters
            var bestCopula = new double[nCopula];
            Array.Copy(solver.BestParameterSet.Values, 0, bestCopula, 0, nCopula);
            copula.SetCopulaParameters(bestCopula);

            // Set marginal 1 parameters
            var par1 = new double[np1];
            Array.Copy(solver.BestParameterSet.Values, nCopula, par1, 0, np1);
            copula.MarginalDistributionX.SetParameters(par1);

            // Set marginal 2 parameters
            var par2 = new double[np2];
            Array.Copy(solver.BestParameterSet.Values, nCopula + np1, par2, 0, np2);
            copula.MarginalDistributionY.SetParameters(par2);
        }

        /// <summary>
        /// Bayesian estimation using MCMC with pseudo log-likelihood.
        /// Uses uniform priors on the parameter constraints by default.
        /// </summary>
        /// <param name="copula">The copula to estimate.</param>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        /// <param name="priors">Optional. The prior distributions for each copula parameter. If null, uniform priors are used.</param>
        /// <returns>The MCMC sampler with posterior samples.</returns>
        private static MCMCSampler BayesianMPL(BivariateCopula copula, IList<double> sampleDataX, IList<double> sampleDataY, List<IUnivariateDistribution> priors = null)
        {
            int nParams = copula.NumberOfCopulaParameters;
            var constraints = copula.ParameterConstraints(sampleDataX, sampleDataY);

            // Build default priors if not provided: Uniform over constraint bounds
            if (priors == null)
            {
                priors = new List<IUnivariateDistribution>();
                for (int i = 0; i < nParams; i++)
                {
                    priors.Add(new Uniform(constraints[i, 0], constraints[i, 1]));
                }
            }

            // First, get a good initial estimate using MPL (fast deterministic optimization)
            var mplCopula = copula.Clone();
            MPL(mplCopula, sampleDataX, sampleDataY);
            copula.SetCopulaParameters(mplCopula.GetCopulaParameters);

            // Get plotting positions for pseudo-likelihood
            var rank1 = Statistics.RanksInPlace(sampleDataX.ToArray());
            var rank2 = Statistics.RanksInPlace(sampleDataY.ToArray());
            for (int i = 0; i < rank1.Length; i++)
            {
                rank1[i] = rank1[i] / (rank1.Length + 1d);
                rank2[i] = rank2[i] / (rank2.Length + 1d);
            }

            // Log-likelihood function (pseudo-likelihood + prior)
            LogLikelihood logLH = (double[] parameters) =>
            {
                var C = copula.Clone();
                C.SetCopulaParameters(parameters);
                double ll = C.PseudoLogLikelihood(rank1, rank2);

                // Add prior log-likelihood
                for (int i = 0; i < nParams; i++)
                    ll += priors[i].LogPDF(parameters[i]);

                return ll;
            };

            // Use ARWMH sampler (self-tuning, robust for low-dimensional problems)
            // Initialize near the MPL estimate rather than running expensive MAP optimization
            var sampler = new ARWMH(priors, logLH);
            sampler.Initialize = MCMCSampler.InitializationType.Randomize;
            sampler.NumberOfChains = 1;
            sampler.Sample();

            // Set the best estimate on the copula
            // Use the MPL estimate if MCMC didn't find something better
            var mplParams = mplCopula.GetCopulaParameters;
            double mplLogLH = logLH(mplParams);
            if (sampler.MAP.Values != null && sampler.MAP.Fitness >= mplLogLH)
            {
                copula.SetCopulaParameters(sampler.MAP.Values);
            }
            else
            {
                copula.SetCopulaParameters(mplParams);
            }

            return sampler;
        }

    }
}
