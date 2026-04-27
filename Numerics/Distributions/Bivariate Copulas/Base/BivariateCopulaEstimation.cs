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
            }
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
            IMaximumLikelihoodEstimation? margin1 = copula.MarginalDistributionX as IMaximumLikelihoodEstimation;
            IMaximumLikelihoodEstimation? margin2 = copula.MarginalDistributionY as IMaximumLikelihoodEstimation;
            if (margin1 == null || margin2 == null) throw new ArgumentOutOfRangeException("marginal distributions", "The marginal distributions must implement the IMaximumLikelihoodEstimation interface to use this method.");

            int nCopula = copula.NumberOfCopulaParameters;
            int np1 = copula.MarginalDistributionX!.NumberOfParameters;
            int np2 = copula.MarginalDistributionY!.NumberOfParameters;
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
            ((IEstimation)copula.MarginalDistributionX!).Estimate(sampleDataX, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY!).Estimate(sampleDataY, ParameterEstimationMethod.MaximumLikelihood);

            var con = margin1.GetParameterConstraints(sampleDataX);
            var parms = copula.MarginalDistributionX!.GetParameters;
            for (int i = 0; i < np1; i++)
            {
                initials[nCopula + i] = parms[i];
                lowers[nCopula + i] = con.Item2[i];
                uppers[nCopula + i] = con.Item3[i];
            }
            con = margin2.GetParameterConstraints(sampleDataY);
            parms = copula.MarginalDistributionY!.GetParameters;
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
                var m1 = ((UnivariateDistributionBase)copula.MarginalDistributionX!).Clone();
                var p1 = new double[np1];
                Array.Copy(x, nCopula, p1, 0, np1);
                m1.SetParameters(p1);

                // Marginal 2
                var m2 = ((UnivariateDistributionBase)copula.MarginalDistributionY!).Clone();
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
            copula.MarginalDistributionX!.SetParameters(par1);

            // Set marginal 2 parameters
            var par2 = new double[np2];
            Array.Copy(solver.BestParameterSet.Values, nCopula + np1, par2, 0, np2);
            copula.MarginalDistributionY!.SetParameters(par2);
        }

    }
}
