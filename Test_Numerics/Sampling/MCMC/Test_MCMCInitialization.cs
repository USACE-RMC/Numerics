using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Regression tests for MCMC initialization behavior.
    /// </summary>
    [TestClass]
    public class Test_MCMCInitialization
    {
        /// <summary>
        /// Verifies that MAP initialization stores the original log-likelihood sign without an extra evaluation.
        /// </summary>
        [TestMethod]
        public void MapInitializationStoresLogLikelihoodWithoutReevaluation()
        {
            int optimizerEvaluations = 0;
            var optimizer = new DifferentialEvolution(
                parameters =>
                {
                    optimizerEvaluations++;
                    return QuadraticLogLikelihood(parameters);
                },
                1,
                new[] { -5d },
                new[] { 5d })
            {
                ReportFailure = false
            };
            optimizer.Maximize();
            Assert.AreEqual(OptimizationStatus.Success, optimizer.Status);

            int samplerEvaluations = 0;
            var sampler = new InitializableRwmh(
                new List<IUnivariateDistribution> { new Uniform(-5d, 5d) },
                parameters =>
                {
                    samplerEvaluations++;
                    return QuadraticLogLikelihood(parameters);
                })
            {
                NumberOfChains = 1,
                InitialIterations = 2,
                Initialize = MCMCSampler.InitializationType.MAP
            };

            sampler.InitializeOnly();

            Assert.IsFalse(sampler.MAPInitializationFailed);
            Assert.AreEqual(optimizerEvaluations + sampler.InitialIterations, samplerEvaluations);
            Assert.AreEqual(QuadraticLogLikelihood(sampler.MAP.Values), sampler.MAP.Fitness, 1E-12);
            Assert.AreEqual(-3d, sampler.MAP.Fitness, 1E-6);
        }

        /// <summary>
        /// Verifies that NUTS uses an analytical gradient during its step-size initialization heuristic.
        /// </summary>
        [TestMethod]
        public void NutsInitializationUsesConfiguredGradientAndReducesLikelihoodWork()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(-5d, 5d) };
            int analyticalLikelihoodEvaluations = 0;
            int analyticalGradientEvaluations = 0;
            var analytical = new InitializableNuts(
                priors,
                parameters =>
                {
                    analyticalLikelihoodEvaluations++;
                    return -parameters[0] * parameters[0];
                },
                parameters =>
                {
                    analyticalGradientEvaluations++;
                    return new Vector(new[] { -2d * parameters[0] });
                });
            ConfigureSingleChainInitialization(analytical);
            analytical.InitializeOnly();

            int numericalLikelihoodEvaluations = 0;
            var numerical = new InitializableNuts(
                priors,
                parameters =>
                {
                    numericalLikelihoodEvaluations++;
                    return -parameters[0] * parameters[0];
                });
            ConfigureSingleChainInitialization(numerical);
            numerical.InitializeOnly();

            Assert.IsGreaterThan(0, analyticalGradientEvaluations);
            Assert.IsGreaterThan(analyticalLikelihoodEvaluations, numericalLikelihoodEvaluations);
        }

        /// <summary>
        /// Verifies that randomized initialization resolves log-likelihood ties by draw order.
        /// </summary>
        /// <remarks>
        /// A wide prior can leave most of the initial population at exactly negative infinity. The chain
        /// starting states are taken from the front of the sorted population, so under an unstable sort
        /// the tied draws that seed the chains would be implementation-defined, and with them the entire
        /// trajectory of every chain.
        /// </remarks>
        [TestMethod]
        public void RandomizedInitializationBreaksLogLikelihoodTiesByDrawOrder()
        {
            var evaluationOrder = new List<double[]>();
            var sampler = new InitializableRwmh(
                new List<IUnivariateDistribution> { new Uniform(-5d, 5d) },
                parameters =>
                {
                    evaluationOrder.Add((double[])parameters.Clone());

                    // Only the twenty-first draw has finite support. The remaining thirty-nine draws tie
                    // at exactly negative infinity, which is more than the insertion-sort threshold of
                    // the framework sort, so an unstable sort is free to permute them.
                    return evaluationOrder.Count == 21 ? -1d : double.NegativeInfinity;
                })
            {
                NumberOfChains = 4,
                InitialIterations = 40,
                PRNGSeed = 12345,
                Initialize = MCMCSampler.InitializationType.Randomize
            };

            var initials = sampler.InitializeOnly();

            Assert.HasCount(40, evaluationOrder);
            Assert.HasCount(4, initials);

            // The single finite draw sorts to the front.
            CollectionAssert.AreEqual(evaluationOrder[20], initials[0].Values);

            // The tied draws follow in the order they were generated.
            CollectionAssert.AreEqual(evaluationOrder[0], initials[1].Values);
            CollectionAssert.AreEqual(evaluationOrder[1], initials[2].Values);
            CollectionAssert.AreEqual(evaluationOrder[2], initials[3].Values);
        }

        /// <summary>
        /// Evaluates the quadratic log-likelihood used to verify MAP fitness and evaluation counts.
        /// </summary>
        /// <param name="parameters">The parameter vector to evaluate.</param>
        /// <returns>The quadratic log-likelihood centered at one with a maximum of negative three.</returns>
        private static double QuadraticLogLikelihood(double[] parameters)
        {
            double difference = parameters[0] - 1d;
            return -3d - difference * difference;
        }

        /// <summary>
        /// Applies deterministic, minimal settings to an initialization-only sampler.
        /// </summary>
        /// <param name="sampler">The sampler to configure.</param>
        private static void ConfigureSingleChainInitialization(MCMCSampler sampler)
        {
            sampler.NumberOfChains = 1;
            sampler.InitialIterations = 1;
            sampler.PRNGSeed = 12345;
        }

        /// <summary>
        /// Exposes random-walk Metropolis-Hastings initialization for focused regression testing.
        /// </summary>
        private sealed class InitializableRwmh : RWMH
        {
            /// <summary>
            /// Initializes a new instance of the <see cref="InitializableRwmh"/> class.
            /// </summary>
            /// <param name="priorDistributions">The prior distributions for the sampled parameters.</param>
            /// <param name="logLikelihood">The log-likelihood function to evaluate.</param>
            internal InitializableRwmh(
                List<IUnivariateDistribution> priorDistributions,
                LogLikelihood logLikelihood)
                : base(priorDistributions, logLikelihood, new Matrix(1))
            {
            }

            /// <summary>
            /// Initializes the chains without running the sampler.
            /// </summary>
            /// <returns>The initialized chain parameter sets.</returns>
            internal ParameterSet[] InitializeOnly()
            {
                return InitializeChains();
            }
        }

        /// <summary>
        /// Exposes No-U-Turn sampler initialization for focused gradient-routing tests.
        /// </summary>
        private sealed class InitializableNuts : NUTS
        {
            /// <summary>
            /// Initializes a new instance of the <see cref="InitializableNuts"/> class.
            /// </summary>
            /// <param name="priorDistributions">The prior distributions for the sampled parameters.</param>
            /// <param name="logLikelihood">The log-likelihood function to evaluate.</param>
            /// <param name="gradient">The optional analytical gradient of the log-likelihood.</param>
            internal InitializableNuts(
                List<IUnivariateDistribution> priorDistributions,
                LogLikelihood logLikelihood,
                HMC.Gradient gradient = null)
                : base(priorDistributions, logLikelihood, maxTreeDepth: 2, gradientFunction: gradient)
            {
            }

            /// <summary>
            /// Initializes the chains and sampler-specific step-size state without drawing samples.
            /// </summary>
            internal void InitializeOnly()
            {
                _chainStates = InitializeChains();
                InitializeCustomSettings();
            }
        }
    }
}
