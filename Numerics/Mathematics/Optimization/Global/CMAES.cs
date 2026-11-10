using Numerics.Distributions;
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The Covariance Matrix Adaptation Evolution Strategy (CMA-ES) with IPOP restarts, which finds a global minimum when no gradient is available.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Each generation samples lambda points x_k = m + sigma B D z_k with z_k ~ N(0, I), where C = B D^2 B^T. The mean moves to
    /// the weighted average of the mu best points. The covariance matrix C learns from the evolution path p_c (rank-one update)
    /// and from all lambda ranked steps (rank-mu update, with the negative weights of active CMA). The step size sigma follows
    /// cumulative step-size adaptation on the conjugate path p_sigma. All strategy parameters are the defaults of Table 1 in
    /// Hansen (2023), including the h_sigma stall of p_c, and B and D are recomputed lazily, after more than
    /// lambda / (c1 + c_mu) / n / 10 evaluations, as in the tutorial's reference code.
    /// </para>
    /// <para>
    /// The search runs in the unit cube: each parameter is scaled by its bounds, the default mean is the box center, and the
    /// default initial step size is 0.3. Bounds are handled as in the pycma BoundTransform: the strategy samples an unbounded
    /// genotype y, and the objective is evaluated at T(y), where T is the periodic linear-quadratic map of pycma's
    /// BoxConstraintsLinQuadTransformation. T is the identity on [0.05, 0.95], quadratic within 0.05 of each bound, and maps
    /// every real number into [0, 1], so every evaluated point lies inside the bounds and the CMA-ES update is unchanged.
    /// As in pycma (option maxstd_boundrange = 1/3), each coordinate standard deviation sigma sqrt(C_ii) is capped at one
    /// third of the bound range by rescaling row and column i of C and component i of p_c. Without the cap the periodic
    /// map makes the fitness insensitive to very long steps, sigma can drift far beyond the box, and the run degrades to
    /// a random search.
    /// </para>
    /// <para>
    /// A run stops on the termination criteria of Hansen (2023, Sect. B.3): TolFun, TolX, EqualFunValues, NoEffectAxis,
    /// NoEffectCoord, ConditionCov, Stagnation and TolXUp. After each run the IPOP scheme of Auger and Hansen (2005) restarts
    /// from a uniformly random mean with the population size doubled, up to <see cref="Restarts"/> times.
    /// </para>
    /// <para>
    ///     <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Hansen, N. (2023). The CMA Evolution Strategy: A Tutorial. arXiv:1604.00772v2. <see href="https://arxiv.org/abs/1604.00772"/>
    /// </description></item>
    /// <item><description>
    /// Hansen, N., and Ostermeier, A. (2001). Completely derandomized self-adaptation in evolution strategies.
    /// Evolutionary Computation, 9(2), 159-195. <see href="https://doi.org/10.1162/106365601750190398"/>
    /// </description></item>
    /// <item><description>
    /// Auger, A., and Hansen, N. (2005). A restart CMA evolution strategy with increasing population size.
    /// 2005 IEEE Congress on Evolutionary Computation, 2, 1769-1776. <see href="https://doi.org/10.1109/CEC.2005.1554902"/>
    /// </description></item>
    /// <item><description>
    /// Hansen, N., Akimoto, Y., Baudis, P., et al. CMA-ES/pycma, cma.transformations.BoxConstraintsLinQuadTransformation.
    /// Zenodo. <see href="https://doi.org/10.5281/zenodo.2559634"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class CMAES : Optimizer
    {

        /// <summary>
        /// Construct a new CMA-ES optimization method that starts from the center of the bounds.
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        /// <exception cref="ArgumentOutOfRangeException">The bounds have the wrong length, or an upper bound is not greater than its lower bound.</exception>
        public CMAES(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
        {
            // Check if the length of the lower and upper bounds equal the number of parameters
            if (lowerBounds.Count != numberOfParameters || upperBounds.Count != numberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(lowerBounds), "The lower and upper bounds must be the same length as the number of parameters.");
            }
            // Check if the lower bounds are less than the upper bounds
            for (int i = 0; i < lowerBounds.Count; i++)
            {
                if (upperBounds[i] <= lowerBounds[i])
                {
                    throw new ArgumentOutOfRangeException(nameof(upperBounds), "The upper bound cannot be less than or equal to the lower bound.");
                }
            }
            LowerBounds = lowerBounds.ToArray();
            UpperBounds = upperBounds.ToArray();
            // The default evaluation budget bounds the restart sequence (Auger and Hansen, 2005).
            MaxFunctionEvaluations = numberOfParameters > int.MaxValue / 10000 ? int.MaxValue : 10000 * numberOfParameters;
            PopulationSize = 4 + (int)Math.Floor(3d * Math.Log(NumberOfParameters));
        }

        /// <summary>
        /// Construct a new CMA-ES optimization method that starts from the given point.
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">The initial mean of the first run. Restarts draw a new mean uniformly within the bounds.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        /// <exception cref="ArgumentOutOfRangeException">An array has the wrong length, an upper bound is not greater than its lower bound, or an initial value is outside the bounds.</exception>
        public CMAES(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds) : this(objectiveFunction, numberOfParameters, lowerBounds, upperBounds)
        {
            if (initialValues.Count != numberOfParameters)
                throw new ArgumentOutOfRangeException(nameof(initialValues), "The initial values must be the same length as the number of parameters.");
            for (int i = 0; i < numberOfParameters; i++)
            {
                if (initialValues[i] < lowerBounds[i] || initialValues[i] > upperBounds[i])
                    throw new ArgumentOutOfRangeException(nameof(initialValues), "The initial values must be between the lower and upper bounds.");
            }
            InitialValues = initialValues.ToArray();
        }

        /// <summary>
        /// An array of lower bounds (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double[] LowerBounds { get; private set; }

        /// <summary>
        /// An array of upper bounds (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double[] UpperBounds { get; private set; }

        /// <summary>
        /// The initial mean of the first run, or null to start from the center of the bounds.
        /// </summary>
        public double[]? InitialValues { get; private set; }

        /// <inheritdoc />
        protected override double[]? ParameterLowerBounds => LowerBounds;

        /// <inheritdoc />
        protected override double[]? ParameterUpperBounds => UpperBounds;

        /// <summary>
        /// The population size lambda of the first run, at least 2. Default = 4 + floor(3 ln n) (Hansen, 2023, Table 1).
        /// Each IPOP restart doubles it.
        /// </summary>
        public int PopulationSize { get; set; }

        /// <summary>
        /// The initial step size sigma in the unit-cube coordinates, where each parameter range has length one. Default = 0.3
        /// (Hansen, 2023, Fig. 6). Restarts use the same value.
        /// </summary>
        public double InitialStepSize { get; set; } = 0.3;

        /// <summary>
        /// The maximum number of IPOP restarts after the first run. Default = 9, as in pycma's IPOP setting.
        /// Set it to 0 for a single CMA-ES run.
        /// </summary>
        /// <remarks>
        /// The restarts share the evaluation budget, <see cref="Optimizer.MaxFunctionEvaluations"/>, which CMA-ES
        /// sets to 10,000 times the number of parameters by default (Auger and Hansen, 2005). Once a run has
        /// converged, reaching the budget during a later restart ends the search with
        /// <see cref="OptimizationStatus.Success"/>. Raise the budget for more restarts on hard multimodal
        /// problems, or lower it for expensive objectives.
        /// </remarks>
        public int Restarts { get; set; } = 9;

        /// <summary>
        /// The number of restarts performed in the last optimization.
        /// </summary>
        public int RestartsPerformed { get; private set; }

        /// <summary>
        /// The termination criterion, population size and evaluations of each run in the last optimization, as
        /// "criterion|lambda|evaluations". A run cut by the evaluation or iteration limit is listed as Budget.
        /// </summary>
        internal List<string> RunTerminations { get; } = new List<string>();

        /// <summary>
        /// Determines whether the best point is refined by a local search when the restarts end. Default = false.
        /// </summary>
        /// <remarks>
        /// The polish runs after the last restart when at least one run has converged, never after the
        /// evaluation budget is spent. Its evaluations count toward <see cref="Optimizer.FunctionEvaluations"/>
        /// and the reported status is the global search's.
        /// </remarks>
        public bool Polish { get; set; }

        /// <summary>
        /// The local method used to polish: BFGS, Nelder-Mead, or Powell. Default = Nelder-Mead.
        /// </summary>
        public LocalMethod PolishMethod { get; set; } = LocalMethod.NelderMead;

        /// <summary>
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// The largest coordinate standard deviation sigma sqrt(C_ii) in unit-cube coordinates: one third of the bound
        /// range, pycma's maxstd_boundrange default.
        /// </summary>
        private const double MaxStd = 1d / 3d;

        /// <summary>
        /// True once a run has ended on a convergence criterion; the budget limits then end the restarts normally.
        /// </summary>
        private bool _converged;

        /// <summary>
        /// The reason a single CMA-ES run stopped.
        /// </summary>
        private enum Termination
        {
            TolFun,
            TolX,
            EqualFunValues,
            NoEffectAxis,
            NoEffectCoord,
            ConditionCov,
            Stagnation,
            TolXUp,
            Budget
        }

        /// <summary>
        /// Orders sample indexes by fitness, with ties broken by index so every sort algorithm gives the same order.
        /// </summary>
        private sealed class FitnessComparer : IComparer<int>
        {
            private readonly double[] _fitness;

            /// <summary>
            /// Creates a comparer over the given fitness array.
            /// </summary>
            /// <param name="fitness">The fitness of each sample.</param>
            public FitnessComparer(double[] fitness) { _fitness = fitness; }

            /// <inheritdoc />
            public int Compare(int a, int b)
            {
                double fa = _fitness[a], fb = _fitness[b];
                if (fa < fb) return -1;
                if (fa > fb) return 1;
                return a.CompareTo(b);
            }
        }

        /// <inheritdoc />
        /// <remarks>
        /// After a run has converged, reaching the iteration or evaluation limit during a later restart ends the
        /// optimization with <see cref="OptimizationStatus.Success"/> and the best point found so far.
        /// </remarks>
        protected override void UpdateStatus(OptimizationStatus status, Exception? exception = null)
        {
            if (_converged && (status == OptimizationStatus.MaximumIterationsReached || status == OptimizationStatus.MaximumFunctionEvaluationsReached))
                status = OptimizationStatus.Success;
            base.UpdateStatus(status, exception);
        }

        /// <inheritdoc/>
        protected override void Optimize()
        {
            if (PopulationSize < 2) throw new ArgumentOutOfRangeException(nameof(PopulationSize), "The population size must be at least 2.");
            if (!(InitialStepSize > 0) || double.IsInfinity(InitialStepSize)) throw new ArgumentOutOfRangeException(nameof(InitialStepSize), "The initial step size must be positive and finite.");
            if (Restarts < 0) throw new ArgumentOutOfRangeException(nameof(Restarts), "The number of restarts cannot be negative.");
            if (Polish) ValidatePolishMethod(PolishMethod);

            int D = NumberOfParameters;
            _converged = false;
            RestartsPerformed = 0;
            RunTerminations.Clear();
            var prng = new MersenneTwister(PRNGSeed);

            // Initial genotype mean in the unit cube
            var mean = new double[D];
            for (int j = 0; j < D; j++)
                mean[j] = InitialValues == null ? 0.5 : InverseTransform((InitialValues[j] - LowerBounds[j]) / (UpperBounds[j] - LowerBounds[j]));

            int lambda = PopulationSize;
            for (int run = 0; ; run++)
            {
                int fe0 = FunctionEvaluations;
                var reason = Run(lambda, mean, prng);
                RunTerminations.Add(reason.ToString() + "|" + lambda + "|" + (FunctionEvaluations - fe0));
                if (reason == Termination.Budget) return;
                if (reason != Termination.TolXUp) _converged = true;
                if (run >= Restarts) break;

                // IPOP restart: double the population and draw a new mean uniformly in the box
                RestartsPerformed += 1;
                if (lambda > int.MaxValue / 2) break;
                lambda *= 2;
                for (int j = 0; j < D; j++)
                    mean[j] = InverseTransform(prng.NextDouble());
            }

            if (_converged && Polish && !PolishBestParameterSet(PolishMethod)) return;
            if (_converged)
                UpdateStatus(OptimizationStatus.Success);
            else
                UpdateStatus(OptimizationStatus.Failure, new InvalidOperationException("Every CMA-ES run ended on TolXUp: the step size diverged."));
        }

        /// <summary>
        /// Performs one CMA-ES run from the given genotype mean.
        /// </summary>
        /// <param name="lambda">The population size.</param>
        /// <param name="initialMean">The initial genotype mean in unit-cube coordinates.</param>
        /// <param name="prng">The random number generator.</param>
        /// <returns>The termination criterion that stopped the run.</returns>
        /// <remarks>
        /// Implements Figure 6 and Table 1 of Hansen (2023), with the termination criteria of Sect. B.3.
        /// The exponent of the h_sigma correction is 2g with g the number of completed generations, as in the tutorial's code.
        /// </remarks>
        private Termination Run(int lambda, double[] initialMean, MersenneTwister prng)
        {
            int n = NumberOfParameters, i, j, k, l;
            bool cancel = false;

            // Strategy parameters (Table 1)
            var w = new double[lambda];
            double sumPos = 0, sumPos2 = 0, sumNeg = 0, sumNeg2 = 0;
            for (i = 0; i < lambda; i++)
            {
                w[i] = Math.Log((lambda + 1) / 2d) - Math.Log(i + 1);
                if (w[i] > 0) { sumPos += w[i]; sumPos2 += w[i] * w[i]; }
                else if (w[i] < 0) { sumNeg -= w[i]; sumNeg2 += w[i] * w[i]; }
            }
            int mu = lambda / 2;
            double mueff = sumPos * sumPos / sumPos2;
            double mueffNeg = sumNeg2 > 0 ? sumNeg * sumNeg / sumNeg2 : 0;
            double c1 = 2d / ((n + 1.3) * (n + 1.3) + mueff);
            double cmu = Math.Min(1 - c1, 2d * (0.25 + mueff + 1 / mueff - 2) / ((n + 2d) * (n + 2d) + mueff));
            double alphaMu = 1 + c1 / cmu;
            double alphaMueff = 1 + 2 * mueffNeg / (mueff + 2);
            double alphaPosDef = (1 - c1 - cmu) / (n * cmu);
            double negScale = Math.Min(alphaMu, Math.Min(alphaMueff, alphaPosDef));
            double sumW = 0;
            for (i = 0; i < lambda; i++)
            {
                w[i] = w[i] >= 0 ? w[i] / sumPos : negScale * w[i] / sumNeg;
                sumW += w[i];
            }
            double cs = (mueff + 2) / (n + mueff + 5);
            double ds = 1 + 2 * Math.Max(0, Math.Sqrt((mueff - 1) / (n + 1)) - 1) + cs;
            double cc = (4 + mueff / n) / (n + 4 + 2 * mueff / n);
            double chiN = Math.Sqrt(n) * (1 - 1d / (4 * n) + 1d / (21d * n * n));
            double hsThreshold = (1.4 + 2d / (n + 1)) * chiN;
            double sqrtCsMueff = Math.Sqrt(cs * (2 - cs) * mueff);
            double sqrtCcMueff = Math.Sqrt(cc * (2 - cc) * mueff);
            double lazyGap = lambda / (c1 + cmu) / n / 10d;

            // Termination settings (Sect. B.3)
            double sigma0 = InitialStepSize;
            double tolX = 1E-12 * sigma0;
            int histLength = 10 + (int)Math.Ceiling(30d * n / lambda);
            double stagMin = 120 + 30d * n / lambda;

            // State
            var m = (double[])initialMean.Clone();
            double sigma = sigma0;
            var C = new double[n, n];
            var B = new double[n, n];
            var Dd = new double[n];
            for (j = 0; j < n; j++) { C[j, j] = 1; B[j, j] = 1; Dd[j] = 1; }
            var ps = new double[n];
            var pc = new double[n];
            double csDecay = 1; // (1 - cs)^(2g), kept as a product

            // Sample buffers
            var Z = new double[lambda][];
            var Y = new double[lambda][];
            for (k = 0; k < lambda; k++) { Z[k] = new double[n]; Y[k] = new double[n]; }
            var x = new double[n];
            var fit = new double[lambda];
            var idx = new int[lambda];
            var comparer = new FitnessComparer(fit);
            var bdz = new double[n];
            var yw = new double[n];
            var zw = new double[n];
            var wc = new double[lambda];
            var scale = new double[n];
            var bestHist = new double[histLength];
            var stagBest = new List<double>();
            var stagMedian = new List<double>();
            double[] scratch = new double[16];
            int gen = 0;
            int eigenEval = FunctionEvaluations;

            while (true)
            {
                if (Iterations >= MaxIterations)
                {
                    UpdateStatus(OptimizationStatus.MaximumIterationsReached);
                    return Termination.Budget;
                }

                // Sample and evaluate the population (Eq. 38-40)
                for (k = 0; k < lambda; k++)
                {
                    var z = Z[k];
                    var y = Y[k];
                    for (j = 0; j < n; j++)
                    {
                        z[j] = Normal.StandardZ(prng.NextDouble());
                        bdz[j] = Dd[j] * z[j];
                    }
                    for (j = 0; j < n; j++)
                    {
                        double s = 0;
                        for (l = 0; l < n; l++) s += B[j, l] * bdz[l];
                        y[j] = s;
                        x[j] = LowerBounds[j] + (UpperBounds[j] - LowerBounds[j]) * Transform(m[j] + sigma * s);
                        if (x[j] < LowerBounds[j]) x[j] = LowerBounds[j];
                        else if (x[j] > UpperBounds[j]) x[j] = UpperBounds[j];
                    }
                    fit[k] = Evaluate(x, ref cancel);
                    if (cancel) return Termination.Budget;
                }
                gen += 1;
                Iterations += 1;

                for (k = 0; k < lambda; k++) idx[k] = k;
                Array.Sort(idx, comparer);

                // Recombination (Eq. 41-42)
                Array.Clear(yw, 0, n);
                Array.Clear(zw, 0, n);
                for (i = 0; i < mu; i++)
                {
                    var y = Y[idx[i]];
                    var z = Z[idx[i]];
                    for (j = 0; j < n; j++) { yw[j] += w[i] * y[j]; zw[j] += w[i] * z[j]; }
                }
                for (j = 0; j < n; j++) m[j] += sigma * yw[j];

                // Step-size path (Eq. 43), using C^(-1/2) <y>_w = B <z>_w
                double normPs2 = 0;
                for (j = 0; j < n; j++)
                {
                    double s = 0;
                    for (l = 0; l < n; l++) s += B[j, l] * zw[l];
                    ps[j] = (1 - cs) * ps[j] + sqrtCsMueff * s;
                    normPs2 += ps[j] * ps[j];
                }
                double normPs = Math.Sqrt(normPs2);
                csDecay *= (1 - cs) * (1 - cs);
                bool hsig = normPs / Math.Sqrt(1 - csDecay) < hsThreshold;

                // Covariance path (Eq. 45)
                for (j = 0; j < n; j++)
                    pc[j] = (1 - cc) * pc[j] + (hsig ? sqrtCcMueff * yw[j] : 0);

                // Covariance matrix update (Eq. 46-47); ||C^(-1/2) y|| = ||z|| for the sampled decomposition
                for (i = 0; i < lambda; i++)
                {
                    wc[i] = w[i];
                    if (w[i] < 0)
                    {
                        var z = Z[idx[i]];
                        double z2 = 0;
                        for (j = 0; j < n; j++) z2 += z[j] * z[j];
                        wc[i] = w[i] * n / z2;
                    }
                }
                double delta = hsig ? 0 : cc * (2 - cc);
                double decay = 1 + c1 * delta - c1 - cmu * sumW;
                for (j = 0; j < n; j++)
                {
                    for (l = j; l < n; l++)
                    {
                        double s = 0;
                        for (i = 0; i < lambda; i++)
                        {
                            var y = Y[idx[i]];
                            s += wc[i] * y[j] * y[l];
                        }
                        double v = decay * C[j, l] + c1 * pc[j] * pc[l] + cmu * s;
                        C[j, l] = v;
                        C[l, j] = v;
                    }
                }

                // Step-size update (Eq. 44)
                sigma *= Math.Exp(cs / ds * (normPs / chiN - 1));

                // Cap the coordinate standard deviations at one third of the bound range (pycma maxstd_boundrange)
                bool capped = false;
                for (j = 0; j < n; j++)
                {
                    double sd = sigma * Math.Sqrt(C[j, j]);
                    scale[j] = sd > MaxStd ? MaxStd / sd : 1;
                    capped |= sd > MaxStd;
                }
                if (capped)
                {
                    for (j = 0; j < n; j++)
                    {
                        pc[j] *= scale[j];
                        for (l = 0; l < n; l++) C[j, l] *= scale[j] * scale[l];
                    }
                }

                // Lazy eigendecomposition (Sect. B.2), and at once after a cap so that sampling uses the capped C
                bool conditionStop = false;
                if (capped || FunctionEvaluations - eigenEval > lazyGap)
                {
                    eigenEval = FunctionEvaluations;
                    conditionStop = !Decompose(C, B, Dd);
                }

                // Termination criteria (Sect. B.3)
                double fBest = fit[idx[0]], fMedian = fit[idx[lambda / 2]];
                bestHist[(gen - 1) % histLength] = fBest;
                stagBest.Add(fBest);
                stagMedian.Add(fMedian);

                if (gen >= histLength)
                {
                    double hMin = double.PositiveInfinity, hMax = double.NegativeInfinity;
                    for (i = 0; i < histLength; i++) { hMin = Math.Min(hMin, bestHist[i]); hMax = Math.Max(hMax, bestHist[i]); }
                    if (hMax - hMin == 0) return Termination.EqualFunValues;
                    hMin = Math.Min(hMin, fit[idx[0]]);
                    hMax = Math.Max(hMax, fit[idx[lambda - 1]]);
                    if (hMax - hMin < AbsoluteTolerance + RelativeTolerance * Math.Abs(fBest)) return Termination.TolFun;
                }

                bool smallX = true;
                for (j = 0; j < n && smallX; j++)
                    smallX = sigma * Math.Sqrt(C[j, j]) < tolX && sigma * Math.Abs(pc[j]) < tolX;
                if (smallX) return Termination.TolX;

                int axis = (gen - 1) % n;
                bool noEffect = true;
                for (j = 0; j < n && noEffect; j++)
                    noEffect = m[j] == m[j] + 0.1 * sigma * Dd[axis] * B[j, axis];
                if (noEffect) return Termination.NoEffectAxis;

                for (j = 0; j < n; j++)
                    if (m[j] == m[j] + 0.2 * sigma * Math.Sqrt(C[j, j])) return Termination.NoEffectCoord;

                double dMin = double.PositiveInfinity, dMax = 0;
                for (j = 0; j < n; j++) { dMin = Math.Min(dMin, Dd[j]); dMax = Math.Max(dMax, Dd[j]); }
                if (conditionStop || dMax * dMax > 1E14 * dMin * dMin) return Termination.ConditionCov;

                if (stagBest.Count >= stagMin)
                {
                    // The 20000 cap is applied last so that it also wins when stagMin exceeds it, which
                    // Tools.Clamp(x, stagMin, 20000) would not do.
                    int window = (int)Math.Min(20000, Math.Max(stagMin, 0.2 * stagBest.Count));
                    int part = Math.Max(1, (int)(0.3 * window));
                    if (scratch.Length < part) scratch = new double[2 * part];
                    int start = stagBest.Count - window;
                    stagBest.CopyTo(stagBest.Count - part, scratch, 0, part);
                    double recentBest = Statistics.MedianInPlace(scratch, part);
                    stagBest.CopyTo(start, scratch, 0, part);
                    if (!(recentBest < Statistics.MedianInPlace(scratch, part)))
                    {
                        stagMedian.CopyTo(stagMedian.Count - part, scratch, 0, part);
                        double recentMedian = Statistics.MedianInPlace(scratch, part);
                        stagMedian.CopyTo(start, scratch, 0, part);
                        if (!(recentMedian < Statistics.MedianInPlace(scratch, part)))
                            return Termination.Stagnation;
                    }
                }

                if (sigma * dMax > 1E4 * sigma0) return Termination.TolXUp;
            }
        }

        /// <summary>
        /// Updates B and D from C by decomposing B^T C B, which is nearly diagonal between lazy updates.
        /// </summary>
        /// <param name="C">The covariance matrix, which is not changed.</param>
        /// <param name="B">On input the previous eigenvectors; on output the new eigenvectors.</param>
        /// <param name="D">On output the square roots of the eigenvalues.</param>
        /// <returns>False if an eigenvalue is not positive, which ends the run as ConditionCov.</returns>
        /// <remarks>
        /// If B^T C B = V L V^T, then C = (B V) L (B V)^T, so the result is the eigensystem of C. Starting from the previous
        /// basis keeps the Jacobi iteration of <see cref="EigenValueDecomposition"/> short.
        /// </remarks>
        private static bool Decompose(double[,] C, double[,] B, double[] D)
        {
            int n = D.Length, i, j, k;
            var CB = new double[n, n];
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                {
                    double s = 0;
                    for (k = 0; k < n; k++) s += C[i, k] * B[k, j];
                    CB[i, j] = s;
                }
            var M = new double[n, n];
            for (i = 0; i < n; i++)
                for (j = i; j < n; j++)
                {
                    double s = 0;
                    for (k = 0; k < n; k++) s += B[k, i] * CB[k, j];
                    M[i, j] = s;
                    M[j, i] = s;
                }
            var evd = new EigenValueDecomposition(new Matrix(M));
            var V = evd.EigenVectors.Array;
            for (i = 0; i < n; i++)
                if (!(evd.EigenValues[i] > 0)) return false;
            var Bn = new double[n, n];
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                {
                    double s = 0;
                    for (k = 0; k < n; k++) s += B[i, k] * V[k, j];
                    Bn[i, j] = s;
                }
            Array.Copy(Bn, B, Bn.Length);
            for (i = 0; i < n; i++) D[i] = Math.Sqrt(evd.EigenValues[i]);
            return true;
        }

        /// <summary>
        /// Maps a genotype coordinate into [0, 1] with pycma's periodic linear-quadratic transformation.
        /// </summary>
        /// <param name="y">The genotype coordinate in unit-cube units.</param>
        /// <returns>The feasible coordinate in [0, 1].</returns>
        /// <remarks>
        /// For bounds [0, 1] the pycma margins are a = 0.05 on both sides. The map has period 2.2, is the identity on
        /// [0.05, 0.95], is (y + 0.05)^2 / 0.2 on [-0.05, 0.05], and mirrors outside [-0.05, 1.05].
        /// </remarks>
        internal static double Transform(double y)
        {
            const double a = 0.05;
            if (y < -2 * a - 0.5 || y > 1 + 2 * a + 0.5)
            {
                const double r = 2 * (1 + 2 * a), s = -2 * a - 0.5;
                y -= r * Math.Floor((y - s) / r);
            }
            if (y > 1 + a) y -= 2 * (y - 1 - a);
            if (y < -a) y += 2 * (-a - y);
            if (y < a) y = (y + a) * (y + a) / (4 * a);
            else if (y > 1 - a) y = 1 - (y - 1 - a) * (y - 1 - a) / (4 * a);
            return y;
        }

        /// <summary>
        /// Inverts <see cref="Transform(double)"/> on [0, 1], returning the genotype in [-0.05, 1.05].
        /// </summary>
        /// <param name="u">The feasible coordinate in [0, 1].</param>
        /// <returns>The genotype coordinate.</returns>
        internal static double InverseTransform(double u)
        {
            const double a = 0.05;
            if (u < a) return Math.Sqrt(4 * a * u) - a;
            if (u > 1 - a) return 1 + a - Math.Sqrt(4 * a * (1 - u));
            return u;
        }

    }
}
