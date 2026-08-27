using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Sampling.MCMC
{
    /// <summary>
    /// This class performs Bayesian MCMC using the Hamiltonian Monte Carlo (HMC) method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b>Authors:</b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description:</b>
    /// </para>
    /// <para>
    ///     The optimal acceptance rate for this sampler is 65%, whereas Metropolis-Hastings samplers have an optimal rate of 23.4%.
    /// </para>
    /// <para>
    ///     A trajectory of L leapfrog steps costs at most L + 1 gradient evaluations: the closing
    ///     half-step of each step is fused with the opening half-step of the next. After a chain's first
    ///     transition the opening evaluation is served from a per-chain memo of the previous transition's
    ///     closing evaluation, reducing the cost to L.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    ///    <see href="https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public class HMC : MCMCSampler
    {

        /// <summary>
        /// The function for evaluating the gradient of the log-likelihood function.
        /// </summary>
        /// <param name="parameters">The list of parameters to evaluate.</param>
        /// <returns>Returns the gradient of the log-likelihood function.</returns>
        /// <remarks>
        /// The sampler hands the delegate a private working array, never the chain state itself, so an
        /// implementation that writes through its argument cannot corrupt the chain. A gradient the
        /// sampler keeps for reuse is copied out of the returned vector, so an implementation may return
        /// one reused buffer.
        /// </remarks>
        public delegate Vector Gradient(IList<double> parameters);

        /// <summary>
        /// Constructs a new HMC sampler.
        /// </summary>
        /// <param name="priorDistributions"></param>
        /// <param name="logLikelihoodFunction"></param>
        /// <param name="mass">Optional. The mass vector for the momentum distribution. Default = Identity.</param>
        /// <param name="stepSize">Optional. The leapfrog step size. Default = 0.1.</param>
        /// <param name="steps">Optional. The number of leapfrog steps. Default = 50.</param>
        /// <param name="gradientFunction">Optional. The function for evaluating the gradient of the log-likelihood. Numerical finite difference will be used by default.</param>
        public HMC(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction, Vector mass = null!, double stepSize = 0.1, int steps = 50, Gradient gradientFunction = null!) : base(priorDistributions, logLikelihoodFunction)
        {

            InitialIterations = 100 * NumberOfParameters;
           
            // Set the mass vector 
            if (mass == null)
            {
                Mass = new Vector(NumberOfParameters, 1d);
            }
            else
            {
                Mass = mass;
            }
       
            // Set the inverse mass vector
            _inverseMass = new Vector(NumberOfParameters);
            for (int i = 0; i < NumberOfParameters; i++)
                _inverseMass[i] = 1d / Mass[i];

            // Set leapfrog inputs
            StepSize = stepSize;
            Steps = steps;

            // Cache prior distribution bounds for the gradient function
            _lowerBounds = new double[NumberOfParameters];
            _upperBounds = new double[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                _lowerBounds[i] = priorDistributions[i].Minimum;
                _upperBounds[i] = priorDistributions[i].Maximum;
            }

            // Set the gradient function with prior bounds so finite-difference probes stay in valid region
            if (gradientFunction == null)
            {
                GradientFunction = (x) => new Vector(NumericalDerivative.Gradient((y) => SafeLogLikelihood(y), x.ToArray(), _lowerBounds, _upperBounds));
            }
            else
            {
                GradientFunction = gradientFunction;
            }

        }

        private Uniform _stepSizeU = new Uniform(0.00, 0.2);
        private UniformDiscrete _stepsU = new UniformDiscrete(1, 100);
        private Vector _inverseMass;
        private double _stepSize = 0.1;
        private int _steps = 50;
        private double[] _lowerBounds;
        private double[] _upperBounds;

        // Per-chain memo of the gradient at the current chain state. A leapfrog trajectory closes with a
        // gradient evaluation at its final position, and the next transition opens by evaluating the
        // gradient at that same position when the proposal was accepted, or at the unchanged state when it
        // was rejected, so the opening evaluation is redundant after the chain's first transition. The
        // stored arrays are owned by the memo and are copies in both directions: positions are recorded
        // before the delegate runs and gradients are copied out of the returned vector, so neither a
        // delegate that writes through its argument nor one that returns a reused buffer can corrupt an
        // entry. Entries are matched on bit patterns, never with ==, so +0 and -0 are distinct positions.
        private double[][] _startGradientPositions = null!;
        private double[][] _startGradientValues = null!;
        private bool[] _startGradientOccupied = null!;

        /// <summary>
        /// The mass vector for the momentum distribution.
        /// </summary>
        public Vector Mass { get; }

        /// <summary>
        /// The leapfrog step size. Default = 0.1.
        /// </summary>
        /// <remarks>
        /// This controls the size of each leapfrog step in the simulation. A smaller step size can lead to more accurate 
        /// simulations of the Hamiltonian dynamics but requires more steps to cover the same distance, increasing 
        /// computational cost. A larger step size reduces computational cost but can lead to inaccurate simulations 
        /// and a higher rejection rate. The step size is often tuned during the warm-up phase of the HMC algorithm to 
        /// achieve an optimal balance.
        /// </remarks>
        public double StepSize
        {
            get { return _stepSize; }
            set 
            { 
                _stepSize = value;
                _stepSizeU = new Uniform(0, 2.0 * _stepSize);
            }
        }

        /// <summary>
        /// The number of leapfrog steps. Default = 50.
        /// </summary>
        /// <remarks>
        /// This refers to the number of leapfrog steps taken in each HMC iteration. 
        /// The total distance covered in the parameter space during each iteration is the product 
        /// of the step size and the number of steps. The number of steps determines how far the 
        /// algorithm moves in parameter space before proposing a new sample. A larger number of steps 
        /// allows the algorithm to explore more of the parameter space but can increase 
        /// computational time. Like the step size, the number of steps can also be tuned, though 
        /// it is often set to a fixed value.
        /// </remarks>
        public int Steps
        {
            get { return _steps; }
            set
            {
                _steps= value;
                _stepsU = new UniformDiscrete(1, 2.0 * _steps);
            }
        }

        /// <summary>
        /// The function for evaluating the gradient of the log-likelihood.
        /// </summary>
        public Gradient GradientFunction { get; }

        /// <inheritdoc/>
        protected override void ValidateCustomSettings()
        {
            if (Mass.Length != NumberOfParameters) throw new ArgumentException("The mass vector must be the same length as the number of parameters.", nameof(Mass));
            if (StepSize < 0) throw new ArgumentException("The leapfrog step size must be positive.", nameof(StepSize));
            if (Steps < 1) throw new ArgumentException("The number of leapfrog steps must be at least one.", nameof(Steps));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Allocates the per-chain start-gradient memo. A resumed simulation skips this, keeping the
        /// entries from the run it continues; those entries still describe the chain states the resume
        /// starts from, so they remain valid.
        /// </remarks>
        protected override void InitializeCustomSettings()
        {
            _startGradientPositions = new double[NumberOfChains][];
            _startGradientValues = new double[NumberOfChains][];
            _startGradientOccupied = new bool[NumberOfChains];
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {

            // Update the sample count
            SampleCount[index] += 1;

            try
            {
                // Jigger the step size and number of steps
                var _stepSize = _stepSizeU.InverseCDF(_chainPRNGs[index].NextDouble());
                var _steps = (int)Math.Ceiling(_stepsU.InverseCDF(_chainPRNGs[index].NextDouble()));

                // Step 1. Sample phi from a N~(0,M)
                var phi = new Vector(NumberOfParameters);
                for (int i = 0; i < NumberOfParameters; i++)
                    phi[i] = Math.Sqrt(Mass[i]) * Normal.StandardZ(_chainPRNGs[index].NextDouble());

                // Get kinetic energy of the current state
                var logKi = -0.5 * QuadraticForm(phi, _inverseMass);

                // Step 2. Perform leapfrog steps to get proposal vector. The trajectory works on a
                // private copy of the chain state, so the gradient delegate never sees the chain's own
                // array.
                var xp = new Vector((double[])state.Values.Clone());
                phi += StartGradient(index, xp.Array) * _stepSize * 0.5;
                double[] proposalPosition = null!;
                Vector proposalGradient = null!;
                for (int i = 0; i < _steps; i++)
                {
                    xp += _inverseMass * phi * _stepSize;

                    // Ensure the parameters are feasible (within the constraints)
                    for (int j = 0; j < NumberOfParameters; j++)
                    {
                        if (xp[j] < PriorDistributions[j].Minimum)
                        {
                            xp[j] = PriorDistributions[j].Minimum + Tools.DoubleMachineEpsilon;
                            phi[j] = -phi[j];
                        }
                        if (xp[j] > PriorDistributions[j].Maximum)
                        {
                            xp[j] = PriorDistributions[j].Maximum - Tools.DoubleMachineEpsilon;
                            phi[j] = -phi[j];
                        }
                    }

                    if (i == _steps - 1)
                    {
                        // The closing half-step evaluates the gradient at the proposal position. Keep the
                        // position, recorded before the delegate runs, and the gradient, so an accepted
                        // proposal seeds the memo for the next transition's opening half-step.
                        proposalPosition = (double[])xp.Array.Clone();
                        proposalGradient = GradientFunction(xp.Array);
                        phi += proposalGradient * _stepSize * 0.5;
                    }
                    else
                    {
                        phi += GradientFunction(xp.Array) * _stepSize;
                    }
                }
                phi *= -1d;

                // Get kinetic energy of the proposal state
                var logKp = -0.5 * QuadraticForm(phi, _inverseMass);

                // Evaluate fitness
                var logLHp = SafeLogLikelihood(xp.Array);
                var logLHi = state.Fitness;

                // Calculate the Metropolis ratio
                var logRatio = logLHp - logLHi + logKp - logKi;

                // Accept the proposal with probability min(1,r)
                // otherwise leave xi unchanged
                var logU = Math.Log(_chainPRNGs[index].NextDouble());
                if (logU <= logRatio)
                {
                    // The proposal is accepted. The closing half-step's gradient becomes the opening
                    // gradient of the next transition, so store it against the accepted position. On the
                    // reject path the memo already holds the entry for the unchanged state.
                    AcceptCount[index] += 1;
                    StoreStartGradient(index, proposalPosition, proposalGradient);
                    return new ParameterSet(xp.Array, logLHp);
                }
                else
                {
                    return state;
                }
            }
            catch (ArithmeticException)
            {
                // Non-finite gradient encountered during leapfrog integration.
                // This occurs when parameters drift into regions where the log-likelihood
                // returns -Infinity. Reject the proposal and return the current state,
                // consistent with Metropolis rejection behavior.
                return state;
            }

        }

        /// <summary>
        /// Returns the gradient at the trajectory start, served from the per-chain memo when the position
        /// matches the stored entry bit for bit, and evaluated through <see cref="GradientFunction"/>
        /// otherwise.
        /// </summary>
        /// <param name="index">The Markov Chain zero-based index.</param>
        /// <param name="position">The trajectory start position. The array is the trajectory's private
        /// working copy, so the delegate may be handed it directly.</param>
        /// <returns>The gradient at <paramref name="position"/>. On a memo hit the returned vector wraps
        /// the memo's stored array and must be treated as read-only; every use here only feeds allocating
        /// vector arithmetic.</returns>
        /// <remarks>
        /// A miss records the position before the delegate runs, so a delegate that writes through its
        /// argument cannot pair a mutated position with another point's gradient; such a mutation only
        /// produces a harmless later miss. Matching is on <see cref="BitConverter.DoubleToInt64Bits(double)"/>
        /// rather than <c>==</c>, so a position differing only in the sign of a zero is a distinct
        /// position, and NaN never matches so a non-finite state can never be served a stale gradient.
        /// </remarks>
        private Vector StartGradient(int index, double[] position)
        {
            if (_startGradientOccupied[index] && MatchesBitwise(_startGradientPositions[index], position))
                return new Vector(_startGradientValues[index]);

            var recorded = (double[])position.Clone();
            var gradient = GradientFunction(position);
            StoreStartGradient(index, recorded, gradient);
            return gradient;
        }

        /// <summary>
        /// Stores a position and its gradient as the chain's start-gradient memo entry, copying the
        /// gradient out of the delegate's vector so a reused buffer cannot rewrite the entry later.
        /// </summary>
        /// <param name="index">The Markov Chain zero-based index.</param>
        /// <param name="position">The position the gradient was evaluated at, recorded before the
        /// delegate ran. The memo takes ownership of this array.</param>
        /// <param name="gradient">The gradient the delegate returned. A gradient of the wrong length
        /// clears the entry instead of storing it.</param>
        private void StoreStartGradient(int index, double[] position, Vector gradient)
        {
            if (gradient.Length != NumberOfParameters)
            {
                _startGradientOccupied[index] = false;
                return;
            }
            _startGradientPositions[index] = position;
            _startGradientValues[index] = (double[])gradient.Array.Clone();
            _startGradientOccupied[index] = true;
        }

        /// <summary>
        /// Compares a stored memo position against a query position on the bit patterns of every
        /// coordinate.
        /// </summary>
        /// <param name="stored">The memo's stored position.</param>
        /// <param name="position">The query position.</param>
        /// <returns>True when every coordinate matches bit for bit.</returns>
        private static bool MatchesBitwise(double[] stored, double[] position)
        {
            if (stored.Length != position.Length) return false;
            for (int i = 0; i < stored.Length; i++)
            {
                if (BitConverter.DoubleToInt64Bits(stored[i]) != BitConverter.DoubleToInt64Bits(position[i]))
                    return false;
            }
            return true;
        }

        /// <summary>
        /// Evaluates the log-likelihood, returning negative infinity if the parameters are out of range.
        /// This prevents ArgumentOutOfRangeException from propagating during leapfrog integration
        /// when the sampler explores parameter values that violate distribution constraints.
        /// </summary>
        private double SafeLogLikelihood(double[] parameters)
        {
            try
            {
                return LogLikelihoodFunction(parameters);
            }
            catch (ArgumentOutOfRangeException)
            {
                return double.NegativeInfinity;
            }
        }

        /// <summary>
        /// Computes the quadratic form φᵀ M⁻¹ φ, which is used to calculate the kinetic energy 
        /// in Hamiltonian Monte Carlo (HMC) sampling. This method avoids allocating intermediate arrays.
        /// </summary>
        /// <param name="phi">
        /// The momentum vector φ (phi), typically drawn from a normal distribution scaled by mass.
        /// </param>
        /// <param name="inverseMass">
        /// The element-wise inverse of the mass vector M. Each component corresponds to 1 / mass[i].
        /// </param>
        /// <returns>
        /// The scalar result of the quadratic form φᵀ M⁻¹ φ, which equals the sum over i of (phi[i]^2 * inverseMass[i]).
        /// </returns>
        /// <exception cref="ArgumentException">
        /// Thrown if <paramref name="phi"/> and <paramref name="inverseMass"/> are not the same length.
        /// </exception>
        public static double QuadraticForm(Vector phi, Vector inverseMass)
        {
            if (phi.Length != inverseMass.Length)
                throw new ArgumentException("Vectors must be the same length to compute the quadratic form.");

            double sum = 0;
            for (int i = 0; i < phi.Length; i++)
                sum += inverseMass[i] * phi[i] * phi[i];
            return sum;
        }

    }
}
