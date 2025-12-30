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

            // Set the gradient function
            if (gradientFunction == null)
            {
                GradientFunction = (x) => new Vector(NumericalDerivative.Gradient((y) => LogLikelihoodFunction(y), x.ToArray()));
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
            if (Mass.Length != NumberOfParameters) throw new ArgumentException(nameof(Mass), "The mass vector must be the same length as the number of parameters.");
            if (StepSize < 0) throw new ArgumentException(nameof(StepSize), "The leapfrog step size must be positive.");
            if (Steps < 1) throw new ArgumentException(nameof(Steps), "The number of leapfrog steps must be at least one.");
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {

            // Update the sample count
            SampleCount[index] += 1;

            // Jigger the step size and number of steps
            var _stepSize = _stepSizeU.InverseCDF(_chainPRNGs[index].NextDouble());
            var _steps = (int)Math.Ceiling(_stepsU.InverseCDF(_chainPRNGs[index].NextDouble()));

            // Step 1. Sample phi from a N~(0,M)
            var phi = new Vector(NumberOfParameters);
            for (int i = 0; i < NumberOfParameters; i++)
                phi[i] = Math.Sqrt(Mass[i]) * Normal.StandardZ(_chainPRNGs[index].NextDouble());

            // Get kinetic energy of the current state
            var logKi = -0.5 * QuadraticForm(phi, _inverseMass);

            // Step 2. Perform leapfrog steps to get proposal vector
            var xp = new Vector(state.Values);
            phi += GradientFunction(xp.Array) * _stepSize * 0.5;
            for (int i = 0; i < _steps; i++)
            {
                xp += _inverseMass * phi * _stepSize;

                // Ensure the parameters are feasible (within the constraints)
                for (int j = 0; j < NumberOfParameters; j++)
                {
                    if (xp[j] < PriorDistributions[j].Minimum) 
                        xp[j] = PriorDistributions[j].Minimum + Tools.DoubleMachineEpsilon;
                    if (xp[j] > PriorDistributions[j].Maximum) 
                        xp[j] = PriorDistributions[j].Maximum - Tools.DoubleMachineEpsilon;
                }

                phi += GradientFunction(xp.Array) * _stepSize * (i == _steps-1 ? 0.5: 1.0);
            }
            phi *= -1d;

            // Get kinetic energy of the proposal state
            var logKp = -0.5 * QuadraticForm(phi, _inverseMass);

            // Evaluate fitness
            var logLHp = LogLikelihoodFunction(xp.Array);
            var logLHi = state.Fitness;

            // Calculate the Metropolis ratio
            var logRatio = logLHp - logLHi + logKp - logKi;

            // Accept the proposal with probability min(1,r)
            // otherwise leave xi unchanged
            var logU = Math.Log(_chainPRNGs[index].NextDouble());
            if (logU <= logRatio)
            {
                // The proposal is accepted
                AcceptCount[index] += 1;
                return new ParameterSet(xp.Array, logLHp);
            }
            else
            {
                return state;
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
