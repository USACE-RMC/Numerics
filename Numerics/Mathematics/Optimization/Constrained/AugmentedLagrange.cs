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
using System.Collections.ObjectModel;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The Augmented Lagrange constrained optimization method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This method replaces the constrained optimization problem with a series of unconstrained problems and 
    /// adds a penalty term to the objective, while also adding a third term designed to mimic a 
    /// Lagrange multiplier. In other words, the unconstrained objective is the Lagrangian of the constrained 
    /// problem, with an additional penalty term.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Augmented_Lagrangian_method#cite_note-Nocedal_2006-6"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public class AugmentedLagrange : Optimizer
    {

        /// <summary>
        /// Constructs new Augmented Lagrange optimizer.
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="constraints">The list of constraints.</param>
        /// <param name="optimizer">The internal optimizer to use in the Augmented Lagrange method.</param>
        public AugmentedLagrange(Func<double[], double> objectiveFunction, Optimizer optimizer, IList<IConstraint> constraints) : base(objectiveFunction, optimizer.NumberOfParameters)
        {        
            // Check if there are constraints
            if (constraints == null || constraints.Count == 0)
            {
                throw new ArgumentOutOfRangeException(nameof(constraints), "There must be at least one constraint.");
            }
            // Constraints
            _constraints = constraints.ToArray();
            _lambda = new double[_constraints.Where((x) => x.Type == ConstraintType.EqualTo).ToArray().Length];
            _mu = new double[_constraints.Where((x) => x.Type == ConstraintType.LesserThanOrEqualTo).ToArray().Length];
            _nu = new double[_constraints.Where((x) => x.Type == ConstraintType.GreaterThanOrEqualTo).ToArray().Length];

            // Set up objective functions and optimizer
            if (optimizer.GetType() == typeof(AugmentedLagrange)) throw new ArgumentException(nameof(optimizer), "The inner optimizer cannot also be an Augmented Lagrange optimizer.");
            _primaryObjectiveFunction = objectiveFunction;
            this.Optimizer = optimizer;
            this.Optimizer.ObjectiveFunction = augmentedLagrangianFunction;

        }

        private double rho = 1;
        private double rhoMax = 1e+1;
        private double rhoMin = 1e-6;
        private double[] _lambda; // equality multipliers
        private double[] _mu;     // "lesser than"  inequality multipliers
        private double[] _nu;     // "greater than" inequality multipliers
        IConstraint[] _constraints;
        private Func<double[], double> _primaryObjectiveFunction;

        /// <summary>
        /// The internal optimizer to use in the Augmented Lagrange method. 
        /// </summary>
        public Optimizer Optimizer { get; }

        /// <summary>
        /// Returns the list of constraints.
        /// </summary>
        public ReadOnlyCollection<IConstraint> Constraints => new ReadOnlyCollection<IConstraint>(_constraints);

        /// <summary>
        /// Returns the Lagrangian equality multipliers.
        /// </summary>
        public ReadOnlyCollection<double> Lambda => new ReadOnlyCollection<double>(_lambda);

        /// <summary>
        /// Returns the "Lesser than" inequality multipliers.
        /// </summary>
        public ReadOnlyCollection<double> Mu => new ReadOnlyCollection<double>(_mu);

        /// <summary>
        /// Returns the "Greater than" inequality multipliers.
        /// </summary>
        public ReadOnlyCollection<double> Nu => new ReadOnlyCollection<double>(_nu);

        /// <summary>
        /// The Augmented Lagrangian objective function. 
        /// </summary>
        private double augmentedLagrangianFunction(double[] x)
        {
            double phi = _primaryObjectiveFunction(x);
            double rho2 = 0.5 * rho;

            // For each equality constraint
            for (int i = 0; i < _constraints.Length; i++)
            {
                double actual = _constraints[i].Function(x);
                double c = 0;
                
                switch (_constraints[i].Type)
                {
                    case ConstraintType.EqualTo:
                        c = actual - _constraints[i].Value;
                        phi += rho2 * Math.Pow(c + _lambda[i] / rho, 2d);
                        break;

                    case ConstraintType.LesserThanOrEqualTo:
                        c = actual - _constraints[i].Value;
                        if (c > 0) phi += rho2 * Math.Pow(c + _mu[i] / rho, 2d);
                        break;

                    case ConstraintType.GreaterThanOrEqualTo:
                        c = _constraints[i].Value - actual;
                        if (c > 0) phi += rho2 * Math.Pow(c + _nu[i] / rho, 2d);
                        break;
                }
            }

            return phi;
        }

        /// <inheritdoc/>
        protected override void Optimize()
        {
            double ICM = Double.PositiveInfinity;
            double minPenalty = Double.PositiveInfinity;
            double penalty;
            double minFitness = Double.PositiveInfinity;
            double currentFitness;
            bool minFeasible = false;
            bool cancel = false;

            // magic parameters from Birgin & Martinez
            const double tau = 0.5, gam = 10;

            this.Optimizer.Minimize();
            var currentValues = this.Optimizer.BestParameterSet.Values.ToArray();

            Array.Clear(_lambda, 0, _lambda.Length);
            Array.Clear(_mu, 0, _mu.Length);
            Array.Clear(_nu, 0, _nu.Length);
            rho = 1;

            // Starting rho suggested by B & M 
            if (_lambda.Length > 0 || _mu.Length > 0 || _nu.Length > 0)
            {
                bool feasible = true;
                double con2 = 0;
                penalty = 0;

                // Evaluate function
                currentFitness = Evaluate(currentValues, ref cancel);

                // For each constraint
                for (int i = 0; i < _constraints.Length; i++)
                {
                    double actual = _constraints[i].Function(currentValues);
                    double c = 0;

                    switch (_constraints[i].Type)
                    {
                        case ConstraintType.EqualTo:
                            c = actual - _constraints[i].Value;
                            penalty += Math.Abs(c);
                            con2 += c * c;
                            feasible = feasible && Math.Abs(c) <= _constraints[i].Tolerance;
                            break;

                        case ConstraintType.GreaterThanOrEqualTo:
                            c = _constraints[i].Value - actual;
                            if (c > 0)
                            {
                                penalty += c;
                                con2 += c * c;
                            }
                            feasible = feasible && c <= _constraints[i].Tolerance;
                            break;

                        case ConstraintType.LesserThanOrEqualTo:
                            c = actual - _constraints[i].Value;
                            if (c > 0)
                            {
                                penalty += c;
                                con2 += c * c;
                            }
                            feasible = feasible && c <= _constraints[i].Tolerance;
                            break;
                    }
                }

                minFitness = currentFitness;
                minPenalty = penalty;
                minFeasible = feasible;
                double num = 2.0 * Math.Abs(minFitness);
                if (num < 1e-300)
                {
                    rho = rhoMin;
                }
                else if (con2 < 1e-300)
                {
                    rho = rhoMax;
                }
                else
                {
                    rho = num / con2;
                    if (rho < rhoMin)
                        rho = rhoMin;
                    else if (rho > rhoMax)
                        rho = rhoMax;
                }

            }

            while (Iterations < MaxIterations)
            {
                double prevICM = ICM;

                this.Optimizer.Minimize();
                currentValues = this.Optimizer.BestParameterSet.Values.ToArray();
                currentFitness = Evaluate(currentValues, ref cancel);
                if (cancel) return;

                ICM = 0;
                penalty = 0;
                bool feasible = true;

                // Update lambdas
                for (int i = 0; i < _constraints.Length; i++)
                {
                    double actual = _constraints[i].Function(currentValues);
                    double c = 0;
                    double newLambda = 0;

                    switch (_constraints[i].Type)
                    {
                        case ConstraintType.EqualTo:
                            c = actual - _constraints[i].Value;
                            newLambda = _lambda[i] + rho * c;
                            penalty += Math.Abs(c);
                            feasible = feasible && Math.Abs(c) <= _constraints[i].Tolerance;
                            ICM = Math.Max(ICM, Math.Abs(c));
                            _lambda[i] = Math.Min(Math.Max(-1e20, newLambda), 1e20);                           
                            break;

                        case ConstraintType.LesserThanOrEqualTo:
                            c = actual - _constraints[i].Value;
                            newLambda = _mu[i] + rho * c;
                            penalty += c > 0 ? c : 0;
                            feasible = feasible && c <= _constraints[i].Tolerance;
                            ICM = Math.Max(ICM, Math.Abs(Math.Max(c, -_mu[i] / rho)));
                            _mu[i] = Math.Min(Math.Max(0.0, newLambda), 1e20);
                            break;

                        case ConstraintType.GreaterThanOrEqualTo:
                            c = _constraints[i].Value - actual;
                            newLambda = _nu[i] + rho * c;
                            penalty += c > 0 ? c : 0;
                            feasible = feasible && c <= _constraints[i].Tolerance;
                            ICM = Math.Max(ICM, Math.Abs(Math.Max(c, -_nu[i] / rho)));
                            _nu[i] = Math.Min(Math.Max(0.0, newLambda), 1e20);
                            break;
                    }
            }
                // Update rho
                if (ICM > tau * prevICM)
                {
                    rho *= gam;
                }

                // Determine if the solution as converged. 
                bool a = !minFeasible || penalty < minPenalty || currentFitness < minFitness;
                bool b = !minFeasible && penalty < minPenalty;
                if ((feasible && a) || b)
                {
                    if (feasible && CheckConvergence(minFitness, currentFitness))
                    {
                        UpdateStatus(OptimizationStatus.Success);
                        return;
                    }

                    BestParameterSet = new ParameterSet(currentValues, currentFitness);
                    minFitness = currentFitness;
                    minPenalty = penalty;
                    minFeasible = feasible;

                }
                else if (ICM == 0)
                {
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }

                Iterations += 1;
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

    }
}
