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
using Numerics.Mathematics.LinearAlgebra;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// A base class for all optimization methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class Optimizer
    {
        /// <summary>
        /// Construct a new optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        protected Optimizer(Func<double[], double> objectiveFunction, int numberOfParameters)
        {
            ObjectiveFunction = objectiveFunction;
            NumberOfParameters = numberOfParameters < 1 ? throw new ArgumentNullException(nameof(numberOfParameters), "There must be at least 1 parameter to evaluate.") : numberOfParameters;
        }

        #region Inputs

        private Func<double[], double> _objectiveFunction;

        /// <summary>
        /// The maximum number of optimization iterations allowed. Default = 10,000.
        /// </summary>
        public int MaxIterations { get; set; } = 10000;

        /// <summary>
        /// The maximum number of function evaluations allowed. Default = int.MaxValue.
        /// </summary>
        public int MaxFunctionEvaluations { get; set; } = int.MaxValue;

        /// <summary>
        /// The desired absolute tolerance for the solution. Default = ~Sqrt(Machine Epsilon), or 1E-8.
        /// </summary>
        public double AbsoluteTolerance { get; set; } = 1E-8;

        /// <summary>
        /// The desired relative tolerance for the solution. Default = ~Sqrt(Machine Epsilon), or 1E-8.
        /// </summary>
        public double RelativeTolerance { get; set; } = 1E-8;

        /// <summary>
        /// Determines if an exception will be thrown if the optimization solver fails to converge.
        /// </summary>
        public bool ReportFailure { get; set; } = true;

        /// <summary>
        /// Determines whether the parameter set traces are saved or not. Default = true.
        /// </summary>
        public bool RecordTraces { get; set; } = true;

        /// <summary>
        /// Determines whether to compute a numerically differentiated Hessian matrix when the optimization was successful. 
        /// </summary>
        public bool ComputeHessian { get; set; } = true;

        /// <summary>
        /// The objective function to evaluate. 
        /// </summary>
        public Func<double[], double> ObjectiveFunction 
        {
            get { return _objectiveFunction; }
            set
            {
                _objectiveFunction = value ?? throw new ArgumentNullException(nameof(ObjectiveFunction), "The objective function cannot be null.");
            }
        }

        /// <summary>
        /// The number of parameters to evaluate in the objective function.
        /// </summary>
        public int NumberOfParameters { get; protected set; }

        /// <summary>
        /// Objective function scaling factor (set to -1 for a maximization problem). By default it is a minimization problem.
        /// </summary>
        protected int functionScale = 1;


        #endregion

        #region Output

        /// <summary>
        /// Returns the number of iterations required to find the solution.
        /// </summary>
        public int Iterations { get; protected set; }

        /// <summary>
        /// Returns the number of function evaluations required to find the solution.
        /// </summary>
        public int FunctionEvaluations { get; protected set; }

        /// <summary>
        /// The optimal point, or parameter set. 
        /// </summary>
        public ParameterSet BestParameterSet { get; protected set; }

        /// <summary>
        /// A trace of the parameter set and fitness evaluated until convergence.
        /// </summary>
        public List<ParameterSet> ParameterSetTrace { get; protected set; }

        /// <summary>
        /// Determines the optimization method status. 
        /// </summary>
        public OptimizationStatus Status { get; protected set; } = OptimizationStatus.None;

        /// <summary>
        /// The numerically differentiated Hessian matrix. This is only computed when the optimization is successful.
        /// </summary>
        public Matrix Hessian { get; protected set; }

        #endregion

        /// <summary>
        /// Clears the results.
        /// </summary>
        public virtual void ClearResults()
        {
            Iterations = 0;
            FunctionEvaluations = 0;
            BestParameterSet = new ParameterSet();
            ParameterSetTrace = new List<ParameterSet>();
            Status = OptimizationStatus.None;
            Hessian = null;
        }

        /// <summary>
        /// Validate inputs. 
        /// </summary>
        protected virtual void Validate()
        {
            if (MaxIterations < 10) throw new ArgumentOutOfRangeException(nameof(MaxIterations), "The maximum number of optimization iterations must be greater than 10.");
            if (MaxFunctionEvaluations < 10) throw new ArgumentOutOfRangeException(nameof(MaxFunctionEvaluations), "The maximum number of function evaluations must be greater than 10.");
            if (RelativeTolerance <= 0 || RelativeTolerance > 1) throw new ArgumentOutOfRangeException(nameof(RelativeTolerance), "The relative tolerance must be between 0 and 1.");
            if (AbsoluteTolerance <= 0 || AbsoluteTolerance > 1) throw new ArgumentOutOfRangeException(nameof(AbsoluteTolerance), "The absolute tolerance must be between 0 and 1.");
        }

        /// <summary>
        /// Implements the actual optimization algorithm. This method should minimize the objective function. 
        /// </summary>
        protected abstract void Optimize();

        /// <summary>
        /// Finds the parameter set that minimizes the objective function.
        /// </summary>
        public virtual void Minimize()
        {
            Validate();
            ClearResults();
            functionScale = 1;
            try
            {               
                Optimize();
                if (Status == OptimizationStatus.Success && ComputeHessian)
                {
                    Hessian = new Matrix(NumericalDerivative.Hessian((x) => { return ObjectiveFunction(x); }, BestParameterSet.Values));
                }
            }
            catch (ArgumentException ex)
            {
                if (ex.ParamName != nameof(MaxIterations) && ex.ParamName != nameof(MaxFunctionEvaluations))
                    UpdateStatus(OptimizationStatus.Failure, ex);
            }
            catch (Exception ex)
            {
                UpdateStatus(OptimizationStatus.Failure, ex);
            }

        }

        /// <summary>
        /// Finds the parameter set that maximizes the objective function.
        /// </summary>
        public virtual void Maximize()
        {
            Validate();
            ClearResults();
            functionScale = -1;
            try
            {
                Optimize();
                if (Status == OptimizationStatus.Success && ComputeHessian)
                {
                    Hessian = new Matrix(NumericalDerivative.Hessian((x) => { return ObjectiveFunction(x); }, BestParameterSet.Values));
                }
            }
            catch (ArgumentException ex)
            {
                if (ex.ParamName != nameof(MaxIterations) && ex.ParamName != nameof(MaxFunctionEvaluations))
                    UpdateStatus(OptimizationStatus.Failure, ex);
            }
            catch (Exception ex)
            {
                UpdateStatus(OptimizationStatus.Failure, ex);
            }
        }

        /// <summary>
        /// Evaluates the objective function and returns the fitness.
        /// </summary>
        /// <param name="values">The parameter values to evaluate.</param>
        /// <param name="cancel">By ref. Determines if the solver should be canceled.</param>
        protected virtual double Evaluate(double[] values, ref bool cancel)
        {        
            double fitness = functionScale * ObjectiveFunction(values);

            // keep track of best fit parameter set
            if (BestParameterSet.Values == null || fitness <= BestParameterSet.Fitness)
                BestParameterSet = new ParameterSet((double[])values.Clone(), fitness);
            
            // Update trace. This is tracked every evaluation
            if (ParameterSetTrace == null) ParameterSetTrace = new List<ParameterSet>();
            if (RecordTraces) ParameterSetTrace.Add(BestParameterSet.Clone());

            // update evaluation counter
            FunctionEvaluations += 1;
            if (FunctionEvaluations >= MaxFunctionEvaluations)
            {
                cancel = true;
                UpdateStatus(OptimizationStatus.MaximumFunctionEvaluationsReached);          
            }

            return fitness;
        }

        /// <summary>
        /// Repair the trial parameter to keep it within bounds.
        /// </summary>
        /// <param name="value">The parameter value to evaluate.</param>
        /// <param name="lowerBound">The lower bound of the interval containing the minima.</param>
        /// <param name="upperBound">The upper bound of the interval containing the minima.</param>
        protected virtual double RepairParameter(double value, double lowerBound, double upperBound)
        {
            // check if the trial parameter set is out of bounds
            if (value < lowerBound) value = lowerBound;
            if (value > upperBound) value = upperBound;
            return value;
        }

        /// <summary>
        /// Update the optimization status. Exceptions will be throw depending on the status. 
        /// </summary>
        /// <param name="status">Optimization status.</param>
        /// <param name="exception">Inner exception.</param>
        protected virtual void UpdateStatus(OptimizationStatus status, Exception exception = null)
        {
            Status = status;
            if (status == OptimizationStatus.MaximumIterationsReached)
            {
                if (ReportFailure) throw new ArgumentException("The maximum number of iterations has been reached.", nameof(MaxIterations), exception);
            }
            else if (status == OptimizationStatus.MaximumFunctionEvaluationsReached)
            {
                if (ReportFailure) throw new ArgumentException("The maximum number of function evaluations has been reached.", nameof(MaxFunctionEvaluations), exception);
            }
            else if (status == OptimizationStatus.Failure)
            {
              if (ReportFailure && exception != null) throw exception;
            }

        }

        /// <summary>
        /// Checks convergence. 
        /// </summary>
        /// <param name="oldValue">The old value.</param>
        /// <param name="newValue">The new value.</param>
        protected virtual bool CheckConvergence(double oldValue, double newValue)
        {
            if (double.IsNaN(oldValue) || double.IsNaN(newValue) || double.IsInfinity(oldValue) || double.IsInfinity(newValue)) return false;
            return 2.0d * Math.Abs(newValue - oldValue) / (Math.Abs(newValue) + Math.Abs(oldValue) + AbsoluteTolerance) < RelativeTolerance;
        }



    }
}
