using System;

namespace Numerics.Mathematics.Integration
{
    /// <summary>
    /// A base class for all integration methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class Integrator
    {

        #region Inputs

        /// <summary>
        /// The minimum number of integration iterations allowed. Default = 1.
        /// </summary>
        public int MinIterations { get; set; } = 1;

        /// <summary>
        /// The maximum number of integration iterations allowed. Default = 1E7.
        /// </summary>
        public int MaxIterations { get; set; } = 10000000;

        /// <summary>
        /// The minimum number of function evaluations allowed. Default = 1.
        /// </summary>
        public int MinFunctionEvaluations { get; set; } = 1;

        /// <summary>
        /// The maximum number of function evaluations allowed. Default = 1E7.
        /// </summary>
        public int MaxFunctionEvaluations { get; set; } = 10000000;

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
        /// Protected result property.
        /// </summary>
        protected double _result;

        /// <summary>
        /// The numerically computed result of the definite integral. 
        /// </summary>
        public double Result { get { return _result; } set { _result = value; } }


        /// <summary>
        /// Determines the integration method status. 
        /// </summary>
        public IntegrationStatus Status { get; protected set; } = IntegrationStatus.None;

        #endregion

        #region Methods

        /// <summary>
        /// Evaluates the integral.
        /// </summary>
        public abstract void Integrate();

        /// <summary>
        /// Clears the results.
        /// </summary>
        protected virtual void ClearResults()
        {
            Iterations = 0;
            FunctionEvaluations = 0;
            Result = double.NaN;
            Status = IntegrationStatus.None;
        }

        /// <summary>
        /// Validate inputs. 
        /// </summary>
        protected virtual void Validate()
        {
            if (MinIterations < 1) throw new ArgumentOutOfRangeException(nameof(MinIterations), "The minimum number of iterations must be greater than or equal to 1.");
            if (MinFunctionEvaluations < 1) throw new ArgumentOutOfRangeException(nameof(MinFunctionEvaluations), "The minimum number of function evaluations must be greater than or equal to 1.");
            if (MaxIterations < 1) throw new ArgumentOutOfRangeException(nameof(MaxIterations), "The maximum number of iterations must be greater than 1.");
            if (MaxFunctionEvaluations < 1) throw new ArgumentOutOfRangeException(nameof(MaxFunctionEvaluations), "The maximum number of function evaluations must be greater than 1.");
            if (RelativeTolerance < 1E-15 || RelativeTolerance > 1) throw new ArgumentOutOfRangeException(nameof(RelativeTolerance), "The relative tolerance must be between 1E-15 and 1.");
            if (AbsoluteTolerance < 1E-15 || AbsoluteTolerance > 1) throw new ArgumentOutOfRangeException(nameof(AbsoluteTolerance), "The absolute tolerance must be between 1E-15 and 1.");
        }

        /// <summary>
        /// Update the optimization status. Exceptions will be throw depending on the status. 
        /// </summary>
        /// <param name="status">Optimization status.</param>
        /// <param name="exception">Inner exception.</param>
        protected virtual void UpdateStatus(IntegrationStatus status, Exception exception = null!)
        {
            Status = status;
            if (status == IntegrationStatus.Failure)
            {
                if (ReportFailure && exception != null) throw exception;
            }

        }

        /// <summary>
        /// Evaluate convergence. 
        /// </summary>
        /// <param name="oldValue">The old value.</param>
        /// <param name="newValue">The new value.</param>
        /// <returns>
        /// Returns true for convergence
        /// </returns>
        protected virtual bool EvaluateConvergence(double oldValue, double newValue)
        {
            if (double.IsNaN(oldValue) || double.IsNaN(newValue) || double.IsInfinity(oldValue) || double.IsInfinity(newValue)) return false;
            return Math.Abs(newValue - oldValue) < AbsoluteTolerance && Math.Abs(newValue - oldValue) / Math.Abs(newValue) < RelativeTolerance;
        }

        #endregion

    }
}
