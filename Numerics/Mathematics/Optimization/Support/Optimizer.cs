using System;
using System.Collections.Generic;
using System.Runtime.ExceptionServices;
using System.Threading;
using System.Threading.Tasks;
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

        private Func<double[], double> _objectiveFunction = null!;

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
        /// Determines whether an exception that stops the solver is rethrown to the caller. Default = true.
        /// </summary>
        /// <remarks>
        /// Reaching <see cref="MaxIterations"/> or <see cref="MaxFunctionEvaluations"/> never throws from
        /// <see cref="Minimize"/> or <see cref="Maximize"/>: it is reported through <see cref="Status"/> only,
        /// whatever this setting. When true, any other exception raised by the objective function or the
        /// algorithm propagates after <see cref="Status"/> is set to <see cref="OptimizationStatus.Failure"/>;
        /// when false, it is recorded as that status and the call returns.
        /// </remarks>
        public bool ReportFailure { get; set; } = true;

        /// <summary>
        /// Determines whether the parameter set traces are saved or not. Default = true.
        /// </summary>
        public bool RecordTraces { get; set; } = true;

        /// <summary>
        /// Determines whether to compute a numerically differentiated Hessian matrix when the optimization was successful
        /// or a line search failed after producing a usable best parameter set.
        /// </summary>
        public bool ComputeHessian { get; set; } = true;

        /// <summary>
        /// A token that stops the run early. Default = none.
        /// </summary>
        /// <remarks>
        /// The token is checked at every function evaluation, including those of a local search run inside a
        /// global one. A cancelled run stops without throwing, keeps the best parameter set found so far, and
        /// ends with <see cref="OptimizationStatus.Cancelled"/>; no Hessian is computed. A token that is never
        /// cancelled does not change the results.
        /// </remarks>
        public CancellationToken CancellationToken { get; set; }

        /// <summary>
        /// Receives a snapshot of the run every <see cref="ProgressInterval"/> function evaluations and once
        /// when the run ends. Default = none.
        /// </summary>
        /// <remarks>
        /// Reports are made on the optimizer's thread; the snapshot holds a copy of the best values. Reporting
        /// does not change the results.
        /// </remarks>
        public IProgress<OptimizationProgress>? Progress { get; set; }

        /// <summary>
        /// The number of function evaluations between <see cref="Progress"/> reports. Default = 100.
        /// </summary>
        public int ProgressInterval { get; set; } = 100;

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

        /// <summary>
        /// The inclusive lower bounds of the parameter space, or null when the search is unbounded.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The base class returns null, which is correct for an unconstrained problem. A derived class that
        /// constrains its search to a box overrides this to expose those bounds, so that the finite-difference
        /// Hessian computed at the end of a successful <see cref="Minimize"/> or <see cref="Maximize"/> keeps
        /// its perturbations inside the feasible region.
        /// </para>
        /// <para>
        /// This is deliberately separate from the public bounds properties declared by the individual solvers.
        /// Those properties are declared independently on each derived class, so a member of the same name on
        /// this base class would be hidden by every one of them. The array returned here must have one entry
        /// per parameter, in the same order as <see cref="ParameterSet.Values"/>.
        /// </para>
        /// </remarks>
        protected virtual double[]? ParameterLowerBounds => null;

        /// <summary>
        /// The inclusive upper bounds of the parameter space, or null when the search is unbounded.
        /// </summary>
        /// <remarks>
        /// See <see cref="ParameterLowerBounds"/> for the rationale and the ordering contract.
        /// </remarks>
        protected virtual double[]? ParameterUpperBounds => null;


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
        /// <remarks>
        /// <see cref="ParameterSet.Fitness"/> is the value the solver minimized. After <see cref="Maximize"/>
        /// that is the negated objective, so the maximum of the objective is <c>-BestParameterSet.Fitness</c>.
        /// The same holds for every entry of <see cref="ParameterSetTrace"/>.
        /// </remarks>
        public ParameterSet BestParameterSet { get; protected set; } = new ParameterSet();

        /// <summary>
        /// A trace of the best-so-far parameter set and fitness at every function evaluation.
        /// </summary>
        /// <remarks>
        /// The trace is read-only: entries recorded between improvements share one values array,
        /// so the trace stores the search history without one array allocation per function
        /// evaluation.
        /// </remarks>
        public List<ParameterSet> ParameterSetTrace { get; protected set; } = new List<ParameterSet>();

        /// <summary>
        /// Determines the optimization method status.
        /// </summary>
        public OptimizationStatus Status { get; protected set; } = OptimizationStatus.None;

        /// <summary>
        /// The numerically differentiated Hessian matrix of the objective function at the best parameter set.
        /// </summary>
        /// <remarks>
        /// It is computed when <see cref="ComputeHessian"/> is true and the run ends with
        /// <see cref="OptimizationStatus.Success"/> or <see cref="OptimizationStatus.LineSearchFailed"/>;
        /// otherwise it is null. It is not negated after <see cref="Maximize"/>.
        /// </remarks>
        public Matrix? Hessian { get; protected set; }

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
            if (ProgressInterval < 1) throw new ArgumentOutOfRangeException(nameof(ProgressInterval), "The progress interval must be at least one evaluation.");
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
                if ((Status == OptimizationStatus.Success || Status == OptimizationStatus.LineSearchFailed) && ComputeHessian)
                {
                    Hessian = new Matrix(NumericalDerivative.Hessian((x) => { return ObjectiveFunction(x); }, BestParameterSet.Values, ParameterLowerBounds!, ParameterUpperBounds!));
                }
            }
            catch (OperationCanceledException) when (CancellationToken.IsCancellationRequested)
            {
                Status = OptimizationStatus.Cancelled;
            }
            catch (ArgumentException ex)
            {
                if (CancellationToken.IsCancellationRequested) Status = OptimizationStatus.Cancelled;
                else if (ex.ParamName != nameof(MaxIterations) && ex.ParamName != nameof(MaxFunctionEvaluations))
                    UpdateStatus(OptimizationStatus.Failure, ex);
            }
            catch (Exception ex)
            {
                // Once cancellation is requested, a later failure is a consequence of the stop, not an error.
                if (CancellationToken.IsCancellationRequested) Status = OptimizationStatus.Cancelled;
                else UpdateStatus(OptimizationStatus.Failure, ex);
            }
            ReportProgress();

        }

        /// <summary>
        /// Finds the parameter set that maximizes the objective function.
        /// </summary>
        /// <remarks>
        /// The solver minimizes the negated objective, so the fitness values it reports are negated; see
        /// <see cref="BestParameterSet"/>. The <see cref="Hessian"/> is that of the objective itself.
        /// </remarks>
        public virtual void Maximize()
        {
            Validate();
            ClearResults();
            functionScale = -1;
            try
            {
                Optimize();
                if ((Status == OptimizationStatus.Success || Status == OptimizationStatus.LineSearchFailed) && ComputeHessian)
                {
                    Hessian = new Matrix(NumericalDerivative.Hessian((x) => { return ObjectiveFunction(x); }, BestParameterSet.Values, ParameterLowerBounds!, ParameterUpperBounds!));
                }
            }
            catch (OperationCanceledException) when (CancellationToken.IsCancellationRequested)
            {
                Status = OptimizationStatus.Cancelled;
            }
            catch (ArgumentException ex)
            {
                if (CancellationToken.IsCancellationRequested) Status = OptimizationStatus.Cancelled;
                else if (ex.ParamName != nameof(MaxIterations) && ex.ParamName != nameof(MaxFunctionEvaluations))
                    UpdateStatus(OptimizationStatus.Failure, ex);
            }
            catch (Exception ex)
            {
                // Once cancellation is requested, a later failure is a consequence of the stop, not an error.
                if (CancellationToken.IsCancellationRequested) Status = OptimizationStatus.Cancelled;
                else UpdateStatus(OptimizationStatus.Failure, ex);
            }
            ReportProgress();
        }

        /// <summary>
        /// Evaluates the objective function and returns the fitness.
        /// </summary>
        /// <param name="values">The parameter values to evaluate.</param>
        /// <param name="cancel">By ref. Determines if the solver should be canceled.</param>
        /// <returns>The scaled fitness. A NaN objective value is returned as positive infinity, the worst fitness.</returns>
        protected virtual double Evaluate(double[] values, ref bool cancel)
        {
            // Stop on cancellation before spending another evaluation. The exception unwinds nested local
            // searches too; Minimize and Maximize turn it into the Cancelled status.
            if (CancellationToken.IsCancellationRequested)
            {
                cancel = true;
                Status = OptimizationStatus.Cancelled;
                throw new OperationCanceledException(CancellationToken);
            }

            return RecordEvaluation(values, ObjectiveFunction(values), ref cancel);
        }

        /// <summary>
        /// Records one objective value: scales it, keeps the best parameter set, updates the trace and the
        /// evaluation counter, and checks the evaluation budget.
        /// </summary>
        /// <param name="values">The evaluated parameter values.</param>
        /// <param name="objective">The unscaled objective value at <paramref name="values"/>.</param>
        /// <param name="cancel">By ref. Set when the evaluation budget is reached.</param>
        /// <returns>The scaled fitness. A NaN objective value is returned as positive infinity, the worst fitness.</returns>
        private double RecordEvaluation(double[] values, double objective, ref bool cancel)
        {
            double fitness = functionScale * objective;
            // A NaN fails every comparison, so it would never be replaced. Rank it as the worst value.
            if (double.IsNaN(fitness)) fitness = double.PositiveInfinity;

            // keep track of best fit parameter set
            if (BestParameterSet.Values == null || fitness <= BestParameterSet.Fitness)
                BestParameterSet = new ParameterSet((double[])values.Clone(), fitness);
            
            // Update trace. This is tracked every evaluation
            if (ParameterSetTrace == null) ParameterSetTrace = new List<ParameterSet>();
            if (RecordTraces) ParameterSetTrace.Add(BestParameterSet.Clone(deep: false));

            // update evaluation counter
            FunctionEvaluations += 1;
            if (Progress != null && FunctionEvaluations % ProgressInterval == 0) ReportProgress();
            if (FunctionEvaluations >= MaxFunctionEvaluations)
            {
                cancel = true;
                UpdateStatus(OptimizationStatus.MaximumFunctionEvaluationsReached);          
            }

            return fitness;
        }

        /// <summary>
        /// Evaluates a batch of points in parallel and records the results in index order.
        /// </summary>
        /// <param name="points">The points to evaluate; the first <paramref name="count"/> are used.</param>
        /// <param name="count">The number of points.</param>
        /// <param name="fitness">Receives the scaled fitness of each evaluated point.</param>
        /// <param name="maxDegreeOfParallelism">The maximum number of concurrent evaluations, or -1 for no limit.</param>
        /// <param name="cancel">By ref. Set when the evaluation budget is reached.</param>
        /// <returns>The number of points evaluated: <paramref name="count"/>, or fewer when the budget runs out.</returns>
        /// <remarks>
        /// The objective is called concurrently and must be thread safe. The results are recorded exactly as
        /// <see cref="Evaluate"/> would record them one after the other, in index order, so the best point,
        /// the trace and the counters do not depend on the number of threads. Only as many points as the
        /// remaining evaluation budget allows are evaluated.
        /// </remarks>
        protected int EvaluateBatch(double[][] points, int count, double[] fitness, int maxDegreeOfParallelism, ref bool cancel)
        {
            if (CancellationToken.IsCancellationRequested)
            {
                cancel = true;
                Status = OptimizationStatus.Cancelled;
                throw new OperationCanceledException(CancellationToken);
            }
            int n = Math.Min(count, MaxFunctionEvaluations - FunctionEvaluations);
            var options = new ParallelOptions { MaxDegreeOfParallelism = maxDegreeOfParallelism, CancellationToken = CancellationToken };
            // Each point's exception is kept, so that the first failing index, not the first thread to fail,
            // decides what is thrown: the earlier points are recorded first, as in a sequential run.
            var errors = new ExceptionDispatchInfo?[n];
            Parallel.For(0, n, options, k =>
            {
                try
                {
                    fitness[k] = ObjectiveFunction(points[k]);
                }
                catch (Exception ex)
                {
                    errors[k] = ExceptionDispatchInfo.Capture(ex);
                }
            });
            for (int k = 0; k < n; k++)
            {
                errors[k]?.Throw();
                fitness[k] = RecordEvaluation(points[k], fitness[k], ref cancel);
                if (cancel) return k + 1;
            }
            return n;
        }

        /// <summary>
        /// Reports the current state to <see cref="Progress"/>, when set.
        /// </summary>
        protected void ReportProgress()
        {
            if (Progress == null || BestParameterSet.Values == null) return;
            Progress.Report(new OptimizationProgress(Iterations, FunctionEvaluations, BestParameterSet.Fitness, (double[])BestParameterSet.Values.Clone()));
        }

        /// <summary>
        /// Throws when a local method cannot be used to polish a global search.
        /// </summary>
        /// <param name="method">The local method.</param>
        /// <exception cref="NotSupportedException"><paramref name="method"/> is not BFGS, Nelder-Mead, or Powell.</exception>
        protected static void ValidatePolishMethod(LocalMethod method)
        {
            if (method != LocalMethod.BFGS && method != LocalMethod.NelderMead && method != LocalMethod.Powell)
                throw new NotSupportedException("The polishing method must be BFGS, Nelder-Mead, or Powell.");
        }

        /// <summary>
        /// Refines the best parameter set with a bounded local search, the polish step of a global search.
        /// </summary>
        /// <param name="method">The local method: BFGS, Nelder-Mead, or Powell.</param>
        /// <returns>False when the evaluation budget is spent or the run is cancelled, so the caller must stop.</returns>
        /// <remarks>
        /// The local search starts at <see cref="BestParameterSet"/>, uses this optimizer's tolerances and its
        /// remaining evaluation budget, and evaluates through <see cref="Evaluate"/>, so every evaluation is
        /// counted, traced, and can improve <see cref="BestParameterSet"/>. It keeps no trace or Hessian of
        /// its own. The caller's status is not changed by the polish. A local search needs at least 10
        /// evaluations, so the polish is skipped when fewer remain in the budget.
        /// </remarks>
        protected bool PolishBestParameterSet(LocalMethod method)
        {
            ValidatePolishMethod(method);
            if (MaxFunctionEvaluations - FunctionEvaluations < 10) return FunctionEvaluations < MaxFunctionEvaluations && !CancellationToken.IsCancellationRequested;
            var lower = ParameterLowerBounds!;
            var upper = ParameterUpperBounds!;
            var start = new double[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
                start[i] = RepairParameter(BestParameterSet.Values[i], lower[i], upper[i]);

            bool localCancel = false;
            Func<double[], double> objective = x => Evaluate(x, ref localCancel);
            Optimizer solver = method == LocalMethod.BFGS ? new BFGS(objective, NumberOfParameters, start, lower, upper)
                : method == LocalMethod.NelderMead ? new NelderMead(objective, NumberOfParameters, start, lower, upper)
                : new Powell(objective, NumberOfParameters, start, lower, upper);
            solver.RelativeTolerance = RelativeTolerance;
            solver.AbsoluteTolerance = AbsoluteTolerance;
            solver.MaxFunctionEvaluations = MaxFunctionEvaluations - FunctionEvaluations;
            solver.RecordTraces = false;
            solver.ComputeHessian = false;
            // An exception from the objective propagates as it would from the global search; reaching the
            // local limits never throws and leaves the global search's status unchanged.
            solver.ReportFailure = ReportFailure;
            solver.Minimize();
            return !localCancel && FunctionEvaluations < MaxFunctionEvaluations && !CancellationToken.IsCancellationRequested;
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
        protected virtual void UpdateStatus(OptimizationStatus status, Exception? exception = null)
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
              if (ReportFailure && exception != null) ExceptionDispatchInfo.Capture(exception).Throw();
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
