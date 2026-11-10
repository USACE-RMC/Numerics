using System;
using System.Collections.Generic;
using System.Runtime.ExceptionServices;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// Writes the outcome distribution of a stochastic stage transition into caller-owned buffers.
    /// </summary>
    /// <param name="stage">The zero-based stage index.</param>
    /// <param name="state">The current state index.</param>
    /// <param name="decision">The decision index.</param>
    /// <param name="nextStates">The buffer that receives the next-state index of each outcome.</param>
    /// <param name="probabilities">The buffer that receives the probability of each outcome.</param>
    /// <returns>
    /// The number of outcomes written, at most the buffer length; zero or a negative value marks the
    /// decision infeasible.
    /// </returns>
    /// <remarks>
    /// The buffers are reused for every call, so the delegate runs with no allocation. Probabilities
    /// must be finite, non-negative and sum to one within 1E-8.
    /// </remarks>
    public delegate int TransitionDistribution(int stage, int state, int decision, int[] nextStates, double[] probabilities);

    /// <summary>
    /// The state and decision path of a forward simulation of a dynamic programming policy.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    public sealed class DynamicProgrammingPath
    {
        /// <summary>
        /// Constructs a path.
        /// </summary>
        /// <param name="states">The visited states, one per stage plus the terminal state.</param>
        /// <param name="decisions">The decision taken at each stage.</param>
        /// <param name="stageValues">The undiscounted stage cost or reward of each decision.</param>
        /// <param name="totalValue">The discounted total value of the path including the terminal value.</param>
        internal DynamicProgrammingPath(int[] states, int[] decisions, double[] stageValues, double totalValue)
        {
            States = states;
            Decisions = decisions;
            StageValues = stageValues;
            TotalValue = totalValue;
        }

        /// <summary>
        /// Gets the visited states: element t is the state at the start of stage t, and the last element is the terminal state.
        /// </summary>
        public int[] States { get; }

        /// <summary>
        /// Gets the decision taken at each stage.
        /// </summary>
        public int[] Decisions { get; }

        /// <summary>
        /// Gets the undiscounted stage cost or reward of each decision.
        /// </summary>
        public double[] StageValues { get; }

        /// <summary>
        /// Gets the discounted total value of the path, accumulated backward as
        /// v = c(t) + γ v from the terminal value, so a deterministic optimal path reproduces
        /// <see cref="DynamicProgramming.Value"/> at stage 0 exactly.
        /// </summary>
        public double TotalValue { get; }
    }

    /// <summary>
    /// Discrete dynamic programming over integer states and decisions: finite-horizon backward
    /// recursion with deterministic or stochastic (Markov) transitions, and value iteration for
    /// stationary infinite-horizon problems.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Stages are t = 0..T-1, states s = 0..S-1 and decisions a = 0..A-1. With stage cost or reward
    /// c(t, s, a), transition probabilities p(s' | t, s, a), discount factor γ in (0, 1] and terminal
    /// value V_T(s), <see cref="Solve"/> computes the Bellman recursion
    /// V_t(s) = opt_a [ c(t, s, a) + γ Σ p(s' | t, s, a) V_{t+1}(s') ], where opt is min or max.
    /// A non-finite stage value or a transition to an infeasible next state marks a decision infeasible;
    /// a state with no feasible decision has value +∞ when minimizing (-∞ when maximizing) and policy -1.
    /// Ties go to the lowest decision index. The work is O(T·S·A·K) for at most K outcomes per transition.
    /// </para>
    /// <para>
    /// For stochastic reservoir operation with Markov inflow classes, encode the state as
    /// s = storage · K + class (the current class is known when the release is chosen) and let the
    /// transition return one outcome per next inflow class. The delegates must be deterministic
    /// functions of their arguments.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Bellman, R. (1957). Dynamic Programming. Princeton University Press.
    /// </description></item>
    /// <item><description>
    /// Puterman, M. L. (1994). Markov Decision Processes: Discrete Stochastic Dynamic Programming. Wiley. https://doi.org/10.1002/9780470316887
    /// </description></item>
    /// <item><description>
    /// Labadie, J. W. (2004). Optimal operation of multireservoir systems: State-of-the-art review. Journal of Water Resources Planning and Management, 130(2), 93-111. https://doi.org/10.1061/(ASCE)0733-9496(2004)130:2(93)
    /// </description></item>
    /// <item><description>
    /// Loucks, D. P., and van Beek, E. (2017). Water Resource Systems Planning and Management. Springer. https://doi.org/10.1007/978-3-319-44234-1
    /// </description></item>
    /// </list>
    /// </remarks>
    public sealed class DynamicProgramming
    {
        /// <summary>The tolerance on the sum of the outcome probabilities of one transition.</summary>
        private const double PROBABILITY_SUM_TOLERANCE = 1E-8;

        private readonly Func<int, int, int, double> _stageValue;
        private readonly Func<int, int, int, int>? _transition;
        private readonly TransitionDistribution? _distribution;
        private double _discountFactor = 1.0;

        /// <summary>
        /// Constructs a deterministic dynamic programming problem.
        /// </summary>
        /// <param name="stages">The number of stages T, at least 1.</param>
        /// <param name="states">The number of states S, at least 1.</param>
        /// <param name="decisions">The number of decisions A, at least 1.</param>
        /// <param name="stageValue">The stage cost or reward c(t, s, a); a non-finite value marks the decision infeasible.</param>
        /// <param name="transition">The next state f(t, s, a); -1 marks the decision infeasible.</param>
        /// <exception cref="ArgumentOutOfRangeException">A dimension is less than 1.</exception>
        /// <exception cref="ArgumentNullException">A delegate is null.</exception>
        public DynamicProgramming(int stages, int states, int decisions, Func<int, int, int, double> stageValue, Func<int, int, int, int> transition)
        {
            ValidateDimensions(stages, states, decisions);
            _stageValue = stageValue ?? throw new ArgumentNullException(nameof(stageValue));
            _transition = transition ?? throw new ArgumentNullException(nameof(transition));
            Stages = stages;
            States = states;
            Decisions = decisions;
            MaxOutcomes = 1;
        }

        /// <summary>
        /// Constructs a stochastic dynamic programming problem.
        /// </summary>
        /// <param name="stages">The number of stages T, at least 1.</param>
        /// <param name="states">The number of states S, at least 1.</param>
        /// <param name="decisions">The number of decisions A, at least 1.</param>
        /// <param name="stageValue">The expected stage cost or reward c(t, s, a); a non-finite value marks the decision infeasible.</param>
        /// <param name="transition">Writes the outcome distribution of a transition into caller-owned buffers.</param>
        /// <param name="maxOutcomes">The buffer length: the largest number of outcomes of any transition, at least 1.</param>
        /// <exception cref="ArgumentOutOfRangeException">A dimension or <paramref name="maxOutcomes"/> is less than 1.</exception>
        /// <exception cref="ArgumentNullException">A delegate is null.</exception>
        public DynamicProgramming(int stages, int states, int decisions, Func<int, int, int, double> stageValue, TransitionDistribution transition, int maxOutcomes)
        {
            ValidateDimensions(stages, states, decisions);
            if (maxOutcomes < 1) throw new ArgumentOutOfRangeException(nameof(maxOutcomes), "The maximum number of outcomes must be at least 1.");
            _stageValue = stageValue ?? throw new ArgumentNullException(nameof(stageValue));
            _distribution = transition ?? throw new ArgumentNullException(nameof(transition));
            Stages = stages;
            States = states;
            Decisions = decisions;
            MaxOutcomes = maxOutcomes;
        }

        /// <summary>Gets the number of stages T.</summary>
        public int Stages { get; }

        /// <summary>Gets the number of states S.</summary>
        public int States { get; }

        /// <summary>Gets the number of decisions A.</summary>
        public int Decisions { get; }

        /// <summary>Gets the largest number of outcomes of a transition (1 for a deterministic problem).</summary>
        public int MaxOutcomes { get; }

        /// <summary>Gets a value indicating whether the transitions are stochastic.</summary>
        public bool IsStochastic => _distribution != null;

        /// <summary>Gets or sets a value indicating whether the stage values are rewards to maximize. The default is false (minimize costs).</summary>
        public bool Maximize { get; set; }

        /// <summary>
        /// Gets or sets the discount factor γ in (0, 1] applied to the value of the next stage. The default is 1.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">The value is not in (0, 1].</exception>
        public double DiscountFactor
        {
            get => _discountFactor;
            set
            {
                if (!(value > 0 && value <= 1)) throw new ArgumentOutOfRangeException(nameof(DiscountFactor), "The discount factor must be in (0, 1].");
                _discountFactor = value;
            }
        }

        /// <summary>
        /// Gets or sets the terminal value V_T(s), one per state; null means zero. A non-finite entry marks the
        /// terminal state infeasible. For <see cref="SolveSteadyState"/> it is the initial value estimate.
        /// </summary>
        public double[]? TerminalValue { get; set; }

        /// <summary>Gets a value indicating whether <see cref="Solve"/> has completed.</summary>
        public bool IsSolved { get; private set; }

        /// <summary>
        /// Gets the optimal value table V[t, s] for t = 0..T, filled by <see cref="Solve"/>; row T is the terminal value.
        /// </summary>
        public double[,] Value { get; private set; } = new double[0, 0];

        /// <summary>
        /// Gets the optimal policy table [t, s] for t = 0..T-1, filled by <see cref="Solve"/>; -1 where no decision is feasible.
        /// </summary>
        public int[,] Policy { get; private set; } = new int[0, 0];

        /// <summary>Gets the stationary value of each state from the last call to <see cref="SolveSteadyState"/>.</summary>
        public double[] SteadyStateValue { get; private set; } = Array.Empty<double>();

        /// <summary>Gets the stationary policy of each state from the last call to <see cref="SolveSteadyState"/>; -1 where no decision is feasible.</summary>
        public int[] SteadyStatePolicy { get; private set; } = Array.Empty<int>();

        /// <summary>Gets the number of value-iteration sweeps of the last call to <see cref="SolveSteadyState"/>.</summary>
        public int SteadyStateIterations { get; private set; }

        /// <summary>Gets a value indicating whether the last call to <see cref="SolveSteadyState"/> converged.</summary>
        public bool SteadyStateConverged { get; private set; }

        /// <summary>
        /// Solves the finite-horizon problem by backward recursion and fills <see cref="Value"/> and <see cref="Policy"/>.
        /// </summary>
        /// <exception cref="ArgumentException"><see cref="TerminalValue"/> does not have one entry per state.</exception>
        /// <exception cref="InvalidOperationException">A transition returned an invalid state, outcome count or probability.</exception>
        public void Solve()
        {
            int T = Stages, S = States, A = Decisions;
            var value = new double[T + 1, S];
            var policy = new int[T, S];
            var next = InitialValue();
            for (int s = 0; s < S; s++) value[T, s] = next[s];

            var current = new double[S];
            for (int t = T - 1; t >= 0; t--)
            {
                int stage = t;
                var nextValue = next;
                var currentValue = current;
                if (UseParallel && S > 1)
                {
                    // Each state's exception is kept, so that the lowest failing state, not the first thread to
                    // fail, decides what is thrown: the same exception as the serial solve.
                    var errors = new ExceptionDispatchInfo?[S];
                    Parallel.For(0, S,
                        () => new Buffers(MaxOutcomes),
                        (s, loop, buffers) =>
                        {
                            try
                            {
                                currentValue[s] = SolveState(stage, s, nextValue, buffers, out int decision);
                                policy[stage, s] = decision;
                            }
                            catch (Exception ex)
                            {
                                errors[s] = ExceptionDispatchInfo.Capture(ex);
                            }
                            return buffers;
                        },
                        buffers => { });
                    for (int s = 0; s < S; s++) errors[s]?.Throw();
                }
                else
                {
                    var buffers = new Buffers(MaxOutcomes);
                    for (int s = 0; s < S; s++)
                    {
                        currentValue[s] = SolveState(stage, s, nextValue, buffers, out int decision);
                        policy[stage, s] = decision;
                    }
                }
                for (int s = 0; s < S; s++) value[t, s] = currentValue[s];
                next = currentValue; current = nextValue;
            }

            Value = value;
            Policy = policy;
            IsSolved = true;
        }

        /// <summary>
        /// Gets or sets a value indicating whether <see cref="Solve"/> evaluates the states of each stage in parallel.
        /// The default is false.
        /// </summary>
        /// <remarks>
        /// States within a stage are independent, so the result is identical to the serial solve. The stage-value and
        /// transition delegates must then be safe to call concurrently; each thread receives its own outcome buffers.
        /// An exception thrown during a parallel stage is rethrown as it is, from the lowest-numbered state that
        /// failed, so the parallel solve throws the same exception as the serial solve.
        /// </remarks>
        public bool UseParallel { get; set; }

        /// <summary>
        /// Per-thread outcome buffers for the stochastic transition.
        /// </summary>
        /// <param name="length">The buffer length.</param>
        private sealed class Buffers(int length)
        {
            /// <summary>The next-state buffer.</summary>
            public readonly int[] NextStates = new int[length];

            /// <summary>The probability buffer.</summary>
            public readonly double[] Probabilities = new double[length];
        }

        /// <summary>
        /// Returns the optimal value of one state at one stage and its best decision (-1 if none is feasible).
        /// </summary>
        /// <param name="t">The stage.</param>
        /// <param name="s">The state.</param>
        /// <param name="next">The next-stage value.</param>
        /// <param name="buffers">The outcome buffers.</param>
        /// <param name="bestDecision">The best decision, or -1.</param>
        private double SolveState(int t, int s, double[] next, Buffers buffers, out int bestDecision)
        {
            bool maximize = Maximize;
            double gamma = _discountFactor;
            double best = Worst;
            bestDecision = -1;
            int A = Decisions;
            var transition = _transition;
            for (int a = 0; a < A; a++)
            {
                double c = _stageValue(t, s, a);
                if (!Tools.IsFinite(c)) continue;
                double expected;
                if (transition != null)
                {
                    int ns = transition(t, s, a);
                    if (ns == -1) continue;
                    if ((uint)ns >= (uint)States) throw new InvalidOperationException($"The transition from state {s} with decision {a} at stage {t} returned state {ns}, outside 0..{States - 1}.");
                    expected = next[ns];
                }
                else
                {
                    expected = Expected(t, s, a, next, buffers.NextStates, buffers.Probabilities);
                }
                double candidate = c + gamma * expected;
                if (!Tools.IsFinite(candidate)) continue;
                if (bestDecision < 0 || (maximize ? candidate > best : candidate < best))
                {
                    best = candidate;
                    bestDecision = a;
                }
            }
            return best;
        }

        /// <summary>
        /// Simulates the optimal policy of a deterministic problem forward from an initial state.
        /// </summary>
        /// <param name="initialState">The state at stage 0.</param>
        /// <returns>The optimal state and decision path.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="initialState"/> is not a valid state.</exception>
        /// <exception cref="InvalidOperationException">The problem is stochastic or unsolved, or the path reaches a state with no feasible decision.</exception>
        public DynamicProgrammingPath Simulate(int initialState)
        {
            if (IsStochastic) throw new InvalidOperationException("A stochastic problem needs a random number generator or an outcome sequence to simulate.");
            return SimulateCore(initialState, null, null);
        }

        /// <summary>
        /// Simulates the optimal policy forward from an initial state, sampling each transition outcome with a random number generator.
        /// </summary>
        /// <param name="initialState">The state at stage 0.</param>
        /// <param name="random">The random number generator; one uniform draw is used per stage.</param>
        /// <returns>The simulated state and decision path.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="random"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="initialState"/> is not a valid state.</exception>
        /// <exception cref="InvalidOperationException">The problem is unsolved, or the path reaches a state with no feasible decision.</exception>
        public DynamicProgrammingPath Simulate(int initialState, Random random)
        {
            if (random == null) throw new ArgumentNullException(nameof(random));
            return SimulateCore(initialState, random, null);
        }

        /// <summary>
        /// Simulates the optimal policy forward from an initial state along a given sequence of outcome indices,
        /// such as an observed sequence of inflow classes.
        /// </summary>
        /// <param name="initialState">The state at stage 0.</param>
        /// <param name="outcomes">The outcome index at each stage, an index into the outcomes the transition writes (0 for a deterministic problem).</param>
        /// <returns>The state and decision path.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="outcomes"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="outcomes"/> does not have one entry per stage.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="initialState"/> or an outcome index is out of range.</exception>
        /// <exception cref="InvalidOperationException">The problem is unsolved, or the path reaches a state with no feasible decision.</exception>
        public DynamicProgrammingPath Simulate(int initialState, IList<int> outcomes)
        {
            if (outcomes == null) throw new ArgumentNullException(nameof(outcomes));
            if (outcomes.Count != Stages) throw new ArgumentException("The outcome sequence must have one entry per stage.", nameof(outcomes));
            return SimulateCore(initialState, null, outcomes);
        }

        /// <summary>
        /// Solves the stationary infinite-horizon problem by value iteration and fills <see cref="SteadyStateValue"/>,
        /// <see cref="SteadyStatePolicy"/>, <see cref="SteadyStateIterations"/> and <see cref="SteadyStateConverged"/>.
        /// </summary>
        /// <param name="tolerance">The convergence tolerance on max_s |V_{k+1}(s) - V_k(s)|. The default is 1E-8.</param>
        /// <param name="maxIterations">The maximum number of sweeps. The default is 10,000.</param>
        /// <returns>True if the iteration converged.</returns>
        /// <remarks>
        /// The delegates are evaluated once at stage 0 and cached. Starting from <see cref="TerminalValue"/> (zero if null),
        /// each sweep applies V_{k+1}(s) = opt_a [ c(0, s, a) + γ Σ p(s' | 0, s, a) V_k(s') ]; the reported policy is the
        /// greedy decision of the last sweep. With γ &lt; 1 the error of the last iterate is at most γ/(1 - γ) times the last
        /// change. With γ = 1 the iteration converges only when values are bounded, as in a stochastic shortest-path problem.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="tolerance"/> is not positive and finite, or <paramref name="maxIterations"/> is less than 1.</exception>
        /// <exception cref="ArgumentException"><see cref="TerminalValue"/> does not have one entry per state.</exception>
        /// <exception cref="InvalidOperationException">A transition returned an invalid state, outcome count or probability.</exception>
        public bool SolveSteadyState(double tolerance = 1E-8, int maxIterations = 10000)
        {
            if (!(tolerance > 0) || double.IsInfinity(tolerance)) throw new ArgumentOutOfRangeException(nameof(tolerance), "The tolerance must be positive and finite.");
            if (maxIterations < 1) throw new ArgumentOutOfRangeException(nameof(maxIterations), "The maximum number of iterations must be at least 1.");

            int S = States, A = Decisions, K = MaxOutcomes;
            // Cache the stationary model: stage values, outcome counts, next states and probabilities.
            var cost = new double[S * A];
            var count = new int[S * A];
            var nextStates = new int[S * A * K];
            var probabilities = new double[S * A * K];
            var bufferStates = new int[K];
            var bufferProbabilities = new double[K];
            for (int s = 0; s < S; s++)
            {
                for (int a = 0; a < A; a++)
                {
                    int sa = s * A + a;
                    double c = _stageValue(0, s, a);
                    cost[sa] = c;
                    if (!Tools.IsFinite(c)) continue;
                    int n = Outcomes(0, s, a, bufferStates, bufferProbabilities);
                    count[sa] = n;
                    Array.Copy(bufferStates, 0, nextStates, sa * K, n);
                    Array.Copy(bufferProbabilities, 0, probabilities, sa * K, n);
                }
            }

            var v = InitialValue();
            var vNew = new double[S];
            var policy = new int[S];
            bool maximize = Maximize;
            double gamma = _discountFactor;
            double worst = Worst;
            bool converged = false;
            int iteration = 0;

            while (iteration < maxIterations)
            {
                iteration++;
                double delta = 0;
                for (int s = 0; s < S; s++)
                {
                    double best = worst;
                    int bestDecision = -1;
                    for (int a = 0; a < A; a++)
                    {
                        int sa = s * A + a;
                        int n = count[sa];
                        if (n <= 0) continue;
                        double expected = 0;
                        int offset = sa * K;
                        for (int k = 0; k < n; k++)
                        {
                            double p = probabilities[offset + k];
                            if (p == 0) continue;
                            expected += p * v[nextStates[offset + k]];
                        }
                        double candidate = cost[sa] + gamma * expected;
                        if (!Tools.IsFinite(candidate)) continue;
                        if (bestDecision < 0 || (maximize ? candidate > best : candidate < best))
                        {
                            best = candidate;
                            bestDecision = a;
                        }
                    }
                    vNew[s] = best;
                    policy[s] = bestDecision;
                    double old = v[s];
                    if (best != old)
                    {
                        double d = Math.Abs(best - old);
                        if (!(d <= delta)) delta = double.IsNaN(d) ? double.PositiveInfinity : d;
                    }
                }
                var swap = v; v = vNew; vNew = swap;
                if (delta <= tolerance) { converged = true; break; }
            }

            SteadyStateValue = v;
            SteadyStatePolicy = policy;
            SteadyStateIterations = iteration;
            SteadyStateConverged = converged;
            return converged;
        }

        /// <summary>Gets the value of a state with no feasible decision.</summary>
        private double Worst => Maximize ? double.NegativeInfinity : double.PositiveInfinity;

        /// <summary>
        /// Validates the problem dimensions.
        /// </summary>
        /// <param name="stages">The number of stages.</param>
        /// <param name="states">The number of states.</param>
        /// <param name="decisions">The number of decisions.</param>
        private static void ValidateDimensions(int stages, int states, int decisions)
        {
            if (stages < 1) throw new ArgumentOutOfRangeException(nameof(stages), "The number of stages must be at least 1.");
            if (states < 1) throw new ArgumentOutOfRangeException(nameof(states), "The number of states must be at least 1.");
            if (decisions < 1) throw new ArgumentOutOfRangeException(nameof(decisions), "The number of decisions must be at least 1.");
        }

        /// <summary>
        /// Returns a copy of the terminal value with non-finite entries replaced by the infeasible value.
        /// </summary>
        private double[] InitialValue()
        {
            var v = new double[States];
            var terminal = TerminalValue;
            if (terminal == null) return v;
            if (terminal.Length != States) throw new ArgumentException("The terminal value must have one entry per state.", nameof(TerminalValue));
            double worst = Worst;
            for (int s = 0; s < States; s++) v[s] = Tools.IsFinite(terminal[s]) ? terminal[s] : worst;
            return v;
        }

        /// <summary>
        /// Writes the validated outcomes of a transition into the buffers.
        /// </summary>
        /// <param name="t">The stage.</param>
        /// <param name="s">The state.</param>
        /// <param name="a">The decision.</param>
        /// <param name="nextStates">The next-state buffer.</param>
        /// <param name="probabilities">The probability buffer.</param>
        /// <returns>The number of outcomes, or 0 if the decision is infeasible.</returns>
        /// <exception cref="InvalidOperationException">The transition returned an invalid state, outcome count or probability.</exception>
        private int Outcomes(int t, int s, int a, int[] nextStates, double[] probabilities)
        {
            if (_transition != null)
            {
                int ns = _transition(t, s, a);
                if (ns == -1) return 0;
                if (ns < 0 || ns >= States) throw new InvalidOperationException($"The transition from state {s} with decision {a} at stage {t} returned state {ns}, outside 0..{States - 1}.");
                nextStates[0] = ns;
                probabilities[0] = 1.0;
                return 1;
            }

            int n = _distribution!(t, s, a, nextStates, probabilities);
            if (n <= 0) return 0;
            if (n > MaxOutcomes) throw new InvalidOperationException($"The transition from state {s} with decision {a} at stage {t} wrote {n} outcomes, more than the maximum of {MaxOutcomes}.");
            double sum = 0;
            for (int k = 0; k < n; k++)
            {
                double p = probabilities[k];
                int ns = nextStates[k];
                if (!(p >= 0) || double.IsInfinity(p)) throw new InvalidOperationException($"The transition from state {s} with decision {a} at stage {t} returned probability {p}.");
                if (ns < 0 || ns >= States) throw new InvalidOperationException($"The transition from state {s} with decision {a} at stage {t} returned state {ns}, outside 0..{States - 1}.");
                sum += p;
            }
            if (Math.Abs(sum - 1.0) > PROBABILITY_SUM_TOLERANCE) throw new InvalidOperationException($"The outcome probabilities of the transition from state {s} with decision {a} at stage {t} sum to {sum}, not 1.");
            return n;
        }

        /// <summary>
        /// Returns the expected next-stage value of a decision, or NaN if the decision is infeasible.
        /// </summary>
        /// <param name="t">The stage.</param>
        /// <param name="s">The state.</param>
        /// <param name="a">The decision.</param>
        /// <param name="next">The next-stage value.</param>
        /// <param name="nextStates">The next-state buffer.</param>
        /// <param name="probabilities">The probability buffer.</param>
        private double Expected(int t, int s, int a, double[] next, int[] nextStates, double[] probabilities)
        {
            int n = Outcomes(t, s, a, nextStates, probabilities);
            if (n == 0) return double.NaN;
            if (_transition != null) return next[nextStates[0]];
            double expected = 0;
            for (int k = 0; k < n; k++)
            {
                double p = probabilities[k];
                if (p == 0) continue;
                expected += p * next[nextStates[k]];
            }
            return expected;
        }

        /// <summary>
        /// Simulates the solved policy forward.
        /// </summary>
        /// <param name="initialState">The state at stage 0.</param>
        /// <param name="random">The generator used to sample outcomes, or null.</param>
        /// <param name="outcomes">The outcome sequence, or null.</param>
        private DynamicProgrammingPath SimulateCore(int initialState, Random? random, IList<int>? outcomes)
        {
            if (!IsSolved) throw new InvalidOperationException("Call Solve before Simulate.");
            if (initialState < 0 || initialState >= States) throw new ArgumentOutOfRangeException(nameof(initialState), $"The initial state must be in 0..{States - 1}.");

            int T = Stages;
            var states = new int[T + 1];
            var decisions = new int[T];
            var stageValues = new double[T];
            var nextStates = new int[MaxOutcomes];
            var probabilities = new double[MaxOutcomes];
            int s = initialState;
            for (int t = 0; t < T; t++)
            {
                states[t] = s;
                int a = Policy[t, s];
                if (a < 0) throw new InvalidOperationException($"State {s} has no feasible decision at stage {t}.");
                decisions[t] = a;
                stageValues[t] = _stageValue(t, s, a);
                int n = Outcomes(t, s, a, nextStates, probabilities);
                int k;
                if (outcomes != null)
                {
                    k = outcomes[t];
                    if (k < 0 || k >= n) throw new ArgumentOutOfRangeException(nameof(outcomes), $"The outcome index {k} at stage {t} is outside 0..{n - 1}.");
                }
                else if (random != null)
                {
                    double u = random.NextDouble();
                    double cumulative = 0;
                    k = -1;
                    for (int j = 0; j < n; j++)
                    {
                        if (probabilities[j] == 0) continue;
                        k = j;
                        cumulative += probabilities[j];
                        if (u < cumulative) break;
                    }
                }
                else
                {
                    k = 0;
                }
                s = nextStates[k];
            }
            states[T] = s;

            double total = Value[T, s];
            for (int t = T - 1; t >= 0; t--) total = stageValues[t] + _discountFactor * total;
            return new DynamicProgrammingPath(states, decisions, stageValues, total);
        }
    }
}
