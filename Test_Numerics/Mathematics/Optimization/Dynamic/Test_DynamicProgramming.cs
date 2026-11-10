using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Tests the discrete dynamic programming solver against brute-force enumeration, a hand-solved
    /// reservoir problem, a closed-form Markov decision problem and an independent numpy oracle.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_DynamicProgramming
    {
        #region Fixtures

        /// <summary>
        /// A random deterministic problem: stage costs (about 15% infeasible) and next-state tables.
        /// </summary>
        private sealed class RandomDeterministic
        {
            public readonly double[,,] Cost;
            public readonly int[,,] Next;
            public readonly double[] Terminal;

            public RandomDeterministic(int seed, int T, int S, int A)
            {
                var rng = new MersenneTwister(seed);
                Cost = new double[T, S, A];
                Next = new int[T, S, A];
                Terminal = new double[S];
                for (int t = 0; t < T; t++)
                    for (int s = 0; s < S; s++)
                        for (int a = 0; a < A; a++)
                        {
                            double u = rng.NextDouble();
                            Cost[t, s, a] = u < 0.1 ? double.NaN : u < 0.15 ? double.PositiveInfinity : rng.NextDouble() * 10 - 3;
                            Next[t, s, a] = rng.NextDouble() < 0.05 ? -1 : rng.Next(0, S);
                        }
                for (int s = 0; s < S; s++) Terminal[s] = rng.NextDouble() < 0.1 ? double.PositiveInfinity : rng.NextDouble() * 5;
            }
        }

        /// <summary>
        /// A random stochastic problem with K outcomes per transition and some zero probabilities.
        /// </summary>
        private sealed class RandomStochastic
        {
            public readonly double[,,] Cost;
            public readonly int[,,,] Next;
            public readonly double[,,,] Prob;
            public readonly double[] Terminal;
            public readonly int K;

            public RandomStochastic(int seed, int T, int S, int A, int K)
            {
                var rng = new MersenneTwister(seed);
                this.K = K;
                Cost = new double[T, S, A];
                Next = new int[T, S, A, K];
                Prob = new double[T, S, A, K];
                Terminal = new double[S];
                for (int t = 0; t < T; t++)
                    for (int s = 0; s < S; s++)
                        for (int a = 0; a < A; a++)
                        {
                            Cost[t, s, a] = rng.NextDouble() < 0.1 ? double.NaN : rng.NextDouble() * 10 - 3;
                            double total = 0;
                            for (int k = 0; k < K; k++)
                            {
                                Next[t, s, a, k] = rng.Next(0, S);
                                double w = rng.NextDouble() < 0.2 ? 0 : rng.NextDouble();
                                Prob[t, s, a, k] = w;
                                total += w;
                            }
                            if (total == 0) { Prob[t, s, a, 0] = 1; total = 1; }
                            for (int k = 0; k < K; k++) Prob[t, s, a, k] /= total;
                        }
                for (int s = 0; s < S; s++) Terminal[s] = rng.NextDouble() * 5;
            }

            public int Transition(int t, int s, int a, int[] next, double[] prob)
            {
                for (int k = 0; k < K; k++) { next[k] = Next[t, s, a, k]; prob[k] = Prob[t, s, a, k]; }
                return K;
            }
        }

        /// <summary>
        /// Brute-force best value from (t, s): enumerates every decision sequence, accumulating backward.
        /// Returns +∞ (or -∞) if no feasible sequence exists.
        /// </summary>
        private static double BruteForce(RandomDeterministic p, int t, int s, int T, int A, double gamma, bool maximize, out int firstDecision)
        {
            firstDecision = -1;
            double worst = maximize ? double.NegativeInfinity : double.PositiveInfinity;
            if (t == T) return double.IsInfinity(p.Terminal[s]) ? worst : p.Terminal[s];
            double best = worst;
            for (int a = 0; a < A; a++)
            {
                double c = p.Cost[t, s, a];
                int ns = p.Next[t, s, a];
                if (double.IsNaN(c) || double.IsInfinity(c) || ns < 0) continue;
                double rest = BruteForce(p, t + 1, ns, T, A, gamma, maximize, out _);
                if (double.IsInfinity(rest)) continue;
                double v = c + gamma * rest;
                if (firstDecision < 0 || (maximize ? v > best : v < best)) { best = v; firstDecision = a; }
            }
            return best;
        }

        /// <summary>
        /// Evaluates a Markov policy of a stochastic problem by forward propagation of the state distribution.
        /// </summary>
        private static double EvaluatePolicy(RandomStochastic p, int[,] policy, int s0, int T, int S, double gamma)
        {
            var mass = new double[S];
            mass[s0] = 1;
            double total = 0, discount = 1;
            for (int t = 0; t < T; t++)
            {
                var nextMass = new double[S];
                for (int s = 0; s < S; s++)
                {
                    if (mass[s] == 0) continue;
                    int a = policy[t, s];
                    double c = p.Cost[t, s, a];
                    if (double.IsNaN(c)) return double.PositiveInfinity;
                    total += discount * mass[s] * c;
                    for (int k = 0; k < p.K; k++) nextMass[p.Next[t, s, a, k]] += mass[s] * p.Prob[t, s, a, k];
                }
                mass = nextMass;
                discount *= gamma;
            }
            for (int s = 0; s < S; s++) total += discount * mass[s] * p.Terminal[s];
            return total;
        }

        // Reservoir model shared with docs/mathematics/oracles/dynamic-programming/dynamic_programming_oracle.py.
        private static int Inflow(int t, int A) => (3 + 5 * t) % A;
        private static int Demand(int t, int A) => (1 + 3 * t) % A;

        private static double ReservoirCost(int storage, int q, int r, int D, int S)
        {
            int x = storage + q - r;
            if (x < 0) return double.NaN;
            int spill = Math.Max(x - (S - 1), 0);
            int nxt = Math.Min(x, S - 1);
            double d = r - D;
            return d * d + 0.5 * spill + 0.1 * (S - 1 - nxt);
        }

        private static DynamicProgramming OracleDeterministic(int T, int S, int A)
        {
            var dp = new DynamicProgramming(T, S, A,
                (t, s, r) => ReservoirCost(s, Inflow(t, A), r, Demand(t, A), S),
                (t, s, r) => { int x = s + Inflow(t, A) - r; return x < 0 ? -1 : Math.Min(x, S - 1); });
            var terminal = new double[S];
            for (int s = 0; s < S; s++) terminal[s] = 0.25 * (S - 1 - s);
            dp.TerminalValue = terminal;
            return dp;
        }

        private static DynamicProgramming OracleStochastic(int T, int S, int A, int K, int qStep)
        {
            var P = new double[K, K];
            for (int j = 0; j < K; j++)
            {
                var w = new double[K];
                double total = 0;
                for (int jn = 0; jn < K; jn++) { w[jn] = 1.0 / (1.0 + Math.Abs(j - jn)); total += w[jn]; }
                for (int jn = 0; jn < K; jn++) P[j, jn] = w[jn] / total;
            }
            var dp = new DynamicProgramming(T, S * K, A,
                (t, s, r) => ReservoirCost(s / K, (s % K) * qStep, r, Demand(t, A), S),
                (t, s, r, next, prob) =>
                {
                    int j = s % K;
                    int x = s / K + j * qStep - r;
                    if (x < 0) return 0;
                    int nxt = Math.Min(x, S - 1);
                    for (int jn = 0; jn < K; jn++) { next[jn] = nxt * K + jn; prob[jn] = P[j, jn]; }
                    return K;
                }, K);
            var terminal = new double[S * K];
            for (int s = 0; s < S * K; s++) terminal[s] = 0.25 * (S - 1 - s / K);
            dp.TerminalValue = terminal;
            return dp;
        }

        #endregion

        /// <summary>
        /// Value and policy of random deterministic problems equal brute-force enumeration exactly, for every
        /// stage and state, minimizing and maximizing, with and without discounting.
        /// </summary>
        [TestMethod]
        public void DeterministicMatchesBruteForce()
        {
            int T = 5, S = 4, A = 3;
            for (int seed = 1; seed <= 8; seed++)
            {
                var p = new RandomDeterministic(seed, T, S, A);
                foreach (bool maximize in new[] { false, true })
                {
                    foreach (double gamma in new[] { 1.0, 0.9 })
                    {
                        var dp = new DynamicProgramming(T, S, A, (t, s, a) => p.Cost[t, s, a], (t, s, a) => p.Next[t, s, a])
                        {
                            Maximize = maximize,
                            DiscountFactor = gamma,
                            TerminalValue = p.Terminal
                        };
                        dp.Solve();
                        for (int t = 0; t < T; t++)
                            for (int s = 0; s < S; s++)
                            {
                                double expected = BruteForce(p, t, s, T, A, gamma, maximize, out int first);
                                Assert.AreEqual(expected, dp.Value[t, s], $"seed {seed}, t {t}, s {s}");
                                Assert.AreEqual(first, dp.Policy[t, s], $"seed {seed}, t {t}, s {s}");
                            }
                    }
                }
            }
        }

        /// <summary>
        /// Value of random stochastic problems equals the best of all Markov policies (each evaluated by forward
        /// propagation) to 1E-12, and the stage-0 decision is that of the best policy.
        /// </summary>
        [TestMethod]
        public void StochasticMatchesPolicyEnumeration()
        {
            int T = 3, S = 3, A = 2, K = 3;
            int cells = T * S;
            for (int seed = 11; seed <= 16; seed++)
            {
                var p = new RandomStochastic(seed, T, S, A, K);
                var dp = new DynamicProgramming(T, S, A, (t, s, a) => p.Cost[t, s, a], p.Transition, K)
                {
                    DiscountFactor = 0.95,
                    TerminalValue = p.Terminal
                };
                dp.Solve();

                var policy = new int[T, S];
                for (int s0 = 0; s0 < S; s0++)
                {
                    double best = double.PositiveInfinity;
                    int bestFirst = -1;
                    int count = 1 << cells;
                    for (int code = 0; code < count; code++)
                    {
                        for (int c = 0; c < cells; c++) policy[c / S, c % S] = (code >> c) & 1;
                        double v = EvaluatePolicy(p, policy, s0, T, S, 0.95);
                        if (v < best - 1E-13) { best = v; bestFirst = policy[0, s0]; }
                    }
                    if (double.IsInfinity(best))
                    {
                        Assert.IsTrue(double.IsPositiveInfinity(dp.Value[0, s0]));
                        Assert.AreEqual(-1, dp.Policy[0, s0]);
                    }
                    else
                    {
                        Assert.AreEqual(best, dp.Value[0, s0], 1E-12, $"seed {seed}, s0 {s0}");
                        Assert.AreEqual(bestFirst, dp.Policy[0, s0], $"seed {seed}, s0 {s0}");
                    }
                }
            }
        }

        /// <summary>
        /// A hand-solved reservoir allocation: capacity 3, initial storage 1, inflows (2, 1, 1), concave benefit
        /// B(r) = (0, 5, 8, 9), no spill, at least one unit left at the end. Total release is at most 4, and the
        /// concave benefit spreads it as (2, 1, 1) in some order, worth 18. The three orders tie, so the lowest
        /// decision at each stage gives releases (1, 1, 2) through storages (1, 2, 2, 1).
        /// </summary>
        [TestMethod]
        public void ReservoirAllocationHandSolved()
        {
            int[] inflow = { 2, 1, 1 };
            double[] benefit = { 0, 5, 8, 9 };
            var dp = new DynamicProgramming(3, 4, 4,
                (t, s, r) => benefit[r],
                (t, s, r) => { int x = s + inflow[t] - r; return x < 0 || x > 3 ? -1 : x; })
            {
                Maximize = true,
                TerminalValue = new[] { double.NegativeInfinity, 0, 0, 0 }
            };
            dp.Solve();

            Assert.AreEqual(18.0, dp.Value[0, 1]);
            var path = dp.Simulate(1);
            CollectionAssert.AreEqual(new[] { 1, 2, 2, 1 }, path.States);
            CollectionAssert.AreEqual(new[] { 1, 1, 2 }, path.Decisions);
            CollectionAssert.AreEqual(new[] { 5.0, 5.0, 8.0 }, path.StageValues);
            Assert.AreEqual(18.0, path.TotalValue);
            // Empty reservoir: water 0 + 4, one unit kept, so 3 released as (1, 1, 1).
            Assert.AreEqual(15.0, dp.Value[0, 0]);
            // A terminal state that is infeasible is reported as -∞ when maximizing.
            Assert.IsTrue(double.IsNegativeInfinity(dp.Value[3, 0]));
        }

        /// <summary>
        /// The deterministic and stochastic reservoir fixtures match the independent numpy oracle
        /// (dynamic-programming-golden.csv) exactly, value at stage 0 and policy at every stage.
        /// </summary>
        [TestMethod]
        public void ReservoirMatchesNumpyOracle()
        {
            double[] detValue = { 10.864375465307155, 9.721750465307153, 8.624250465307153, 7.56988061871677, 6.515510772126383, 5.461140925535999, 4.406771078945614, 3.356771078945614, 2.356771078945614, 2.856771078945614, 3.356771078945614 };
            int[,] detPolicy =
            {
                { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 },
                { 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 },
                { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 },
                { 0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4 },
                { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                { 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4 },
                { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 },
                { 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4 },
                { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 },
                { 0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4 },
                { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 },
            };
            double[] stoValue = { 30.451012215610056, 19.080895874554486, 13.391419537017438, 25.62737328121764, 16.830889743843976, 11.91549076610738, 21.991014200643924, 14.9130810245639, 10.615751431735973, 19.333379765966335, 13.214065326737154, 9.50362660250607, 17.02673411388698, 11.766984869974474, 8.50362660250607, 14.975443675304094, 10.521788981820139, 7.731656635383914, 13.297173037878935, 9.507944966993776, 7.058195224957174, 11.834978596232848, 8.507944966993776, 6.643554810864019, 10.612129611132856, 7.584690530575265, 7.143554810864019, 9.549355423119959, 6.916592322789108, 7.643554810864019, 8.549355423119959, 7.416592322789108, 8.143554810864018 };
            int[,] stoPolicy =
            {
                { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1 },
                { 0, 2, 3, 1, 3, 3, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4 },
                { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1 },
                { 0, 2, 3, 1, 3, 3, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4 },
                { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1 },
                { 0, 2, 3, 1, 3, 3, 2, 3, 3, 2, 3, 3, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4 },
                { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                { 0, 2, 3, 1, 3, 3, 2, 3, 3, 2, 3, 3, 3, 3, 4, 3, 3, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4 },
                { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                { 0, 2, 3, 1, 3, 3, 2, 3, 4, 3, 3, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 },
                { 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                { 0, 2, 4, 1, 3, 4, 2, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 },
            };

            var det = OracleDeterministic(12, 11, 6);
            det.DiscountFactor = 0.95;
            det.Solve();
            for (int s = 0; s < 11; s++) Assert.AreEqual(detValue[s], det.Value[0, s], 1E-12, $"det s {s}");
            for (int t = 0; t < 12; t++)
                for (int s = 0; s < 11; s++) Assert.AreEqual(detPolicy[t, s], det.Policy[t, s], $"det t {t}, s {s}");

            var sto = OracleStochastic(12, 11, 6, 3, 2);
            sto.DiscountFactor = 0.95;
            sto.Solve();
            for (int s = 0; s < 33; s++) Assert.AreEqual(stoValue[s], sto.Value[0, s], 1E-12, $"sto s {s}");
            for (int t = 0; t < 12; t++)
                for (int s = 0; s < 33; s++) Assert.AreEqual(stoPolicy[t, s], sto.Policy[t, s], $"sto t {t}, s {s}");
        }

        /// <summary>
        /// The parallel stage sweep gives value and policy tables identical to the serial solve.
        /// </summary>
        [TestMethod]
        public void ParallelSolveIsIdentical()
        {
            var serial = OracleStochastic(30, 25, 8, 4, 3);
            serial.DiscountFactor = 0.97;
            serial.Solve();
            var parallel = OracleStochastic(30, 25, 8, 4, 3);
            parallel.DiscountFactor = 0.97;
            parallel.UseParallel = true;
            parallel.Solve();
            CollectionAssert.AreEqual(serial.Value, parallel.Value);
            CollectionAssert.AreEqual(serial.Policy, parallel.Policy);

            var det = OracleDeterministic(30, 25, 8);
            det.Solve();
            var detParallel = OracleDeterministic(30, 25, 8);
            detParallel.UseParallel = true;
            detParallel.Solve();
            CollectionAssert.AreEqual(det.Value, detParallel.Value);
            CollectionAssert.AreEqual(det.Policy, detParallel.Policy);
        }

        /// <summary>
        /// A parallel solve throws the same exception as the serial solve: the one from the lowest-numbered
        /// failing state, not an <see cref="AggregateException"/>.
        /// </summary>
        [TestMethod]
        public void ParallelSolveThrowsSerialException()
        {
            Func<int, int, int, double> cost = (t, s, a) => 1.0;
            Func<int, int, int, int> badTransition = (t, s, a) => s >= 2 ? 99 : 0;
            Func<int, int, int, double> badCost = (t, s, a) => s >= 3 ? throw new ArithmeticException($"state {s}") : 1.0;
            Func<int, int, int, int> stay = (t, s, a) => s;
            foreach (var (stageValue, transition) in new[] { (cost, badTransition), (badCost, stay) })
            {
                var serial = new DynamicProgramming(3, 16, 2, stageValue, transition);
                var parallel = new DynamicProgramming(3, 16, 2, stageValue, transition) { UseParallel = true };
                Exception expected = null, actual = null;
                try { serial.Solve(); } catch (Exception ex) { expected = ex; }
                try { parallel.Solve(); } catch (Exception ex) { actual = ex; }
                Assert.IsNotNull(expected);
                Assert.IsNotNull(actual);
                Assert.AreEqual(expected.GetType(), actual.GetType());
                Assert.AreEqual(expected.Message, actual.Message);
                Assert.IsFalse(parallel.IsSolved);
            }
        }

        /// <summary>
        /// Value iteration on a 2-state, 2-decision discounted MDP matches the exact solution of the linear
        /// Bellman system (I - γ P_π) v = c_π of the best of the four stationary policies.
        /// </summary>
        [TestMethod]
        public void ValueIterationMatchesClosedForm()
        {
            double gamma = 0.9;
            double[,] c = { { 2.0, 1.0 }, { 3.0, 0.5 } };
            double[,,] P = { { { 0.9, 0.1 }, { 0.2, 0.8 } }, { { 0.7, 0.3 }, { 0.05, 0.95 } } };
            var dp = new DynamicProgramming(1, 2, 2, (t, s, a) => c[s, a],
                (t, s, a, next, prob) => { next[0] = 0; next[1] = 1; prob[0] = P[s, a, 0]; prob[1] = P[s, a, 1]; return 2; }, 2)
            {
                DiscountFactor = gamma
            };
            Assert.IsTrue(dp.SolveSteadyState(1E-13, 100000));
            Assert.IsTrue(dp.SteadyStateConverged);
            Assert.IsGreaterThan(1, dp.SteadyStateIterations);

            double best0 = double.PositiveInfinity, best1 = double.PositiveInfinity;
            int pol0 = -1, pol1 = -1;
            for (int a0 = 0; a0 < 2; a0++)
                for (int a1 = 0; a1 < 2; a1++)
                {
                    double m00 = 1 - gamma * P[0, a0, 0], m01 = -gamma * P[0, a0, 1];
                    double m10 = -gamma * P[1, a1, 0], m11 = 1 - gamma * P[1, a1, 1];
                    double det = m00 * m11 - m01 * m10;
                    double v0 = (c[0, a0] * m11 - m01 * c[1, a1]) / det;
                    double v1 = (m00 * c[1, a1] - m10 * c[0, a0]) / det;
                    if (v0 < best0) { best0 = v0; pol0 = a0; }
                    if (v1 < best1) { best1 = v1; pol1 = a1; }
                }
            Assert.AreEqual(best0, dp.SteadyStateValue[0], 1E-9);
            Assert.AreEqual(best1, dp.SteadyStateValue[1], 1E-9);
            Assert.AreEqual(pol0, dp.SteadyStatePolicy[0]);
            Assert.AreEqual(pol1, dp.SteadyStatePolicy[1]);
        }

        /// <summary>
        /// Value iteration reports non-convergence when the sweep budget is too small, and converges for an
        /// undiscounted stochastic shortest-path problem with an absorbing zero-cost state.
        /// </summary>
        [TestMethod]
        public void ValueIterationConvergenceReporting()
        {
            var dp = new DynamicProgramming(1, 2, 1, (t, s, a) => 1.0, (t, s, a) => s);
            dp.DiscountFactor = 0.99;
            Assert.IsFalse(dp.SolveSteadyState(1E-10, 5));
            Assert.IsFalse(dp.SteadyStateConverged);
            Assert.AreEqual(5, dp.SteadyStateIterations);

            // Chain 2 -> 1 -> 0 with unit costs and an absorbing free state 0: V = (0, 1, 2).
            var ssp = new DynamicProgramming(1, 3, 1, (t, s, a) => s == 0 ? 0.0 : 1.0, (t, s, a) => Math.Max(s - 1, 0));
            Assert.IsTrue(ssp.SolveSteadyState());
            CollectionAssert.AreEqual(new[] { 0.0, 1.0, 2.0 }, ssp.SteadyStateValue);
            CollectionAssert.AreEqual(new[] { 0, 0, 0 }, ssp.SteadyStatePolicy);
        }

        /// <summary>
        /// A stationary state with no feasible decision has value +∞ and policy -1 in value iteration, and states
        /// that can only reach it become infeasible too.
        /// </summary>
        [TestMethod]
        public void ValueIterationInfeasibleStates()
        {
            // State 2 has no feasible decision; state 1 can only move to 2; state 0 stays at 0.
            var dp = new DynamicProgramming(1, 3, 1, (t, s, a) => s == 2 ? double.NaN : 1.0, (t, s, a) => s == 0 ? 0 : 2)
            {
                DiscountFactor = 0.5
            };
            Assert.IsTrue(dp.SolveSteadyState());
            Assert.AreEqual(2.0, dp.SteadyStateValue[0], 1E-7);
            Assert.IsTrue(double.IsPositiveInfinity(dp.SteadyStateValue[1]));
            Assert.IsTrue(double.IsPositiveInfinity(dp.SteadyStateValue[2]));
            CollectionAssert.AreEqual(new[] { 0, -1, -1 }, dp.SteadyStatePolicy);
        }

        /// <summary>
        /// Discounting: a single state paying 1 per stage for 3 stages with γ = 0.5 and terminal value 8 is worth
        /// 1 + 0.5 (1 + 0.5 (1 + 0.5 · 8)) = 2.75.
        /// </summary>
        [TestMethod]
        public void DiscountFactorKnownValue()
        {
            var dp = new DynamicProgramming(3, 1, 1, (t, s, a) => 1.0, (t, s, a) => 0)
            {
                DiscountFactor = 0.5,
                TerminalValue = new[] { 8.0 }
            };
            dp.Solve();
            Assert.AreEqual(2.75, dp.Value[0, 0]);
            Assert.AreEqual(3.5, dp.Value[1, 0]);
            Assert.AreEqual(5.0, dp.Value[2, 0]);
            Assert.AreEqual(8.0, dp.Value[3, 0]);
        }

        /// <summary>
        /// Minimize and maximize pick opposite decisions, ties go to the lowest index, infeasible decisions are
        /// skipped, and a state with no feasible decision has value ±∞ and policy -1.
        /// </summary>
        [TestMethod]
        public void MinimizeMaximizeTiesAndInfeasible()
        {
            // State 0: costs (3, NaN, 1, 1, 5); state 1: all infeasible (+∞ or transition -1).
            double[] costs = { 3, double.NaN, 1, 1, 5 };
            Func<int, int, int, double> cost = (t, s, a) => s == 0 ? costs[a] : (a % 2 == 0 ? double.PositiveInfinity : 1.0);
            Func<int, int, int, int> next = (t, s, a) => s == 0 ? 0 : -1;

            var min = new DynamicProgramming(1, 2, 5, cost, next);
            min.Solve();
            Assert.AreEqual(1.0, min.Value[0, 0]);
            Assert.AreEqual(2, min.Policy[0, 0]);
            Assert.IsTrue(double.IsPositiveInfinity(min.Value[0, 1]));
            Assert.AreEqual(-1, min.Policy[0, 1]);

            var max = new DynamicProgramming(1, 2, 5, cost, next) { Maximize = true };
            max.Solve();
            Assert.AreEqual(5.0, max.Value[0, 0]);
            Assert.AreEqual(4, max.Policy[0, 0]);
            Assert.IsTrue(double.IsNegativeInfinity(max.Value[0, 1]));
            Assert.AreEqual(-1, max.Policy[0, 1]);
            Assert.Throws<InvalidOperationException>(() => max.Simulate(1));
        }

        /// <summary>
        /// A decision whose outcome has positive probability of reaching an infeasible state is infeasible, while
        /// a zero-probability outcome to that state is ignored.
        /// </summary>
        [TestMethod]
        public void StochasticInfeasibleOutcome()
        {
            // State 1 is infeasible at the terminal stage. Decision 0 reaches it with probability 0.5, decision 1 with 0.
            var dp = new DynamicProgramming(1, 2, 2, (t, s, a) => a == 0 ? 0.0 : 10.0,
                (t, s, a, next, prob) => { next[0] = 0; next[1] = 1; prob[0] = a == 0 ? 0.5 : 1.0; prob[1] = a == 0 ? 0.5 : 0.0; return 2; }, 2)
            {
                TerminalValue = new[] { 0.0, double.NaN }
            };
            dp.Solve();
            Assert.AreEqual(1, dp.Policy[0, 0]);
            Assert.AreEqual(10.0, dp.Value[0, 0]);
            Assert.IsTrue(double.IsPositiveInfinity(dp.Value[1, 1]));
        }

        /// <summary>
        /// Forward simulation of a deterministic optimal policy reproduces the stage-0 value exactly.
        /// </summary>
        [TestMethod]
        public void SimulateDeterministicReproducesValue()
        {
            var dp = OracleDeterministic(24, 11, 6);
            dp.DiscountFactor = 0.97;
            dp.Solve();
            for (int s0 = 0; s0 < 11; s0++)
            {
                var path = dp.Simulate(s0);
                Assert.AreEqual(dp.Value[0, s0], path.TotalValue);
                Assert.HasCount(25, path.States);
                Assert.AreEqual(s0, path.States[0]);
                for (int t = 0; t < 24; t++) Assert.AreEqual(dp.Policy[t, path.States[t]], path.Decisions[t]);
            }
        }

        /// <summary>
        /// Stochastic simulation follows a given outcome (inflow-class) sequence, and the mean simulated return
        /// with a seeded generator agrees with the expected value within four standard errors.
        /// </summary>
        [TestMethod]
        public void SimulateStochastic()
        {
            int K = 3;
            var dp = OracleStochastic(12, 11, 6, K, 2);
            dp.Solve();
            int s0 = 5 * K + 1;

            var classes = new[] { 2, 2, 1, 0, 0, 0, 1, 2, 2, 1, 1, 0 };
            var path = dp.Simulate(s0, classes);
            for (int t = 0; t < 12; t++) Assert.AreEqual(classes[t], path.States[t + 1] % K);
            Assert.AreEqual(path.TotalValue, dp.Simulate(s0, classes).TotalValue);

            var rng = new MersenneTwister(4242);
            int n = 20000;
            double sum = 0, sumSq = 0;
            for (int i = 0; i < n; i++)
            {
                double v = dp.Simulate(s0, rng).TotalValue;
                sum += v; sumSq += v * v;
            }
            double mean = sum / n;
            double se = Math.Sqrt((sumSq / n - mean * mean) / n);
            Assert.AreEqual(dp.Value[0, s0], mean, 4 * se);

            Assert.Throws<InvalidOperationException>(() => dp.Simulate(s0));
            Assert.Throws<ArgumentOutOfRangeException>(() => dp.Simulate(s0, new[] { 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }));
            Assert.Throws<ArgumentException>(() => dp.Simulate(s0, new[] { 0, 1 }));
        }

        /// <summary>
        /// Invalid dimensions, delegates, discount factors, terminal values and transition outputs are rejected
        /// with clear exceptions, and simulating before solving is an error.
        /// </summary>
        [TestMethod]
        public void InvalidInputs()
        {
            Func<int, int, int, double> cost = (t, s, a) => 0.0;
            Func<int, int, int, int> next = (t, s, a) => 0;
            TransitionDistribution dist = (t, s, a, ns, p) => { ns[0] = 0; p[0] = 1; return 1; };

            Assert.Throws<ArgumentOutOfRangeException>(() => new DynamicProgramming(0, 1, 1, cost, next));
            Assert.Throws<ArgumentOutOfRangeException>(() => new DynamicProgramming(1, 0, 1, cost, next));
            Assert.Throws<ArgumentOutOfRangeException>(() => new DynamicProgramming(1, 1, 0, cost, next));
            Assert.Throws<ArgumentOutOfRangeException>(() => new DynamicProgramming(1, 1, 1, cost, dist, 0));
            Assert.Throws<ArgumentNullException>(() => new DynamicProgramming(1, 1, 1, null!, next));
            Assert.Throws<ArgumentNullException>(() => new DynamicProgramming(1, 1, 1, cost, (Func<int, int, int, int>)null!));
            Assert.Throws<ArgumentNullException>(() => new DynamicProgramming(1, 1, 1, cost, (TransitionDistribution)null!, 1));

            var dp = new DynamicProgramming(2, 2, 1, cost, next);
            Assert.Throws<ArgumentOutOfRangeException>(() => dp.DiscountFactor = 0);
            Assert.Throws<ArgumentOutOfRangeException>(() => dp.DiscountFactor = 1.5);
            Assert.Throws<ArgumentOutOfRangeException>(() => dp.DiscountFactor = double.NaN);
            Assert.Throws<InvalidOperationException>(() => dp.Simulate(0));
            dp.TerminalValue = new double[3];
            Assert.Throws<ArgumentException>(() => dp.Solve());
            dp.TerminalValue = null;
            dp.Solve();
            Assert.IsTrue(dp.IsSolved);
            Assert.Throws<ArgumentOutOfRangeException>(() => dp.Simulate(2));
            Assert.Throws<ArgumentNullException>(() => dp.Simulate(0, (Random)null!));
            Assert.Throws<ArgumentNullException>(() => dp.Simulate(0, (IList<int>)null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => dp.SolveSteadyState(0));
            Assert.Throws<ArgumentOutOfRangeException>(() => dp.SolveSteadyState(double.PositiveInfinity));
            Assert.Throws<ArgumentOutOfRangeException>(() => dp.SolveSteadyState(1E-8, 0));

            Assert.Throws<InvalidOperationException>(() => new DynamicProgramming(1, 2, 1, cost, (t, s, a) => 2).Solve());
            Assert.Throws<InvalidOperationException>(() => new DynamicProgramming(1, 2, 1, cost, (t, s, a) => -2).Solve());
            Assert.Throws<InvalidOperationException>(() => new DynamicProgramming(1, 2, 1, cost,
                (t, s, a, ns, p) => { ns[0] = 0; p[0] = 0.9; return 1; }, 1).Solve());
            Assert.Throws<InvalidOperationException>(() => new DynamicProgramming(1, 2, 1, cost,
                (t, s, a, ns, p) => { ns[0] = 0; ns[1] = 1; p[0] = 1.5; p[1] = -0.5; return 2; }, 2).Solve());
            Assert.Throws<InvalidOperationException>(() => new DynamicProgramming(1, 2, 1, cost,
                (t, s, a, ns, p) => { ns[0] = 5; p[0] = 1; return 1; }, 1).Solve());
            Assert.Throws<InvalidOperationException>(() => new DynamicProgramming(1, 2, 1, cost,
                (t, s, a, ns, p) => 2, 1).Solve());
        }
    }
}
