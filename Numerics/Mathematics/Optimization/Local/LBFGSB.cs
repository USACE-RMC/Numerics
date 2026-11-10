using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// Contains the limited-memory BFGS algorithm for bound-constrained optimization (L-BFGS-B).
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// L-BFGS-B minimizes a smooth function subject to simple bounds <c>l ≤ x ≤ u</c>, where each bound may be
    /// infinite. Each iteration builds a limited-memory BFGS model of the objective from the last
    /// <see cref="MemorySize"/> steps in the compact representation, finds the generalized Cauchy point of the model
    /// along the projected steepest-descent path, which fixes the active bounds, minimizes the model over the free
    /// variables, and takes a More-Thuente line search along the resulting direction. The run stops when the largest
    /// projected gradient component is at most <see cref="Optimizer.AbsoluteTolerance"/> or the relative reduction of
    /// the objective, <c>(f_k - f_(k+1)) / max(|f_k|, |f_(k+1)|, 1)</c>, is at most
    /// <see cref="Optimizer.RelativeTolerance"/>. These are the <c>pgtol</c> and <c>factr·ε</c> tests of the reference
    /// code; SciPy's <c>gtol</c> and <c>ftol</c>.
    /// </para>
    /// <para>
    /// The subspace step follows L-BFGS-B 3.0: the unconstrained subspace minimizer is projected onto the box and kept
    /// when it is a descent direction, otherwise it is truncated at the first bound. The memory is reset when a
    /// triangular system is singular, a factorization fails, or a line search fails with a non-empty memory; an update
    /// is skipped when the curvature <c>sᵀy</c> is not sufficiently positive.
    /// </para>
    /// <para>
    /// Beyond the reference code: a start outside the box is projected onto it, and a trial point whose objective or
    /// gradient is not finite restarts the line search with half the step as its maximum step.
    /// </para>
    /// <para>
    /// <b> Third-party notice: </b>
    /// This is a C# translation of the SciPy C translation (scipy/optimize/__lbfgsb.c, SciPy 1.17.1) of the L-BFGS-B
    /// 3.0 Fortran subroutines by Ciyou Zhu, Richard Byrd, Peihuang Lu, Jorge Nocedal and Jose Luis Morales, which are
    /// distributed under the New BSD License. The SciPy translation is Copyright (C) 2024 SciPy developers and is
    /// distributed under the 3-clause BSD license reproduced in the source file.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, "A limited memory algorithm for bound constrained optimization",
    /// SIAM Journal on Scientific Computing 16(5), 1190-1208, 1995.
    /// </description></item>
    /// <item><description>
    /// C. Zhu, R. H. Byrd, P. Lu and J. Nocedal, "Algorithm 778: L-BFGS-B: Fortran subroutines for large-scale
    /// bound-constrained optimization", ACM Transactions on Mathematical Software 23(4), 550-560, 1997.
    /// </description></item>
    /// <item><description>
    /// J. L. Morales and J. Nocedal, "Remark on 'Algorithm 778: L-BFGS-B: Fortran subroutines for large-scale bound
    /// constrained optimization'", ACM Transactions on Mathematical Software 38(1), article 7, 2011.
    /// </description></item>
    /// <item><description>
    /// J. J. More and D. J. Thuente, "Line search algorithms with guaranteed sufficient decrease", ACM Transactions on
    /// Mathematical Software 20(3), 286-307, 1994.
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class LBFGSB : Optimizer
    {

        // Copyright (C) 2024 SciPy developers
        //
        // Redistribution and use in source and binary forms, with or without
        // modification, are permitted provided that the following conditions are met:
        //
        // a. Redistributions of source code must retain the above copyright notice,
        //    this list of conditions and the following disclaimer.
        // b. Redistributions in binary form must reproduce the above copyright
        //    notice, this list of conditions and the following disclaimer in the
        //    documentation and/or other materials provided with the distribution.
        // c. Names of the SciPy Developers may not be used to endorse or promote
        //    products derived from this software without specific prior written
        //    permission.
        //
        // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
        // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
        // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
        // ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
        // BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
        // OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
        // SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
        // INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
        // CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
        // ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
        // THE POSSIBILITY OF SUCH DAMAGE.
        //
        // L-BFGS-B (version 3.0, April 25, 2011) by Ciyou Zhu, Richard Byrd, Peihuang Lu, Jorge Nocedal and
        // Jose Luis Morales is released under the "New BSD License" (the 3-clause license).

        /// <summary>
        /// Construct a new L-BFGS-B optimization method.
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values; a value outside its bounds is projected onto them.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive). Use negative infinity for no lower bound.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive). Use positive infinity for no upper bound.</param>
        /// <param name="gradient">Optional. Function to evaluate the gradient. Default uses finite differences.</param>
        /// <exception cref="ArgumentOutOfRangeException">The arrays have the wrong length, a bound is NaN, an upper bound
        /// is less than its lower bound, or an initial value is not finite.</exception>
        public LBFGSB(Func<double[], double> objectiveFunction, int numberOfParameters,
                      IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds,
                      Func<double[], double[]>? gradient = null) : base(objectiveFunction, numberOfParameters)
        {
            if (initialValues.Count != numberOfParameters || lowerBounds.Count != numberOfParameters || upperBounds.Count != numberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(lowerBounds), "The initial values and lower and upper bounds must be the same length as the number of parameters.");
            }
            for (int j = 0; j < numberOfParameters; j++)
            {
                if (double.IsNaN(lowerBounds[j]) || double.IsNaN(upperBounds[j]))
                    throw new ArgumentOutOfRangeException(nameof(lowerBounds), "The bounds cannot be NaN.");
                if (upperBounds[j] < lowerBounds[j])
                    throw new ArgumentOutOfRangeException(nameof(upperBounds), "The upper bound cannot be less than the lower bound.");
                if (!Tools.IsFinite(initialValues[j]))
                    throw new ArgumentOutOfRangeException(nameof(initialValues), "The initial values must be finite.");
            }
            InitialValues = initialValues.ToArray();
            LowerBounds = lowerBounds.ToArray();
            UpperBounds = upperBounds.ToArray();
            Gradient = gradient;
        }

        /// <summary>
        /// An array of initial values to evaluate.
        /// </summary>
        public double[] InitialValues { get; private set; }

        /// <summary>
        /// An array of lower bounds (inclusive). Negative infinity means no lower bound.
        /// </summary>
        public double[] LowerBounds { get; private set; }

        /// <summary>
        /// An array of upper bounds (inclusive). Positive infinity means no upper bound.
        /// </summary>
        public double[] UpperBounds { get; private set; }

        /// <inheritdoc />
        protected override double[]? ParameterLowerBounds => LowerBounds;

        /// <inheritdoc />
        protected override double[]? ParameterUpperBounds => UpperBounds;

        /// <summary>
        /// The function for evaluating the gradient of the objective function.
        /// </summary>
        public Func<double[], double[]>? Gradient;

        /// <summary>
        /// Gets or sets the number of correction pairs kept in the limited-memory matrix. Default = 10.
        /// </summary>
        public int MemorySize { get; set; } = 10;

        /// <summary>
        /// Gets or sets the maximum number of function evaluations in one line search. Default = 20.
        /// </summary>
        public int MaxLineSearchSteps { get; set; } = 20;

        /// <inheritdoc/>
        protected override void Validate()
        {
            base.Validate();
            if (MemorySize < 1) throw new ArgumentOutOfRangeException(nameof(MemorySize), "The memory size must be at least one.");
            if (MaxLineSearchSteps < 1) throw new ArgumentOutOfRangeException(nameof(MaxLineSearchSteps), "The maximum number of line search steps must be at least one.");
        }

        /// <summary>The machine epsilon used by the reference code.</summary>
        private const double Epsilon = 2.220446049250313e-016;

        /// <summary>The sufficient decrease, curvature and interval tolerances of the line search.</summary>
        private const double Ftol = 1e-3, Gtol = 0.9, Xtol = 0.1;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            int n = NumberOfParameters, m = MemorySize;
            double pgtol = AbsoluteTolerance, tol = RelativeTolerance;
            var l = LowerBounds;
            var u = UpperBounds;
            var nbd = new int[n];
            for (int i = 0; i < n; i++)
            {
                bool lower = !double.IsNegativeInfinity(l[i]), upper = !double.IsPositiveInfinity(u[i]);
                nbd[i] = lower && upper ? 2 : lower ? 1 : upper ? 3 : 0;
            }

            var ws = new double[m * n];
            var wy = new double[m * n];
            var sy = new double[m * m];
            var ss = new double[m * m];
            var wt = new double[m * m];
            var wn = new double[4 * m * m];
            var snd = new double[4 * m * m];
            var z = new double[n];
            var r = new double[n];
            var d = new double[n];
            var t = new double[n];
            var xp = new double[n];
            var wa = new double[8 * m];
            var index = new int[n];
            var iwhere = new int[n];
            var indx2 = new int[n];

            var x = (double[])InitialValues.Clone();
            Active(n, l, u, nbd, x, iwhere, out bool cnstnd, out bool boxed);

            int col = 0, head = 0, iupdat = 0, itail = 0, iter = 0, nfree = n, nenter = 0, ileave = 0;
            bool updatd = false;
            double theta = 1.0;

            bool cancel = false;
            double f = EvaluateObjective(x, ref cancel);
            if (cancel) return;
            if (!Tools.IsFinite(f))
                throw new ArgumentException("The initial objective value must be finite.", nameof(ObjectiveFunction));
            var g = EvaluateGradient(x, ref cancel);
            if (cancel) return;
            if (g == null) throw new ArgumentException("The gradient must contain only finite values.", nameof(Gradient));

            double sbgnrm = ProjectedGradientNorm(n, l, u, nbd, x, g);
            if (sbgnrm <= pgtol)
            {
                Finish(x, f, OptimizationStatus.Success);
                return;
            }

            var search = new MoreThuenteSearch();
            while (true)
            {
                bool wrk;
                if (!cnstnd && col > 0)
                {
                    Array.Copy(x, z, n);
                    wrk = updatd;
                }
                else
                {
                    // Compute the generalized Cauchy point.
                    if (Cauchy(n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m, wy, ws, sy, wt, theta, col, head, wa, sbgnrm) != 0)
                    {
                        // A singular triangular system; refresh the memory.
                        ResetMemory(ref col, ref head, ref theta, ref iupdat, ref updatd);
                        continue;
                    }
                    FreeVariables(n, ref nfree, index, out nenter, out ileave, indx2, iwhere, out wrk, updatd, cnstnd, iter);
                }

                if (nfree != 0 && col != 0)
                {
                    // Subspace minimization over the free variables.
                    int info = 0;
                    if (wrk) info = FormK(n, nfree, index, nenter, ileave, indx2, iupdat, updatd, wn, snd, m, ws, wy, sy, theta, col, head);
                    if (info == 0) info = ReducedGradient(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, col, head, nfree, cnstnd);
                    if (info == 0) info = SubspaceMinimization(n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy, theta, x, g, col, head, wa, wn);
                    if (info != 0)
                    {
                        ResetMemory(ref col, ref head, ref theta, ref iupdat, ref updatd);
                        continue;
                    }
                }

                for (int i = 0; i < n; i++) d[i] = z[i] - x[i];

                bool accepted = LineSearch(n, l, u, nbd, x, ref f, out double fold, out double gd, out double gdold, g, d, r, t, z,
                    out double stp, out double dtd, iter, boxed, cnstnd, search, ref cancel);
                if (cancel) return;
                if (!accepted)
                {
                    Array.Copy(t, x, n);
                    Array.Copy(r, g, n);
                    f = fold;
                    if (col == 0)
                    {
                        iter++;
                        Iterations = iter;
                        Finish(x, f, OptimizationStatus.LineSearchFailed);
                        return;
                    }
                    ResetMemory(ref col, ref head, ref theta, ref iupdat, ref updatd);
                    continue;
                }

                iter++;
                Iterations = iter;
                sbgnrm = ProjectedGradientNorm(n, l, u, nbd, x, g);
                if (iter >= MaxIterations)
                {
                    BestParameterSet = new ParameterSet((double[])x.Clone(), f);
                    UpdateStatus(OptimizationStatus.MaximumIterationsReached);
                    return;
                }

                // Termination tests.
                if (sbgnrm <= pgtol || fold - f <= tol * Math.Max(Math.Max(Math.Abs(fold), Math.Abs(f)), 1.0))
                {
                    Finish(x, f, OptimizationStatus.Success);
                    return;
                }

                // Form the correction pair s = stp·d and y = g_new - g_old.
                for (int i = 0; i < n; i++) r[i] = g[i] - r[i];
                double rr = CrossProductKernel.Dot(r, 0, r, 0, n);
                double dr, ddum;
                if (stp == 1.0)
                {
                    dr = gd - gdold;
                    ddum = -gdold;
                }
                else
                {
                    dr = (gd - gdold) * stp;
                    for (int i = 0; i < n; i++) d[i] *= stp;
                    ddum = -gdold * stp;
                }
                if (dr <= Epsilon * ddum)
                {
                    // Skip the update when the curvature is not sufficiently positive.
                    updatd = false;
                    continue;
                }

                updatd = true;
                iupdat++;
                UpdateCorrections(n, m, ws, wy, sy, ss, d, r, ref itail, iupdat, ref col, ref head, ref theta, rr, dr, stp, dtd);
                if (FormT(m, wt, sy, ss, col, theta) != 0)
                    ResetMemory(ref col, ref head, ref theta, ref iupdat, ref updatd);
            }
        }

        /// <summary>Records the accepted point as the result and sets the status.</summary>
        /// <param name="x">The accepted point.</param>
        /// <param name="f">Its scaled objective.</param>
        /// <param name="status">The final status.</param>
        private void Finish(double[] x, double f, OptimizationStatus status)
        {
            BestParameterSet = new ParameterSet((double[])x.Clone(), f);
            UpdateStatus(status);
        }

        /// <summary>Discards the limited-memory corrections.</summary>
        private static void ResetMemory(ref int col, ref int head, ref double theta, ref int iupdat, ref bool updatd)
        {
            col = 0;
            head = 0;
            theta = 1.0;
            iupdat = 0;
            updatd = false;
        }

        /// <summary>Evaluates an objective trial while retaining only finite incumbents.</summary>
        /// <param name="x">The trial point.</param>
        /// <param name="cancel">The evaluation-budget cancellation flag.</param>
        /// <returns>The scaled objective value, positive infinity for a NaN.</returns>
        private double EvaluateObjective(double[] x, ref bool cancel)
        {
            var incumbent = BestParameterSet;
            double f = Evaluate(x, ref cancel);
            if (!Tools.IsFinite(f)) BestParameterSet = incumbent;
            return f;
        }

        /// <summary>Evaluates the gradient in minimization coordinates.</summary>
        /// <param name="x">The point at which to differentiate.</param>
        /// <param name="cancel">The evaluation-budget cancellation flag.</param>
        /// <returns>A private copy of the scaled gradient, or null when a component is not finite.</returns>
        /// <exception cref="ArgumentException">The supplied gradient has the wrong length.</exception>
        private double[]? EvaluateGradient(double[] x, ref bool cancel)
        {
            double[] g;
            if (Gradient != null)
            {
                var supplied = Gradient((double[])x.Clone());
                if (supplied == null || supplied.Length != NumberOfParameters)
                    throw new ArgumentException("The gradient must contain one value per parameter.", nameof(Gradient));
                g = (double[])supplied.Clone();
                for (int i = 0; i < g.Length; i++) g[i] *= functionScale;
            }
            else
            {
                bool stopped = cancel;
                // Finite differences may request more probes after cancellation; do not spend beyond the budget.
                g = NumericalDerivative.Gradient(p => stopped ? double.NaN : EvaluateObjective(p, ref stopped),
                    x, LowerBounds, UpperBounds);
                cancel = stopped;
                if (cancel) return g;
            }
            for (int i = 0; i < g.Length; i++)
                if (!Tools.IsFinite(g[i])) return null;
            return g;
        }

        /// <summary>
        /// Projects the start onto the box and classifies the variables (reference routine <c>active</c>).
        /// </summary>
        private static void Active(int n, double[] l, double[] u, int[] nbd, double[] x, int[] iwhere, out bool cnstnd, out bool boxed)
        {
            cnstnd = false;
            boxed = true;
            for (int i = 0; i < n; i++)
            {
                if (nbd[i] > 0)
                {
                    if (nbd[i] <= 2 && x[i] <= l[i])
                    {
                        if (x[i] < l[i]) x[i] = l[i];
                    }
                    else if (nbd[i] >= 2 && x[i] >= u[i])
                    {
                        if (x[i] > u[i]) x[i] = u[i];
                    }
                }
            }
            for (int i = 0; i < n; i++)
            {
                if (nbd[i] != 2) boxed = false;
                if (nbd[i] == 0)
                {
                    iwhere[i] = -1;
                }
                else
                {
                    cnstnd = true;
                    iwhere[i] = nbd[i] == 2 && u[i] - l[i] <= 0.0 ? 3 : 0;
                }
            }
        }

        /// <summary>
        /// Returns the infinity norm of the projected gradient (reference routine <c>projgr</c>).
        /// </summary>
        private static double ProjectedGradientNorm(int n, double[] l, double[] u, int[] nbd, double[] x, double[] g)
        {
            double norm = 0.0;
            for (int i = 0; i < n; i++)
            {
                double gi = g[i];
                if (nbd[i] != 0)
                {
                    if (gi < 0.0)
                    {
                        if (nbd[i] >= 2) gi = Math.Max(x[i] - u[i], gi);
                    }
                    else
                    {
                        if (nbd[i] <= 2) gi = Math.Min(x[i] - l[i], gi);
                    }
                }
                norm = Math.Max(norm, Math.Abs(gi));
            }
            return norm;
        }

        /// <summary>
        /// Computes the product of the middle matrix of the compact L-BFGS formula with a 2·col vector
        /// (reference routine <c>bmv</c>).
        /// </summary>
        /// <returns>Zero, or a nonzero value when a triangular system is singular.</returns>
        private static int Bmv(int m, double[] sy, double[] wt, int col, double[] v, int vo, double[] p, int po)
        {
            if (col == 0) return 0;
            // Solve [ D^(1/2)  O ; -L D^(-1/2)  J ] [p1; p2] = [v1; v2].
            p[po + col] = v[vo + col];
            for (int i = 1; i < col; i++)
            {
                double sum = 0.0;
                for (int k = 0; k <= i - 1; k++) sum += sy[i + m * k] * v[vo + k] / sy[k + m * k];
                p[po + col + i] = v[vo + col + i] + sum;
            }
            if (TriangularSolve(wt, 0, m, col, p, po + col, true) != 0) return 1;
            for (int i = 0; i < col; i++) p[po + i] = v[vo + i] / Math.Sqrt(sy[i + m * i]);

            // Solve [ -D^(1/2)  D^(-1/2) L' ; 0  J' ] [p1; p2] = [p1; p2].
            if (TriangularSolve(wt, 0, m, col, p, po + col, false) != 0) return 1;
            for (int i = 0; i < col; i++) p[po + i] = -p[po + i] / Math.Sqrt(sy[i + m * i]);
            for (int i = 0; i < col; i++)
            {
                double sum = 0.0;
                for (int k = i + 1; k < col; k++) sum += sy[k + m * i] * p[po + col + k] / sy[i + m * i];
                p[po + i] += sum;
            }
            return 0;
        }

        /// <summary>
        /// Finds the generalized Cauchy point along the projected steepest-descent path (reference routine <c>cauchy</c>).
        /// </summary>
        /// <returns>Zero, or a nonzero value when a triangular system is singular.</returns>
        /// <remarks>Byrd, Lu, Nocedal and Zhu (1995), section 4; the breakpoints are ordered with a heap.</remarks>
        private static int Cauchy(int n, double[] x, double[] l, double[] u, int[] nbd, double[] g, int[] iorder,
            int[] iwhere, double[] t, double[] d, double[] xcp, int m, double[] wy, double[] ws, double[] sy,
            double[] wt, double theta, int col, int head, double[] wa, double sbgnrm)
        {
            // Work vectors in wa: p = [0, 2m), c = [2m, 4m), wbp = [4m, 6m), v = [6m, 8m).
            int pOff = 0, cOff = 2 * m, wbpOff = 4 * m, vOff = 6 * m;
            if (sbgnrm <= 0.0)
            {
                Array.Copy(x, xcp, n);
                return 0;
            }
            bool bnded = true;
            int nfree = n, nbreak = -1, ibkmin = 0, col2 = 2 * col;
            double bkmin = 0.0, f1 = 0.0, tl = 0.0, tu = 0.0;
            for (int i = 0; i < col2; i++) wa[pOff + i] = 0.0;

            for (int i = 0; i < n; i++)
            {
                double neggi = -g[i];
                if (iwhere[i] != 3 && iwhere[i] != -1)
                {
                    // Distances to the bounds; classify the variable at the start of the path.
                    if (nbd[i] <= 2) tl = x[i] - l[i];
                    if (nbd[i] >= 2) tu = u[i] - x[i];
                    bool xlower = nbd[i] <= 2 && tl <= 0.0;
                    bool xupper = nbd[i] >= 2 && tu <= 0.0;
                    iwhere[i] = 0;
                    if (xlower)
                    {
                        if (neggi <= 0.0) iwhere[i] = 1;
                    }
                    else if (xupper)
                    {
                        if (neggi >= 0.0) iwhere[i] = 2;
                    }
                    else
                    {
                        if (Math.Abs(neggi) <= 0.0) iwhere[i] = -3;
                    }
                }
                int pointr = head;
                if (iwhere[i] != 0 && iwhere[i] != -1)
                {
                    d[i] = 0.0;
                }
                else
                {
                    d[i] = neggi;
                    f1 -= neggi * neggi;
                    // p = W'd.
                    for (int j = 0; j < col; j++)
                    {
                        wa[pOff + j] += wy[i + n * pointr] * neggi;
                        wa[pOff + col + j] += ws[i + n * pointr] * neggi;
                        pointr = (pointr + 1) % m;
                    }
                    if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.0)
                    {
                        // x(i) + d(i) is bounded; record the breakpoint.
                        nbreak++;
                        iorder[nbreak] = i;
                        t[nbreak] = tl / (-neggi);
                        if (nbreak == 0 || t[nbreak] < bkmin)
                        {
                            bkmin = t[nbreak];
                            ibkmin = nbreak;
                        }
                    }
                    else if (nbd[i] >= 2 && neggi > 0.0)
                    {
                        nbreak++;
                        iorder[nbreak] = i;
                        t[nbreak] = tu / neggi;
                        if (nbreak == 0 || t[nbreak] < bkmin)
                        {
                            bkmin = t[nbreak];
                            ibkmin = nbreak;
                        }
                    }
                    else
                    {
                        // x(i) + d(i) is not bounded.
                        nfree--;
                        iorder[nfree] = i;
                        if (Math.Abs(neggi) > 0.0) bnded = false;
                    }
                }
            }

            // The indices of the nonzero components of d are in iorder[0..nbreak] and iorder[nfree..n).
            if (theta != 1.0)
                for (int j = 0; j < col; j++) wa[pOff + col + j] *= theta;
            Array.Copy(x, xcp, n);
            if (nbreak == -1 && nfree == n) return 0;

            for (int j = 0; j < col2; j++) wa[cOff + j] = 0.0;
            double f2 = -theta * f1;
            double f2Org = f2;
            if (col > 0)
            {
                if (Bmv(m, sy, wt, col, wa, pOff, wa, vOff) != 0) return 1;
                f2 -= CrossProductKernel.Dot(wa, vOff, wa, pOff, col2);
            }
            double dtm = -f1 / f2;
            double tsum = 0.0;

            if (nbreak != -1)
            {
                int nleft = nbreak, iter = 0;
                double tj = 0.0;
                while (true)
                {
                    // Find the next smallest breakpoint.
                    double tj0 = tj;
                    int ibp;
                    if (iter == 0)
                    {
                        tj = bkmin;
                        ibp = iorder[ibkmin];
                    }
                    else
                    {
                        if (iter == 1 && ibkmin != nbreak)
                        {
                            // Replace the already used smallest breakpoint with the last one, then heap sort.
                            t[ibkmin] = t[nbreak];
                            iorder[ibkmin] = iorder[nbreak];
                        }
                        HeapSort(nleft, t, iorder, iter - 1);
                        tj = t[nleft];
                        ibp = iorder[nleft];
                    }
                    double dt = tj - tj0;

                    // The minimizer lies within this interval.
                    if (dtm < dt) break;

                    // Otherwise fix one variable and reset the coefficients of the quadratic.
                    tsum += dt;
                    nleft--;
                    iter++;
                    double dibp = d[ibp];
                    d[ibp] = 0.0;
                    double zibp;
                    if (dibp > 0.0)
                    {
                        zibp = u[ibp] - x[ibp];
                        xcp[ibp] = u[ibp];
                        iwhere[ibp] = 2;
                    }
                    else
                    {
                        zibp = l[ibp] - x[ibp];
                        xcp[ibp] = l[ibp];
                        iwhere[ibp] = 1;
                    }
                    if (nleft == -1 && nbreak == n)
                    {
                        // All n variables are fixed.
                        dtm = dt;
                        if (col > 0) CrossProductKernel.Axpy(dtm, wa, pOff, wa, cOff, col2);
                        return 0;
                    }

                    double dibp2 = dibp * dibp;
                    f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
                    f2 -= theta * dibp2;
                    if (col > 0)
                    {
                        CrossProductKernel.Axpy(dt, wa, pOff, wa, cOff, col2);
                        // wbp is the row of W corresponding to the breakpoint.
                        int pointr = head;
                        for (int j = 0; j < col; j++)
                        {
                            wa[wbpOff + j] = wy[ibp + n * pointr];
                            wa[wbpOff + col + j] = theta * ws[ibp + n * pointr];
                            pointr = (pointr + 1) % m;
                        }
                        if (Bmv(m, sy, wt, col, wa, wbpOff, wa, vOff) != 0) return 1;
                        double wmc = CrossProductKernel.Dot(wa, cOff, wa, vOff, col2);
                        double wmp = CrossProductKernel.Dot(wa, pOff, wa, vOff, col2);
                        double wmw = CrossProductKernel.Dot(wa, wbpOff, wa, vOff, col2);
                        CrossProductKernel.Axpy(-dibp, wa, wbpOff, wa, pOff, col2);
                        f1 += dibp * wmc;
                        f2 = f2 + dibp * 2.0 * wmp - dibp2 * wmw;
                    }
                    f2 = Math.Max(Epsilon * f2Org, f2);
                    if (nleft >= 0)
                    {
                        dtm = -f1 / f2;
                    }
                    else if (bnded)
                    {
                        dtm = 0.0;
                        break;
                    }
                    else
                    {
                        dtm = -f1 / f2;
                        break;
                    }
                }
            }

            if (dtm <= 0.0) dtm = 0.0;
            tsum += dtm;
            // Move free variables and those whose breakpoints were not reached.
            CrossProductKernel.Axpy(tsum, d, 0, xcp, 0, n);
            if (col > 0) CrossProductKernel.Axpy(dtm, wa, pOff, wa, cOff, col2);
            return 0;
        }

        /// <summary>
        /// Computes the reduced gradient <c>r = -Z'B(xcp - x) - Z'g</c> (reference routine <c>cmprlb</c>).
        /// </summary>
        /// <returns>Zero, or a nonzero value when a triangular system is singular.</returns>
        private static int ReducedGradient(int n, int m, double[] x, double[] g, double[] ws, double[] wy, double[] sy,
            double[] wt, double[] z, double[] r, double[] wa, int[] index, double theta, int col, int head, int nfree, bool cnstnd)
        {
            if (!cnstnd && col > 0)
            {
                for (int i = 0; i < n; i++) r[i] = -g[i];
                return 0;
            }
            for (int i = 0; i < nfree; i++)
            {
                int k = index[i];
                r[i] = -theta * (z[k] - x[k]) - g[k];
            }
            if (Bmv(m, sy, wt, col, wa, 2 * m, wa, 0) != 0) return -8;
            int pointr = head;
            for (int j = 0; j < col; j++)
            {
                double a1 = wa[j], a2 = theta * wa[col + j];
                for (int i = 0; i < nfree; i++)
                {
                    int k = index[i];
                    r[i] += wy[k + n * pointr] * a1 + ws[k + n * pointr] * a2;
                }
                pointr = (pointr + 1) % m;
            }
            return 0;
        }

        /// <summary>
        /// Forms and factors the matrix of the subspace problem (reference routine <c>formk</c>).
        /// </summary>
        /// <returns>Zero, or a nonzero value when a Cholesky factorization fails.</returns>
        /// <remarks>
        /// Forms <c>K = [-D - Y'ZZ'Y/θ, L_a' - R_z'; L_a - R_z, θS'AA'S]</c> in the upper triangle of wn as the
        /// factorization <c>K = LEL'</c> with <c>E = [-I 0; 0 I]</c>. The inner products of the free and active
        /// columns are updated from the previous iteration in wn1 using the entering and leaving variables.
        /// </remarks>
        private static int FormK(int n, int nsub, int[] ind, int nenter, int ileave, int[] indx2, int iupdat, bool updatd,
            double[] wn, double[] wn1, int m, double[] ws, double[] wy, double[] sy, double theta, int col, int head)
        {
            int m2 = 2 * m, upcl;
            if (updatd)
            {
                if (iupdat > m)
                {
                    // Shift the old part of wn1.
                    for (int jy = 0; jy < m - 1; jy++)
                    {
                        int js = m + jy;
                        int len = m - (jy + 1);
                        Array.Copy(wn1, (jy + 1) + m2 * (jy + 1), wn1, jy + m2 * jy, len);
                        Array.Copy(wn1, (js + 1) + m2 * (js + 1), wn1, js + m2 * js, len);
                        Array.Copy(wn1, (m + 1) + m2 * (jy + 1), wn1, m + m2 * jy, m - 1);
                    }
                }

                // Put the new rows in blocks (1,1), (2,1) and (2,2).
                int iy = col - 1, isr = m + col - 1;
                int ipntr = head + col - 1;
                if (ipntr >= m) ipntr -= m;
                int jpntr = head;
                for (int jy = 0; jy < col; jy++)
                {
                    int js = m + jy;
                    double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0;
                    for (int k = 0; k < nsub; k++)
                    {
                        int k1 = ind[k];
                        temp1 += wy[k1 + n * ipntr] * wy[k1 + n * jpntr];
                    }
                    for (int k = nsub; k < n; k++)
                    {
                        int k1 = ind[k];
                        temp2 += ws[k1 + n * ipntr] * ws[k1 + n * jpntr];
                        temp3 += ws[k1 + n * ipntr] * wy[k1 + n * jpntr];
                    }
                    wn1[iy + m2 * jy] = temp1;
                    wn1[isr + m2 * js] = temp2;
                    wn1[isr + m2 * jy] = temp3;
                    jpntr = (jpntr + 1) % m;
                }

                // Put the new column in block (2,1).
                int jyc = col - 1;
                jpntr = head + col - 1;
                if (jpntr >= m) jpntr -= m;
                ipntr = head;
                for (int i = 0; i < col; i++)
                {
                    int ist = m + i;
                    double temp3 = 0.0;
                    for (int k = 0; k < nsub; k++)
                    {
                        int k1 = ind[k];
                        temp3 += ws[k1 + n * ipntr] * wy[k1 + n * jpntr];
                    }
                    ipntr = (ipntr + 1) % m;
                    wn1[ist + m2 * jyc] = temp3;
                }
                upcl = col - 1;
            }
            else
            {
                upcl = col;
            }

            // Modify the old parts in blocks (1,1) and (2,2) for the changes in the free set.
            {
                int ipntr = head;
                for (int iy = 0; iy < upcl; iy++)
                {
                    int isr = m + iy;
                    int jpntr = head;
                    for (int jy = 0; jy <= iy; jy++)
                    {
                        int js = m + jy;
                        double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0;
                        for (int k = 0; k < nenter; k++)
                        {
                            int k1 = indx2[k];
                            temp1 += wy[k1 + n * ipntr] * wy[k1 + n * jpntr];
                            temp2 += ws[k1 + n * ipntr] * ws[k1 + n * jpntr];
                        }
                        for (int k = ileave; k < n; k++)
                        {
                            int k1 = indx2[k];
                            temp3 += wy[k1 + n * ipntr] * wy[k1 + n * jpntr];
                            temp4 += ws[k1 + n * ipntr] * ws[k1 + n * jpntr];
                        }
                        wn1[iy + m2 * jy] = wn1[iy + m2 * jy] + temp1 - temp3;
                        wn1[isr + m2 * js] = wn1[isr + m2 * js] - temp2 + temp4;
                        jpntr = (jpntr + 1) % m;
                    }
                    ipntr = (ipntr + 1) % m;
                }

                // Modify the old parts in block (2,1).
                ipntr = head;
                for (int isr = m; isr < m + upcl; isr++)
                {
                    int jpntr = head;
                    for (int jy = 0; jy < upcl; jy++)
                    {
                        double temp1 = 0.0, temp3 = 0.0;
                        for (int k = 0; k < nenter; k++)
                        {
                            int k1 = indx2[k];
                            temp1 += ws[k1 + n * ipntr] * wy[k1 + n * jpntr];
                        }
                        for (int k = ileave; k < n; k++)
                        {
                            int k1 = indx2[k];
                            temp3 += ws[k1 + n * ipntr] * wy[k1 + n * jpntr];
                        }
                        if (isr <= jy + m)
                            wn1[isr + m2 * jy] = wn1[isr + m2 * jy] + temp1 - temp3;
                        else
                            wn1[isr + m2 * jy] = wn1[isr + m2 * jy] - temp1 + temp3;
                        jpntr = (jpntr + 1) % m;
                    }
                    ipntr = (ipntr + 1) % m;
                }
            }

            // Form the upper triangle of wn = [D + Y'ZZ'Y/θ, -L_a' + R_z'; -L_a + R_z, θS'AA'S].
            for (int iy = 0; iy < col; iy++)
            {
                int isr = col + iy, is1 = m + iy;
                for (int jy = 0; jy <= iy; jy++)
                {
                    int js = col + jy, js1 = m + jy;
                    wn[jy + m2 * iy] = wn1[iy + m2 * jy] / theta;
                    wn[js + m2 * isr] = wn1[is1 + m2 * js1] * theta;
                }
                for (int jy = 0; jy < iy; jy++) wn[jy + m2 * isr] = -wn1[is1 + m2 * jy];
                for (int jy = iy; jy < col; jy++) wn[jy + m2 * isr] = wn1[is1 + m2 * jy];
                wn[iy + m2 * iy] += sy[iy + m * iy];
            }

            // Form the upper triangle of wn = [LL', L^-1(-L_a' + R_z'); L^-1(-L_a + R_z), S'AA'Sθ + ...].
            if (Cholesky(wn, 0, m2, col) != 0) return -1;
            int col2 = 2 * col;
            for (int js = col; js < col2; js++) TriangularSolve(wn, 0, m2, col, wn, m2 * js, true);
            for (int isr = col; isr < col2; isr++)
                for (int js = isr; js < col2; js++)
                    wn[isr + m2 * js] += CrossProductKernel.Dot(wn, m2 * isr, wn, m2 * js, col);
            if (Cholesky(wn, col + m2 * col, m2, col) != 0) return -2;
            return 0;
        }

        /// <summary>
        /// Forms and factors <c>T = θS'S + LD^(-1)L'</c> (reference routine <c>formt</c>).
        /// </summary>
        /// <returns>Zero, or a nonzero value when the Cholesky factorization fails.</returns>
        private static int FormT(int m, double[] wt, double[] sy, double[] ss, int col, double theta)
        {
            for (int j = 0; j < col; j++) wt[m * j] = theta * ss[m * j];
            for (int i = 1; i < col; i++)
            {
                for (int j = i; j < col; j++)
                {
                    int k1 = Math.Min(i, j);
                    double sum = 0.0;
                    for (int k = 0; k < k1; k++) sum += sy[i + m * k] * sy[j + m * k] / sy[k + m * k];
                    wt[i + m * j] = sum + theta * ss[i + m * j];
                }
            }
            return Cholesky(wt, 0, m, col) != 0 ? -3 : 0;
        }

        /// <summary>
        /// Builds the sets of free and active variables at the Cauchy point and counts the variables that enter and
        /// leave the free set (reference routine <c>freev</c>).
        /// </summary>
        private static void FreeVariables(int n, ref int nfree, int[] idx, out int nenter, out int ileave, int[] idx2,
            int[] iwhere, out bool wrk, bool updatd, bool cnstnd, int iter)
        {
            nenter = 0;
            ileave = n;
            if (iter > 0 && cnstnd)
            {
                for (int i = 0; i < nfree; i++)
                {
                    int k = idx[i];
                    if (iwhere[k] > 0)
                    {
                        ileave--;
                        idx2[ileave] = k;
                    }
                }
                for (int i = nfree; i < n; i++)
                {
                    int k = idx[i];
                    if (iwhere[k] <= 0)
                    {
                        idx2[nenter] = k;
                        nenter++;
                    }
                }
            }
            wrk = ileave < n || nenter > 0 || updatd;

            // Free variables first, then active ones in reverse order.
            nfree = 0;
            int iact = n;
            for (int i = 0; i < n; i++)
            {
                if (iwhere[i] <= 0)
                {
                    idx[nfree] = i;
                    nfree++;
                }
                else
                {
                    iact--;
                    idx[iact] = i;
                }
            }
        }

        /// <summary>
        /// Moves the least breakpoint to position n of a binary heap (reference routine <c>hpsolb</c>).
        /// </summary>
        /// <param name="n">The index of the last element of the heap.</param>
        /// <param name="t">The breakpoints.</param>
        /// <param name="iorder">The variable of each breakpoint.</param>
        /// <param name="iheap">Zero when t[0..n] is not yet a heap.</param>
        private static void HeapSort(int n, double[] t, int[] iorder, int iheap)
        {
            if (iheap == 0)
            {
                // Rearrange t[0..n] into a heap.
                for (int k = 2; k <= n + 1; k++)
                {
                    double value = t[k - 1];
                    int indxin = iorder[k - 1];
                    int i = k;
                    while (i > 1)
                    {
                        int j = i / 2;
                        if (value < t[j - 1])
                        {
                            t[i - 1] = t[j - 1];
                            iorder[i - 1] = iorder[j - 1];
                            i = j;
                        }
                        else break;
                    }
                    t[i - 1] = value;
                    iorder[i - 1] = indxin;
                }
            }
            if (n > 0)
            {
                // Move the least element to t[n] and restore the heap on t[0..n-1].
                int i = 1;
                double outValue = t[0];
                int indxout = iorder[0];
                double value = t[n];
                int indxin = iorder[n];
                while (true)
                {
                    int j = i + i;
                    if (j > n) break;
                    if (t[j] < t[j - 1]) j++;
                    if (t[j - 1] < value)
                    {
                        t[i - 1] = t[j - 1];
                        iorder[i - 1] = iorder[j - 1];
                        i = j;
                    }
                    else break;
                }
                t[i - 1] = value;
                iorder[i - 1] = indxin;
                t[n] = outValue;
                iorder[n] = indxout;
            }
        }

        /// <summary>
        /// Searches along d with the More-Thuente line search, keeping the trial points in the box
        /// (reference routine <c>lnsrlb</c>).
        /// </summary>
        /// <returns>Whether a step was accepted; on failure t and r hold the starting point and gradient.</returns>
        private bool LineSearch(int n, double[] l, double[] u, int[] nbd, double[] x, ref double f, out double fold,
            out double gd, out double gdold, double[] g, double[] d, double[] r, double[] t, double[] z, out double stp,
            out double dtd, int iter, bool boxed, bool cnstnd, MoreThuenteSearch search, ref bool cancel)
        {
            double dnorm = Math.Sqrt(CrossProductKernel.Dot(d, 0, d, 0, n));
            dtd = dnorm * dnorm;

            // Determine the maximum step length.
            double stpmx = 1.0e+10;
            if (cnstnd)
            {
                if (iter == 0)
                {
                    stpmx = 1.0;
                }
                else
                {
                    for (int i = 0; i < n; i++)
                    {
                        double a1 = d[i];
                        if (nbd[i] == 0) continue;
                        if (a1 < 0.0 && nbd[i] <= 2)
                        {
                            double a2 = l[i] - x[i];
                            if (a2 >= 0.0) stpmx = 0.0;
                            else if (a1 * stpmx < a2) stpmx = a2 / a1;
                        }
                        else if (a1 > 0.0 && nbd[i] >= 2)
                        {
                            double a2 = u[i] - x[i];
                            if (a2 <= 0.0) stpmx = 0.0;
                            else if (a1 * stpmx > a2) stpmx = a2 / a1;
                        }
                    }
                }
            }
            stp = iter == 0 && !boxed ? Math.Min(1.0 / dnorm, stpmx) : 1.0;

            Array.Copy(x, t, n);
            Array.Copy(g, r, n);
            fold = f;
            gd = CrossProductKernel.Dot(g, 0, d, 0, n);
            gdold = gd;
            // The direction is not a descent direction, which can only arise from rounding.
            if (!(gd < 0.0)) return false;

            var status = search.Start(f, gd, stp, Ftol, Gtol, Xtol, 0.0, stpmx);
            int ifun = 0;
            var gNew = g;
            while (true)
            {
                if (status == MoreThuenteStatus.Converged || status == MoreThuenteStatus.Warning)
                {
                    if (!ReferenceEquals(gNew, g)) Array.Copy(gNew, g, n);
                    return true;
                }
                if (status == MoreThuenteStatus.Error) return false;

                ifun++;
                if (ifun - 1 >= MaxLineSearchSteps) return false;
                if (stp == 1.0)
                {
                    Array.Copy(z, x, n);
                }
                else
                {
                    for (int i = 0; i < n; i++)
                    {
                        x[i] = stp * d[i] + t[i];
                        if (nbd[i] == 1 || nbd[i] == 2) x[i] = Math.Max(x[i], l[i]);
                        if (nbd[i] == 2 || nbd[i] == 3) x[i] = Math.Min(x[i], u[i]);
                    }
                }

                f = EvaluateObjective(x, ref cancel);
                if (cancel) return false;
                double[]? trial = Tools.IsFinite(f) ? EvaluateGradient(x, ref cancel) : null;
                if (cancel) return false;
                if (trial == null)
                {
                    // A non-finite objective or gradient: restart the search with half the step as its maximum.
                    stp *= 0.5;
                    status = search.Start(fold, gdold, stp, Ftol, Gtol, Xtol, 0.0, stp);
                    continue;
                }
                gNew = trial;
                gd = CrossProductKernel.Dot(gNew, 0, d, 0, n);
                status = search.Next(f, gd, ref stp);
            }
        }

        /// <summary>
        /// Adds the newest correction pair to the limited-memory matrices (reference routine <c>matupd</c>).
        /// </summary>
        private static void UpdateCorrections(int n, int m, double[] ws, double[] wy, double[] sy, double[] ss, double[] d,
            double[] r, ref int itail, int iupdat, ref int col, ref int head, ref double theta, double rr, double dr,
            double stp, double dtd)
        {
            // Set the pointers for the matrices ws and wy.
            if (iupdat <= m)
            {
                col = iupdat;
                itail = (head + iupdat - 1) % m;
            }
            else
            {
                itail = (itail + 1) % m;
                head = (head + 1) % m;
            }
            Array.Copy(d, 0, ws, itail * n, n);
            Array.Copy(r, 0, wy, itail * n, n);
            theta = rr / dr;

            // Move the old information when the memory is full.
            if (iupdat > m)
            {
                for (int j = 1; j < col; j++)
                {
                    // The upper triangle of ss and the lower triangle of sy move one column to the left.
                    Array.Copy(ss, 1 + m * j, ss, m * (j - 1), j);
                    Array.Copy(sy, j + m * j, sy, (j - 1) + m * (j - 1), col - j);
                }
            }

            // Add the new information: the last row of sy and the last column of ss.
            int pointr = head;
            for (int j = 0; j < col - 1; j++)
            {
                sy[col - 1 + m * j] = CrossProductKernel.Dot(d, 0, wy, pointr * n, n);
                ss[j + m * (col - 1)] = CrossProductKernel.Dot(ws, pointr * n, d, 0, n);
                pointr = (pointr + 1) % m;
            }
            ss[(col - 1) + m * (col - 1)] = stp == 1.0 ? dtd : stp * stp * dtd;
            sy[(col - 1) + m * (col - 1)] = dr;
        }

        /// <summary>
        /// Minimizes the model over the free variables and projects the result onto the box
        /// (reference routine <c>subsm</c>, with the L-BFGS-B 3.0 projection).
        /// </summary>
        /// <remarks>On entry x is the Cauchy point and d the reduced gradient; on exit x is the end point of the
        /// subspace step and d the subspace Newton direction. xx and gg are the current iterate and gradient.</remarks>
        /// <returns>Zero, or a nonzero value when a triangular system is singular.</returns>
        private static int SubspaceMinimization(int n, int m, int nsub, int[] ind, double[] l, double[] u, int[] nbd,
            double[] x, double[] d, double[] xp, double[] ws, double[] wy, double theta, double[] xx, double[] gg,
            int col, int head, double[] wv, double[] wn)
        {
            if (nsub <= 0) return 0;

            // Compute wv = W'Zd.
            int pointr = head;
            for (int i = 0; i < col; i++)
            {
                double temp1 = 0.0, temp2 = 0.0;
                for (int j = 0; j < nsub; j++)
                {
                    int k = ind[j];
                    temp1 += wy[k + n * pointr] * d[j];
                    temp2 += ws[k + n * pointr] * d[j];
                }
                wv[i] = temp1;
                wv[col + i] = theta * temp2;
                pointr = (pointr + 1) % m;
            }

            // Compute wv := K^(-1) wv.
            int m2 = 2 * m, col2 = 2 * col;
            if (TriangularSolve(wn, 0, m2, col2, wv, 0, true) != 0) return 1;
            for (int i = 0; i < col; i++) wv[i] = -wv[i];
            if (TriangularSolve(wn, 0, m2, col2, wv, 0, false) != 0) return 1;

            // Compute d = (1/θ)d + (1/θ²)Z'W wv.
            pointr = head;
            for (int jy = 0; jy < col; jy++)
            {
                int js = col + jy;
                for (int i = 0; i < nsub; i++)
                {
                    int k = ind[i];
                    d[i] += wy[k + n * pointr] * wv[jy] / theta + ws[k + n * pointr] * wv[js];
                }
                pointr = (pointr + 1) % m;
            }
            double scale = 1.0 / theta;
            for (int i = 0; i < nsub; i++) d[i] *= scale;

            // Project the subspace minimizer onto the box.
            bool projected = false;
            Array.Copy(x, xp, n);
            for (int i = 0; i < nsub; i++)
            {
                int k = ind[i];
                double dk = d[i], xk = x[k];
                if (nbd[k] == 0)
                {
                    x[k] = xk + dk;
                }
                else if (nbd[k] == 1)
                {
                    x[k] = Math.Max(l[k], xk + dk);
                    if (x[k] == l[k]) projected = true;
                }
                else if (nbd[k] == 2)
                {
                    xk = Math.Max(l[k], xk + dk);
                    x[k] = Math.Min(u[k], xk);
                    if (x[k] == l[k] || x[k] == u[k]) projected = true;
                }
                else
                {
                    x[k] = Math.Min(u[k], xk + dk);
                    if (x[k] == u[k]) projected = true;
                }
            }
            if (!projected) return 0;

            // Keep the projected point when it gives a descent direction.
            double ddp = 0.0;
            for (int i = 0; i < n; i++) ddp += (x[i] - xx[i]) * gg[i];
            if (!(ddp > 0.0)) return 0;

            // Otherwise truncate the subspace step at the first bound.
            Array.Copy(xp, x, n);
            double alpha = 1.0, temp = 1.0;
            int ibd = 0;
            for (int i = 0; i < nsub; i++)
            {
                int k = ind[i];
                double dk = d[i];
                if (nbd[k] == 0) continue;
                if (dk < 0.0 && nbd[k] <= 2)
                {
                    double t2 = l[k] - x[k];
                    if (t2 >= 0.0) temp = 0.0;
                    else if (dk * alpha < t2) temp = t2 / dk;
                }
                else if (dk > 0.0 && nbd[k] >= 2)
                {
                    double t2 = u[k] - x[k];
                    if (t2 <= 0.0) temp = 0.0;
                    else if (dk * alpha > t2) temp = t2 / dk;
                }
                if (temp < alpha)
                {
                    alpha = temp;
                    ibd = i;
                }
            }
            if (alpha < 1.0)
            {
                double dk = d[ibd];
                int k = ind[ibd];
                if (dk > 0.0)
                {
                    x[k] = u[k];
                    d[ibd] = 0.0;
                }
                else if (dk < 0.0)
                {
                    x[k] = l[k];
                    d[ibd] = 0.0;
                }
            }
            for (int i = 0; i < nsub; i++) x[ind[i]] += alpha * d[i];
            return 0;
        }

        /// <summary>
        /// Factors the leading n-by-n block of a column-major symmetric matrix as U'U, in its upper triangle.
        /// </summary>
        /// <param name="a">The matrix storage.</param>
        /// <param name="offset">The index of the block's first element.</param>
        /// <param name="lda">The leading dimension.</param>
        /// <param name="n">The order of the block.</param>
        /// <returns>Zero, or j + 1 when the leading minor of order j + 1 is not positive.</returns>
        private static int Cholesky(double[] a, int offset, int lda, int n)
        {
            for (int j = 0; j < n; j++)
            {
                int cj = offset + lda * j;
                double ajj = a[cj + j] - CrossProductKernel.Dot(a, cj, a, cj, j);
                if (!(ajj > 0.0))
                {
                    a[cj + j] = ajj;
                    return j + 1;
                }
                ajj = Math.Sqrt(ajj);
                a[cj + j] = ajj;
                for (int i = j + 1; i < n; i++)
                {
                    int ci = offset + lda * i;
                    a[ci + j] = (a[ci + j] - CrossProductKernel.Dot(a, cj, a, ci, j)) / ajj;
                }
            }
            return 0;
        }

        /// <summary>
        /// Solves U'x = b or Ux = b in place for the upper triangular leading n-by-n block of a column-major matrix.
        /// </summary>
        /// <param name="a">The matrix storage, starting at its first element.</param>
        /// <param name="offset">The index of the block's first element.</param>
        /// <param name="lda">The leading dimension.</param>
        /// <param name="n">The order of the block.</param>
        /// <param name="b">The right-hand side, overwritten by the solution.</param>
        /// <param name="bo">The index of the first element of b.</param>
        /// <param name="transpose">True to solve U'x = b, false to solve Ux = b.</param>
        /// <returns>Zero, or i + 1 when the diagonal element i is zero.</returns>
        private static int TriangularSolve(double[] a, int offset, int lda, int n, double[] b, int bo, bool transpose)
        {
            for (int i = 0; i < n; i++)
                if (a[offset + i + lda * i] == 0.0) return i + 1;
            if (transpose)
            {
                for (int i = 0; i < n; i++)
                {
                    int ci = offset + lda * i;
                    b[bo + i] = (b[bo + i] - CrossProductKernel.Dot(a, ci, b, bo, i)) / a[ci + i];
                }
            }
            else
            {
                for (int i = n - 1; i >= 0; i--)
                {
                    double sum = b[bo + i];
                    for (int k = i + 1; k < n; k++) sum -= a[offset + i + lda * k] * b[bo + k];
                    b[bo + i] = sum / a[offset + i + lda * i];
                }
            }
            return 0;
        }

        /// <summary>The state of the More-Thuente search after a call.</summary>
        private enum MoreThuenteStatus
        {
            /// <summary>Evaluate the function and derivative at the returned step.</summary>
            Evaluate,
            /// <summary>The step satisfies the strong Wolfe conditions.</summary>
            Converged,
            /// <summary>Rounding, the interval tolerance or a step bound stops the search; the step is accepted.</summary>
            Warning,
            /// <summary>The input is invalid.</summary>
            Error
        }

        /// <summary>
        /// The More-Thuente line search with safeguarded cubic and quadratic steps (reference routines <c>dcsrch</c>
        /// and <c>dcstep</c>, MINPACK-2).
        /// </summary>
        [Serializable]
        private sealed class MoreThuenteSearch
        {
            private bool _brackt;
            private int _stage;
            private double _ginit, _gtest, _gx, _gy, _finit, _fx, _fy, _stx, _sty, _stmin, _stmax, _width, _width1;
            private double _ftol, _gtol, _xtol, _stpmin, _stpmax;

            /// <summary>Starts a search at step zero with the given value and slope.</summary>
            /// <returns>Evaluate, or Error for invalid input.</returns>
            public MoreThuenteStatus Start(double f, double g, double stp, double ftol, double gtol, double xtol, double stpmin, double stpmax)
            {
                if (stp < stpmin || stp > stpmax || g >= 0.0 || ftol < 0.0 || gtol < 0.0 || xtol < 0.0 || stpmin < 0.0 || stpmax < stpmin)
                    return MoreThuenteStatus.Error;
                _ftol = ftol; _gtol = gtol; _xtol = xtol; _stpmin = stpmin; _stpmax = stpmax;
                _brackt = false;
                _stage = 1;
                _finit = f;
                _ginit = g;
                _gtest = ftol * _ginit;
                _width = stpmax - stpmin;
                _width1 = _width / 0.5;
                _stx = 0.0; _fx = _finit; _gx = _ginit;
                _sty = 0.0; _fy = _finit; _gy = _ginit;
                _stmin = 0.0;
                _stmax = stp + 4.0 * stp;
                return MoreThuenteStatus.Evaluate;
            }

            /// <summary>Takes the value and slope at the current step and returns the next step or a verdict.</summary>
            public MoreThuenteStatus Next(double f, double g, ref double stp)
            {
                const double xtrapl = 1.1, xtrapu = 4.0;
                double ftest = _finit + stp * _gtest;
                if (_stage == 1 && f <= ftest && g >= 0.0) _stage = 2;

                var status = MoreThuenteStatus.Evaluate;
                if (_brackt && (stp <= _stmin || stp >= _stmax)) status = MoreThuenteStatus.Warning;
                if (_brackt && _stmax - _stmin <= _xtol * _stmax) status = MoreThuenteStatus.Warning;
                if (stp == _stpmax && f <= ftest && g <= _gtest) status = MoreThuenteStatus.Warning;
                if (stp == _stpmin && (f > ftest || g >= _gtest)) status = MoreThuenteStatus.Warning;
                if (f <= ftest && Math.Abs(g) <= _gtol * (-_ginit)) status = MoreThuenteStatus.Converged;
                if (status != MoreThuenteStatus.Evaluate) return status;

                if (_stage == 1 && f <= _fx && f > ftest)
                {
                    // Use the modified function to compute the step while the function is above the Armijo line.
                    double fm = f - stp * _gtest, fxm = _fx - _stx * _gtest, fym = _fy - _sty * _gtest;
                    double gm = g - _gtest, gxm = _gx - _gtest, gym = _gy - _gtest;
                    Step(ref _stx, ref fxm, ref gxm, ref _sty, ref fym, ref gym, ref stp, fm, gm, ref _brackt, _stmin, _stmax);
                    _fx = fxm + _stx * _gtest;
                    _fy = fym + _sty * _gtest;
                    _gx = gxm + _gtest;
                    _gy = gym + _gtest;
                }
                else
                {
                    Step(ref _stx, ref _fx, ref _gx, ref _sty, ref _fy, ref _gy, ref stp, f, g, ref _brackt, _stmin, _stmax);
                }

                // Decide whether a bisection step is needed.
                if (_brackt)
                {
                    if (Math.Abs(_sty - _stx) >= 0.66 * _width1) stp = _stx + 0.5 * (_sty - _stx);
                    _width1 = _width;
                    _width = Math.Abs(_sty - _stx);
                }

                // Set the minimum and maximum steps allowed.
                if (_brackt)
                {
                    _stmin = Math.Min(_stx, _sty);
                    _stmax = Math.Max(_stx, _sty);
                }
                else
                {
                    _stmin = stp + xtrapl * (stp - _stx);
                    _stmax = stp + xtrapu * (stp - _stx);
                }
                stp = Math.Max(stp, _stpmin);
                stp = Math.Min(stp, _stpmax);

                // With no further progress possible, use the best step obtained.
                if (_brackt && (stp <= _stmin || stp >= _stmax || _stmax - _stmin <= _xtol * _stmax)) stp = _stx;
                return MoreThuenteStatus.Evaluate;
            }

            /// <summary>Computes a safeguarded step and updates the interval that contains a step satisfying the
            /// Wolfe conditions (reference routine <c>dcstep</c>).</summary>
            private static void Step(ref double stx, ref double fx, ref double dx, ref double sty, ref double fy, ref double dy,
                ref double stp, double fp, double dp, ref bool brackt, double stpmin, double stpmax)
            {
                double sgnd = dp * (dx / Math.Abs(dx));
                double stpf;
                if (fp > fx)
                {
                    // A higher function value: the minimum is bracketed.
                    double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
                    double s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
                    double gamma = s * Math.Sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
                    if (stp < stx) gamma = -gamma;
                    double p = (gamma - dx) + theta;
                    double q = ((gamma - dx) + gamma) + dp;
                    double r = p / q;
                    double stpc = stx + r * (stp - stx);
                    double stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.0) * (stp - stx);
                    stpf = Math.Abs(stpc - stx) < Math.Abs(stpq - stx) ? stpc : stpc + (stpq - stpc) / 2.0;
                    brackt = true;
                }
                else if (sgnd < 0.0)
                {
                    // Derivatives of opposite sign: the minimum is bracketed.
                    double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
                    double s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
                    double gamma = s * Math.Sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
                    if (stp > stx) gamma = -gamma;
                    double p = (gamma - dp) + theta;
                    double q = ((gamma - dp) + gamma) + dx;
                    double r = p / q;
                    double stpc = stp + r * (stx - stp);
                    double stpq = stp + (dp / (dp - dx)) * (stx - stp);
                    stpf = Math.Abs(stpc - stp) > Math.Abs(stpq - stp) ? stpc : stpq;
                    brackt = true;
                }
                else if (Math.Abs(dp) < Math.Abs(dx))
                {
                    // Same sign and decreasing magnitude of the derivative.
                    double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
                    double s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dx), Math.Abs(dp)));
                    double gamma = s * Math.Sqrt(Math.Max(0.0, (theta / s) * (theta / s) - (dx / s) * (dp / s)));
                    if (stp > stx) gamma = -gamma;
                    double p = (gamma - dp) + theta;
                    double q = (gamma + (dx - dp)) + gamma;
                    double r = p / q;
                    double stpc;
                    if (r < 0.0 && gamma != 0.0) stpc = stp + r * (stx - stp);
                    else if (stp > stx) stpc = stpmax;
                    else stpc = stpmin;
                    double stpq = stp + (dp / (dp - dx)) * (stx - stp);
                    if (brackt)
                    {
                        stpf = Math.Abs(stpc - stp) < Math.Abs(stpq - stp) ? stpc : stpq;
                        if (stp > stx) stpf = Math.Min(stp + 0.66 * (sty - stp), stpf);
                        else stpf = Math.Max(stp + 0.66 * (sty - stp), stpf);
                    }
                    else
                    {
                        stpf = Math.Abs(stpc - stp) > Math.Abs(stpq - stp) ? stpc : stpq;
                        stpf = Math.Min(stpmax, stpf);
                        stpf = Math.Max(stpmin, stpf);
                    }
                }
                else
                {
                    // Same sign and non-decreasing magnitude of the derivative.
                    if (brackt)
                    {
                        double theta = 3.0 * (fp - fy) / (sty - stp) + dy + dp;
                        double s = Math.Max(Math.Abs(theta), Math.Max(Math.Abs(dy), Math.Abs(dp)));
                        double gamma = s * Math.Sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
                        if (stp > sty) gamma = -gamma;
                        double p = (gamma - dp) + theta;
                        double q = ((gamma - dp) + gamma) + dy;
                        double r = p / q;
                        stpf = stp + r * (sty - stp);
                    }
                    else if (stp > stx) stpf = stpmax;
                    else stpf = stpmin;
                }

                // Update the interval that contains a minimizer.
                if (fp > fx)
                {
                    sty = stp; fy = fp; dy = dp;
                }
                else
                {
                    if (sgnd < 0.0)
                    {
                        sty = stx; fy = fx; dy = dx;
                    }
                    stx = stp; fx = fp; dx = dp;
                }
                stp = stpf;
            }
        }

    }
}
