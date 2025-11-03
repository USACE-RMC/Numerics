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

namespace Numerics.Mathematics.Integration
{
    /// <summary>
    /// A class that performs adaptive Simpson's integration in two dimensions. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Adaptive Simpson's rule is extended to 2D by integrating over a rectangular region by subdividing the domain
    /// recursively in both the x and y directions.
    /// </para>
    /// </remarks>
    [Serializable]
    public class AdaptiveSimpsonsRule2D : Integrator
    {

        /// <summary>
        /// Constructs a new adaptive Simpson's rule for 2D integration.
        /// </summary>
        /// <param name="function">The function to integrate.</param>
        /// <param name="minX">The minimum x-value under which the integral must be computed.</param>
        /// <param name="maxX">The maximum x-value under which the integral must be computed.</param>
        /// <param name="minY">The minimum y-value under which the integral must be computed.</param>
        /// <param name="maxY">The maximum y-value under which the integral must be computed.</param>
        public AdaptiveSimpsonsRule2D(Func<double, double, double> function, double minX, double maxX, double minY, double maxY)
        {
            if (maxX <= minX) throw new ArgumentNullException(nameof(maxX), "The maximum x-value cannot be less than or equal to the minimum x-value.");
            if (maxY <= minY) throw new ArgumentNullException(nameof(maxY), "The maximum y-value cannot be less than or equal to the minimum y-value.");
            Function = function ?? throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            ax = minX;
            bx = maxX;
            ay = minY;
            by = maxY;
        }

        // 2D bounds
        private double ax, bx, ay, by;
        private double _squaredError;

        /// <summary>
        /// The 2D function to integrate.
        /// </summary>
        public Func<double, double, double> Function { get; }

        /// <summary>
        /// The minimum x-value under which the integral must be computed.
        /// </summary>
        public double MinX => ax;

        /// <summary>
        /// The maximum x-value under which the integral must be computed.
        /// </summary>
        public double MaxX => bx;

        /// <summary>
        /// The minimum y-value under which the integral must be computed.
        /// </summary>
        public double MinY => ay;

        /// <summary>
        /// The maximum y-value under which the integral must be computed.
        /// </summary>
        public double MaxY => by;

        /// <summary>
        /// The minimum recursion depth. Default = 0.
        /// </summary>
        public int MinDepth { get; set; } = 0;

        /// <summary>
        /// The maximum recursion depth. Default = 100.
        /// </summary>
        public int MaxDepth { get; set; } = 100;

        /// <summary>
        /// Returns an approximate measure of the standard error of the integration. 
        /// </summary>
        public double StandardError { get; private set; }

        /// <inheritdoc/>
        public override void Integrate()
        {
            _squaredError = 0;
            StandardError = 0;
            ClearResults();
            Validate();

            try
            {
                // Initial evaluation of Simpson's Rule on the entire 2D domain
                double mx = (ax + bx) / 2d;
                double my = (ay + by) / 2d;

                // Evaluate at the four corners
                double f_ax_ay = Function(ax, ay);
                double f_bx_ay = Function(bx, ay);
                double f_ax_by = Function(ax, by);
                double f_bx_by = Function(bx, by);

                // Evaluate at edge midpoints
                double f_mx_ay = Function(mx, ay);
                double f_mx_by = Function(mx, by);
                double f_ax_my = Function(ax, my);
                double f_bx_my = Function(bx, my);

                // Evaluate at the center
                double f_mx_my = Function(mx, my);

                FunctionEvaluations += 9; // Count all 9 evaluations

                // 2D Simpson's rule: (hx/3) * (hy/3) * weighted sum
                // Apply Simpson's in x for each y, then Simpson's in y
                double hx = (bx - ax) / 2d;
                double hy = (by - ay) / 2d;
                double whole = (hx / 3d) * (hy / 3d) * (
                    (f_ax_ay + 4 * f_mx_ay + f_bx_ay) +
                    4 * (f_ax_my + 4 * f_mx_my + f_bx_my) +
                    (f_ax_by + 4 * f_mx_by + f_bx_by)
                );

                // Recursively sub-divide
                Result = AdaptiveSimpsons2D(Function, ax, bx, ay, by,
                    f_ax_ay, f_bx_ay, f_ax_by, f_bx_by,
                    f_mx_ay, f_mx_by, f_ax_my, f_bx_my, f_mx_my,
                    MaxDepth, whole, mx, my, ax, bx, ay, by);

                // Standard error calculated after recursion completes
                StandardError = Math.Sqrt(_squaredError);

                if (FunctionEvaluations >= MaxFunctionEvaluations)
                {
                    Status = IntegrationStatus.MaximumFunctionEvaluationsReached;
                }
                else
                {
                    Status = IntegrationStatus.Success;
                }
            }
            catch (Exception ex)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw ex;
            }
        }

        /// <summary>
        /// A helper function to perform adaptive Simpson's rule in two dimensions.
        /// Recursively subdivides the 2D integration domain and applies Simpson's rule.
        /// </summary>
        /// <param name="f">The 2D function to integrate.</param>
        /// <param name="ax">The minimum X value of the current rectangular interval.</param>
        /// <param name="bx">The maximum X value of the current rectangular interval.</param>
        /// <param name="ay">The minimum Y value of the current rectangular interval.</param>
        /// <param name="by">The maximum Y value of the current rectangular interval.</param>
        /// <param name="f_ax_ay">The function evaluated at the lower-left corner (ax, ay).</param>
        /// <param name="f_bx_ay">The function evaluated at the lower-right corner (bx, ay).</param>
        /// <param name="f_ax_by">The function evaluated at the upper-left corner (ax, by).</param>
        /// <param name="f_bx_by">The function evaluated at the upper-right corner (bx, by).</param>
        /// <param name="f_mx_ay">The function evaluated at the bottom edge midpoint (mx, ay).</param>
        /// <param name="f_mx_by">The function evaluated at the top edge midpoint (mx, by).</param>
        /// <param name="f_ax_my">The function evaluated at the left edge midpoint (ax, my).</param>
        /// <param name="f_bx_my">The function evaluated at the right edge midpoint (bx, my).</param>
        /// <param name="f_mx_my">The function evaluated at the center point (mx, my).</param>
        /// <param name="depth">The current recursion depth remaining.</param>
        /// <param name="whole">The Simpson's rule evaluation over the current 2D rectangular interval.</param>
        /// <param name="mx">The midpoint of the current X interval.</param>
        /// <param name="my">The midpoint of the current Y interval.</param>
        /// <param name="ax0">The original lower bound of the X domain (used for Richardson scaling).</param>
        /// <param name="bx0">The original upper bound of the X domain (used for Richardson scaling).</param>
        /// <param name="ay0">The original lower bound of the Y domain (used for Richardson scaling).</param>
        /// <param name="by0">The original upper bound of the Y domain (used for Richardson scaling).</param>
        /// <returns>
        /// An evaluation of the 2D integral using adaptive Simpson's rule with error less than the specified tolerance. 
        /// This is accomplished by subdividing the rectangular interval into four quadrants until the error between 
        /// the previous evaluation and the current evaluation is sufficiently small.
        /// </returns>
        private double AdaptiveSimpsons2D(
            Func<double, double, double> f,
            double ax, double bx, double ay, double by,
            double f_ax_ay, double f_bx_ay, double f_ax_by, double f_bx_by,
            double f_mx_ay, double f_mx_by, double f_ax_my, double f_bx_my, double f_mx_my,
            int depth, double whole, double mx, double my,
            double ax0, double bx0, double ay0, double by0)
        {
            // Midpoints for subdividing into four quadrants
            double lmx = (ax + mx) / 2d;  // left midpoint in X
            double rmx = (mx + bx) / 2d;  // right midpoint in X
            double lmy = (ay + my) / 2d;  // lower midpoint in Y
            double rmy = (my + by) / 2d;  // upper midpoint in Y

            // Evaluate at new edge midpoints and quadrant centers
            double f_lmx_ay = f(lmx, ay);
            double f_rmx_ay = f(rmx, ay);
            double f_lmx_by = f(lmx, by);
            double f_rmx_by = f(rmx, by);

            double f_ax_lmy = f(ax, lmy);
            double f_bx_lmy = f(bx, lmy);
            double f_ax_rmy = f(ax, rmy);
            double f_bx_rmy = f(bx, rmy);

            double f_mx_lmy = f(mx, lmy);
            double f_mx_rmy = f(mx, rmy);
            double f_lmx_my = f(lmx, my);
            double f_rmx_my = f(rmx, my);

            double f_lmx_lmy = f(lmx, lmy);
            double f_rmx_lmy = f(rmx, lmy);
            double f_lmx_rmy = f(lmx, rmy);
            double f_rmx_rmy = f(rmx, rmy);

            FunctionEvaluations += 16;  // Count the 16 new function evaluations

            // Calculate Simpson's rule for each of the four quadrants
            double hx_left = (mx - ax) / 2d;
            double hx_right = (bx - mx) / 2d;
            double hy_lower = (my - ay) / 2d;
            double hy_upper = (by - my) / 2d;

            // Bottom-left quadrant: (ax, mx) x (ay, my)
            double q1 = (hx_left / 3d) * (hy_lower / 3d) * (
                (f_ax_ay + 4 * f_lmx_ay + f_mx_ay) +
                4 * (f_ax_lmy + 4 * f_lmx_lmy + f_mx_lmy) +
                (f_ax_my + 4 * f_lmx_my + f_mx_my)
            );

            // Bottom-right quadrant: (mx, bx) x (ay, my)
            double q2 = (hx_right / 3d) * (hy_lower / 3d) * (
                (f_mx_ay + 4 * f_rmx_ay + f_bx_ay) +
                4 * (f_mx_lmy + 4 * f_rmx_lmy + f_bx_lmy) +
                (f_mx_my + 4 * f_rmx_my + f_bx_my)
            );

            // Top-left quadrant: (ax, mx) x (my, by)
            double q3 = (hx_left / 3d) * (hy_upper / 3d) * (
                (f_ax_my + 4 * f_lmx_my + f_mx_my) +
                4 * (f_ax_rmy + 4 * f_lmx_rmy + f_mx_rmy) +
                (f_ax_by + 4 * f_lmx_by + f_mx_by)
            );

            // Top-right quadrant: (mx, bx) x (my, by)
            double q4 = (hx_right / 3d) * (hy_upper / 3d) * (
                (f_mx_my + 4 * f_rmx_my + f_bx_my) +
                4 * (f_mx_rmy + 4 * f_rmx_rmy + f_bx_rmy) +
                (f_mx_by + 4 * f_rmx_by + f_bx_by)
            );

            double sum = q1 + q2 + q3 + q4;

            // Calculate the error
            double error = sum - whole;
            double delta = error / 15d; // Richardson correction

            // Richardson-based convergence tolerance (scaled by domain size ratio)
            double areaRatio = ((bx - ax) * (by - ay)) / ((bx0 - ax0) * (by0 - ay0));
            double toleranceScaled = 15d * RelativeTolerance * areaRatio;

            // Absolute and Relative tolerance checks
            bool absoluteToleranceReached = Math.Abs(error) <= AbsoluteTolerance;
            bool relativeToleranceReached = Math.Abs(error) <= toleranceScaled;

            // Check if convergence criteria are met
            if (depth <= 0 ||
                Math.Abs(ax - bx) <= Tools.DoubleMachineEpsilon ||
                Math.Abs(ay - by) <= Tools.DoubleMachineEpsilon ||
                FunctionEvaluations >= MaxFunctionEvaluations ||
                (FunctionEvaluations >= MinFunctionEvaluations &&
                 depth <= MaxDepth - MinDepth &&
                 (absoluteToleranceReached || relativeToleranceReached)))
            {
                // Convergence is reached, accumulate the error
                _squaredError += delta * delta;
                return sum + delta;
            }
            else
            {
                // Recursively subdivide the four quadrants
                var q1_result = AdaptiveSimpsons2D(f, ax, mx, ay, my,
                    f_ax_ay, f_mx_ay, f_ax_my, f_mx_my,
                    f_lmx_ay, f_lmx_my, f_ax_lmy, f_mx_lmy, f_lmx_lmy,
                    depth - 1, q1, lmx, lmy, ax0, bx0, ay0, by0);

                var q2_result = AdaptiveSimpsons2D(f, mx, bx, ay, my,
                    f_mx_ay, f_bx_ay, f_mx_my, f_bx_my,
                    f_rmx_ay, f_rmx_my, f_mx_lmy, f_bx_lmy, f_rmx_lmy,
                    depth - 1, q2, rmx, lmy, ax0, bx0, ay0, by0);

                var q3_result = AdaptiveSimpsons2D(f, ax, mx, my, by,
                    f_ax_my, f_mx_my, f_ax_by, f_mx_by,
                    f_lmx_my, f_lmx_by, f_ax_rmy, f_mx_rmy, f_lmx_rmy,
                    depth - 1, q3, lmx, rmy, ax0, bx0, ay0, by0);

                var q4_result = AdaptiveSimpsons2D(f, mx, bx, my, by,
                    f_mx_my, f_bx_my, f_mx_by, f_bx_by,
                    f_rmx_my, f_rmx_by, f_mx_rmy, f_bx_rmy, f_rmx_rmy,
                    depth - 1, q4, rmx, rmy, ax0, bx0, ay0, by0);

                return q1_result + q2_result + q3_result + q4_result;
            }
        }
    }
}
