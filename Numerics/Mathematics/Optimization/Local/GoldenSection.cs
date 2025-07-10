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

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The Golden-Section optimization algorithm. The function need not be differentiable, and no derivatives are taken.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This class contains methods for finding the minimum or maximum of a function using the Golden-Section method.
    /// This method operates by successively narrowing the range of values on the specified interval. When searching for a
    /// minimum, the central point is always less than or equal to the outer points, assuring that a minimum is contained 
    /// between the outer points.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </description></item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Golden-section_search"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class GoldenSection : Optimizer
    {

        /// <summary>
        /// Construct a new Golden-Section optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="lowerBound">The lower bound (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBound">The upper bound (inclusive) of the interval containing the optimal point.</param>
        public GoldenSection(Func<double, double> objectiveFunction, double lowerBound, double upperBound) : base((x) => objectiveFunction(x[0]), 1)
        {
            // validate inputs
            if (upperBound < lowerBound)
            {
                throw new ArgumentOutOfRangeException("upperBound", "The upper bound cannot be less than the lower bound.");
            }
            LowerBound = lowerBound;
            UpperBound = upperBound;
        }

        /// <summary>
        /// The lower bound (inclusive) of the interval containing the optimal point. 
        /// </summary>
        public double LowerBound { get; private set; }

        /// <summary>
        /// The upper bound (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double UpperBound { get; private set; }

        /// <inheritdoc/>
        protected override void Optimize()
        {
            // Define variables
            bool cancel = false;
            // Bracket problem
            double ax = LowerBound, bx = 0.5 * (UpperBound + LowerBound), cx = UpperBound;
            double R = 0.61803399, C = 1.0 - R;
            double x1, x2;
            double x0 = ax;
            double x3 = cx;
            if (Math.Abs(cx - bx) > Math.Abs(bx - ax))
            {
                x1 = bx;
                x2 = bx + C * (cx - bx);
            }
            else
            {
                x2 = bx;
                x1 = bx - C * (bx - ax);
            }
            double f1 = Evaluate(new double[] { x1 }, ref cancel);
            double f2 = Evaluate(new double[] { x2 }, ref cancel);
            while (Math.Abs(x3 - x0) > AbsoluteTolerance || Math.Abs(x3 - x0) > RelativeTolerance * (Math.Abs(x1) + Math.Abs(x2)))
            {
                if (f2 < f1)
                {
                    x0 = x1;
                    x1 = x2;
                    x2 = R * x2 + C * x3;
                    f1 = f2;
                    f2 = Evaluate(new double[] { x2 }, ref cancel);                 
                }
                else
                {
                    x3 = x2;
                    x2 = x1;
                    x1 = R * x1 + C * x0;
                    f2 = f1;
                    f1 = Evaluate(new double[] { x1 }, ref cancel);
                }
                if (cancel == true) return;
                Iterations += 1;
                if (Iterations >= MaxIterations)
                {
                    UpdateStatus(OptimizationStatus.MaximumIterationsReached);
                    return;
                }
            }
            UpdateStatus(OptimizationStatus.Success);
        }

    }
}