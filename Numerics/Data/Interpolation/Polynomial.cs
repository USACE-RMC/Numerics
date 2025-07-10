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
using System.Linq;

namespace Numerics.Data
{
    /// <summary>
    /// A class for polynomial interpolation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This method is the interpolation of a given bivariate data set by the polynomial of lowest possible degree
    /// that passes through the points of the dataset.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </description></item>
    /// <item><description>
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Polynomial_interpolation"/>
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class Polynomial : Interpolater
    {
        /// <summary>
        /// Construct new polynomial interpolation. 
        /// </summary>
        /// <param name="order">The polynomial order. There are order + 1 terms for each polynomial function.</param>
        /// <param name="xValues">List of x-values.</param>
        /// <param name="yValues">List of y-values.</param>
        /// <param name="sortOrder">The sort order of the x-values, either ascending or descending. Default = Ascending. </param>
        public Polynomial(int order, IList<double> xValues, IList<double> yValues, SortOrder sortOrder = SortOrder.Ascending) : base(xValues, yValues, sortOrder)
        {
            if (order >= Count) throw new ArgumentException(nameof(order), "The order must be less than the length of the x value list.");
            Order = order;
            
        }

        /// <summary>
        /// The error estimate for the most recent call to the interpolation function.
        /// </summary>
        public double Error { get; private set; }

        /// <summary>
        /// The polynomial order. There are order + 1 terms for each polynomial function. 
        /// </summary>
        public int Order { get; set; }

        /// <inheritdoc/>
        public override double BaseInterpolate(double x, int index)
        {
            // Given a value x, this routine returns an interpolate value y, and stores an error estimate. 
            // The return value is obtained by Degree-point polynomial interpolation on the subrange x[start..start + Order]
            if (index < 0 || index >= Count) index = 0;
            int mm = Order + 1, ns = 0;
            double y, den, dif, dift, ho, hp, w;
            int jl = index + mm > Count ? 0 : index;
            var xa = XValues.ToArray().Subset(jl);
            var ya = YValues.ToArray().Subset(jl);
            var c = new double[mm];
            var d = new double[mm];
            dif = Math.Abs(x - xa[0]);
            // Here we find the index ns of the closest table entry
            for (int i = 0; i < mm; i++)
            {
                if ((dift = Math.Abs(x - xa[i])) < dif)
                {
                    ns = i;
                    dif = dift;
                }
                // and initialize the tableau of c's and d's.
                c[i] = ya[i];
                d[i] = ya[i];
            }
            // This is the initial approximation to y.
            y = ya[ns--];
            // For each column in the tableau, 
            // we loop over the current c's and d's and update them.
            for (int m = 1; m < mm; m++)
            {
                for (int i = 0; i < mm - m; i++)
                {
                    ho = xa[i] - x;
                    hp = xa[i + m] - x;
                    w = c[i + 1] - d[i];
                    den = ho - hp;
                    den = w / den;
                    // Here the c's and d's are updated. 
                    d[i] = hp * den;
                    c[i] = ho * den;
                }
                y += (Error = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]));
                // After each column in the tableau is completed, we decide which correction, c or d, we
                // want to add to our accumulating value of y, i.e., which path to take through the tableau
                // - forking up or down. TWe do this in such a way as to take the most "straight line"
                // route through the tableau to its apex, updating ns accordingly to keep track of where
                // we are. This route keeps the partial approximations centered (insofar as possible) on 
                // the target x. The last dy added is thus the error indication. 
            }
            return y;
        }

    }
}
