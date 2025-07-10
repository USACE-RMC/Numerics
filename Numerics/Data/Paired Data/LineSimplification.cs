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
    /// Class to perform line simplification using the RDP algorithm.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// <para>
    /// <b> Description: </b>
    /// The Ramer-Douglas-Peucker algorithm decimates a curve composed of line segments to a similar curve with fewer points. The purpose is, given a curve
    /// composed of line segments,, to find a similar curve with fewer points. "Dissimilar" in this case is defined based on the maximum distance
    /// between the original curve and the simplified curve (in this class, that maximum distance is user defined). The simplified curve consists of
    /// a subset of the points that defined the original curve.
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm"/>
    /// </remarks>
    public class LineSimplification
    {
        /// <summary>
        /// Implements the Ramer-Douglas-Peucker algorithm
        /// </summary>
        /// <param name="ordinates"> The list of input ordinates that make up the original curve </param>
        /// <param name="epsilon"> The max allowable distance between the original curve and outputted simplified curve </param>
        /// <param name="output"> The ordinates that make up simplified curve produced by the algorithm </param>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        /// <exception cref="Exception"></exception>
        public static void RamerDouglasPeucker(List<Ordinate> ordinates, double epsilon, ref List<Ordinate> output)
        {
            if (ordinates.Count < 2)
                throw new ArgumentOutOfRangeException("Not enough points to simplify");

            // Find the point with the maximum distance from line between the start and end
            double dmax = 0.0;
            int index = 0;
            int end = ordinates.Count - 1;
            for (int i = 1; i < end; ++i)
            {
                double d = PerpendicularDistance(ordinates[i], ordinates[0], ordinates[end]);
                if (d > dmax)
                {
                    index = i;
                    dmax = d;
                }
            }

            // If max distance is greater than epsilon, recursively simplify
            if (dmax > epsilon)
            {
                List<Ordinate> recResults1 = new List<Ordinate>();
                List<Ordinate> recResults2 = new List<Ordinate>();
                List<Ordinate> firstLine = ordinates.Take(index + 1).ToList();
                List<Ordinate> lastLine = ordinates.Skip(index).ToList();
                RamerDouglasPeucker(firstLine, epsilon, ref recResults1);
                RamerDouglasPeucker(lastLine, epsilon, ref recResults2);

                // build the result list
                output.AddRange(recResults1.Take(recResults1.Count - 1));
                output.AddRange(recResults2);
                if (output.Count < 2) throw new Exception("Problem assembling output");
            }
            else
            {
                // Just return start and end points
                output.Clear();
                output.Add(ordinates[0]);
                output.Add(ordinates[ordinates.Count - 1]);
            }
            
        }

        /// <summary>
        /// Calculates the distance of a point from a given line
        /// </summary>
        /// <param name="pt"> The point who's distance is being measured </param>
        /// <param name="lineStart"> The starting point of line </param>
        /// <param name="lineEnd"> The ending point of the line</param>
        /// <returns></returns>
        public static double PerpendicularDistance(Ordinate pt, Ordinate lineStart, Ordinate lineEnd)
        {
            double dx = lineEnd.X - lineStart.X;
            double dy = lineEnd.Y - lineStart.Y;

            // Normalize
            double mag = Math.Sqrt(dx * dx + dy * dy);
            if (mag > 0.0)
            {
                dx /= mag;
                dy /= mag;
            }
            double pvx = pt.X - lineStart.X;
            double pvy = pt.Y - lineStart.Y;

            // Get dot product (project pv onto normalized direction)
            double pvdot = dx * pvx + dy * pvy;

            // Scale line direction vector and subtract it from pv
            double ax = pvx - pvdot * dx;
            double ay = pvy - pvdot * dy;

            return Math.Sqrt(ax * ax + ay * ay);
        }
    }
}
