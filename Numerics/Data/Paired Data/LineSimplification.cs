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
