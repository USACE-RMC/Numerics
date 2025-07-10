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

using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Integration
{
    /// <summary>
    /// A class for Miser, the recursive stratified sampling algorithm for multidimensional integration.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The Miser method aims to reduce overall error by concentrating integration points in the regions of highest variance. The smallest error estimate is
    /// obtained by allocating sample points in proportion to the standard deviation of the function in each sub-region. This algorithm bisects the integration
    /// region along one coordinate axis at each step. The direction is chose by examining all possible bisections and selecting whichever will minimize the combined variance.
    /// This procedure is repeated recursively for each of the two half spaces down to a user-specified depth where each sub-region is integrated using a plain Monte
    /// Carlo estimate.
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
    /// <see href="https://en.wikipedia.org/wiki/Monte_Carlo_integration"/>
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class Miser : Integrator
    {
        /// <summary>
        /// Creates a new Miser class, the recursive stratified sampling algorithm for multidimensional integration.
        /// </summary>
        /// <param name="function">The multidimensional function to integrate.</param>
        /// <param name="dimensions">The number of dimensions in the function to evaluate.</param>
        /// <param name="min">The minimum values under which the integral must be computed.</param>
        /// <param name="max">The maximum values under which the integral must be computed.</param>
        public Miser(Func<double[], double> function, int dimensions, IList<double> min, IList<double> max)
        {
            if (dimensions < 1) throw new ArgumentOutOfRangeException(nameof(dimensions), "There must be at least 1 dimension to evaluate.");

            // Check if the length of the min and max values equal the number of dimensions
            if (min.Count != dimensions || max.Count != dimensions)
            {
                throw new ArgumentOutOfRangeException(nameof(min), "The minimum and maximum values must be the same length as the number of dimensions.");
            }

            // Check if the minimum values are less than the maximum values
            for (int i = 0; i < min.Count; i++)
            {
                if (max[i] <= min[i])
                {
                    throw new ArgumentOutOfRangeException(nameof(max), "The maximum values cannot be less than or equal to the minimum values.");
                }
            }

            Function = function ?? throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            Dimensions = dimensions;
            Min = min.ToArray();
            Max = max.ToArray();
            Random = new MersenneTwister();
            _sobol = new SobolSequence(Dimensions);

        }

        private SobolSequence _sobol;

        /// <summary>
        /// The multidimensional function to integrate.
        /// </summary>
        public Func<double[], double> Function { get; }

        /// <summary>
        /// The number of dimensions in the function to evaluate.
        /// </summary>
        public int Dimensions { get; }

        /// <summary>
        /// The minimum values under which the integral must be computed.
        /// </summary>
        public double[] Min { get; }

        /// <summary>
        /// The maximum values under which the integral must be computed. 
        /// </summary>
        public double[] Max { get; }

        /// <summary>
        /// The random number generator to be used within the Monte Carlo integration.
        /// </summary>
        public Random Random { get; set; }

        /// <summary>
        /// The integration standard error. 
        /// </summary>
        public double StandardError { get; protected set; }

        /// <summary>
        /// The fraction of remaining function evaluations used at each stage to explore the variance of the function. Default = 0.1.
        /// </summary>
        public double Fraction { get; set; } = 0.1;

        /// <summary>
        /// The minimum number of points and function evaluations performed in each subregion. Default = 15. 
        /// </summary>
        public int MinimumNumberOfSubregionPoints { get; set; } = 15;

        /// <summary>
        /// A subregion is further bisected only if this number of function evaluations are available. Default = 4 * 15 = 60. 
        /// </summary>
        public int MinimumNumberOfBisections { get; set; } = 60;

        /// <summary>
        /// Dither should normally be set to 0.0, but can be set to 0.1 if the integrands active region fall on the boundary of a power-of-two subdivision of a region.
        /// </summary>
        public double Dither { get; set; } = 0.0;

        /// <summary>
        /// Determines whether to use a Sobol sequence or a pseudo-Random number generator. 
        /// </summary>
        public bool UseSobolSequence { get; set; } = true;

        /// <inheritdoc/>
        public override void Integrate()
        {
            ClearResults();
            Validate();

            try
            {
                // Create region array
                var regn = new double[2 * Dimensions];
                for (int i = 0; i < Dimensions; i++)
                    regn[i] = Min[i];
                for (int i = Dimensions; i < 2 * Dimensions; i++)
                    regn[i] = Max[i - Dimensions];

                double volume = 1;
                for (int i = 0; i < Dimensions; i++)
                    volume *= (Max[i] - Min[i]);

                // Compute
                double result = 0;
                double error = 0;

                miser(Function, regn, MaxFunctionEvaluations, Dither, ref result, ref error);
                Result = result * volume;
                StandardError = Math.Sqrt(error) * volume;

            }
            catch (Exception ex)
            {
                UpdateStatus(IntegrationStatus.Failure, ex);
            }

        }
        
        /// <summary>
        /// Helper function for Integrate() function that includes the Monte Carlo sampler for the method
        /// </summary>
        /// <param name="function">n dimensional function being evaluated</param>
        /// <param name="regn"> A vector consisting of ndim "lower-left" coordinates of the region followed by ndim "upper - right"
        /// coordinates. Specifies the rectangular volume by regn[0...2 * ndim -1]</param>
        /// <param name="npts"> Number of times tbe function is sampled </param>
        /// <param name="dith"> Should be normally set to 0, but can be est to 0.1 if the function's active region falls on the boundary of
        /// a power-of-two subdivision region. </param>
        /// <param name="ave"> The mean value of the function in the region </param>
        /// <param name="var"> An estimate of the statistical uncertainty of ave (square of standard deviation), i.e. variance</param>
        private void miser(Func<double[], double> function, double[] regn, int npts, double dith, ref double ave, ref double var)
        {
            // Monte Carlo samples a user-supplied ndim-dimensional function func in a rectangular volume
            // specified by regn[0..2 * ndim - 1], a vector consisting of ndim “lower - left” coordinates of the
            // region followed by ndim “upper - right” coordinates. The function is sampled a total of npts
            // times, at locations determined by the method of recursive stratified sampling. The mean value
            // of the function in the region is returned as ave; an estimate of the statistical uncertainty of ave
            // (square of standard deviation) is returned as var. The input parameter dith should normally
            // be set to zero, but can be set to(e.g.) 0.1 if func’s active region falls on the boundary of a
            // power-of-two subdivision of region.

            int MNPT = MinimumNumberOfSubregionPoints, MNBS = MinimumNumberOfBisections;
            double PFAC = Fraction, TINY = 1.0e-30, BIG = 1.0e30;

            // Here PFAC is the fraction of remaining function evaluations used at each stage to explore
            // the variance of func. At least MNPT function evaluations are performed in any terminal
            // subregion; a subregion is further bisected only if at least MNBS function evaluations are
            // available. We take MNBS = 4 * MNPT.

            int iran = 0;
            int j, jb, n, ndim, npre, nptl, nptr;
            double avel = 0, varl = 0, fracl, fval, rgl, rgm, rgr, s, sigl, siglb, sigr, sigrb;
            double sum, sumb, summ, summ2;

            ndim = Dimensions;
            var pt = new double[ndim];
            if (npts < MNBS)
            {
                // Too few points to bisect, do straight Monte Carlo
                summ = summ2 = 0.0;
                for (n = 0; n < npts; n++)
                {
                    ranpt(pt, regn);
                    fval = function(pt);
                    FunctionEvaluations++;
                    summ += fval;
                    summ2 += fval * fval;
                }
                ave = summ / npts;
                var = Math.Max(TINY, (summ2 - summ * summ / npts) / (npts * npts));
            }
            else
            {
                // Do the preliminary (uniform sampling)
                var rmid = new double[ndim];
                npre = Math.Max((int)(npts * PFAC), MNPT);
                var fmaxl = new double[ndim];
                var fmaxr = new double[ndim];
                var fminl = new double[ndim];
                var fminr = new double[ndim];
                // Initialize the left and right bounds for each dimension
                for (j = 0; j < ndim; j++)
                {              
                    iran = (iran * 2661 + 36979) % 175000;
                    s = Tools.Sign(dith, (double)(iran - 87500));
                    rmid[j] = (0.5 + s) * regn[j] + (0.5 - s) * regn[ndim + j];
                    fminl[j] = fminr[j] = BIG;
                    fmaxl[j] = fmaxr[j] = (-BIG);
                }
                // Loop over the points in the sample
                for (n = 0; n < npre; n++)
                {
                    ranpt(pt, regn);
                    fval = function(pt);
                    FunctionEvaluations++;
                    // Find the left and right bounds for each dimension
                    for (j = 0; j < ndim; j++)
                    {
                        if (pt[j] <= rmid[j])
                        {
                            fminl[j] = Math.Min(fminl[j], fval);
                            fmaxl[j] = Math.Max(fmaxl[j], fval);
                        }
                        else
                        {
                            fminr[j] = Math.Min(fminr[j], fval);
                            fmaxr[j] = Math.Max(fmaxr[j], fval);
                        }
                    }
                }
                // Choose which dimension jb to bisect
                sumb = BIG;
                jb = -1;
                siglb = sigrb = 1.0;
                for (j = 0; j < ndim; j++)
                {
                    if (fmaxl[j] > fminl[j] && fmaxr[j] > fminr[j])
                    {
                        sigl = Math.Max(TINY, Math.Pow(fmaxl[j] - fminl[j], 2.0 / 3.0));
                        sigr = Math.Max(TINY, Math.Pow(fmaxr[j] - fminr[j], 2.0 / 3.0));
                        sum = sigl + sigr; // Equation (7.9.24), see text.
                        if (sum <= sumb)
                        {
                            sumb = sum;
                            jb = j;
                            siglb = sigl;
                            sigrb = sigr;
                        }
                    }
                }
                if (jb == -1) jb = (ndim * iran) / 175000; // MNPT may be too small
                // Apportion the remain points between left and right
                rgl = regn[jb];
                rgm = rmid[jb];
                rgr = regn[ndim + jb];
                fracl = Math.Abs((rgm - rgl) / (rgr - rgl));
                nptl = (int)(MNPT + (npts - npre - 2 * MNPT) * fracl * siglb / (fracl * siglb + (1.0 - fracl) * sigrb)); // Equation (7.9.23)
                nptr = npts - npre - nptl;
                // Now allocate and integrate the two sub-regions
                var regn_temp = new double[2 * ndim];
                for (j = 0; j < ndim; j++)
                {
                    regn_temp[j] = regn[j];
                    regn_temp[ndim + j] = regn[ndim + j];
                }
                regn_temp[ndim + jb] = rmid[jb];
                miser(function, regn_temp, nptl, dith, ref avel, ref varl);
                regn_temp[jb] = rmid[jb];
                regn_temp[ndim + jb] = regn[ndim + jb];
                miser(function, regn_temp, nptr, dith, ref ave, ref var);
                ave = fracl * avel + (1 - fracl) * ave;
                var = fracl * fracl * varl + (1 - fracl) * (1 - fracl) * var;
                // Combine left and right regions by equation (7.9.11)
            }
        }

        /// <summary>
        /// Returns a uniformly random point pt in an n-dimensional rectangular region. Used by miser.
        /// </summary>
        /// <param name="pt"> An array the length of ndim </param>
        /// <param name="regn">  A vector consisting of ndim "lower-left" coordinates of the region followed by ndim "upper - right"
        /// coordinates. Specifies the rectangular volume by regn[0...2 * ndim -1] </param>
        private void ranpt(double[] pt, double[] regn)
        {
            int j, n = pt.Length;
            double[] rnd = null;
            if (UseSobolSequence)
                rnd = _sobol.NextDouble();

            for (j = 0; j < n; j++)
                pt[j] = regn[j] + (regn[n + j] - regn[j]) * (UseSobolSequence ? rnd[j] : Random.NextDouble());
        }
    }
}
