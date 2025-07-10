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
using Numerics.Data;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Bivariate Empirical distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class BivariateEmpirical : MultivariateDistribution
    {

        /// <summary>
        /// Constructs an empty Bivariate Empirical CDF.
        /// </summary>
        /// <param name="x1Transform">Optional. Linear, logarithmic, or normal z. Default = None. </param>
        /// <param name="x2Transform">Optional. Linear, logarithmic, or normal z. Default = None. </param>
        /// <param name="probabilityTransform">Optional. Linear, logarithmic, or normal z. Default = None. </param>
        public BivariateEmpirical(Transform x1Transform = Transform.None, Transform x2Transform = Transform.None, Transform probabilityTransform = Transform.None)
        {
            X1Transform = x1Transform;
            X2Transform = x2Transform;
            ProbabilityTransform = probabilityTransform;
        }

        /// <summary>
        /// Constructs a Bivariate Empirical CDF with specified parameters.
        /// </summary>
        /// <param name="x1Values">Array of X1 values. The X1-values represent the primary values. There are rows in the table of probability values.</param>
        /// <param name="x2Values">Array of X2 values. The X2-values represent the secondary values. There are columns in the table of probability values.</param>
        /// <param name="pValues">Array of probability values. Range 0 ≤ p ≤ 1.</param>
        /// <param name="x1Transform">Optional. Linear, logarithmic, or normal z. Default = None. </param>
        /// <param name="x2Transform">Optional. Linear, logarithmic, or normal z. Default = None. </param>
        /// <param name="probabilityTransform">Optional. Linear, logarithmic, or normal z. Default = None. </param>
        public BivariateEmpirical(IList<double> x1Values, IList<double> x2Values, double[,] pValues, Transform x1Transform = Transform.None, Transform x2Transform = Transform.None, Transform probabilityTransform = Transform.None)
        {
            X1Transform = x1Transform;
            X2Transform = x2Transform;
            ProbabilityTransform = probabilityTransform;
            SetParameters(x1Values, x2Values, pValues);
        }

        //    X2 ... X2n
        // X1 P(1,1)   P(1,n)
        // ...
        // X1n P(n,1)   P(n,n)        

        private Bilinear bilinear = null;
        private bool _parametersValid = true;

        /// <summary>
        /// Return the array of X1 values (distribution 1). Points On the cumulative curve are specified
        /// with increasing value and increasing probability.
        /// </summary>
        public double[] X1Values { get; private set; }

        /// <summary>
        /// Return the array of X2 values (distribution 2). Points on the cumulative curve are specified
        /// with increasing value and increasing probability.
        /// </summary>
        public double[] X2Values { get; private set; }

        /// <summary>
        /// Return the array of probability values. Points on the cumulative curve are specified
        /// with increasing value and increasing probability.
        /// </summary>
        public double[,] ProbabilityValues { get; private set; }

        /// <summary>
        /// Determines the interpolation transform for the X1-values.
        /// </summary>
        public Transform X1Transform { get; set; }

        /// <summary>
        /// Determines the interpolation transform for the X2-values.
        /// </summary>
        public Transform X2Transform { get; set; }

        /// <summary>
        /// Determines the interpolation transform for the Probability-values.
        /// </summary>
        public Transform ProbabilityTransform { get; set; }

        /// <summary>
        /// Gets the number of variables for the distribution.
        /// </summary>
        public override int Dimension
        {
            get { return 2; }
        }

        /// <summary>
        /// Returns the multivariate distribution type.
        /// </summary>
        public override MultivariateDistributionType Type
        {
            get { return MultivariateDistributionType.BivariateEmpiricalDistribution; }
        }

        /// <summary>
        /// Returns the name of the distribution type as a string.
        /// </summary>
        public override string DisplayName
        {
            get { return "Bivariate Empirical"; }
        }

        /// <summary>
        /// Returns the short display name of the distribution as a string.
        /// </summary>
        public override string ShortDisplayName
        {
            get { return "Bi. Emp"; }
        }

        /// <summary>
        /// Determines whether the parameters are valid or not.
        /// </summary>
        public override bool ParametersValid
        {
            get { return _parametersValid; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="x1Values">Array of X1 values. The X1-values represent the primary values. There are rows in the table of probability values.</param>
        /// <param name="x2Values">Array of X2 values. The X2-values represent the secondary values. There are columns in the table of probability values.</param>
        /// <param name="pValues">Array of probability values. Range 0 ≤ p ≤ 1.</param>
        public void SetParameters(IList<double> x1Values, IList<double> x2Values, double[,] pValues)
        {
            // validate parameters
            _parametersValid = ValidateParameters(x1Values, x2Values, pValues, false) is null;
            X1Values = x1Values.ToArray();
            X2Values = x2Values.ToArray();
            ProbabilityValues = pValues;
            
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="x1Values">Array of X1 values. The X1-values represent the primary values. There are rows in the table of probability values.</param>
        /// <param name="x2Values">Array of X2 values. The X2-values represent the secondary values. There are columns in the table of probability values.</param>
        /// <param name="pValues">Array of probability values. Range 0 ≤ p ≤ 1.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(IList<double> x1Values, IList<double> x2Values, double[,] pValues, bool throwException)
        {

            if (x1Values.Count < 2)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException("Primary Values", "There must be at least 2 primary values.");
                return new ArgumentOutOfRangeException("Primary Values", "There must be at least 2 primary values.");
            }
            if (x2Values.Count < 2)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException("Secondary Values", "There must be at least 2 secondary values.");
                return new ArgumentOutOfRangeException("Secondary Values", "There must be at least 2 secondary values.");
            }

            // Check if X1 and X2 are in ascending order
            for (int i = 1; i < x1Values.Count; i++)
            {
                if (x1Values[i] <= x1Values[i - 1])
                {
                    if (throwException)
                        throw new ArgumentOutOfRangeException("Primary Values", "Primary values must be in ascending order.");
                    return new ArgumentOutOfRangeException("Primary Values", "Primary values must be in ascending order.");
                }
            }
            for (int i = 1; i < x2Values.Count; i++)
            {
                if (x2Values[i] <= x2Values[i - 1])
                {
                    if (throwException)
                        throw new ArgumentOutOfRangeException("Secondary Values", "Secondary values must be in ascending order.");
                    return new ArgumentOutOfRangeException("Secondary Values", "Secondary values must be in ascending order.");
                }
            }

            // Check if probabilities are between 0 and 1 and in ascending order
            if (pValues.GetLength(0) != x1Values.Count)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(ProbabilityValues), "The number of rows in the probability array must be the same length as the X1 array.");
                return new ArgumentOutOfRangeException(nameof(ProbabilityValues), "The number of rows in the probability array must be the same length as the X1 array.");
            }
            if (pValues.GetLength(1) != x2Values.Count)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(ProbabilityValues), "The number of columns in the probability array must be the same length as the X2 array.");
                return new ArgumentOutOfRangeException(nameof(ProbabilityValues), "The number of columns in the probability array must be the same length as the X2 array.");
            }

            for (int i = 0; i < pValues.GetLength(0); i++) // Rows
            {
                for (int j = 0; j < pValues.GetLength(1); j++) // Columns
                {
                    if (pValues[i, j] < 0.0d || pValues[i, j] > 1.0d)
                    {
                        if (throwException) throw new ArgumentOutOfRangeException(nameof(ProbabilityValues), "Probability values must be equal to or between 0 and 1.");
                        return new ArgumentOutOfRangeException(nameof(ProbabilityValues), "Probability values must be equal to or between 0 and 1.");
                    }
                }
            }

            return null;
        }

        /// <summary>
        /// The Probability Density Function (PDF) of the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        public override double PDF(double[] x)
        {
            // Validate parameters
            //if (_parametersValid == false) ValidateParameters(X1Values, X2Values, ProbabilityValues, true);
            return PDF(x[0], x[1]);
        }

        /// <summary>
        /// The Probability Density Function (PDF) of the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x1">The x1-value.</param>
        /// <param name="x2">The x2-value.</param>
        public double PDF(double x1, double x2)
        {
            // The PDF is estimated using numerical differentiation of the CDF.
            // This approach is not ideal, and is a temporary place holder, 
            // until I learn how to do this better. 
            // double h = 0.0001d;
            // return (CDF(x1 + h, x2 + h) - CDF(x1 - h, x2 - h)) / (2d * h);
            return double.NaN;
        }

        /// <summary>
        /// The Cumulative Distribution Function (CDF) for the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        public override double CDF(double[] x)
        {
            return CDF(x[0], x[1]);
        }

        /// <summary>
        /// The Cumulative Distribution Function (CDF) for the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x1">The x1-value.</param>
        /// <param name="x2">The x2-value.</param>
        public double CDF(double x1, double x2)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters(X1Values, X2Values, ProbabilityValues, true);
            // Make sure the transforms are up-to-date.
            if (bilinear == null) bilinear = new Bilinear(X1Values, X2Values, ProbabilityValues);
            bilinear.X1Transform = X1Transform;
            bilinear.X2Transform = X2Transform;
            bilinear.YTransform = ProbabilityTransform;
            return bilinear.Interpolate(x1, x2);
        }

        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public override MultivariateDistribution Clone()
        {
            return new BivariateEmpirical(X1Values, X2Values, ProbabilityValues) { X1Transform = X1Transform, X2Transform = X2Transform, ProbabilityTransform = ProbabilityTransform };
        }

    }
}