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

namespace Numerics.Distributions.Copulas
{
    /// <summary>
    /// Declares common functionality of all Bivariate Copulas.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class BivariateCopula : IBivariateCopula
    {

        #region Members

        /// <summary>
        /// Protected dependency property.
        /// </summary>
        protected double _theta;

        /// <summary>
        /// Protected parameter is valid property.
        /// </summary>
        protected bool _parametersValid = true;

        /// <inheritdoc/>
        public abstract CopulaType Type { get; }

        /// <inheritdoc/>
        public double Theta
        {
            get { return _theta; }
            set
            {
                _parametersValid = ValidateParameter(value, false) is null;
                _theta = value;
            }
        }

        /// <inheritdoc/>
        public abstract double ThetaMinimum { get; }

        /// <inheritdoc/>
        public abstract double ThetaMaximum { get; }

        /// <summary>
        /// Returns the name and value of the theta parameter in 2-column array of string.
        /// </summary>
        public abstract string[,] ParameterToString { get; }


        /// <summary>
        /// Returns the distribution parameter name in short form (e.g. θ).
        /// </summary>
        public abstract string ParameterNameShortForm { get; }

        /// <summary>
        /// Returns the distribution parameter property name.
        /// </summary>
        public string GetParameterPropertyName
        {
            get { return nameof(Theta); }
        }

        /// <summary>
        /// Returns a boolean value describing if the current parameters are valid or not.
        /// If not, an ArgumentOutOfRange exception will be thrown when trying to use statistical functions (e.g. CDF())
        /// </summary>
        public bool ParametersValid
        {
            get { return _parametersValid; }
        }

        /// <inheritdoc/>
        public virtual IUnivariateDistribution MarginalDistributionX { get; set; }

        /// <inheritdoc/>
        public virtual IUnivariateDistribution MarginalDistributionY { get; set; }

        /// <inheritdoc/>
        public abstract string DisplayName { get; }

        /// <inheritdoc/>
        public abstract string ShortDisplayName { get; }

        #endregion

        #region Methods

        /// <inheritdoc/>
        public abstract double PDF(double u, double v);

        /// <summary>
        /// Returns the natural log of the PDF.
        /// </summary>
        /// <param name="u">The reduced variate between 0 and 1.</param>
        /// <param name="v">The reduced variate between 0 and 1.</param>
        public double LogPDF(double u, double v)
        {
            double f = PDF(u, v);
            return Math.Log(f);
        }

        /// <inheritdoc/>
        public abstract double CDF(double u, double v);

        /// <inheritdoc/>
        public abstract double[] InverseCDF(double u, double v);

        /// <inheritdoc/>
        public abstract double[] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY);

        /// <inheritdoc/>
        public abstract ArgumentOutOfRangeException ValidateParameter(double parameter, bool throwException);

        /// <summary>
        /// Create a deep copy of the copula.
        /// </summary>
        public abstract BivariateCopula Clone();

        /// <summary>
        /// Returns the OR joint exceedance probability. When either of the variables exceeds a particular threshold value
        /// </summary>
        /// <param name="u">The reduced variate between 0 and 1.</param>
        /// <param name="v">The reduced variate between 0 and 1.</param>
        public double ORJointExceedanceProbability(double u, double v)
        {
            return 1 - CDF(u, v);
        }

        /// <summary>
        /// Returns the AND joint exceedance probability. When both variables exceed a particular threshold value simultaneously.
        /// </summary>
        /// <param name="u">The reduced variate between 0 and 1.</param>
        /// <param name="v">The reduced variate between 0 and 1.</param>
        public double ANDJointExceedanceProbability(double u, double v)
        {
            return 1 - u - v + CDF(u, v);
        }

        /// <inheritdoc/>
        public double[,] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            var rand = LatinHypercube.Random(sampleSize, 2, seed);
            var sample = new double[sampleSize, 2];
            for (int i = 0; i < sampleSize; i++)
            {
                var vals = InverseCDF(rand[i, 0], rand[i, 1]);
                if (MarginalDistributionX != null && MarginalDistributionY != null)
                {
                    sample[i, 0] = MarginalDistributionX.InverseCDF(vals[0]);
                    sample[i, 1] = MarginalDistributionY.InverseCDF(vals[1]);
                }
                else
                {
                    sample[i, 0] = vals[0];
                    sample[i, 1] = vals[1];
                }
            }
            return sample;
        }

        /// <summary>
        /// The pseudo log-likelihood function.
        /// </summary>
        /// <param name="sampleDataX">The sample data for the X variable. When estimating with pseudo likelihood, this should be the plotting positions of the data.</param>
        /// <param name="sampleDataY">The sample data for the Y variable. When estimating with pseudo likelihood, this should be the plotting positions of the data.</param>
        public double PseudoLogLikelihood(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            if (sampleDataX.Count < 2 || sampleDataY.Count < 2) throw new ArgumentOutOfRangeException("sample data", "There must be at least two items in each sample.");
            if (sampleDataX.Count != sampleDataY.Count) throw new ArgumentOutOfRangeException("sample data", "The sample data arrays must be the same length.");
           
            double LogLH = 0;
            for (int i = 0; i < sampleDataX.Count; i++)
                LogLH += LogPDF(sampleDataX[i], sampleDataY[i]);
            if (double.IsNaN(LogLH) || double.IsInfinity(LogLH)) return double.MinValue;
            return LogLH;
        }

        /// <summary>
        /// The inference from margins (IFM) log-likelihood function. The marginal distribution are assumed to have already been estimated independently. 
        /// </summary>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        public double IFMLogLikelihood(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            if (sampleDataX.Count < 2 || sampleDataY.Count < 2) throw new ArgumentOutOfRangeException("sample data", "There must be at least two items in each sample.");
            if (sampleDataX.Count != sampleDataY.Count) throw new ArgumentOutOfRangeException("sample data", "The sample data arrays must be the same length.");
            if (MarginalDistributionX == null || MarginalDistributionY == null) throw new ArgumentOutOfRangeException("marginal distributions", "There must be 2 marginal distributions to evaluate.");

            double LogLH = 0;
            for (int i = 0; i < sampleDataX.Count; i++)
                LogLH += LogPDF(MarginalDistributionX.CDF(sampleDataX[i]), MarginalDistributionY.CDF(sampleDataY[i]));
            if (double.IsNaN(LogLH) || double.IsInfinity(LogLH)) return double.MinValue;
            return LogLH;
        }

        /// <summary>
        /// The full log-likelihood function. The marginal distributions are estimated simultaneously with the copula. 
        /// </summary>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        public double LogLikelihood(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            if (sampleDataX.Count < 2 || sampleDataY.Count < 2) throw new ArgumentOutOfRangeException("sample data", "There must be at least two items in each sample.");
            if (sampleDataX.Count != sampleDataY.Count) throw new ArgumentOutOfRangeException("sample data", "The sample data arrays must be the same length.");
            if (MarginalDistributionX == null || MarginalDistributionY == null) throw new ArgumentOutOfRangeException("marginal distributions", "There must be 2 marginal distributions to evaluate.");

            double LogLH = 0;
            for (int i = 0; i < sampleDataX.Count; i++)
                LogLH += LogPDF(MarginalDistributionX.CDF(sampleDataX[i]), MarginalDistributionY.CDF(sampleDataY[i])) + MarginalDistributionX.LogPDF(sampleDataX[i]) + MarginalDistributionY.LogPDF(sampleDataY[i]);
            if (double.IsNaN(LogLH) || double.IsInfinity(LogLH)) return double.MinValue;
            return LogLH;
        }


        #endregion
    }
}
