using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;

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
        public virtual IUnivariateDistribution? MarginalDistributionX { get; set; }

        /// <inheritdoc/>
        public virtual IUnivariateDistribution? MarginalDistributionY { get; set; }

        /// <inheritdoc/>
        public abstract string DisplayName { get; }

        /// <inheritdoc/>
        public abstract string ShortDisplayName { get; }

        /// <inheritdoc/>
        public abstract double UpperTailDependence { get; }

        /// <inheritdoc/>
        public abstract double LowerTailDependence { get; }

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
            if (double.IsNaN(f) || double.IsInfinity(f) || f <= 0d)
                return double.NegativeInfinity;
            return Math.Log(f);
        }

        /// <inheritdoc/>
        public abstract double CDF(double u, double v);

        /// <summary>
        /// Returns the forward conditional cumulative distribution function (h-function),
        /// h(v|u) = ∂C(u,v)/∂u = P(V ≤ v | U = u).
        /// </summary>
        /// <param name="u">The conditioning variate, a non-exceedance probability between 0 and 1.</param>
        /// <param name="v">The dependent variate, a non-exceedance probability between 0 and 1.</param>
        /// <returns>The conditional non-exceedance probability P(V ≤ v | U = u).</returns>
        /// <remarks>
        /// Both arguments follow the copula convention of non-exceedance probabilities.
        /// The base implementation approximates the partial derivative with a central finite
        /// difference over <see cref="CDF"/> using ε = 1E-6, clamping u ± ε to [0, 1] and dividing
        /// by the realized step width. Every copula shipped with the library overrides this method
        /// with its exact analytic conditional; the finite-difference fallback exists so an external
        /// subclass that only implements <see cref="CDF"/> still gets a correct approximation.
        /// <see cref="InverseConditionalCDF"/> inverts this function in v.
        /// </remarks>
        public virtual double ConditionalCDF(double u, double v)
        {
            double uPlus = Math.Min(u + 1E-6, 1d);
            double uMinus = Math.Max(u - 1E-6, 0d);
            return (CDF(uPlus, v) - CDF(uMinus, v)) / (uPlus - uMinus);
        }

        /// <summary>
        /// Returns the inverse of the forward conditional CDF with respect to v: the value v
        /// such that h(v|u) = t, where h is <see cref="ConditionalCDF"/>.
        /// </summary>
        /// <param name="u">The conditioning variate, a non-exceedance probability between 0 and 1.</param>
        /// <param name="t">The conditional non-exceedance probability between 0 and 1.</param>
        /// <returns>The dependent variate v such that P(V ≤ v | U = u) = t.</returns>
        /// <remarks>
        /// This is the scalar, non-allocating form of the conditional simulation that
        /// <see cref="InverseCDF(double, double)"/> returns as its second element. The base
        /// implementation delegates to that array form; every copula shipped with the library
        /// overrides this method with the scalar computation and recomposes the array form on
        /// top of it, so hot loops can invert conditional probabilities without per-call
        /// allocation. Subclass authors must not implement <see cref="InverseCDF(double, double)"/>
        /// by delegating its second coordinate back to this inherited method, because the two
        /// base calls would recurse. Override this method when composing the array form from a
        /// scalar conditional inverse.
        /// </remarks>
        public virtual double InverseConditionalCDF(double u, double t)
        {
            return InverseCDF(u, t)[1];
        }

        /// <inheritdoc/>
        public abstract double[] InverseCDF(double u, double v);

        /// <inheritdoc/>
        public abstract int NumberOfCopulaParameters { get; }

        /// <inheritdoc/>
        public abstract double[] GetCopulaParameters { get; }

        /// <inheritdoc/>
        public abstract void SetCopulaParameters(double[] parameters);

        /// <inheritdoc/>
        public abstract double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY);

        /// <inheritdoc/>
        public abstract ArgumentOutOfRangeException? ValidateParameter(double parameter, bool throwException);

        /// <summary>
        /// Create a deep copy of the copula, including independently cloned marginal
        /// distributions.
        /// </summary>
        /// <remarks>
        /// Marginals are deep-copied through <see cref="CloneMarginal"/> because
        /// distributions memoize internal state lazily — clones sharing marginal
        /// references are not safe to use concurrently (e.g., one clone per thread in
        /// Monte Carlo frameworks or parallel likelihood evaluation).
        /// </remarks>
        public abstract BivariateCopula Clone();

        /// <summary>
        /// Deep-copies an attached marginal distribution for <see cref="Clone"/>.
        /// </summary>
        /// <param name="marginal">The marginal distribution, or null when unattached.</param>
        /// <returns>An independent copy of the marginal; null when unattached. A marginal
        /// implementation outside <see cref="UnivariateDistributionBase"/> (which carries the
        /// clone support) is returned by reference, preserving the previous behavior.</returns>
        protected static IUnivariateDistribution? CloneMarginal(IUnivariateDistribution? marginal)
        {
            if (marginal is UnivariateDistributionBase cloneable)
            {
                return cloneable.Clone();
            }
            return marginal;
        }

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
            if (double.IsNaN(LogLH) || double.IsInfinity(LogLH)) return double.NegativeInfinity;
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
            if (double.IsNaN(LogLH) || double.IsInfinity(LogLH)) return double.NegativeInfinity;
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
            if (double.IsNaN(LogLH) || double.IsInfinity(LogLH)) return double.NegativeInfinity;
            return LogLH;
        }

        /// <summary>
        /// Returns an XElement of the copula's type and parameters, which can be used for
        /// serialization.
        /// </summary>
        /// <returns>
        /// An XElement named "Copula" carrying the <see cref="Type"/> by enumeration name and
        /// the copula parameters in a "Parameters" attribute as a pipe-delimited,
        /// invariant-culture, round-trip ("G17") string in <see cref="SetCopulaParameters"/>
        /// order. A zero-parameter copula writes an empty "Parameters" attribute.
        /// </returns>
        /// <remarks>
        /// Marginal distributions are deliberately not serialized: consumers own their
        /// marginals and attach them separately, so the element captures the dependence
        /// structure alone. Serializing the type by name (never by numeric value) keeps the
        /// payload valid as long as <see cref="CopulaType"/> members are only ever appended.
        /// <see cref="CopulaFactory.CreateCopula(XElement)"/> reconstructs the copula.
        /// </remarks>
        public virtual XElement ToXElement()
        {
            var result = new XElement("Copula");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            var parameters = GetCopulaParameters;
            var values = new string[parameters.Length];
            for (int i = 0; i < parameters.Length; i++)
                values[i] = parameters[i].ToString("G17", CultureInfo.InvariantCulture);
            result.SetAttributeValue("Parameters", string.Join("|", values));
            return result;
        }

        #endregion
    }
}
