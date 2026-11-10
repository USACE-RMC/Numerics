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
        /// <returns>The log density; negative infinity where the density is zero, infinite or NaN.</returns>
        /// <remarks>
        /// The base implementation takes the log of <see cref="PDF(double, double)"/>. Families whose density
        /// underflows or overflows in direct form override it with a log-space evaluation, which the
        /// likelihood methods use.
        /// </remarks>
        public virtual double LogPDF(double u, double v)
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
        /// difference over <see cref="CDF(double, double)"/> using ε = 1E-6, clamping u ± ε to [0, 1] and dividing
        /// by the realized step width. Every copula shipped with the library overrides this method
        /// with its exact analytic conditional; the finite-difference fallback exists so an external
        /// subclass that only implements <see cref="CDF(double, double)"/> still gets a correct approximation.
        /// <see cref="InverseConditionalCDF(double, double)"/> inverts this function in v.
        /// </remarks>
        public virtual double ConditionalCDF(double u, double v)
        {
            double uPlus = Math.Min(u + 1E-6, 1d);
            double uMinus = Math.Max(u - 1E-6, 0d);
            return (CDF(uPlus, v) - CDF(uMinus, v)) / (uPlus - uMinus);
        }

        /// <summary>
        /// Returns the inverse of the forward conditional CDF with respect to v: the value v
        /// such that h(v|u) = t, where h is <see cref="ConditionalCDF(double, double)"/>.
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

        /// <summary>
        /// The PDF evaluated over paired lists of reduced variates.
        /// </summary>
        /// <param name="u">The first reduced variates, each between 0 and 1.</param>
        /// <param name="v">The second reduced variates, each between 0 and 1, one per entry of <paramref name="u"/>.</param>
        /// <returns>One density per pair, each equal to <see cref="PDF(double, double)"/> at that pair.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="u"/> or <paramref name="v"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists differ in length.</exception>
        /// <remarks>
        /// One call evaluates every pair, in the same way as the list forms of the univariate distributions
        /// (<see cref="UnivariateDistributionBase.PDF(IList{double})"/>); a caller crossing a language boundary (for
        /// example from Python) pays the crossing once instead of once per pair.
        /// </remarks>
        public double[] PDF(IList<double> u, IList<double> v)
        {
            ValidatePairs(u, v);
            var result = new double[u.Count];
            for (int i = 0; i < u.Count; i++)
                result[i] = PDF(u[i], v[i]);
            return result;
        }

        /// <summary>
        /// The CDF evaluated over paired lists of reduced variates.
        /// </summary>
        /// <param name="u">The first reduced variates, each between 0 and 1.</param>
        /// <param name="v">The second reduced variates, each between 0 and 1, one per entry of <paramref name="u"/>.</param>
        /// <returns>One joint non-exceedance probability per pair, each equal to <see cref="CDF(double, double)"/> at that pair.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="u"/> or <paramref name="v"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists differ in length.</exception>
        /// <remarks>
        /// One call evaluates every pair, in the same way as <see cref="PDF(IList{double}, IList{double})"/>.
        /// </remarks>
        public double[] CDF(IList<double> u, IList<double> v)
        {
            ValidatePairs(u, v);
            var result = new double[u.Count];
            for (int i = 0; i < u.Count; i++)
                result[i] = CDF(u[i], v[i]);
            return result;
        }

        /// <summary>
        /// The forward conditional CDF (h-function) evaluated over paired lists of variates.
        /// </summary>
        /// <param name="u">The conditioning variates, each a non-exceedance probability between 0 and 1.</param>
        /// <param name="v">The dependent variates, each a non-exceedance probability between 0 and 1, one per entry of <paramref name="u"/>.</param>
        /// <returns>One conditional non-exceedance probability per pair, each equal to <see cref="ConditionalCDF(double, double)"/> at that pair.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="u"/> or <paramref name="v"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists differ in length.</exception>
        /// <remarks>
        /// One call evaluates every pair, in the same way as <see cref="PDF(IList{double}, IList{double})"/>.
        /// </remarks>
        public double[] ConditionalCDF(IList<double> u, IList<double> v)
        {
            ValidatePairs(u, v);
            var result = new double[u.Count];
            for (int i = 0; i < u.Count; i++)
                result[i] = ConditionalCDF(u[i], v[i]);
            return result;
        }

        /// <summary>
        /// The inverse of the forward conditional CDF evaluated over paired lists of variates.
        /// </summary>
        /// <param name="u">The conditioning variates, each a non-exceedance probability between 0 and 1.</param>
        /// <param name="t">The conditional non-exceedance probabilities, each between 0 and 1, one per entry of <paramref name="u"/>.</param>
        /// <returns>One dependent variate per pair, each equal to <see cref="InverseConditionalCDF(double, double)"/> at that pair.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="u"/> or <paramref name="t"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists differ in length.</exception>
        /// <remarks>
        /// One call evaluates every pair, in the same way as <see cref="PDF(IList{double}, IList{double})"/>: the
        /// conditional simulation of many events in one call.
        /// </remarks>
        public double[] InverseConditionalCDF(IList<double> u, IList<double> t)
        {
            ValidatePairs(u, t);
            var result = new double[u.Count];
            for (int i = 0; i < u.Count; i++)
                result[i] = InverseConditionalCDF(u[i], t[i]);
            return result;
        }

        /// <summary>
        /// Checks that two lists of variates pair up.
        /// </summary>
        /// <param name="u">The first list.</param>
        /// <param name="v">The second list.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="u"/> or <paramref name="v"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists differ in length.</exception>
        private static void ValidatePairs(IList<double> u, IList<double> v)
        {
            if (u == null) throw new ArgumentNullException(nameof(u));
            if (v == null) throw new ArgumentNullException(nameof(v));
            if (u.Count != v.Count) throw new ArgumentException("The lists of variates must have the same length.", nameof(v));
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
        /// distributions memoize internal state lazily; clones sharing marginal
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
        /// clone support) is returned by reference.</returns>
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
            return TransformDesign(LatinHypercube.Random(sampleSize, 2, seed));
        }

        /// <summary>
        /// Generates a random sample from the copula, drawing the Latin hypercube design from a supplied generator.
        /// </summary>
        /// <param name="sampleSize">The number of samples to generate.</param>
        /// <param name="random">The master generator passed to <see cref="LatinHypercube.Random(int, int, Random)"/>.</param>
        /// <returns>A <paramref name="sampleSize"/> by 2 array of samples, in the units of the marginal distributions
        /// when both are attached and as non-exceedance probabilities otherwise.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="random"/> is null.</exception>
        /// <remarks>
        /// The sample is a deterministic function of the generator state, with no clock fallback. Passing
        /// <c>new MersenneTwister(seed)</c> with a positive seed reproduces
        /// <see cref="GenerateRandomValues(int, int)"/> with that seed exactly.
        /// </remarks>
        public double[,] GenerateRandomValues(int sampleSize, Random random)
        {
            return TransformDesign(LatinHypercube.Random(sampleSize, 2, random));
        }

        /// <summary>
        /// Maps a two-column uniform design to copula samples by conditional inversion, then through the
        /// marginal distributions when both are attached.
        /// </summary>
        /// <param name="rand">The uniform design, one row per sample.</param>
        /// <returns>The sample, one row per design row.</returns>
        private double[,] TransformDesign(double[,] rand)
        {
            int sampleSize = rand.GetLength(0);
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
        /// Marginal distributions are not serialized by this overload: the element captures the
        /// dependence structure alone, for consumers that own their marginals and attach them
        /// separately. <see cref="ToXElement(bool)"/> includes them. Serializing the type by name
        /// (never by numeric value) keeps the payload valid as long as <see cref="CopulaType"/>
        /// members are only ever appended. <see cref="CopulaFactory.CreateCopula(XElement)"/>
        /// reconstructs the copula.
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

        /// <summary>
        /// Returns an XElement of the copula's type and parameters, optionally with its attached
        /// marginal distributions.
        /// </summary>
        /// <param name="includeMarginals">True to append each attached marginal distribution.</param>
        /// <returns>
        /// The element written by <see cref="ToXElement()"/>. When <paramref name="includeMarginals"/>
        /// is true, each attached marginal is appended as a child element named
        /// "MarginalDistributionX" or "MarginalDistributionY" that holds the element written by
        /// <see cref="UnivariateDistributionBase.ToXElement"/>; an unattached marginal is omitted.
        /// </returns>
        /// <exception cref="ArgumentException">Thrown when <paramref name="includeMarginals"/> is true and an
        /// attached marginal does not derive from <see cref="UnivariateDistributionBase"/>, which carries the
        /// serialization support.</exception>
        /// <remarks>
        /// <see cref="CopulaFactory.CreateCopula(XElement, bool)"/> reconstructs the copula and its marginals.
        /// </remarks>
        public XElement ToXElement(bool includeMarginals)
        {
            var result = ToXElement();
            if (!includeMarginals) return result;
            AppendMarginal(result, nameof(MarginalDistributionX), MarginalDistributionX);
            AppendMarginal(result, nameof(MarginalDistributionY), MarginalDistributionY);
            return result;
        }

        /// <summary>
        /// Appends an attached marginal distribution to a serialized copula.
        /// </summary>
        /// <param name="element">The copula element.</param>
        /// <param name="name">The name of the child element.</param>
        /// <param name="marginal">The marginal distribution, or null when unattached.</param>
        /// <exception cref="ArgumentException">Thrown when the marginal does not derive from
        /// <see cref="UnivariateDistributionBase"/>.</exception>
        private static void AppendMarginal(XElement element, string name, IUnivariateDistribution? marginal)
        {
            if (marginal == null) return;
            if (marginal is not UnivariateDistributionBase serializable)
                throw new ArgumentException("The marginal distribution " + name + " does not derive from UnivariateDistributionBase and cannot be serialized.", "includeMarginals");
            element.Add(new XElement(name, serializable.ToXElement()));
        }

        #endregion
    }
}
