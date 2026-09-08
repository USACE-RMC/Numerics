using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Kappa-4 distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// <list type="bullet">
    /// <item><description>
    /// If h = -1 , then the Kappa-4 is the Generalized Logistic distribution.
    /// </description></item>
    /// <item><description>
    /// If h = 0 , then the Kappa-4 is the Generalized Extreme Value distribution.
    /// </description></item>
    /// <item><description>
    /// If h = 1 , then the Kappa-4 is the Generalized Pareto distribution.
    /// </description></item>
    /// </list>
    /// </para>
    /// <para>
    /// Finite shapes and a positive finite scale define the distribution. Moment existence and
    /// estimation restrictions are separate: an absolute moment of order r requires rκ &gt; -1,
    /// and additionally rκh &gt; -1 when h &lt; 0. Moment properties use standardized adaptive
    /// quantile integration with relative tolerance 1e-10 and absolute tolerance 1e-12, returning
    /// NaN when the required moment does not exist or quadrature fails.
    /// </para>
    /// <para>
    /// The unrestricted likelihood can be unbounded at a sample-dependent support endpoint.
    /// A successful numerical fit does not establish the existence of a finite global maximum.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> <see href = "https://ieeexplore.ieee.org/document/5389569" /> </item>
    /// <item> <see href = "https://rdrr.io/cran/nsRFA/src/R/KAPPA.R" /> </item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class KappaFour : UnivariateDistributionBase, IStandardError, IEstimation, IMaximumLikelihoodEstimation, ILinearMomentEstimation, IBootstrappable
    {

        /// <summary>
        /// Constructs a Kappa-4 distribution.
        /// </summary>
        public KappaFour()
        {
            SetParameters([100d, 10d, 0d, 0d]);
        }

        /// <summary>
        /// Constructs a Kappa-4 distribution with the given parameters ξ, α, κ, and h.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        /// <param name="shape2">The shape parameter h (hondo).</param>
        public KappaFour(double location, double scale, double shape, double shape2)
        {
            SetParameters([location, scale, shape, shape2]);
        }
    
        private double _xi; // location
        private double _alpha; // scale
        private double _kappa; // shape
        private double _hondo; // shape 2
        private bool _momentsComputed = false;
        private double[] u = [double.NaN, double.NaN, double.NaN, double.NaN];
        [NonSerialized] private bool _compensatedMinimumComputed;
        [NonSerialized] private double _compensatedMinimum;

        /// <summary>
        /// Gets and sets the location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return _xi; }
            set
            {
                _parametersValid = ValidateParameters([value, Alpha, Kappa, Hondo], false) is null;
                _xi = value;
                _momentsComputed = false;
                _compensatedMinimumComputed = false;
            }
        }

        /// <summary>
        /// Gets and sets the scale parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters([Xi, value, Kappa, Hondo], false) is null;
                _alpha = value;
                _momentsComputed = false;
                _compensatedMinimumComputed = false;
            }
        }

        /// <summary>
        /// Gets and sets the shape parameter κ (kappa).
        /// </summary>
        public double Kappa
        {
            get { return _kappa; }
            set
            {
                _parametersValid = ValidateParameters([Xi, Alpha, value, Hondo], false) is null;
                _kappa = value;
                _momentsComputed = false;
                _compensatedMinimumComputed = false;
            }
        }

        /// <summary>
        /// Gets and sets the shape parameter h (hondo).
        /// </summary>
        public double Hondo
        {
            get { return _hondo; }
            set
            {
                _parametersValid = ValidateParameters([Xi, Alpha, Kappa, value], false) is null;
                _hondo = value;
                _momentsComputed = false;
                _compensatedMinimumComputed = false;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 4; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.KappaFour; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Kappa-4"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "K4"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[4, 2];
                parmString[0, 0] = "Location (ξ)";
                parmString[1, 0] = "Scale (α)";
                parmString[2, 0] = "Shape (κ)";
                parmString[3, 0] = "Shape (h)";
                parmString[0, 1] = Xi.ToString();
                parmString[1, 1] = Alpha.ToString();
                parmString[2, 1] = Kappa.ToString();
                parmString[3, 1] = Hondo.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ξ", "α", "κ", "h"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Xi), nameof(Alpha), nameof(Kappa), nameof(Hondo)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Xi, Alpha, Kappa, Hondo]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = ComputeMoments();
                    _momentsComputed = true;
                }
                return u[0];
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        /// <remarks>Returns the unique density mode, including one-sided endpoint maxima, or NaN when no unique mode exists.</remarks>
        public override double Mode
        {
            get
            {
                EnsureValidParameters();
                double lower = Minimum, upper = Maximum;
                double lowerDensity = Tools.IsFinite(lower) ? LowerEndpointDensity() : 0d;
                double upperDensity = Tools.IsFinite(upper) ? PDF(upper) : 0d;
                if (double.IsPositiveInfinity(lowerDensity) && double.IsPositiveInfinity(upperDensity)) return double.NaN;
                if (double.IsPositiveInfinity(lowerDensity)) return lower;
                if (double.IsPositiveInfinity(upperDensity)) return upper;
                double denominator = 1d - Kappa * Hondo;
                double t = (1d - Kappa) / denominator;
                if (denominator > 0d && t > 0d && (Hondo <= 0d || Hondo * t < 1d))
                {
                    double logT = Tools.Log1p(-Kappa) - Tools.Log1p(-Kappa * Hondo);
                    return AffineQuantile(Xi, Alpha, QuantileFromLogT(logT, Kappa));
                }
                if (lowerDensity > upperDensity) return lower;
                if (upperDensity > lowerDensity) return upper;
                return double.NaN;
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = ComputeMoments();
                    _momentsComputed = true;
                }
                return u[1];
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = ComputeMoments();
                    _momentsComputed = true;
                }
                return u[2];
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = ComputeMoments();
                    _momentsComputed = true;
                }
                return u[3];
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get
            {
                if (Hondo > 0)
                {
                    if (!_compensatedMinimumComputed)
                        _compensatedMinimumComputed = KappaFourBoundary.TryLowerEndpoint(Xi, Alpha, Kappa, Hondo, out _compensatedMinimum);
                    if (_compensatedMinimumComputed) return _compensatedMinimum;
                }
                if (Hondo <= 0d && Kappa < 0d)
                {
                    return LocationPlusScaleOverShape();
                }
                else if (Hondo > 0d && Kappa != 0d)
                {
                    double logH = Math.Log(Hondo);
                    return AffineQuantile(Xi, Alpha, QuantileFromLogT(-logH, Kappa));
                }
                else if (Hondo > 0d && Kappa == 0d)
                {
                    return AffineQuantile(Xi, Alpha, Math.Log(Hondo));
                }
                else if (Hondo <= 0d && Kappa >= 0d)
                {
                    return double.NegativeInfinity;
                }
                return double.NaN;
            }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get
            {
                if (Kappa <= 0d)
                {
                    return double.PositiveInfinity;
                }
                else
                {
                    return LocationPlusScaleOverShape();
                }
            }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, 0.0d, double.NegativeInfinity, double.NegativeInfinity]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        /// <remarks>Parameters are installed only after estimation returns successfully.</remarks>
        /// <exception cref="InvalidOperationException">The selected estimator cannot produce a successful valid fit.</exception>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            ValidateFittingSample(sample);
            if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
            {
                SetParameters(ParametersFromLinearMoments(Statistics.LinearMoments(sample)));
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                SetParameters(MLE(sample));
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <inheritdoc/>
        /// <remarks>Uses the existing seeded random generator and propagates estimation failure without returning a fitted distribution.</remarks>
        /// <exception cref="InvalidOperationException">The resampled data cannot be fitted successfully.</exception>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            if (sampleSize < 4) throw new ArgumentOutOfRangeException(nameof(sampleSize), "At least four observations are required for Kappa Four estimation.");
            var newDistribution = new KappaFour(Xi, Alpha, Kappa, Hondo);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            // Set parameters
            Xi = parameters[0];
            Alpha = parameters[1];
            Kappa = parameters[2];
            Hondo = parameters[3];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
            }
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
            }
            if (double.IsNaN(parameters[2]) || double.IsInfinity(parameters[2]))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be a number.");
            }
            if (double.IsNaN(parameters[3]) || double.IsInfinity(parameters[3]))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Hondo), "The shape parameter h (hondo) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Hondo), "The shape parameter h (hondo) must be a number.");
            }
            return null!;
        }

        /// <summary>
        /// Estimates ξ, α, κ, and h from the first two L-moments and the L-skewness and L-kurtosis.
        /// </summary>
        /// <param name="moments">The values L1, L2, τ3, and τ4, in that order.</param>
        /// <returns>The fitted location, scale, kappa, and hondo parameters.</returns>
        /// <exception cref="ArgumentNullException">The L-moments are null.</exception>
        /// <exception cref="ArgumentException">There are not four L-moments, or their ratios are outside the estimation region.</exception>
        /// <exception cref="ArgumentOutOfRangeException">An L-moment is nonfinite.</exception>
        /// <exception cref="InvalidOperationException">The Newton iteration fails, or the resulting location and scale are not finite and valid.</exception>
        /// <remarks>
        /// Uses Hosking's Newton iteration with tolerance 1e-6, at most 20 iterations, and at most ten
        /// step reductions. This estimator retains its h &gt; -1 estimation region. Continuous
        /// quantile limits and their derivatives are integrated when gamma-ratio differences would
        /// lose accuracy near a zero shape; the actual shape values are retained.
        /// </remarks>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            // This routine is taken and converted directly from Fortran code
            //***********************************************************************
            //*                                                                     *
            //* FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525,  *
            //* 'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3'  *
            //*                                                                     *
            //* J.R.M.HOSKING                                                       *
            //* IBM RESEARCH DIVISION                                               *
            //* T.J.WATSON RESEARCH CENTER                                          *
            //* YORKTOWN HEIGHTS                                                    *
            //* NEW YORK 10598, U.S.A.                                              *
            //*                                                                     *
            //* VERSION 3     AUGUST 1996                                           *
            //*                                                                     *
            //***********************************************************************

            if (moments is null) throw new ArgumentNullException(nameof(moments));
            if (moments.Count != 4) throw new ArgumentException("Exactly four L-moments are required.", nameof(moments));
            for (int index = 0; index < 4; index++)
                if (!Tools.IsFinite(moments[index]))
                    throw new ArgumentOutOfRangeException(nameof(moments), "L-moments must be finite.");

            double L1 = moments[0];
            double L2 = moments[1];
            double T3 = moments[2];
            double T4 = moments[3];
            double eps = 1E-6;
            int maxit = 20, maxsr = 10;

            //  TEST FOR FEASIBILITY

            if (L2 <= 0d) throw new ArgumentException("L-moments invalid.");
            if (Math.Abs(T3) >= 1d || Math.Abs(T4) >= 1d) throw new ArgumentException("L-moments invalid.");
            if (T4 <= (5d * T3 * T3 - 1d) / 4d) throw new ArgumentException("L-moments invalid.");
            if (T4 >= (5d * T3 * T3 + 1d) / 6d) throw new ArgumentException("(TAU-3, TAU-4) lies above the Generalized Logistic (suggests that L-moments are not constant with any Kappa distribution with hondo > -1).");

            //  SET STARTING VALUES FOR N-R ITERATION:
            //  G IS CHOSEN TO GIVE THE CORRECT VALUE OF TAU - 3 ON THE
            //  ASSUMPTION THAT H = 1(I.E.A GENERALIZED PARETO FIT) -
            //  BUT H IS ACTUALLY SET TO 1.001 TO AVOID NUMERICAL
            //  DIFFICULTIES WHICH CAN SOMETIMES ARISE WHEN H = 1 EXACTLY

            double G = (1d - 3d * T3) / (1d + T3);
            double H = 1.001d;
            double Z = G + H * 0.725d;
            double XDIST = 10d, DIST = 0;
            double U1 = 0, U2 = 0, U3 = 0, U4 = 0;
            double ALAM2 = 0, ALAM3, ALAM4;
            double TAU3 = 0, TAU4 = 0;
            double E1 = 0, E2 = 0;
            double DEL1 = 0, DEL2 = 0;
            double XG = 0, XH = 0, XZ, RHH;
            double U1G, U2G, U3G, U4G, U1H, U2H, U3H, U4H;
            double DL2G, DL2H, DL3G, DL3H, DL4G, DL4H;
            double D11, D12, D21, D22, DET, H11, H12, H21, H22;
            double FACTOR;

            //  START OF NEWTON-RAPHSON ITERATION

            bool converged = false;
            for (int i = 1; i <= maxit; i++)
            {

                //  REDUCE STEPLENGTH UNTIL WE ARE NEARER TO THE REQUIRED
                //  VALUES OF TAU-3 AND TAU-4 THAN WE WERE AT THE PREVIOUS STEP

                for (int j = 1; j <= maxsr; j++)
                {
                    if (!Tools.IsFinite(G) || !Tools.IsFinite(H) || G > 53d)
                        throw new InvalidOperationException("L-moment iteration encountered a nonfinite shape or an unsupported gamma calculation.");

                    if (KappaLinearMomentsNeedIntegration(G, H))
                    {
                        double[] linearMoments = KappaStandardLinearMoments(G, H);
                        ALAM2 = linearMoments[1];
                        ALAM3 = linearMoments[2];
                        ALAM4 = linearMoments[3];
                    }
                    else
                    {
                        if (H < 0)
                        {
                            U1 = Math.Exp(Gamma.LogGamma(-1d / H - G) - Gamma.LogGamma(-1d / H + 1d));
                            U2 = Math.Exp(Gamma.LogGamma(-2d / H - G) - Gamma.LogGamma(-2d / H + 1d));
                            U3 = Math.Exp(Gamma.LogGamma(-3d / H - G) - Gamma.LogGamma(-3d / H + 1d));
                            U4 = Math.Exp(Gamma.LogGamma(-4d / H - G) - Gamma.LogGamma(-4d / H + 1d));
                        }
                        else
                        {
                            U1 = Math.Exp(Gamma.LogGamma(1d / H) - Gamma.LogGamma(1d / H + 1d + G));
                            U2 = Math.Exp(Gamma.LogGamma(2d / H) - Gamma.LogGamma(2d / H + 1d + G));
                            U3 = Math.Exp(Gamma.LogGamma(3d / H) - Gamma.LogGamma(3d / H + 1d + G));
                            U4 = Math.Exp(Gamma.LogGamma(4d / H) - Gamma.LogGamma(4d / H + 1d + G));
                        }
                        ALAM2 = U1 - 2d * U2;
                        ALAM3 = -U1 + 6d * U2 - 6d * U3;
                        ALAM4 = U1 - 12d * U2 + 30d * U3 - 20d * U4;
                    }
                    if (ALAM2 == 0d || !Tools.IsFinite(ALAM2) || !Tools.IsFinite(ALAM3) || !Tools.IsFinite(ALAM4))
                        throw new InvalidOperationException("L-moment iteration could not evaluate finite L-moment ratios.");
                    TAU3 = ALAM3 / ALAM2;
                    TAU4 = ALAM4 / ALAM2;
                    E1 = TAU3 - T3;
                    E2 = TAU4 - T4;

                    // IF NEARER THAN BEFORE, EXIT THIS LOOP
                    DIST = Math.Max(Math.Abs(E1), Math.Abs(E2));
                    if (DIST < XDIST) break;

                    // OTHERWISE, HALVE THE STEPLENGTH AND TRY AGAIN
                    DEL1 *= 0.5;
                    DEL2 *= 0.5;
                    G = XG - DEL1;
                    H = XH - DEL2;
                    Z = G + H * 0.725d;

                    // TOO MANY STEPLENGTH REDUCTIONS
                    if (j == maxsr) throw new InvalidOperationException("L-moment iteration failed after ten step reductions.");
                }

                //  TEST FOR CONVERGENCE
                if (DIST < eps)
                {
                    converged = true;
                    break;
                }
                if (i == maxit) break;

                //  NOT CONVERGED: CALCULATE NEXT STEP
                //  NOTATION:
                //  U1G  - DERIVATIVE OF U1 W.R.T.G
                //  DL2G - DERIVATIVE OF ALAM2 W.R.T.G
                //  D..  - MATRIX OF DERIVATIVES OF TAU-3 AND TAU-4 W.R.T.G AND H
                //  H..  - INVERSE OF DERIVATIVE MATRIX
                //  DEL. - STEPLENGTH

                XG = G;
                XH = H;
                XZ = Z;
                XDIST = DIST;
                if (KappaLinearMomentsNeedIntegration(G, H))
                {
                    DL2G = KappaLinearMomentDerivative(G, H, 1, 2);
                    DL2H = KappaLinearMomentDerivative(G, H, 1, 3);
                    DL3G = KappaLinearMomentDerivative(G, H, 2, 2);
                    DL3H = KappaLinearMomentDerivative(G, H, 2, 3);
                    DL4G = KappaLinearMomentDerivative(G, H, 3, 2);
                    DL4H = KappaLinearMomentDerivative(G, H, 3, 3);
                }
                else
                {
                    RHH = 1d / (H * H);

                    if (H > 0)
                    {
                        U1G = -U1 * Gamma.Digamma(1d / H + 1d + G);
                        U2G = -U2 * Gamma.Digamma(2d / H + 1d + G);
                        U3G = -U3 * Gamma.Digamma(3d / H + 1d + G);
                        U4G = -U4 * Gamma.Digamma(4d / H + 1d + G);
                        U1H = RHH * (-U1G - U1 * Gamma.Digamma(1d / H));
                        U2H = 2d * RHH * (-U2G - U2 * Gamma.Digamma(2d / H));
                        U3H = 3d * RHH * (-U3G - U3 * Gamma.Digamma(3d / H));
                        U4H = 4d * RHH * (-U4G - U4 * Gamma.Digamma(4d / H));
                    }
                    else
                    {
                        U1G = -U1 * Gamma.Digamma(-1d / H - G);
                        U2G = -U2 * Gamma.Digamma(-2d / H - G);
                        U3G = -U3 * Gamma.Digamma(-3d / H - G);
                        U4G = -U4 * Gamma.Digamma(-4d / H - G);
                        U1H = RHH * (-U1G - U1 * Gamma.Digamma(-1d / H + 1d));
                        U2H = 2d * RHH * (-U2G - U2 * Gamma.Digamma(-2d / H + 1d));
                        U3H = 3d * RHH * (-U3G - U3 * Gamma.Digamma(-3d / H + 1d));
                        U4H = 4d * RHH * (-U4G - U4 * Gamma.Digamma(-4d / H + 1d));
                    }

                    DL2G = U1G - 2d * U2G;
                    DL2H = U1H - 2d * U2H;
                    DL3G = -U1G + 6d * U2G - 6d * U3G;
                    DL3H = -U1H + 6d * U2H - 6d * U3H;
                    DL4G = U1G - 12d * U2G + 30d * U3G - 20d * U4G;
                    DL4H = U1H - 12d * U2H + 30d * U3H - 20d * U4H;
                }
                D11 = (DL3G - TAU3 * DL2G) / ALAM2;
                D12 = (DL3H - TAU3 * DL2H) / ALAM2;
                D21 = (DL4G - TAU4 * DL2G) / ALAM2;
                D22 = (DL4H - TAU4 * DL2H) / ALAM2;
                DET = D11 * D22 - D12 * D21;
                if (DET == 0d || !Tools.IsFinite(DET))
                    throw new InvalidOperationException("The L-moment derivative matrix is singular or nonfinite.");
                H11 = D22 / DET;
                H12 = -D12 / DET;
                H21 = -D21 / DET;
                H22 = D11 / DET;
                DEL1 = E1 * H11 + E2 * H12;
                DEL2 = E1 * H21 + E2 * H22;

                //  TAKE NEXT N-R STEP

                G = XG - DEL1;
                H = XH - DEL2;
                Z = G + H * 0.725;

                //  REDUCE STEP IF G AND H ARE OUTSIDE THE PARAMETER SPACE
                FACTOR = 1d;
                if (G <= -1d) FACTOR = 0.8 * (XG + 1d) / DEL1;
                if (H <= -1d) FACTOR = Math.Min(FACTOR, 0.8 * (XH + 1d) / DEL2);
                if (Z <= -1d) FACTOR = Math.Min(FACTOR, 0.8 * (XZ + 1d) / (XZ - Z));
                if (H <= 0 && G * H <= -1d) FACTOR = Math.Min(FACTOR, 0.8 * (XG * XH + 1d) / (XG * XH - G * H));
                if (FACTOR != 1d)
                {
                    DEL1 = DEL1 * FACTOR;
                    DEL2 = DEL2 * FACTOR;
                    G = XG - DEL1;
                    H = XH - DEL2;
                    Z = G + H * 0.725;
                }

            }

            if (!converged)
                throw new InvalidOperationException("L-moment iterations failed to converge after 20 iterations.");

            // Reevaluate at the accepted shapes: step reductions must not leave stale scale factors.
            double[] standardized = KappaStandardLinearMoments(G, H);
            double alpha = L2 / standardized[1];
            double xi = L1 - alpha * standardized[0];
            if (!Tools.IsFinite(xi) || !Tools.IsFinite(alpha) || alpha <= 0d)
                throw new InvalidOperationException("L-moment shapes converged, but finite valid location and scale could not be recovered.");
            return [xi, alpha, G, H];

        }

        /// <summary>
        /// Calculates the first two L-moments, L-skewness, and L-kurtosis for the given parameters.
        /// </summary>
        /// <param name="parameters">The location, scale, kappa, and hondo parameters.</param>
        /// <returns>L1, L2, τ3, and τ4, in that order.</returns>
        /// <exception cref="ArgumentNullException">The parameters are null.</exception>
        /// <exception cref="ArgumentException">There are not four parameters.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or the first absolute moment does not exist.</exception>
        /// <exception cref="InvalidOperationException">Finite L-moments could not be evaluated numerically.</exception>
        /// <remarks>
        /// L-moments require κ &gt; -1 and, when h &lt; 0, κh &gt; -1. Gamma products are
        /// evaluated as log-gamma ratios. Near zero shapes, adaptive quantile integration avoids
        /// cancellation without modifying the parameters; the exact Gumbel limit is explicit.
        /// </remarks>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            if (parameters is null) throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count != 4) throw new ArgumentException("Exactly four parameters are required.", nameof(parameters));
            ValidateParameters(parameters, true);
            double kappa = parameters[2], hondo = parameters[3];
            if (kappa <= -1d || (hondo < 0d && kappa >= -1d / hondo))
                throw new ArgumentOutOfRangeException(nameof(parameters), "L-moments require kappa > -1 and, for hondo < 0, kappa*hondo > -1.");

            double[] standardized = KappaStandardLinearMoments(kappa, hondo);
            if (!Tools.IsFinite(standardized[0]) || !Tools.IsFinite(standardized[1]) || standardized[1] <= 0d ||
                !Tools.IsFinite(standardized[2]) || !Tools.IsFinite(standardized[3]))
                throw new InvalidOperationException("Finite standardized L-moments could not be evaluated.");
            double l1 = parameters[0] + parameters[1] * standardized[0];
            double l2 = parameters[1] * standardized[1];
            if (!Tools.IsFinite(l1) || !Tools.IsFinite(l2) || l2 <= 0d)
                throw new InvalidOperationException("The L-moments are outside the finite numerical range.");
            return [l1, l2, standardized[2] / standardized[1], standardized[3] / standardized[1]];
        }

        /// <summary>
        /// Identifies shape neighborhoods where gamma-ratio differences or their derivatives
        /// cancel; this selects an equivalent numerical representation, not a limiting shape.
        /// </summary>
        private static bool KappaLinearMomentsNeedIntegration(double kappa, double hondo)
        {
            return Math.Abs(kappa) < 0.001d || Math.Abs(hondo) < 0.001d;
        }

        /// <summary>
        /// Evaluates standardized L1 through L4 using log-gamma probability-weighted moments
        /// or full-interval adaptive integration of the actual quantile near zero shapes.
        /// </summary>
        private static double[] KappaStandardLinearMoments(double kappa, double hondo)
        {
            if (kappa == 0d && hondo == 0d)
            {
                double log2 = Math.Log(2d), log3 = Math.Log(3d);
                return [0.577215664901532860606512090082402431d, log2, 2d * log3 - 3d * log2, 16d * log2 - 10d * log3];
            }

            // At h=0 the gamma ratio has an exact simpler form; k=0 is evaluated by
            // the continuous quantile itself instead of a small artificial replacement.
            if (Math.Abs(kappa) < 0.001d || (hondo != 0d && Math.Abs(hondo) < 0.001d))
            {
                var integrated = new double[4];
                for (int order = 0; order < 4; order++)
                {
                    int polynomialOrder = order;
                    integrated[order] = IntegrateProbability(logProbability =>
                        StandardQuantile(logProbability, kappa, hondo) * KappaLinearMomentPolynomial(logProbability, polynomialOrder));
                }
                return integrated;
            }

            var beta = new double[4];
            double logGamma = Gamma.LogGamma(1d + kappa);
            for (int r = 1; r <= 4; r++)
            {
                double logG;
                if (hondo == 0d)
                    logG = logGamma - kappa * Math.Log(r);
                else if (hondo > 0d)
                    logG = Math.Log(r) + logGamma + Gamma.LogGamma(r / hondo) -
                        (1d + kappa) * Math.Log(hondo) - Gamma.LogGamma(1d + kappa + r / hondo);
                else
                    logG = Math.Log(r) + logGamma + Gamma.LogGamma(-kappa - r / hondo) -
                        (1d + kappa) * Math.Log(-hondo) - Gamma.LogGamma(1d - r / hondo);
                beta[r - 1] = -Tools.Expm1(logG) / (kappa * r);
            }
            return [beta[0], 2d * beta[1] - beta[0],
                6d * beta[2] - 6d * beta[1] + beta[0],
                20d * beta[3] - 30d * beta[2] + 12d * beta[1] - beta[0]];
        }

        /// <summary>
        /// Integrates a continuous shape derivative against the shifted Legendre polynomial
        /// needed by Hosking's Newton derivative matrix.
        /// </summary>
        private static double KappaLinearMomentDerivative(double kappa, double hondo, int order, int component)
        {
            return IntegrateProbability(logProbability =>
                StandardQuantileGradient(logProbability, kappa, hondo)[component] * KappaLinearMomentPolynomial(logProbability, order));
        }

        /// <summary>Evaluates the first four shifted Legendre polynomials on the probability interval.</summary>
        private static double KappaLinearMomentPolynomial(double logProbability, int order)
        {
            double probability = Math.Exp(logProbability);
            if (order == 0) return 1d;
            if (order == 1) return 2d * probability - 1d;
            if (order == 2) return (6d * probability - 6d) * probability + 1d;
            return ((20d * probability - 30d) * probability + 12d) * probability - 1d;
        }


        /// <inheritdoc/>
        /// <remarks>
        /// Retains the existing L-moment bounds. An initializer must have finite sample likelihood;
        /// the existing GEV initializer is tried within those bounds when the Kappa initializer fails.
        /// If L-moment construction fails before bounds are available, the existing GEV bounds are used.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than four observations or contains a nonfinite value.</exception>
        /// <exception cref="InvalidOperationException">Neither initializer has finite likelihood inside the fitting bounds.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            ValidateFittingSample(sample);
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            bool haveBounds = false;
            Exception? initializationFailure = null;

            // Get initial values
            try
            {
                initialVals = ParametersFromLinearMoments(Statistics.LinearMoments(sample));

                // Get bounds of location
                if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
                lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
                upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));

                // Get bounds of scale
                lowerVals[1] = Tools.DoubleMachineEpsilon;
                upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[1])) + 1d));

                // Get bounds of shape
                lowerVals[2] = -10d;
                upperVals[2] = 10d;

                // Get bounds of shape 2
                lowerVals[3] = -2d;
                upperVals[3] = 2d;
                haveBounds = true;

                // Correct initial value of kappa and hondo if necessary
                if (initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
                {
                    initialVals[2] = 0d;
                }

                if (initialVals[3] <= lowerVals[3] || initialVals[3] >= upperVals[3])
                {
                    initialVals[3] = 0d;
                }

                if (IsUsableInitializer(initialVals, lowerVals, upperVals, sample))
                    return Tuple.Create(initialVals, lowerVals, upperVals);
                initializationFailure = new InvalidOperationException("The Kappa L-moment initializer does not have finite likelihood within the fitting bounds.");
            }
            catch (Exception exception) when (exception is ArgumentException || exception is InvalidOperationException || exception is ArithmeticException)
            {
                initializationFailure = exception;
            }

            try
            {
                // Use the existing GEV initializer, preserving already established Kappa bounds.
                var gev = new GeneralizedExtremeValue();
                var parms = gev.GetParameterConstraints(sample);
                for (int i = 0; i < 3; i++)
                {
                    initialVals[i] = parms.Item1[i];
                    if (!haveBounds)
                    {
                        lowerVals[i] = parms.Item2[i];
                        upperVals[i] = parms.Item3[i];
                    }
                }
                initialVals[3] = 0d;
                lowerVals[3] = -2d;
                upperVals[3] = 2d;
                if (IsUsableInitializer(initialVals, lowerVals, upperVals, sample))
                    return Tuple.Create(initialVals, lowerVals, upperVals);
            }
            catch (Exception exception) when (exception is ArgumentException || exception is InvalidOperationException || exception is ArithmeticException)
            {
                throw new InvalidOperationException("Neither Kappa nor GEV initialization produced a valid Kappa Four fitting candidate.",
                    new AggregateException(initializationFailure!, exception));
            }
            throw new InvalidOperationException("Neither Kappa nor GEV initialization has finite sample likelihood within the fitting bounds.", initializationFailure);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Uses the existing bounded Nelder-Mead optimizer. Returns only after successful solver
        /// termination with valid parameters and finite sample likelihood. The unrestricted Kappa
        /// likelihood can be unbounded, so numerical success is not proof of a finite global MLE.
        /// </remarks>
        /// <exception cref="InvalidOperationException">Initialization fails, the solver does not succeed, or the final fit is invalid or has nonfinite likelihood.</exception>
        public double[] MLE(IList<double> sample)
        {
            // Set constraints
            var tuple = GetParameterConstraints(sample);
            var Initials = tuple.Item1;
            var Lowers = tuple.Item2;
            var Uppers = tuple.Item3;

            // Solve using the existing Nelder-Mead configuration.
            double logLH(double[] x)
            {
                var K4 = new KappaFour();
                K4.SetParameters(x);
                return K4.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            if (solver.Status != OptimizationStatus.Success)
                throw new InvalidOperationException($"Kappa Four maximum likelihood estimation failed with optimizer status {solver.Status}.");
            if (!IsUsableInitializer(solver.BestParameterSet.Values, Lowers, Uppers, sample))
                throw new InvalidOperationException("Kappa Four maximum likelihood estimation returned invalid parameters or nonfinite sample likelihood.");
            return solver.BestParameterSet.Values;

        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            EnsureValidParameters();
            if (double.IsNaN(x)) return double.NaN;
            if (double.IsInfinity(x) || x < Minimum || x > Maximum) return 0d;
            if (x == Minimum) return LowerEndpointDensity();
            if (x == Maximum) return Kappa < 1d ? 0d : Kappa == 1d ? 1d / Alpha : double.PositiveInfinity;
            return Math.Exp(InteriorLogDensity(x));
        }

        /// <summary>Validates the observations needed by the four-parameter fitting initializers.</summary>
        private static void ValidateFittingSample(IList<double> sample)
        {
            if (sample is null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 4) throw new ArgumentOutOfRangeException(nameof(sample), "At least four observations are required for Kappa Four estimation.");
            for (int i = 0; i < sample.Count; i++)
                if (!Tools.IsFinite(sample[i])) throw new ArgumentOutOfRangeException(nameof(sample), "Observations must be finite.");
        }

        /// <summary>Checks parameter validity, unchanged fitting bounds, and finite sample likelihood.</summary>
        private bool IsUsableInitializer(double[] parameters, double[] lower, double[] upper, IList<double> sample)
        {
            if (ValidateParameters(parameters, false) != null) return false;
            for (int i = 0; i < NumberOfParameters; i++)
                if (!Tools.IsFinite(lower[i]) || !Tools.IsFinite(upper[i]) || lower[i] >= upper[i]
                    || parameters[i] < lower[i] || parameters[i] > upper[i]) return false;
            var candidate = new KappaFour(parameters[0], parameters[1], parameters[2], parameters[3]);
            return Tools.IsFinite(candidate.LogLikelihood(sample));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluates interior densities directly in log space, including finite log densities whose
        /// ordinary density underflows or overflows. Infinite endpoint density limits remain
        /// positive infinity; aggregate likelihood evaluation applies its own contribution convention.
        /// </remarks>
        public override double LogPDF(double x)
        {
            EnsureValidParameters();
            if (double.IsNaN(x) || double.IsInfinity(x) || x < Minimum || x > Maximum) return double.NegativeInfinity;
            if (x == Minimum || x == Maximum)
            {
                double density = PDF(x);
                return density > 0d ? Math.Log(density) : double.NegativeInfinity;
            }
            return InteriorLogDensity(x);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            return Math.Exp(LogCDF(x));
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            EnsureValidParameters();
            if (x <= Minimum) return double.NegativeInfinity;
            if (x >= Maximum) return 0d;
            return InteriorLogProbability(x, LogT(x));
        }

        /// <inheritdoc/>
        public override double CCDF(double x)
        {
            return -Tools.Expm1(LogCDF(x));
        }

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            EnsureValidParameters();
            if (x <= Minimum) return 0d;
            if (x >= Maximum) return double.NegativeInfinity;
            double logT = LogT(x);
            double logF = InteriorLogProbability(x, logT);
            // 1-F ~ t when t is too small to represent; retain its logarithm.
            if (logF == 0d) return logT;
            return logF < -Math.Log(2d) ? Tools.Log1p(-Math.Exp(logF)) : Math.Log(-Tools.Expm1(logF));
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            // Validate parameters
            if (_parametersValid == false) ValidateParameters([Xi, Alpha, Kappa, Hondo], true);


            return AffineQuantile(Xi, Alpha, StandardQuantile(Math.Log(probability), Kappa, Hondo));
        }

        /// <summary>Validates the current parameters before numerical evaluation.</summary>
        private void EnsureValidParameters()
        {
            if (!_parametersValid) ValidateParameters(GetParameters, true);
        }

        /// <summary>Returns the continuous exponential divided difference (exp(x)-1)/x.</summary>
        private static double ExponentialRelative(double x)
        {
            return x == 0d ? 1d : Tools.Expm1(x) / x;
        }

        /// <summary>Returns log((exp(x)-1)/x) without intermediate overflow or cancellation.</summary>
        private static double LogExponentialRelative(double x)
        {
            if (Math.Abs(x) < 1E-4)
                return Tools.Log1p(x * (0.5d + x * (1d / 6d + x * (1d / 24d + x * (1d / 120d + x / 720d)))));
            return x > 0d ? x + Math.Log(-Tools.Expm1(-x)) - Math.Log(x)
                : Math.Log(-Tools.Expm1(x)) - Math.Log(-x);
        }

        /// <summary>Evaluates the standardized quantile from log probability, retaining nonzero shapes.</summary>
        private static double StandardQuantile(double logProbability, double k, double h)
        {
            return QuantileFromLogT(QuantileLogT(logProbability, h), k);
        }

        /// <summary>Returns the quantile latent log(t), including products outside the finite double range.</summary>
        private static double QuantileLogT(double logProbability, double h)
        {
            double s = h * logProbability;
            if (double.IsPositiveInfinity(s)) return double.PositiveInfinity;
            if (h > 0d && s < -0.5d) return Math.Log(-Tools.Expm1(s)) - Math.Log(h);
            return Math.Log(-logProbability) + LogExponentialRelative(s);
        }

        /// <summary>Converts log(t) to a standardized quantile without indeterminate products at limits.</summary>
        private static double QuantileFromLogT(double w, double k)
        {
            if (k == 0d) return -w;
            double v = k * w;
            return Math.Abs(v) < 0.5d ? -w * ExponentialRelative(v) : -Tools.Expm1(v) / k;
        }

        /// <summary>Applies location and scale, recovering finite cancellation after product overflow.</summary>
        private static double AffineQuantile(double location, double scale, double value)
        {
            double product = scale * value;
            if (double.IsInfinity(product) && Tools.IsFinite(value) && Math.Sign(location) != Math.Sign(product))
                return scale * (value + location / scale);
            return location + product;
        }

        /// <summary>Evaluates the finite-support affine endpoint without overflowing an intermediate quotient.</summary>
        private double LocationPlusScaleOverShape()
        {
            double shift = Alpha / Kappa;
            if (double.IsInfinity(shift) && Math.Sign(Xi) != Math.Sign(shift)) return (Xi * Kappa + Alpha) / Kappa;
            return Xi + shift;
        }

        /// <summary>
        /// Integrates a function of log probability over the full unit probability interval.
        /// Each half uses p=u^8/2 (or its survival counterpart), regularizing endpoint singularities
        /// and avoiding the loss of upper-tail probabilities to subtraction from one.
        /// </summary>
        /// <param name="integrand">The function evaluated at log(p).</param>
        /// <returns>The integral, or NaN if quadrature fails its status or error checks.</returns>
        private static double IntegrateProbability(Func<double, double> integrand)
        {
            var integration = new AdaptiveGaussKronrod(u =>
            {
                double logU = Math.Log(u);
                double logTail = 8d * logU - Math.Log(2d);
                double jacobian = Math.Exp(Math.Log(4d) + 7d * logU);
                // Pair the two tails before estimating the error so relative error is measured
                // against the complete integral, including cancellation in a near-zero mean.
                double value = integrand(logTail) * jacobian + integrand(Tools.Log1p(-Math.Exp(logTail))) * jacobian;
                if (!Tools.IsFinite(value)) throw new ArithmeticException("The Kappa Four quantile integrand is nonfinite.");
                return value;
            }, 0d, 1d)
            {
                RelativeTolerance = 1E-10,
                AbsoluteTolerance = 1E-12,
                ReportFailure = false
            };
            integration.Integrate();
            if (integration.Status != IntegrationStatus.Success || !Tools.IsFinite(integration.Result)
                || !Tools.IsFinite(integration.StandardError)) return double.NaN;
            return integration.StandardError <= Math.Max(1E-12, 1E-10 * Math.Abs(integration.Result)) ? integration.Result : double.NaN;
        }

        /// <summary>Tests absolute existence of an ordinary moment of the given positive order.</summary>
        private bool MomentExists(int order)
        {
            return order * Kappa > -1d && (Hondo >= 0d || order * Kappa * Hondo > -1d);
        }

        /// <summary>
        /// Computes moments in standardized coordinates, centering before higher powers and then
        /// applying location and scale. Nonexistent or numerically unresolved moments are NaN.
        /// </summary>
        private double[] ComputeMoments()
        {
            EnsureValidParameters();
            double[] result = [double.NaN, double.NaN, double.NaN, double.NaN];
            if (!MomentExists(1)) return result;
            double mean = IntegrateProbability(logP => StandardQuantile(logP, Kappa, Hondo));
            if (!Tools.IsFinite(mean)) return result;
            result[0] = AffineQuantile(Xi, Alpha, mean);
            if (!MomentExists(2)) return result;
            double variance = IntegrateProbability(logP => Math.Pow(StandardQuantile(logP, Kappa, Hondo) - mean, 2d));
            if (!(variance > 0d) || !Tools.IsFinite(variance)) return result;
            double sd = Math.Sqrt(variance);
            result[1] = Alpha * sd;
            if (MomentExists(3)) result[2] = IntegrateProbability(logP => Math.Pow((StandardQuantile(logP, Kappa, Hondo) - mean) / sd, 3d));
            if (MomentExists(4)) result[3] = IntegrateProbability(logP => Math.Pow((StandardQuantile(logP, Kappa, Hondo) - mean) / sd, 4d));
            return result;
        }

        /// <summary>Returns the derivative of the exponential divided difference.</summary>
        private static double ExponentialRelativeDerivative(double x)
        {
            if (Math.Abs(x) < 1E-3)
                return 0.5d + x * (1d / 3d + x * (1d / 8d + x * (1d / 30d + x * (1d / 144d + x / 840d))));
            if (x > 50d) return Math.Exp(x + Math.Log(x - 1d) - 2d * Math.Log(x));
            return (x * Math.Exp(x) - Tools.Expm1(x)) / (x * x);
        }

        /// <summary>Returns the derivative of log(exprel(x)), including its removable singularity.</summary>
        private static double LogExponentialRelativeDerivative(double x)
        {
            if (Math.Abs(x) < 1E-3) return 0.5d + x * (1d / 12d + x * x * (-1d / 720d + x * x / 30240d));
            return x > 0d ? 1d / -Tools.Expm1(-x) - 1d / x : Math.Exp(x) / Tools.Expm1(x) - 1d / x;
        }

        /// <summary>Returns the unit-scale quantile gradient in the order xi, alpha, kappa, hondo.</summary>
        private static double[] StandardQuantileGradient(double logProbability, double k, double h)
        {
            double s = h * logProbability;
            double w = QuantileLogT(logProbability, h);
            double v = k * w;
            double dh = h > 0d && s < -0.5d ? logProbability * (Math.Exp(s) / Tools.Expm1(s)) - 1d / h
                : logProbability * LogExponentialRelativeDerivative(s);
            double dk = v < -50d ? -1d / k / k
                : v > 50d ? -Math.Exp(v + Math.Log(v - 1d) - 2d * Math.Log(Math.Abs(k)))
                : -w * w * ExponentialRelativeDerivative(v);
            return [1d, QuantileFromLogT(w, k), dk, -Math.Exp(v) * dh];
        }

        /// <summary>Returns log(t), where t=(1-k*y)^(1/k), with its k=0 limit.</summary>
        private double LogT(double x)
        {
            double y = (x - Xi) / Alpha;
            if (double.IsInfinity(y) && Tools.IsFinite(x)) y = x / Alpha - Xi / Alpha;
            if (Kappa == 0d) return -y;
            double product = -Kappa * y;
            if (product < -0.9d) return KappaFourBoundary.LogT(x, Xi, Alpha, Kappa);
            double logBase = double.IsPositiveInfinity(product) ? Math.Log(Math.Abs(Kappa)) + Math.Log(Math.Abs(y)) : Tools.Log1p(product);
            return logBase / Kappa;
        }

        /// <summary>Evaluates log F from log(t) without exponentiating a large t.</summary>
        private double LogProbabilityFromLogT(double w)
        {
            if (Hondo == 0d) return -Math.Exp(w);
            double z = Math.Log(Math.Abs(Hondo)) + w;
            if (z < -36d) return -Math.Exp(w); // correction is smaller than a double rounding unit
            if (Hondo < 0d)
                return (z > 0d ? z + Tools.Log1p(Math.Exp(-z)) : Tools.Log1p(Math.Exp(z))) / Hondo;
            return (z < -Math.Log(2d) ? Tools.Log1p(-Math.Exp(z)) : Math.Log(-Tools.Expm1(z))) / Hondo;
        }

        /// <summary>Retains the small lower-support residual when ordinary log terms nearly cancel.</summary>
        private double InteriorLogProbability(double x, double w)
        {
            if (Hondo > 0d && Math.Abs(Math.Log(Hondo) + w) < 1E-5)
                return KappaFourBoundary.LogProbability(x, Xi, Alpha, Kappa, Hondo);
            return LogProbabilityFromLogT(w);
        }

        /// <summary>Evaluates the interior log density using the latent transform.</summary>
        private double InteriorLogDensity(double x)
        {
            double w = LogT(x);
            if (Hondo < 0d && Math.Log(-Hondo) + w > 0d)
            {
                // Combine the linear terms before multiplication when t is large.
                double z = Math.Log(-Hondo) + w;
                double correction = Math.Log(-Hondo) + Tools.Log1p(Math.Exp(-z));
                double result = (1d / Hondo - Kappa) * w + (1d / Hondo - 1d) * correction - Math.Log(Alpha);
                return double.IsNaN(result) ? (w + correction) / Hondo - Kappa * w - correction - Math.Log(Alpha) : result;
            }
            return (1d - Kappa) * w + (1d - Hondo) * InteriorLogProbability(x, w) - Math.Log(Alpha);
        }

        /// <summary>Returns the one-sided density limit at the finite lower support endpoint.</summary>
        private double LowerEndpointDensity()
        {
            if (Hondo > 0d) return Hondo < 1d ? 0d : Hondo == 1d ? 1d / Alpha : double.PositiveInfinity;
            if (Hondo == 0d) return 0d;
            double exponent = 1d / Hondo - Kappa;
            return exponent < 0d ? 0d : exponent > 0d ? double.PositiveInfinity
                : Math.Exp((1d / Hondo - 1d) * Math.Log(-Hondo) - Math.Log(Alpha));
        }

        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public override UnivariateDistributionBase Clone()
        {
            return new KappaFour(Xi, Alpha, Kappa, Hondo);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Local asymptotic MLE covariance in xi, alpha, kappa, hondo order. Finite regular
        /// information requires kappa &lt; 1/2, hondo &lt; 1/2 and kappa*hondo &lt; 1/2.
        /// These conditions do not narrow distribution validity or assert global-MLE existence.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Sample size, parameters or information regularity is invalid.</exception>
        /// <exception cref="InvalidOperationException">Numerical information or its inversion cannot be resolved.</exception>
        /// <exception cref="NotImplementedException">The requested estimator is not maximum likelihood.</exception>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            EnsureValidParameters();
            if (estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
                throw new NotImplementedException("Kappa Four covariance is implemented only for local maximum-likelihood uncertainty.");
            return KappaExpectedInformation.ParameterCovariance(Alpha, Kappa, Hondo, sampleSize, 4);
        }

        /// <inheritdoc/>
        /// <remarks>Applies the local-MLE delta method in common physical quantile coordinates,
        /// avoiding underflow or overflow from forming physical covariance entries first.</remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            EnsureValidParameters();
            var unit = new KappaFour(0, 1, Kappa, Hondo);
            double[,] covariance = unit.ParameterCovariance(sampleSize, estimationMethod);
            double logProbability = Math.Log(probability), w = QuantileLogT(logProbability, Hondo);
            double argument = Kappa * w, s = Hondo * logProbability;
            double dh = Hondo > 0 && s < -.5 ? logProbability * (Math.Exp(s) / Tools.Expm1(s)) - 1 / Hondo
                : logProbability * LogExponentialRelativeDerivative(s);
            double scaleGradient = DistributionNumerics.ScaledExprelProduct(Alpha, -w, argument);
            double shapeGradient = -DistributionNumerics.ScaledExprelDerivativeProduct(Alpha, w, argument);
            double hondoGradient = dh == 0 ? 0 : -Math.Sign(dh) * Math.Exp(Math.Log(Alpha) + argument + Math.Log(Math.Abs(dh)));
            return DistributionNumerics.ScaledQuantileVariance(covariance, [Alpha, scaleGradient, shapeGradient, hondoGradient]);
        }

        /// <inheritdoc/>
        /// <remarks>Uses continuous zero-shape limits in the parameter order xi, alpha, kappa, hondo.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">The probability is not strictly between zero and one, or distribution parameters are invalid.</exception>
        public double[] QuantileGradient(double probability)
        {
            EnsureValidParameters();
            if (double.IsNaN(probability) || probability <= 0d || probability >= 1d)
                throw new ArgumentOutOfRangeException(nameof(probability), "Quantile gradients require a probability strictly between zero and one.");
            double[] gradient = StandardQuantileGradient(Math.Log(probability), Kappa, Hondo);
            gradient[2] *= Alpha;
            gradient[3] *= Alpha;
            return gradient;
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }
    }
}
