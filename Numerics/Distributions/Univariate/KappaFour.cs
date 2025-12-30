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
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;
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
                    u = CentralMoments(1000);
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
        public override double Mode
        {
            get
            {
                var brent = new BrentSearch(PDF, InverseCDF(0.001), InverseCDF(0.999));
                brent.Maximize();
                return brent.BestParameterSet.Values[0];
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
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
                    u = CentralMoments(1000);
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
                    u = CentralMoments(1000);
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
                if (Hondo <= 0d && Kappa < 0d)
                {
                    return Xi + Alpha / Kappa;
                }
                else if (Hondo > 0d && Kappa != 0d)
                {
                    return Xi + Alpha / Kappa * (1d - Math.Pow(Hondo, -Kappa));
                }
                else if (Hondo > 0d && Kappa == 0d)
                {
                    return Xi + Alpha * Math.Log(Hondo);
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
                    return Xi + Alpha / Kappa;
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
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
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
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
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
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
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

        /// <inheritdoc/>
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

            double L1 = moments[0];
            double L2 = moments[1];
            double T3 = moments[2];
            double T4 = moments[3];
            double eps = 1E-6;
            int maxit = 20, maxsr = 10;
            double xi = 0, alpha = 0, kappa = 0, hondo = 0;

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

            for (int i = 1; i < maxit; i++)
            {

                //  REDUCE STEPLENGTH UNTIL WE ARE NEARER TO THE REQUIRED
                //  VALUES OF TAU-3 AND TAU-4 THAN WE WERE AT THE PREVIOUS STEP

                for (int j = 1; j <= maxsr; j++)
                {
                    if (G > 53d) throw new Exception("Iteration encountered numerical difficulties - overflow would have been likely to occur.");

                    if (H <= 0)
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
                    if (ALAM2 == 0d) throw new Exception("Iteration encountered numerical difficulties - overflow would have been likely to occur.");
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

                    // TOO MANY STEPLENGTH REDUCTIONS
                    if (j == maxsr) throw new Exception("Iteration encountered numerical difficulties - overflow would have been likely to occur.");
                }

                //  TEST FOR CONVERGENCE
                if (DIST < eps) break;

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
                    U4G = -U4 * Gamma.Digamma(-4d/ H - G);
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
                D11 = (DL3G - TAU3 * DL2G) / ALAM2;
                D12 = (DL3H - TAU3 * DL2H) / ALAM2;
                D21 = (DL4G - TAU4 * DL2G) / ALAM2;
                D22 = (DL4H - TAU4 * DL2H) / ALAM2;
                DET = D11 * D22 - D12 * D21;
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

                // NOT CONVERGED
                if (i == maxit) throw new Exception("Iterations failed to converge.");
            }

            hondo = H;
            kappa = G;
            var TEMP = Gamma.LogGamma(1d + G);
            if (TEMP > 170d) throw new Exception("Iteration for hondo and kappa converged, but overflow would have occurred when calculating xi and alpha.");
            var GAM = Math.Exp(TEMP);
            TEMP = (1d + G) * Math.Log(Math.Abs(H));
            if (TEMP > 170d) throw new Exception("Iteration for hondo and kappa converged, but overflow would have occurred when calculating xi and alpha.");
            var HH = Math.Exp(TEMP);
            alpha = L2 * G * HH / (ALAM2 * GAM);
            xi = L1 - alpha / G * (1d - GAM * U1 / HH);

            return [xi, alpha, kappa, hondo];

        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double xi = parameters[0];
            double alpha = parameters[1];
            double kappa = parameters[2];
            double hondo = parameters[3];
            if (((kappa < -1) && (hondo >= 0)) || ((hondo < 0) && ((kappa <= -1) || (kappa >= -1 / hondo))))
            {
                throw new ArgumentOutOfRangeException(nameof(parameters), "L-moments can only be defined for hondo (h) >= 0 and kappa (k) > -1, or if h < 0 and -1 < k < -1/h.");
            }

            double L1;
            double L2;
            double T3;
            double T4;
            if (kappa == 0.0d) kappa = Math.Pow(10d, -100);
            if (hondo == 0.0d)
            {
                L1 = xi + alpha * (1.0d - Gamma.Function(1.0d + kappa)) / kappa;
                L2 = alpha * (1.0d - Math.Pow(2.0d, -kappa)) * Gamma.Function(1.0d + kappa) / kappa;
                T3 = 2.0d * (1.0d - Math.Pow(3.0d, -kappa)) / (1.0d - Math.Pow(2.0d, -kappa)) - 3.0d;
                T4 = (5.0d * (1.0d - Math.Pow(4.0d, -kappa)) - 10.0d * (1.0d - Math.Pow(3.0d, -kappa)) + 6.0d * (1.0d - Math.Pow(2.0d, -kappa))) / (1.0d - Math.Pow(2.0d, -kappa));
            }
            else
            {
                var g = new double[4];
                for (int r = 1; r <= 4; r++)
                {
                    if (hondo > 0.0d)
                    {
                        g[r - 1] = r * Gamma.Function(1d + kappa) * Gamma.Function(r / hondo) / (Math.Pow(hondo, 1d + kappa) * Gamma.Function(1d + kappa + r / hondo));
                    }
                    else
                    {
                        g[r - 1] = r * Gamma.Function(1d + kappa) * Gamma.Function(-kappa - r / hondo) / (Math.Pow(-hondo, 1d + kappa) * Gamma.Function(1d - r / hondo));
                    }
                }

                L1 = xi + alpha * (1d - g[0]) / kappa;
                L2 = alpha * (g[0] - g[1]) / kappa;
                T3 = (-g[0] + 3d * g[1] - 2d * g[2]) / (g[0] - g[1]);
                T4 = -(-g[0] + 6d * g[1] - 10d * g[2] + 5d * g[3]) / (g[0] - g[1]);
            }
            // 
            return [L1, L2, T3, T4];
        }


        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];

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

                // Correct initial value of kappa and hondo if necessary
                if (initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
                {
                    initialVals[2] = 0d;
                }

                if (initialVals[3] <= lowerVals[3] || initialVals[3] >= upperVals[3])
                {
                    initialVals[3] = 0d;
                }

            }
            catch
            {

                // Get constraints from GEV
                var gev = new GeneralizedExtremeValue();
                var parms = gev.GetParameterConstraints(sample);
                for (int i = 0; i < 3; i++)
                {
                    initialVals[i] = parms.Item1[i];
                    lowerVals[i] = parms.Item2[i];
                    upperVals[i] = parms.Item3[i];

                    // Get bounds of shape 2
                    initialVals[3] = 0;
                    lowerVals[3] = -2d;
                    upperVals[3] = 2d;
                }

            }


            // 
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        public double[] MLE(IList<double> sample)
        {
            // Set constraints
            var tuple = GetParameterConstraints(sample);
            var Initials = tuple.Item1;
            var Lowers = tuple.Item2;
            var Uppers = tuple.Item3;

            // Solve using Powell
            double logLH(double[] x)
            {
                var K4 = new KappaFour();
                K4.SetParameters(x);
                return K4.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;

        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters([Xi, Alpha, Kappa, Hondo], true);
            if (x < Minimum || x > Maximum) return 0.0d;

            double F = CDF(x);
            double y = (x - Xi) / Alpha;
            double yy = 1 - Kappa * y;
            return (1 / Alpha) * Math.Pow(yy, 1 / Kappa - 1) * Math.Pow(F, 1 - Hondo);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameters([Xi, Alpha, this.Kappa, this.Hondo], true);
            if (x <= Minimum) return 0d;
            if (x >= Maximum) return 1d;

            double y = (x - Xi) / Alpha;
            double yy = 1 - Kappa * y;
            if (Kappa !=0 &&  Hondo != 0)
            {
                return Math.Pow(1 - Hondo * Math.Pow(yy, 1 / Kappa), 1 / Hondo);
            }
            else if (Kappa != 0 && Hondo == 0)
            {
                return Math.Exp(-Math.Pow(yy, 1 / Kappa));
            }
            else if (Kappa == 0 && Hondo != 0)
            {
                return Math.Pow(1 - Hondo * Math.Exp(-y), 1 / Hondo);
            }
            else
            {
                return Math.Exp(-Math.Exp(-y));
            }
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


            if (Kappa != 0 && Hondo != 0)
            {
                return Xi + Alpha / Kappa * (1 - Math.Pow((1 - Math.Pow(probability, Hondo)) / Hondo, Kappa));
            }
            else if (Kappa != 0 && Hondo == 0)
            {
                return Xi + (Alpha / Kappa) * (1 - Math.Pow(-Math.Log(probability), Kappa));
            }
            else if (Kappa == 0 && Hondo != 0)
            {
                return Xi - Alpha * Math.Log(1 - Math.Pow(probability, Hondo) / Hondo);
            }
            else
            {
                return Xi - Alpha * Math.Log(-Math.Log(probability));
            }
        }

        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public override UnivariateDistributionBase Clone()
        {
            return new KappaFour(Xi, Alpha, Kappa, Hondo);
        }

        /// <inheritdoc/>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            throw new NotImplementedException();
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            throw new NotImplementedException();
        }

        /// <inheritdoc/>
        public double[] QuantileGradient(double probability)
        {
            double a = Alpha, k = Kappa, h = Hondo, F = probability;
            double dxi = 1d;
            double da = (1d - Math.Pow((1d - Math.Pow(F, h)) / h, k)) / k;
            double dk = a * (-(1d - Math.Pow((1d - Math.Pow(F, h)) / h, k)) / (k * k) - Math.Pow((1d - Math.Pow(F, h)) / h, k) * Math.Log((1d - Math.Pow(F, h)) / h) / k);
            double x = 1d - Math.Pow(F, h);
            double dh = -(a * (Math.Pow(F, h) - h * Math.Log(F) * Math.Pow(F, h) - 1d) * Math.Sign(x) * Math.Pow(Math.Abs(x), k - 1d)) / (Math.Sign(h) * Math.Pow(Math.Abs(h), k + 1d));
            return [dxi, da, dk, dh];
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            if (probabilities.Count != NumberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(probabilities), "The number of probabilities must be the same length as the number of distribution parameters.");
            }
            // |a b c d|
            // |e f g h|
            // |i j k l|
            // |m n o p| 
            var jacobian = new Matrix(4);
            for (int i = 0; i < 4; i++)
            {
                // Get the gradient
                var dFdx = QuantileGradient(probabilities[i]);
                // Populate the Jacobian matrix
                for (int j = 0; j < 4; j++)
                    jacobian[i, j] = dFdx[j];
            }
            // Solve determinant with LU decomposition
            var LU = new LUDecomposition(jacobian);
            determinant = LU.Determinant();
            // Return Jacobian
            return jacobian.ToArray();
        }
    }
}