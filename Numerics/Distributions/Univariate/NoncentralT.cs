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
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Noncentral t probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Noncentral_t-distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class NoncentralT : UnivariateDistributionBase
    {
 
        /// <summary>
        /// Constructs a Noncentral t distribution with 10 degrees of freedom and noncentrality = 0.
        /// </summary>
        public NoncentralT()
        {
            SetParameters(10d, 0d);
        }

        /// <summary>
        /// Constructs a Noncentral t distribution with given degrees of freedom and noncentrality.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom.</param>
        /// <param name="noncentrality">The noncentrality parameter.</param>
        public NoncentralT(double degreesOfFreedom, double noncentrality)
        {
            SetParameters(degreesOfFreedom, noncentrality);
        }

        private int _degreesOfFreedom;
        private double _noncentrality;
       
        /// <summary>
        /// Gets and sets the degrees of freedom ν (nu) of the distribution.
        /// </summary>
        public int DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
            set
            {
                _parametersValid = ValidateParameters(value, Noncentrality, false) is null;
                _degreesOfFreedom = value;
            }
        }

        /// <summary>
        /// Gets and sets the noncentrality parameter μ (mu) of the distribution.
        /// </summary>
        public double Noncentrality
        {
            get { return _noncentrality; }
            set
            {
                _parametersValid = ValidateParameters(DegreesOfFreedom, value, false) is null;
                _noncentrality = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.NoncentralT; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Noncentral t"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "NCT"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Degrees of Freedom (ν)";
                parmString[1, 0] = "Noncentrality (μ)";
                parmString[0, 1] = DegreesOfFreedom.ToString();
                parmString[1, 1] = Noncentrality.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ν", "μ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(DegreesOfFreedom), nameof(Noncentrality)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [DegreesOfFreedom, Noncentrality]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (_degreesOfFreedom > 1)
                {
                    return Noncentrality * Math.Sqrt(DegreesOfFreedom / 2d) * Gamma.Function((DegreesOfFreedom - 1) / 2d) / Gamma.Function(DegreesOfFreedom / 2d);
                }
                else
                {
                    return double.NaN;
                }
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
                if (DegreesOfFreedom > 2)
                {
                    double a = DegreesOfFreedom * (1d + Math.Pow(Noncentrality, 2d)) / (DegreesOfFreedom - 2);
                    double b = DegreesOfFreedom * Math.Pow(Noncentrality, 2d) / 2d;
                    double c = Gamma.Function((DegreesOfFreedom - 1) / 2d) / Gamma.Function(DegreesOfFreedom / 2d);
                    return Math.Sqrt(a - b * c * c);
                }
                else
                {
                    return double.NaN;
                }
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {                
                return CentralMoments(1E-8)[2];
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                return CentralMoments(1E-8)[3];
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return double.NegativeInfinity; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [1.0d, double.NegativeInfinity]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="v">The degrees of freedom ν (nu). Range: ν > 0.</param>
        /// <param name="mu">The noncentrality parameter μ (mu).</param>
        public void SetParameters(double v, double mu)
        {
            DegreesOfFreedom = (int)v;
            Noncentrality = mu;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="v">The degrees of freedom ν (nu). Range: ν > 0.</param>
        /// <param name="mu">The noncentrality parameter μ (mu).</param>
        /// <param name="throwException"></param>
        public ArgumentOutOfRangeException ValidateParameters(double v, double mu, bool throwException)
        {
            if (v < 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
                return new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
            }
            if (double.IsNaN(mu) || double.IsInfinity(mu))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Noncentrality), "The noncentrality parameter μ (mu) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Noncentrality), "The noncentrality parameter μ (mu) must be a number.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            double u = Noncentrality; 
            double v = DegreesOfFreedom;
            if (x != 0)
            {
                double A = NCT_CDF(x * Math.Sqrt(1 + 2 / v), v + 2, u);
                double B = NCT_CDF(x, v, u);
                double C = v / x;
                return C * (A - B);
            }
            else
            {
                double A = Gamma.Function((v + 1) / 2);
                double B = Math.Sqrt(Math.PI * v) * Gamma.Function(v / 2);
                double C = Math.Exp(-(u * u) / 2);
                return (A / B) * C;
            }
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            return NCTDist(x, DegreesOfFreedom, Noncentrality);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            // 
            return NCT_INV(probability, DegreesOfFreedom, Noncentrality);
        }

        
        /// <summary>
        /// Cumulative probability at T of the non-central t-distribution
        /// with DF degrees of freedom (may be fractional) and non-centrality
        /// parameter DELTA.
        /// </summary>
        /// <param name="t">A single point in the distribution range.</param>
        /// <param name="df">The degrees of freedom.</param>
        /// <param name="delta">The noncentrality parameter.</param>
        private double NCT_CDF(double t, double df, double delta)
        {
            double Z;
            double ANS = 0d;
            try
            {
                ANS = NCTDist(t, df, delta);
            }
            catch (Exception ex)
            {
                // If the solver fails, then use approximation
                if (ex.ToString() == "Maximum number of iterations reached.")
                {
                    Z = (t * (1.0d - 1.0d / (4.0d * df)) - delta) / Math.Sqrt(1.0d + Math.Pow(t, 2d) / (2.0d * df));
                    ANS = Normal.StandardCDF(Z);
                }

                if (ANS < 0d)
                    ANS = 0d;
                if (ANS > 1d)
                    ANS = 1d;
            }

            return ANS;
        }

        /// <summary>
        /// Cumulative probability at T of the non-central t-distribution
        /// with DF degrees of freedom (may be fractional) and non-centrality
        /// parameter DELTA.
        /// </summary>
        /// <param name="t">A single point in the distribution range.</param>
        /// <param name="df">The degrees of freedom.</param>
        /// <param name="delta">The noncentrality parameter.</param>
        /// <remarks>
        /// The function is based on ALGORITHM AS 243  APPL. STATIST. (1989), VOL.38, NO. 1.
        /// Original FORTRAN code can be found at:
        /// http://people.sc.fsu.edu/~jburkardt/f77_src/asa243/asa243.html
        /// </remarks>
        private double NCTDist(double t, double df, double delta)
        {

            // REAL FUNCTION TNC(T, DF, DELTA, IFAULT)
            // 
            // ALGORITHM AS 243  APPL. STATIST. (1989), VOL.38, NO. 1
            // 
            // Cumulative probability at T of the non-central t-distribution
            // with DF degrees of freedom (may be fractional) and non-centrality
            // parameter DELTA.
            // 
            // Note - requires the following auxiliary routines
            // ALOGAM (X)                         - ACM 291 or AS 245
            // BETAIN (X, A, B, ALBETA, IFAULT)   - AS 63 (updated in ASR 19)
            // ALNORM (X, UPPER)                  - AS 66
            // 
            // Translated by William A. Huber.  www.quantdec.com
            // 
            double a;
            double ALBETA;
            double b;
            double DEL;
            var N = default(int);
            double ERRBD;
            double GEVEN;
            double GODD;
            double LAMBDA;
            double P;
            double q;
            double RXB;
            double s;
            double TT;
            double x;
            double XEVEN;
            double XODD;
            double TNC;
            bool NEGDEL;
            // 
            // Note - ITRMAX and ERRMAX may be changed to suit one's needs.
            // 
            const int ITRMAX = 1000;
            const double Errmax = 0.0000001d;

            // DATA ITRMAX/100.1/, ERRMAX/1.E-06/
            // 
            // Constants - R2PI = 1/ {GAMMA(1.5) * SQRT(2)} = SQRT(2 / PI)
            // ALNRPI = Ln(SQRT(Pi))
            // 
            const double zero = 0d;
            const double half = 0.5d;
            const double one = 1.0d;
            const double two = 2.0d;
            const double r2pi = 0.797884560802865d;
            const double alnrpi = 0.5723649429247d;
            TNC = zero;
            TT = t;
            DEL = delta;
            NEGDEL = false;
            if (t < zero)
            {
                NEGDEL = true;
                TT = -TT;
                DEL = -DEL;
            }
            // 
            // Initialize twin series (Guenther, J. Statist. Computn. Simuln.
            // vol.6, 199, 1978).
            // 
            x = t * t / (t * t + df);
            if (x <= zero)
                goto Twenty;
            LAMBDA = DEL * DEL;
            P = half * Math.Exp(-half * LAMBDA);
            q = r2pi * P * DEL;
            s = half - P;
            a = half;
            b = half * df;
            RXB = Math.Pow(one - x, b);
            ALBETA = alnrpi + Gamma.LogGamma(b) - Gamma.LogGamma(a + b);
            XODD = Beta.IncompleteRatio(x, a, b, ALBETA);
            GODD = two * RXB * Math.Exp(a * Math.Log(x) - ALBETA);
            XEVEN = one - RXB;
            GEVEN = b * x * RXB;
            TNC = P * XODD + q * XEVEN;
            // 
            // Repeat until convergence
            // 
            N = 1;
            do
            {
                a = a + one;
                XODD = XODD - GODD;
                XEVEN = XEVEN - GEVEN;
                GODD = GODD * x * (a + b - one) / a;
                GEVEN = GEVEN * x * (a + b - half) / (a + half);
                P = P * LAMBDA / (two * N);
                q = q * LAMBDA / (two * N + one);
                s = s - P;
                N = N + 1;
                TNC = TNC + P * XODD + q * XEVEN;
                ERRBD = two * s * (XODD - GODD);
            }
            while (ERRBD > Errmax & N <= ITRMAX);
            // 
            Twenty:
            ;
            if (N > ITRMAX)
            {
                throw new ArgumentException("Max number of iterations were exceeded.");
            }

            if (NEGDEL)
            {
                // TNC = one - TNC
                TNC = Normal.StandardCDF(DEL) - TNC;
            }
            else
            {
                TNC = TNC + (1.0d - Normal.StandardCDF(DEL));
            } // Upper tail area of N(0,1)
            return TNC;
        }

        /// <summary>
        /// The inverse of the non-central t distribution
        /// </summary>
        /// <param name="p">Probability between 0 and 1.</param>
        /// <param name="df">The degrees of freedom.</param>
        /// <param name="delta">The noncentrality parameter.</param>
        private double NCT_INV(double p, double df, double delta)
        {
            double t0, t1, t2, y0, y1, tInc, Slope;
            int iter;
            //double tEstimate;
            const double ytol = 0.0000001d;         // Y-tolerance
            const double xtol = 0.0000001d;         // X-tolerance
            //const double w = 1.5d;                 // Over-shooting parameter (>= 1)
            const int iterMax = 50;
            t0 = NCTInv0(p, df, delta);
            //tEstimate = t0;                          // Default if we run into trouble
            y0 = NCTDist(t0, df, delta) - p;
            if (y0 > ytol)
            {
                tInc = -1.0d;
            }
            else if (y0 < -ytol)
            {
                tInc = 1.0d;
            }
            else
            {
                return t0;
            }
            // 
            // Find a bracket through overshooting.
            // 
            t1 = t0 + tInc;
            y1 = NCTDist(t1, df, delta) - p;
            iter = 0;
            while (y0 < 0d != y1 > 0d & Math.Abs(t1 - t0) > xtol & iter < iterMax)
            {
                // 
                // Use secant method to extrapolate a zero, but overshoot by w >= 1.
                // 
                Slope = (y1 - y0) / (t1 - t0);
                if (Slope == 0d)
                {
                    return (t1 + t0) / 2.0d;
                }

                t2 = t0 - 2.0d * y0 / Slope;
                t0 = t1;
                t1 = t2;
                y0 = y1;
                y1 = NCTDist(t1, df, delta) - p;
                iter = iter + 1;
            }
            // Solve for T using Brent
            double ANS = ZBrent(t0, t1, y0, y1, xtol, df, delta, p);
            return ANS;
        }

        private double NCTInv0(double P, double N, double D)
        {
            // 
            // Approximates percentage points of the non-central t distribution.
            // P is the percentage, N is the degrees of freedom, D is the non-centrality parameter.
            // 
            // Non-central t is the ratio (U + D)/(Chi/Sqrt(n)) where U and Chi are independent
            // random variables distributed as Normal(0, 1) and Chi(n), respectively.
            // D is the non-centrality parameter.  When equal to 0, TInvNC is Student's t.
            // 
            // Source: Johnson & Kotz, Continuous Univariate Distributions, Volume 2.
            // 
            // NB: The number of iterations appears to go quadratically in D.  Thus, for large
            // D, we will be in trouble.
            // 
            // VB version (c) 2001 Quantitative Decisions.  All rights reserved.
            // Contact William A. Huber.  www.quantdec.com
            // 
            double b;
            double z;
            double b2;
            double u2;
            double T;
            // 
            // Establish entry conditions.
            // 
            // If N < 1.0# Or P <= 0.000001 Or P >= 0.999999 Then
            // Return Double.NaN
            // End If

            var StudentT = new StudentT(N);
            z = Normal.StandardZ(P);
            // 
            // Jennett & Welch approximation, formula (14.1).
            // Intended for large values of D^2, such as are used for most tolerance interval calculations.
            // 
            b = Math.Exp(Gamma.LogGamma((N + 1d) / 2d) - Gamma.LogGamma(N / 2d)) * Math.Sqrt(2d / N);
            u2 = z * z;
            b2 = b * b;
            T = b2 + (1d - b2) * (Math.Pow(D, 2d) - u2);
            if (T > 0d)
            {
                T = (D * b + z * Math.Sqrt(T)) / (b2 - u2 * (1d - b2));
            }
            else
            {
                T = 0d;
            }

            return T;
        }

        private double ZBrent(double X1, double X2, double y1, double y2, double Tol, double N, double Dnc, double Perc)
        {
            double ZBrentRet = default;
            // 
            // Finds the zero of NCTDist(x, n, d) - p given that [x1, x2] brackets the zero and
            // Y1 = value at X1, Y2 = value at X2.
            // 
            // Translated from Numerical Recipes (1986).
            // William A. Huber, 24 March 2001.
            // 
            double a;
            double b;
            var c = default(double);
            double fc;
            var D = default(double);
            var e = default(double);
            double tol1;
            double xm;
            double s;
            double P;
            double q;
            double r;
            double fa;
            double fb;
            int iter;
            const int itmax = 100;
            const double eps = 0.00000003d;
            a = X1;
            b = X2;
            fa = y1;
            fb = y2;
            if (fb * fa > 0d)
            {
                throw new ArgumentException("Brent's method failed because the root is not bracketed.");
            }

            fc = fb;
            for (iter = 1; iter <= itmax; iter++)
            {
                if (fb * fc > 0d)
                {
                    c = a;
                    fc = fa;
                    D = b - a;
                    e = D;
                }

                if (Math.Abs(fc) < Math.Abs(fb))
                {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }

                tol1 = 2.0d * eps * Math.Abs(b) + 0.5d * Tol;
                xm = 0.5d * (c - b);
                if (Math.Abs(xm) <= tol1 | fb == 0d)
                {
                    return b;
                }

                if (Math.Abs(e) >= tol1 & Math.Abs(fa) > Math.Abs(fb))
                {
                    s = fb / fa;
                    if (a == c)
                    {
                        P = 2.0d * xm * s;
                        q = 1.0d - s;
                    }
                    else
                    {
                        q = fa / fc;
                        r = fb / fc;
                        P = s * (2.0d * xm * q * (q - r) - (b - a) * (r - 1.0d));
                        q = (q - 1.0d) * (r - 1.0d) * (s - 1.0d);
                    }

                    if (P > 0d)
                        q = -q;
                    P = Math.Abs(P);
                    if (2.0d * P < 3.0d * xm * q - Math.Abs(tol1 * q) & 2.0d * P < Math.Abs(e * q))
                    {
                        e = D;
                        D = P / q;
                    }
                    else
                    {
                        D = xm;
                        e = D;
                    }
                }
                else
                {
                    D = xm;
                    e = D;
                }

                a = b;
                fa = fb;
                if (Math.Abs(D) > tol1)
                {
                    b = b + D;
                }
                else if (xm > 0d)  // b = b + SIGN(tol1, xm)
                {
                    b = b + Math.Abs(tol1);
                }
                else
                {
                    b = b - Math.Abs(tol1);
                }

                fb = NCTDist(b, N, Dnc) - Perc;
            }

            ZBrentRet = b;
            return ZBrentRet;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new NoncentralT(DegreesOfFreedom, Noncentrality);
        }

    }
}