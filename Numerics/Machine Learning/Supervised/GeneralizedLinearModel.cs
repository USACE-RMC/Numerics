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

using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;


namespace Numerics.MachineLearning
{
    /// <summary>
    /// A class for performing generalized linear regression. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class GeneralizedLinearModel
    {
        #region Construction

        /// <summary>
        /// Constructs a generalized linear model. 
        /// </summary>
        /// <param name="x">The matrix of predictor values.</param>
        /// <param name="y">The response vector.</param>
        /// <param name="hasIntercept">Determines if an intercept should be estimate. Default = true.</param>
        /// <param name="linkType">The link function type. Default = identity.</param>
        public GeneralizedLinearModel(Matrix x, Vector y, bool hasIntercept = true, LinkFunctionType linkType = LinkFunctionType.Identity) 
        {
            if (y.Length != x.NumberOfRows) throw new ArgumentException("X and Y must have the same number of rows.");
            if (y.Length <= 2) throw new ArithmeticException("There must be at least three data points.");
            if (x.NumberOfColumns > y.Length) throw new ArithmeticException($"A regression of the requested order requires at least {x.NumberOfColumns} data points. Only {y.Length} data points have been provided.");

            // Set inputs
            Y = y;
            X = hasIntercept ? AddInterceptColumn(x) : x;
            HasIntercept = hasIntercept;
            DegreesOfFreedom = X.NumberOfRows - X.NumberOfColumns;

            // Set model name
            if (Y.Header == null || Y.Header.Length == 0)
            {
                Y.Header = "Y Data";
            }

            // Set parameter names for summary report.
            ParameterNames = new List<string>();
            if (hasIntercept) ParameterNames.Add("Intercept");
            if (X.Header != null && X.Header.Length == x.NumberOfColumns)
            {
                ParameterNames.AddRange(X.Header);
            }
            else
            {
                for (int i = 1; i <= x.NumberOfColumns; i++)
                {
                    ParameterNames.Add("β" + i);
                }
            }

            LinkType = linkType;
            (_inverseLink, _inverseLinkDerivative, _logLikelihoodTerm, _logLikelihoodGradient) = GetLinkFunctions(linkType);
            SetOptimizer();

        }

        #endregion

        #region Members

        /// <summary>
        /// Determines if the linear model has an intercept. 
        /// </summary>
        public bool HasIntercept { get; private set; }

        /// <summary>
        /// The vector of response values.
        /// </summary>
        public Vector Y { get; private set; }

        /// <summary>
        /// The matrix of predictor values. 
        /// </summary>
        public Matrix X { get; private set; }

        /// <summary>
        /// The list of estimated parameter values.
        /// </summary>
        public double[] Parameters { get; private set; } = Array.Empty<double>();

        /// <summary>
        /// The list of the estimated parameter names. 
        /// </summary>
        public List<string> ParameterNames { get; private set; }

        /// <summary>
        /// The list of the estimated parameter standard errors. 
        /// </summary>
        public double[] ParameterStandardErrors { get; private set; } = Array.Empty<double>();

        /// <summary>
        /// The list of the estimated parameter z-scores.
        /// </summary>
        public double[] ParameterZScores { get; private set; } = Array.Empty<double>();

        /// <summary>
        /// The list of the estimated parameter p-values.
        /// </summary>
        public double[] ParameterPValues { get; private set; } = Array.Empty<double>();

        /// <summary>
        /// The estimate parameter covariance matrix. 
        /// </summary>
        public Matrix Covariance { get; private set; } = new Matrix(1, 1);

        /// <summary>
        /// The residuals of the fitted linear model. 
        /// </summary>
        public double[] Residuals { get; private set; } = Array.Empty<double>();

        /// <summary>
        /// The model standard error.
        /// </summary>
        public double StandardError { get; private set; }

        /// <summary>
        /// The data sample size. 
        /// </summary>
        public int SampleSize => Y.Length;

        /// <summary>
        /// The model degrees of freedom.
        /// </summary>
        public int DegreesOfFreedom { get; private set; }

        /// <summary>
        /// the Akaike Information Criterion (AIC). 
        /// </summary>
        public double AIC { get; private set; }

        /// <summary>
        /// the Akaike Information Criterion corrected for small sample sizes (AICc)
        /// </summary>
        public double AICc { get; private set; }

        /// <summary>
        /// The Bayesian information criterion (BIC). 
        /// </summary>
        public double BIC { get; private set; }

        /// <summary>
        /// Estimate robust standard errors.
        /// </summary>
        public bool UseRobustSE { get; set; } = false;

        /// <summary>
        /// Gets the optimizer used to train the model. Default = Nelder-Mead.
        /// </summary>
        public Optimizer Optimizer { get; private set; } = null!;

        /// <summary>
        /// Gets the link function type. 
        /// </summary>
        public LinkFunctionType LinkType { get; private set; }

        /// <summary>
        /// Enumeration of link function types.
        /// </summary>
        public enum LinkFunctionType
        {
            Identity,
            Log,
            Logit,
            Probit,
            ComplementaryLogLog
        }

        private Func<double, double> _inverseLink;
        private Func<double, double> _inverseLinkDerivative;
        private Func<(double mu, double y), double> _logLikelihoodTerm;
        private Func<(double eta, double y), double> _logLikelihoodGradient;

        #endregion

        #region Methods

        /// <summary>
        /// Helper method to add an intercept column to the covariate matrix.
        /// </summary>
        /// <param name="x">The matrix of predictor values.</param>
        private static Matrix AddInterceptColumn(Matrix x)
        {
            Matrix result = new Matrix(x.NumberOfRows, x.NumberOfColumns + 1);
            for (int i = 0; i < x.NumberOfRows; i++)
            {
                result[i, 0] = 1.0;
                for (int j = 0; j < x.NumberOfColumns; j++)
                    result[i, j + 1] = x[i, j];
            }
            return result;
        }

        /// <summary>
        /// Returns the necessary link functions.
        /// </summary>
        /// <param name="type">The link function type.</param>
        private static (Func<double, double> InverseLink,
                   Func<double, double> InverseLinkDerivative,
                   Func<(double mu, double y), double> LogLikelihoodTerm,
                   Func<(double eta, double y), double> LogLikelihoodGradientTerm)
            GetLinkFunctions(LinkFunctionType type)
        {
            switch (type)
            {
                // ------------------ Identity Link (Normal) ------------------
                case LinkFunctionType.Identity:
                    return (
                        eta => eta,
                        eta => 1.0,
                        (pair) =>
                        {
                            double resid = pair.y - pair.mu;
                            return -0.5 * resid * resid;
                        },
                        (pair) => (pair.y - pair.eta) * 1.0 // d/dη of SSE
                    );
                // ------------------ Log Link (Poisson) ------------------
                case LinkFunctionType.Log:
                    return (
                        eta => Math.Exp(eta),
                        eta => Math.Exp(eta),
                        (pair) =>
                        {
                            double mu = pair.mu;
                            return pair.y * Math.Log(mu) - mu;
                        },
                        (pair) =>
                        {
                            double mu = Math.Exp(pair.eta);
                            return (pair.y - mu) * mu; // mu * dμ/dη
                        }
                    );
                // ------------------ Logit Link (Binomial) ------------------
                case LinkFunctionType.Logit:
                    return (
                        eta => 1.0 / (1.0 + Math.Exp(-eta)),
                        eta =>
                        {
                            double ex = Math.Exp(-eta);
                            return ex / Math.Pow(1 + ex, 2);
                        },
                        (pair) =>
                        {
                            double mu = pair.mu;
                            return pair.y * Math.Log(mu) + (1 - pair.y) * Math.Log(1 - mu);
                        },
                        (pair) =>
                        {
                            double mu = 1.0 / (1.0 + Math.Exp(-pair.eta));
                            double dMu_dEta = mu * (1 - mu);
                            return (pair.y - mu) * dMu_dEta;
                        }
                    );
                // ------------------ Probit Link (Binomial) ------------------
                case LinkFunctionType.Probit:
                    return (
                        eta => Normal.StandardCDF(eta),
                        eta => Normal.StandardPDF(eta),
                        (pair) =>
                        {
                            double mu = pair.mu;
                            return pair.y * Math.Log(mu) + (1 - pair.y) * Math.Log(1 - mu);
                        },
                        (pair) =>
                        {
                            double mu = Normal.StandardCDF(pair.eta);
                            double dMu = Normal.StandardPDF(pair.eta);
                            return (pair.y - mu) * dMu;
                        }
                    );
                // ------------------ Complementary Log-Log Link (Binomial) ------------------
                case LinkFunctionType.ComplementaryLogLog:
                    return (
                        eta => 1.0 - Math.Exp(-Math.Exp(eta)),
                        eta => Math.Exp(eta - Math.Exp(eta)),
                        (pair) =>
                        {
                            double mu = pair.mu;
                            return pair.y * Math.Log(mu) + (1 - pair.y) * Math.Log(1 - mu);
                        },
                        (pair) =>
                        {
                            double mu = 1.0 - Math.Exp(-Math.Exp(pair.eta));
                            double dMu = Math.Exp(pair.eta - Math.Exp(pair.eta));
                            return (pair.y - mu) * dMu;
                        }
                    );

                default:
                    throw new ArgumentOutOfRangeException(nameof(type));
            }
        }


        /// <summary>
        /// Set up the local optimizer.
        /// </summary>
        /// <param name="method">The optimization method. Default = Nelder-Mead.</param>
        public void SetOptimizer(LocalMethod method = LocalMethod.NelderMead)
        {
            int n = X.NumberOfRows;
            int p = X.NumberOfColumns;

            double logLikelihood(double[] beta)
            {
                int n = X.NumberOfRows;
                double logLH = 0.0;
                for (int i = 0; i < n; i++)
                {
                    logLH += _logLikelihoodTerm((_inverseLink(X.GetRow(i).DotProduct(beta)), Y[i]));
                }
                if (double.IsNaN(logLH) || double.IsInfinity(logLH)) return double.MaxValue;
                return -logLH;
            }
            double[] gradient(double[] beta)
            {
                int n = X.NumberOfRows;
                int p = X.NumberOfColumns;
                Vector g = new Vector(p);
                for (int i = 0; i < n; i++)
                {
                    Vector xi = new Vector(X.Row(i));
                    double eta = xi.Array.DotProduct(beta);
                    g -= xi * _logLikelihoodGradient((eta, Y[i]));
                }
                return g.Array;
            }

            // Set the parameter constraints
            var initial = new double[p];
            var lower = new double[p];
            var upper = new double[p];

            if (LinkType == LinkFunctionType.Identity)
            {
                double min = Statistics.Minimum(Y.Array);
                double max = Statistics.Maximum(Y.Array);
                double range = Math.Max(max - min, 1e-6);
                double slope = range / Math.Max(n, 1);
                for (int i = 0; i < p; i++)
                {
                    initial[i] = 0.0;
                    lower[i] = -slope * 100;
                    upper[i] = slope * 100;
                }
                if (HasIntercept)
                {
                    initial[0] = (min + max) / 2.0;
                    lower[0] = initial[0] / 100;
                    upper[0] = initial[0] * 100;
                }
            }
            else if (LinkType == LinkFunctionType.Log)
            {
                var logY = Y.Array.Map(y => Math.Log(Math.Max(y, 1e-6)));
                double min = Statistics.Minimum(logY);
                double max = Statistics.Maximum(logY);
                double range = Math.Max(max - min, 1e-6);
                double slope = range / Math.Max(n, 1);
                for (int i = 0; i < p; i++)
                {
                    initial[i] = 0.0;
                    lower[i] = -slope * 100;
                    upper[i] = slope * 100;
                }
                if (HasIntercept)
                {
                    initial[0] = (min + max) / 2.0;
                    lower[0] = initial[0] / 100;
                    upper[0] = initial[0] * 100;
                }
            }
            else if (LinkType == LinkFunctionType.Logit)
            {
                lower.Fill(-10);
                upper.Fill(10);
                if (HasIntercept)
                {
                    double sum = Y.Array.Sum();
                    initial[0] = Math.Log((sum + 0.5) / (Y.Length - sum + 0.5)); // Log-odds of observed ratio
                }
            }
            else if (LinkType == LinkFunctionType.Probit)
            {
                lower.Fill(-6);
                upper.Fill(6);
                if (HasIntercept)
                {
                    double sum = Y.Array.Sum();
                    double rate = (sum + 0.5) / (Y.Length + 1.0); // Additive smoothing
                    initial[0] = Normal.StandardZ(rate); // Roughly centered on observed success rate
                }
            }
            else if (LinkType == LinkFunctionType.ComplementaryLogLog)
            {
                lower.Fill(-10);
                upper.Fill(10);
                if (HasIntercept)
                {
                    double prob = Statistics.Mean(Y.Array);
                    initial[0] = Math.Log(-Math.Log(1 - prob));

                }
            }

            // Set the optimizer
            switch (method)
            {
                case LocalMethod.ADAM:
                    Optimizer = new ADAM(logLikelihood, p, initial, lower, upper, gradient: gradient);
                    break;
                case LocalMethod.BFGS:
                    Optimizer = new BFGS(logLikelihood, p, initial, lower, upper, gradient: gradient);
                    break;
                case LocalMethod.GradientDescent:
                    Optimizer = new GradientDescent(logLikelihood, p, initial, lower, upper, gradient: gradient);
                    break;
                case LocalMethod.NelderMead:
                    Optimizer = new NelderMead(logLikelihood, p, initial, lower, upper);
                    break;
                case LocalMethod.Powell:
                    Optimizer = new Powell(logLikelihood, p, initial, lower, upper);
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(method));
            }

        }

        /// <summary>
        /// Train the Generalized Linear Model.
        /// </summary>
        public void Train()
        {
           if (Optimizer == null) throw new ArgumentNullException(nameof(Optimizer));
            Optimizer.Minimize();
            Parameters = Optimizer.BestParameterSet.Values;
            ComputeDiagnostics();
        }

        /// <summary>
        /// Compute GLM diagnostics
        /// </summary>
        private void ComputeDiagnostics()
        {
            int n = X.NumberOfRows;
            int p = X.NumberOfColumns;

            Residuals = new double[n];
            var mu = Predict(X);
            var diag = new Vector(n);
            double sse = 0.0;
            for (int i = 0; i < n; i++)
            {
                Residuals[i] = Y[i] - mu[i];
                diag[i] = Tools.Sqr(Residuals[i]);
                sse += diag[i];
            }
            StandardError = Math.Sqrt(sse / (n - p));

            // --- Compute log-likelihood with constant ---
            double logLik = 0.0;
            switch (LinkType)
            {
                case LinkFunctionType.Identity:
                    // Assume residual variance = SE^2
                    double sigma = StandardError;
                    logLik = -n * 0.5 * Math.Log(2 * Math.PI) - n * Math.Log(sigma) - sse / (2 * sigma * sigma);
                    break;

                case LinkFunctionType.Log:
                    for (int i = 0; i < n; i++)
                    {
                        double y = Y[i];
                        double m = mu[i];
                        logLik += y * Math.Log(m) - m - Factorial.LogFactorial((int)Math.Round(y));
                    }
                    break;

                case LinkFunctionType.Logit:
                case LinkFunctionType.Probit:
                case LinkFunctionType.ComplementaryLogLog:
                    for (int i = 0; i < n; i++)
                    {
                        double y = Y[i];
                        double m = mu[i];
                        logLik += y * Math.Log(m) + (1 - y) * Math.Log(1 - m);
                    }
                    break;

                default:
                    throw new InvalidOperationException("Unknown link function.");
            }

            AIC = GoodnessOfFit.AIC(p, logLik);
            AICc = GoodnessOfFit.AICc(n, p, logLik);
            BIC = GoodnessOfFit.BIC(n, p, logLik);

            // Jacobian matrix for the delta method
            Matrix J = new Matrix(n, p);
            for (int i = 0; i < n; i++)
            {
                var xi = X.GetRow(i);
                double eta = xi.DotProduct(Parameters);
                double dMu = _inverseLinkDerivative(eta);
                for (int j = 0; j < p; j++)
                    J[i, j] = dMu * xi[j];
            }

            Matrix JTJ = J.Transpose().Multiply(J);
            Matrix JTJinv = JTJ.Inverse();

            if (!UseRobustSE)
            {
                Covariance = JTJinv;
            }
            else
            {
                Matrix Omega = Matrix.Diagonal(diag);
                Matrix meat = J.Transpose().Multiply(Omega).Multiply(J);
                Covariance = JTJinv.Multiply(meat).Multiply(JTJinv);
            }

            ParameterStandardErrors = LinkType != LinkFunctionType.Log ? Covariance.Diagonal().Map(Math.Sqrt).Multiply(StandardError) : Covariance.Diagonal().Map(Math.Sqrt);
            ParameterZScores = new double[p];
            ParameterPValues = new double[p];
            for (int i = 0; i < p; i++)
            {
                double z = Parameters[i] / ParameterStandardErrors[i];
                ParameterZScores[i] = z;
                ParameterPValues[i] = 2 * (1 - Normal.StandardCDF(Math.Abs(z)));
            }
        }

        /// <summary>
        /// Provides a standard summary output table in a list of strings. 
        /// </summary>
        public List<string> Summary()
        {
            var text = new List<string>
            {
                "",
                $"Generalized Linear Model ({SampleSize} obs, {ParameterNames.Count} parameters):",
                $"Model for predicting {Y.Header}:",
                "Parameters:",
                $"{"",-15}{"Estimate",12}{"Std. Error",12}{"z value",12}{"Pr(>|z|)",12}"
            };

            for (int i = 0; i < Parameters.Length; i++)
            {
                string name = ParameterNames[i];
                double tval = ParameterZScores[i];
                double pval = ParameterPValues[i];
                string pvalStrg = pval > 1E-4 ? pval.ToString("N4") : pval < 1E-15 ? "< 1E-15" : pval.ToString("E2");
                string sig = pval < 1E-3 ? " ***" : pval < 1E-2 ? " **" : pval < 0.05 ? " *" : pval < 0.1 ? " ." : "  ";
                text.Add($"{name,-15}{Parameters[i],12:N5}{ParameterStandardErrors[i],12:N5}{tval,12:N3}{pvalStrg,12}{sig,-2}");
            }

            text.Add("---");
            text.Add("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1");
            text.Add("");
            text.Add($"Residual standard error: {StandardError:N4} on {DegreesOfFreedom} degrees of freedom");
            text.Add($"AIC: {AIC:N4}  AICc: {AICc:N4}  BIC: {BIC:N4}");
            text.Add("");

            // Residuals
            text.Add("Residuals:");
            text.Add($"{"Min",10}" + $"{"1Q",10}" + $"{"Median",10}" + $"{"3Q",10}" + $"{"Max",10}");
            var res = Statistics.FiveNumberSummary(Residuals);
            text.Add($"{res[0].ToString("N4"),10}" + $"{res[1].ToString("N4"),10}" + $"{res[2].ToString("N4"),10}" + $"{res[3].ToString("N4"),10}" + $"{res[4].ToString("N4"),10}");
            text.Add("");

            return text;
        }

        /// <summary>
        /// Returns the mean prediction. 
        /// </summary>
        /// <param name="x">The matrix of predictors.</param>
        public double[] Predict(Matrix x)
        {
            int n = x.NumberOfRows;
            var result = new double[n];
            for (int i = 0; i < n; i++)
                result[i] = _inverseLink(x.GetRow(i).DotProduct(Parameters.ToArray()));
            return result;
        }

        /// <summary>
        /// Returns the prediction with confidence intervals in a 2D array with columns: lower, mean, upper. 
        /// </summary>
        /// <param name="x">The matrix of predictors.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,] Predict(Matrix x, double alpha = 0.1)
        {
            var z = Normal.StandardZ(1 - alpha / 2);
            var result = new double[x.NumberOfRows, 3];
            for (int i = 0; i < x.NumberOfRows; i++)
            {
                double[] xi;
                if (HasIntercept)
                {
                    int p = x.NumberOfColumns;
                    xi = new double[x.NumberOfColumns + 1];
                    xi[0] = 1;
                    for (int j = 1; j < p; j++)
                    {
                        xi[i] = x[i, j];
                    }
                }
                else
                {
                    xi = x.GetRow(i);
                }

                double mu = _inverseLink(xi.DotProduct(Parameters.ToArray()));
                double se = Math.Sqrt(xi.DotProduct(Covariance.Multiply(new Vector(xi)).Array));
                result[i, 0] = mu - z * se;
                result[i, 1] = mu;
                result[i, 2] = mu + z * se;
            }
            return result;
        }


        #endregion

    }
}
