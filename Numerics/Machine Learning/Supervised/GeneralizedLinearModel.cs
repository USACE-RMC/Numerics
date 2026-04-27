using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Functions;
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
            : this(x, y, LinkFunctionFactory.Create(linkType), hasIntercept, linkType)
        {
        }

        /// <summary>
        /// Constructs a generalized linear model with a custom link function.
        /// </summary>
        /// <param name="x">The matrix of predictor values.</param>
        /// <param name="y">The response vector.</param>
        /// <param name="linkFunction">The link function instance.</param>
        /// <param name="hasIntercept">Determines if an intercept should be estimated. Default = true.</param>
        /// <param name="linkType">The link function type for optimizer initialization. Default = identity.</param>
        public GeneralizedLinearModel(Matrix x, Vector y, ILinkFunction linkFunction, bool hasIntercept = true, LinkFunctionType linkType = LinkFunctionType.Identity)
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
            LinkFunction = linkFunction;
            _logLikelihoodTerm = GetFamilyLogLikelihood(linkType);
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
        public double[] Parameters { get; private set; } = null!;

        /// <summary>
        /// The list of the estimated parameter names.
        /// </summary>
        public List<string> ParameterNames { get; private set; } = null!;

        /// <summary>
        /// The list of the estimated parameter standard errors.
        /// </summary>
        public double[] ParameterStandardErrors { get; private set; } = null!;

        /// <summary>
        /// The list of the estimated parameter z-scores.
        /// </summary>
        public double[] ParameterZScores { get; private set; } = null!;

        /// <summary>
        /// The list of the estimated parameter p-values.
        /// </summary>
        public double[] ParameterPValues { get; private set; } = null!;

        /// <summary>
        /// The estimate parameter covariance matrix.
        /// </summary>
        public Matrix Covariance { get; private set; } = null!;

        /// <summary>
        /// The residuals of the fitted linear model.
        /// </summary>
        public double[] Residuals { get; private set; } = null!;

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
        /// Gets the link function instance used by this model.
        /// </summary>
        public ILinkFunction LinkFunction { get; private set; }

        /// <summary>
        /// Per-observation log-likelihood contribution as a function of (mu, y).
        /// This is distribution-family-specific (Normal, Poisson, Binomial).
        /// </summary>
        private Func<(double mu, double y), double> _logLikelihoodTerm;

        #endregion

        #region Methods

        /// <summary>
        /// Prepares the design matrix for prediction by adding an intercept column if needed.
        /// If the matrix already has the expected number of columns (e.g. internal X), it passes through.
        /// If it has one fewer column and HasIntercept is true, the intercept column is added.
        /// </summary>
        /// <param name="x">The matrix of predictor values.</param>
        private Matrix PrepareDesignMatrix(Matrix x)
        {
            int expected = Parameters.Length;
            if (x.NumberOfColumns == expected)
                return x;
            if (HasIntercept && x.NumberOfColumns == expected - 1)
                return AddInterceptColumn(x);
            throw new ArgumentException(
                $"Expected {expected} columns{(HasIntercept ? $" (or {expected - 1} without intercept)" : "")}, but got {x.NumberOfColumns}.",
                nameof(x));
        }

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
        /// Returns the distribution-family-specific per-observation log-likelihood function.
        /// </summary>
        /// <param name="type">The link function type, which determines the assumed distribution family.</param>
        /// <returns>A function that computes the per-observation log-likelihood given (mu, y).</returns>
        private static Func<(double mu, double y), double> GetFamilyLogLikelihood(LinkFunctionType type)
        {
            switch (type)
            {
                // Normal family (Identity link)
                case LinkFunctionType.Identity:
                    return (pair) =>
                    {
                        double resid = pair.y - pair.mu;
                        return -0.5 * resid * resid;
                    };

                // Poisson family (Log link)
                case LinkFunctionType.Log:
                    return (pair) =>
                    {
                        double mu = pair.mu;
                        return pair.y * Math.Log(mu) - mu;
                    };

                // Binomial family (Logit, Probit, or Complementary Log-Log link)
                case LinkFunctionType.Logit:
                case LinkFunctionType.Probit:
                case LinkFunctionType.ComplementaryLogLog:
                    return (pair) =>
                    {
                        double mu = pair.mu;
                        return pair.y * Math.Log(mu) + (1 - pair.y) * Math.Log(1 - mu);
                    };

                default:
                    throw new ArgumentOutOfRangeException(nameof(type));
            }
        }

        /// <summary>
        /// Computes the inverse link derivative d&#956;/d&#951; at the given link-space value.
        /// </summary>
        /// <param name="eta">The link-space value.</param>
        /// <returns>The derivative d&#956;/d&#951; = 1 / h&#8242;(&#956;).</returns>
        private double InverseLinkDerivative(double eta)
        {
            double mu = LinkFunction.InverseLink(eta);
            double dEta_dMu = LinkFunction.DLink(mu);
            return 1.0 / dEta_dMu;
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
                    double mu = LinkFunction.InverseLink(X.GetRow(i).DotProduct(beta));
                    logLH += _logLikelihoodTerm((mu, Y[i]));
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
                    double mu = LinkFunction.InverseLink(eta);
                    double dMu_dEta = InverseLinkDerivative(eta);
                    double gradTerm = (Y[i] - mu) * dMu_dEta;
                    g -= xi * gradTerm;
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
                    lower[0] = Math.Min(initial[0] / 100, initial[0] * 100);
                    upper[0] = Math.Max(initial[0] / 100, initial[0] * 100);
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
                    lower[0] = Math.Min(initial[0] / 100, initial[0] * 100);
                    upper[0] = Math.Max(initial[0] / 100, initial[0] * 100);
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
                double dMu = InverseLinkDerivative(eta);
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
            text.Add("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1");
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
            var xp = PrepareDesignMatrix(x);
            int n = xp.NumberOfRows;
            var result = new double[n];
            for (int i = 0; i < n; i++)
                result[i] = LinkFunction.InverseLink(xp.GetRow(i).DotProduct(Parameters.ToArray()));
            return result;
        }

        /// <summary>
        /// Returns the prediction with confidence intervals in a 2D array with columns: lower, mean, upper.
        /// </summary>
        /// <param name="x">The matrix of predictors.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,] Predict(Matrix x, double alpha = 0.1)
        {
            var xp = PrepareDesignMatrix(x);
            var z = Normal.StandardZ(1 - alpha / 2);
            var result = new double[xp.NumberOfRows, 3];
            for (int i = 0; i < xp.NumberOfRows; i++)
            {
                var xi = xp.GetRow(i);
                double mu = LinkFunction.InverseLink(xi.DotProduct(Parameters.ToArray()));
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
