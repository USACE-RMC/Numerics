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
    /// <para>
    ///     <b> Description: </b>
    /// </para>
    /// <para>
    ///     The link type selects the family: identity for Normal, log for Poisson, and logit, probit or
    ///     complementary log-log for Bernoulli responses. The maximum likelihood estimate is computed by
    ///     iteratively reweighted least squares (Fisher scoring), which also centers the optimizer's start
    ///     and bounds; the optimizer then refines it, and the fit with the higher likelihood is kept. When
    ///     Fisher scoring does not converge, for example under complete separation in a binary model, the
    ///     optimizer starts from a family-specific initial value inside fixed bounds.
    /// </para>
    /// <para>
    ///     The parameter covariance is the inverse Fisher information (X'WX)^-1 with weights
    ///     W = (d&#956;/d&#951;)^2 / V(&#956;), multiplied by the dispersion: the residual variance for the Normal
    ///     family and one for the Poisson and Bernoulli families. The robust option uses the HC0 sandwich
    ///     estimator.
    /// </para>
    /// <para>
    ///     <b> References: </b>
    /// </para>
    /// <para>
    ///     McCullagh, P. and Nelder, J. A. (1989). Generalized Linear Models, 2nd ed. Chapman and Hall.
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
            _rows = new double[X.NumberOfRows][];
            for (int i = 0; i < _rows.Length; i++)
                _rows[i] = X.Row(i);

            // Set model name
            if (Y.Header == null || Y.Header.Length == 0)
            {
                Y.Header = "Y Data";
            }

            // Set parameter names for summary report.
            ParameterNames = new List<string>();
            if (hasIntercept) ParameterNames.Add("Intercept");
            if (x.Header != null && x.Header.Length == x.NumberOfColumns)
            {
                ParameterNames.AddRange(x.Header);
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
        /// <remarks>
        /// "Intercept" comes first when the model has an intercept, followed by the predictor matrix header names
        /// when the header has one name per predictor column, or β1, β2, ... otherwise.
        /// </remarks>
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
        /// <remarks>
        /// AIC = 2k - 2 ln L with k the number of regression coefficients. For the Normal family (identity link),
        /// ln L is evaluated at the maximum likelihood variance σ² = SSE / n, which matches statsmodels GLM; the
        /// residual standard error <see cref="StandardError"/> keeps the divisor n - k. R's glm also counts the
        /// variance as a parameter, so its Gaussian AIC is larger by 2.
        /// </remarks>
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

        /// <summary>
        /// The rows of the design matrix <see cref="X"/>, read by every likelihood and gradient evaluation.
        /// </summary>
        private readonly double[][] _rows;

        /// <summary>
        /// The Fisher scoring estimate that centers the optimizer, or null when Fisher scoring did not converge.
        /// </summary>
        private double[]? _fisherScoringEstimate;

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
        /// <returns>
        /// The matrix with a leading column of ones. When <paramref name="x"/> has one header name per column, the
        /// header is "Intercept" followed by those names.
        /// </returns>
        private static Matrix AddInterceptColumn(Matrix x)
        {
            Matrix result = new Matrix(x.NumberOfRows, x.NumberOfColumns + 1);
            for (int i = 0; i < x.NumberOfRows; i++)
            {
                result[i, 0] = 1.0;
                for (int j = 0; j < x.NumberOfColumns; j++)
                    result[i, j + 1] = x[i, j];
            }
            if (x.Header != null && x.Header.Length == x.NumberOfColumns)
            {
                result.Header = new string[x.NumberOfColumns + 1];
                result.Header[0] = "Intercept";
                Array.Copy(x.Header, 0, result.Header, 1, x.NumberOfColumns);
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
        /// Returns the family variance function V(&#956;) for the model's link type.
        /// </summary>
        /// <param name="mu">The mean response.</param>
        /// <returns>1 for the Normal family, &#956; for Poisson and &#956;(1 - &#956;) for Bernoulli.</returns>
        private double FamilyVariance(double mu)
        {
            switch (LinkType)
            {
                case LinkFunctionType.Identity:
                    return 1.0;
                case LinkFunctionType.Log:
                    return mu;
                default:
                    return mu * (1.0 - mu);
            }
        }

        /// <summary>
        /// Returns the Fisher information X'WX at the given parameters.
        /// </summary>
        /// <param name="beta">The parameter vector.</param>
        /// <returns>The p by p information matrix.</returns>
        private Matrix FisherInformation(double[] beta)
        {
            int n = X.NumberOfRows, p = X.NumberOfColumns;
            var info = new Matrix(p, p);
            for (int i = 0; i < n; i++)
            {
                var xi = _rows[i];
                double eta = xi.DotProduct(beta);
                double dMu = InverseLinkDerivative(eta);
                double w = dMu * dMu / FamilyVariance(LinkFunction.InverseLink(eta));
                for (int a = 0; a < p; a++)
                    for (int b = 0; b < p; b++)
                        info[a, b] += w * xi[a] * xi[b];
            }
            return info;
        }

        /// <summary>
        /// Computes the maximum likelihood estimate by iteratively reweighted least squares.
        /// </summary>
        /// <param name="beta">The estimate when the iteration converges.</param>
        /// <param name="standardErrors">The square roots of the diagonal of the inverse Fisher information at the estimate.</param>
        /// <returns>True when the iteration converged to finite values.</returns>
        /// <remarks>
        /// Starts from the mean values used by R's glm (y for Normal, y + 0.1 for Poisson and (y + 0.5) / 2
        /// for Bernoulli) and stops when the relative change in deviance falls below 1E-12, or after 100
        /// iterations. Each step solves the weighted normal equations (X'WX) b = X'Wz with working
        /// response z = &#951; + (y - &#956;) d&#951;/d&#956;.
        /// </remarks>
        private bool TryFisherScoring(out double[] beta, out double[] standardErrors)
        {
            int n = X.NumberOfRows, p = X.NumberOfColumns;
            beta = new double[p];
            standardErrors = new double[p];
            var eta = new double[n];
            for (int i = 0; i < n; i++)
            {
                double mu0 = LinkType == LinkFunctionType.Identity ? Y[i] :
                             LinkType == LinkFunctionType.Log ? Y[i] + 0.1 : (Y[i] + 0.5) / 2.0;
                eta[i] = LinkFunction.Link(mu0);
            }
            double deviance = double.MaxValue;
            bool converged = false;
            try
            {
                for (int iteration = 0; iteration < 100; iteration++)
                {
                    var xtwx = new Matrix(p, p);
                    var xtwz = new Vector(p);
                    for (int i = 0; i < n; i++)
                    {
                        double mu = LinkFunction.InverseLink(eta[i]);
                        double dMu = InverseLinkDerivative(eta[i]);
                        double w = dMu * dMu / FamilyVariance(mu);
                        double z = eta[i] + (Y[i] - mu) / dMu;
                        for (int a = 0; a < p; a++)
                        {
                            double wx = w * X[i, a];
                            xtwz[a] += wx * z;
                            for (int b = 0; b < p; b++)
                                xtwx[a, b] += wx * X[i, b];
                        }
                    }
                    var next = new CholeskyDecomposition(xtwx).Solve(xtwz).Array;
                    for (int i = 0; i < n; i++)
                        eta[i] = _rows[i].DotProduct(next);
                    double newDeviance = -2.0 * LogLikelihoodKernel(next);
                    beta = next;
                    if (!Tools.IsFinite(newDeviance)) return false;
                    if (Math.Abs(newDeviance - deviance) <= 1E-12 * (Math.Abs(newDeviance) + 0.1))
                    {
                        converged = true;
                        break;
                    }
                    deviance = newDeviance;
                }
                if (!converged || beta.Any(b => !Tools.IsFinite(b))) return false;
                var covariance = FisherInformation(beta).Inverse();
                for (int a = 0; a < p; a++)
                    standardErrors[a] = Math.Sqrt(covariance[a, a]);
                return standardErrors.All(Tools.IsFinite);
            }
            catch (Exception)
            {
                // A singular weighted cross-product (the Cholesky factorization throws) or a saturated
                // link ends Fisher scoring; the optimizer then starts from the family-specific values.
                return false;
            }
        }

        /// <summary>
        /// Returns the log-likelihood kernel (without constants) at the given parameters.
        /// </summary>
        /// <param name="beta">The parameter vector.</param>
        /// <returns>The sum of the per-observation log-likelihood terms.</returns>
        private double LogLikelihoodKernel(double[] beta)
        {
            double logLH = 0.0;
            for (int i = 0; i < X.NumberOfRows; i++)
            {
                double mu = LinkFunction.InverseLink(_rows[i].DotProduct(beta));
                logLH += _logLikelihoodTerm((mu, Y[i]));
            }
            return logLH;
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
                double logLH = LogLikelihoodKernel(beta);
                if (!Tools.IsFinite(logLH)) return double.MaxValue;
                return -logLH;
            }
            double[] gradient(double[] beta)
            {
                int n = X.NumberOfRows;
                int p = X.NumberOfColumns;
                Vector g = new Vector(p);
                for (int i = 0; i < n; i++)
                {
                    var xi = _rows[i];
                    double eta = xi.DotProduct(beta);
                    double mu = LinkFunction.InverseLink(eta);
                    double dMu_dEta = InverseLinkDerivative(eta);
                    double gradTerm = (Y[i] - mu) * dMu_dEta;
                    for (int j = 0; j < p; j++)
                        g[j] = g[j] - xi[j] * gradTerm;
                }
                return g.Array;
            }

            // Set the parameter constraints
            var initial = new double[p];
            var lower = new double[p];
            var upper = new double[p];

            _fisherScoringEstimate = null;
            if (TryFisherScoring(out var start, out var startSE))
            {
                // Center the search on the maximum likelihood estimate with bounds far wider than its
                // sampling uncertainty, so the bounds never bind at the optimum.
                _fisherScoringEstimate = start;
                for (int i = 0; i < p; i++)
                {
                    double width = Math.Max(1.0, 10.0 * Math.Abs(start[i])) + 100.0 * startSE[i];
                    initial[i] = start[i];
                    lower[i] = start[i] - width;
                    upper[i] = start[i] + width;
                }
            }
            else if (LinkType == LinkFunctionType.Identity)
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
        /// <remarks>
        /// Runs the optimizer and keeps its result unless the Fisher scoring estimate has a higher
        /// log-likelihood.
        /// </remarks>
        public void Train()
        {
           if (Optimizer == null) throw new ArgumentNullException(nameof(Optimizer));
            Optimizer.Minimize();
            Parameters = Optimizer.BestParameterSet.Values;
            if (_fisherScoringEstimate != null && LogLikelihoodKernel(_fisherScoringEstimate) >= LogLikelihoodKernel(Parameters))
                Parameters = (double[])_fisherScoringEstimate.Clone();
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
                    // Normal log-likelihood at the maximum likelihood variance sse / n, as statsmodels GLM uses
                    double sigma = Math.Sqrt(sse / n);
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

            // Inverse Fisher information (X'WX)^-1 with W = (dmu/deta)^2 / V(mu)
            Matrix bread = FisherInformation(Parameters).Inverse();
            if (!UseRobustSE)
            {
                // Dispersion: the residual variance for the Normal family, one otherwise
                double dispersion = LinkType == LinkFunctionType.Identity ? StandardError * StandardError : 1.0;
                Covariance = bread.Multiply(dispersion);
            }
            else
            {
                // HC0 sandwich: bread * sum(x x' u^2) * bread with score u = (y - mu) (dmu/deta) / V(mu)
                Matrix meat = new Matrix(p, p);
                for (int i = 0; i < n; i++)
                {
                    var xi = _rows[i];
                    double eta = xi.DotProduct(Parameters);
                    double u = Residuals[i] * InverseLinkDerivative(eta) / FamilyVariance(mu[i]);
                    for (int a = 0; a < p; a++)
                        for (int b = 0; b < p; b++)
                            meat[a, b] += u * u * xi[a] * xi[b];
                }
                Covariance = bread.Multiply(meat).Multiply(bread);
            }

            ParameterStandardErrors = Covariance.Diagonal().Map(Math.Sqrt);
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
        /// Returns the mean prediction with a confidence interval for the mean response, in a 2D array
        /// with columns: lower, mean, upper.
        /// </summary>
        /// <param name="x">The matrix of predictors.</param>
        /// <param name="alpha">The significance level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <returns>The lower bound, mean and upper bound for each row of <paramref name="x"/>.</returns>
        /// <remarks>
        /// The interval is formed on the link scale, &#951; &#177; z se(&#951;) with se(&#951;)^2 = x'Cx, and mapped
        /// through the inverse link, so it stays inside the response range. It describes uncertainty in
        /// the mean response, not the spread of new observations.
        /// </remarks>
        public double[,] Predict(Matrix x, double alpha = 0.1)
        {
            var xp = PrepareDesignMatrix(x);
            var z = Normal.StandardZ(1 - alpha / 2);
            var result = new double[xp.NumberOfRows, 3];
            for (int i = 0; i < xp.NumberOfRows; i++)
            {
                var xi = xp.GetRow(i);
                double eta = xi.DotProduct(Parameters.ToArray());
                double se = Math.Sqrt(xi.DotProduct(Covariance.Multiply(new Vector(xi)).Array));
                double a = LinkFunction.InverseLink(eta - z * se);
                double b = LinkFunction.InverseLink(eta + z * se);
                result[i, 0] = Math.Min(a, b);
                result[i, 1] = LinkFunction.InverseLink(eta);
                result[i, 2] = Math.Max(a, b);
            }
            return result;
        }


        #endregion

    }
}
