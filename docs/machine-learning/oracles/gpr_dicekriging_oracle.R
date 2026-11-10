# DiceKriging oracle values for GaussianProcessRegression with a constant mean.
#
# Uses the training and test points of gpr-sklearn.json (passed as the first argument) and writes
# gpr-dicekriging.json with, for a constant trend estimated by generalized least squares:
#   * covtype "gauss" in two dimensions with fixed coef.cov = (0.7, 1.4), coef.var = 1.3 and
#     nugget = 0.05: the trend coefficient, the log likelihood, predict.km(type = "UK") mean and sd,
#     and leaveOneOut.km(type = "UK", trend.reestim = TRUE) mean and sd;
#   * covtype "matern5_2" and "exp" in one dimension (first feature) with coef.cov = 0.7: the same.
#     DiceKriging forms multidimensional Matern and exponential kernels as products of one-dimensional
#     kernels, so only the Gaussian kernel matches the Euclidean form in more than one dimension.
#   * a maximum likelihood fit of the Gaussian kernel with an estimated nugget: the log likelihood.
#
# Usage: Rscript gpr_dicekriging_oracle.R gpr-sklearn.json gpr-dicekriging.json
suppressMessages(library(DiceKriging))
# jsonlite reads the fixture.
args <- commandArgs(trailingOnly = TRUE)
fixture <- jsonlite::fromJSON(args[1])
x <- fixture$x; y <- fixture$y; xt <- fixture$xTest

fmt <- function(v) paste0("[", paste(sprintf("%.17g", v), collapse = ", "), "]")
case <- function(design, test, covtype, theta) {
  colnames(design) <- paste0("x", seq_len(ncol(design)))
  test <- data.frame(test); colnames(test) <- colnames(design)
  m <- km(~1, design = data.frame(design), response = y, covtype = covtype,
          coef.cov = theta, coef.var = 1.3, nugget = 0.05)
  p <- predict(m, newdata = test, type = "UK", checkNames = FALSE)
  loo <- leaveOneOut.km(m, type = "UK", trend.reestim = TRUE)
  # logLik.km needs an estimation method; with every covariance parameter fixed, the Gaussian log
  # likelihood is evaluated from the model's stored factor T (C = T'T) and whitened residuals z.
  ll <- -0.5 * length(y) * log(2 * pi) - sum(log(diag(m@T))) - 0.5 * sum(m@z^2)
  sprintf('{"covtype": "%s", "beta": %s, "logLikelihood": %.17g, "mean": %s, "sd": %s, "looMean": %s, "looSd": %s}',
          covtype, fmt(m@trend.coef), ll, fmt(p$mean), fmt(p$sd), fmt(loo$mean), fmt(loo$sd))
}
cases <- c(case(x, xt, "gauss", c(0.7, 1.4)),
           case(x[, 1, drop = FALSE], xt[, 1, drop = FALSE], "matern5_2", 0.7),
           case(x[, 1, drop = FALSE], xt[, 1, drop = FALSE], "exp", 0.7))
d <- data.frame(x); colnames(d) <- c("x1", "x2")
set.seed(20261105)
fit <- km(~1, design = d, response = y, covtype = "gauss", nugget.estim = TRUE, multistart = 20,
          control = list(trace = FALSE))
out <- sprintf('{\n "DiceKriging": "%s",\n "R": "%s",\n "fixed": [\n  %s\n ],\n "fitted": {"covtype": "gauss", "logLikelihood": %.17g, "theta": %s, "variance": %.17g, "nugget": %.17g}\n}\n',
               packageVersion("DiceKriging"), getRversion(), paste(cases, collapse = ",\n  "),
               logLik(fit), fmt(fit@covariance@range.val), fit@covariance@sd2, fit@covariance@nugget)
con <- file(args[2], "wb"); writeLines(out, con, sep = "\n"); close(con)
