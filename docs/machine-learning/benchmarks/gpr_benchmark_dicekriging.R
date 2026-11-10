# Timing of DiceKriging::km for the GaussianProcessRegression benchmark.
#
# Uses the deterministic data of gpr_benchmark_sklearn.py (training rows 0..n-1, m = 1000 prediction
# rows n..n+m-1). Fit: km(~1, covtype = "matern5_2", nugget.estim = TRUE) with one BFGS start
# (multistart = 1). Predict: predict.km(type = "UK") mean and standard deviation at the m points.
# The median wall time over the repetitions is recorded for each size.
#
# Usage: Rscript gpr_benchmark_dicekriging.R <d> <out.json> n:reps [n:reps ...]
suppressMessages(library(DiceKriging))
args <- commandArgs(trailingOnly = TRUE)
d <- as.integer(args[1])
primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71)
make <- function(rows, offset = 0) {
  i <- (offset + 1):(offset + rows)
  v <- outer(i, sqrt(primes[1:d]))
  x <- v - floor(v)
  e <- i * sqrt(97)
  y <- as.vector((sin(2 * pi * x)) %*% (1 / (1:d))) + 0.1 * (e - floor(e) - 0.5)
  colnames(x) <- paste0("x", 1:d)
  list(x = data.frame(x), y = y)
}
sizes <- lapply(strsplit(args[-(1:2)], ":"), as.integer)
rows <- list()
for (s in sizes) {
  n <- s[1]; reps <- s[2]
  tr <- make(n); te <- make(1000, n)
  fit <- numeric(reps); pred <- numeric(reps); ll <- NA
  for (r in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    m <- km(~1, design = tr$x, response = tr$y, covtype = "matern5_2", nugget.estim = TRUE,
            multistart = 1, control = list(trace = FALSE))
    fit[r] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    p <- predict(m, newdata = te$x, type = "UK", checkNames = FALSE)
    pred[r] <- proc.time()[["elapsed"]] - t0
    ll <- logLik(m)
  }
  rows[[length(rows) + 1]] <- list(n = n, reps = reps, fitSeconds = median(fit), predictSeconds = median(pred),
                                   logLikelihood = ll)
  message(sprintf("n=%d fit=%.3f predict=%.3f", n, median(fit), median(pred)))
}
blas <- extSoftVersion()[["BLAS"]]
out <- list(library = "DiceKriging", DiceKriging = as.character(packageVersion("DiceKriging")),
            R = as.character(getRversion()), blas = if (nzchar(blas)) blas else "R reference BLAS",
            d = d, m = 1000, sizes = rows)
writeLines(jsonlite::toJSON(out, digits = NA, auto_unbox = TRUE, pretty = TRUE), args[2], useBytes = TRUE)
