# Independent expected-information oracle for Hosking GNO/GLO, Kappa Four and GEV.
# No Numerics assembly, production numerical helper, fitted sample, or package is used.
# Run from the repository root with Rscript --vanilla docs/distributions/oracles/generalized-fisher.R.
# Public parameter order is xi, alpha, kappa[, hondo]. GEV fixes hondo=0;
# GLO fixes hondo=-1. Covariance in a restricted family is the inverse of its
# principal INFORMATION block, not a principal block of the larger covariance.
options(digits = 17, warn = 2)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
directory <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
target <- file.path(directory, "generalized-fisher.csv")
notes <- file.path(directory, "generalized-fisher.md")
relative_tolerance <- 2e-12
absolute_tolerance <- 2e-13

# expm1(x)-x evaluated from its convergent series near the removable singularity.
# This independently implements the derivative of the defining density, not a
# finite difference of a quantile with its observation moving with the parameter.
expm1_minus_x <- function(x) {
  out <- expm1(x) - x
  small <- abs(x) < 0.05
  if (any(small)) {
    y <- x[small]
    term <- y * y / 2
    total <- term
    for (j in 3:16) { term <- term * y / j; total <- total + term }
    out[small] <- total
  }
  out
}

# Scores at a fixed observation, written in terms of its true latent variate.
normal_scores <- function(z, k) {
  A <- k - z
  r <- exp(k * z)
  d1 <- if (k == 0) z else expm1(k * z) / k
  d2 <- if (k == 0) z * z / 2 else expm1_minus_x(k * z) / k^2
  cbind(-A * r, -1 - A * d1, z + A * d2)
}

# Each row below is score * sqrt(phi(z)); combining the Gaussian exponent first
# avoids overflow before multiplication by the integration weight in remote tails.
normal_weighted_scores <- function(z, k) {
  A <- k - z
  root <- exp(-z * z / 4) / (2 * pi)^0.25
  rroot <- exp(k * z - z * z / 4) / (2 * pi)^0.25
  if (k == 0) {
    d1root <- z * root
    d2root <- z * z * root / 2
  } else {
    d1root <- (rroot - root) / k
    d2root <- (rroot - root - k * z * root) / k^2
    small <- abs(k * z) < 0.05
    d1root[small] <- expm1(k * z[small]) / k * root[small]
    d2root[small] <- expm1_minus_x(k * z[small]) / k^2 * root[small]
  }
  cbind(-A * rroot, -root - A * d1root, z * root + A * d2root)
}

kappa_scores <- function(logp, logq, k, h) {
  if (h == 0) {
    B <- -logp
    w <- log(-logp)
  } else {
    B <- expm1(-h * logp) / h
    w <- numeric(length(logp))
    large <- h * logp > 36
    if (any(large)) w[large] <- h * logp[large] + log1p(-exp(-h * logp[large])) - log(-h)
    w[!large] <- log(-expm1(h * logp[!large]) / h)
  }
  # If the survival probability is below rounding resolution, log(t)=log(q)
  # to binary64 relative precision. Retain its logarithm even if q underflows.
  upper <- logq < -36
  w[upper] <- logq[upper]
  A <- 1 - k - (1 - h) * B
  R <- exp(-k * w)
  d1 <- if (k == 0) -w else expm1(-k * w) / k
  d2 <- if (k == 0) -w * w / 2 else -expm1_minus_x(-k * w) / k^2
  if (h == 0) {
    sh <- -logp - logp * logp / 2
  } else {
    # Rearranged to preserve the exact nonzero shape when h is small.
    sh <- B - expm1_minus_x(-h * logp) / h^2
  }
  cbind(A * R, -1 + A * d1, -w + A * d2, sh)
}

# Score times sqrt(Jacobian), combining exponents before exponentiation. This
# matters for large negative h: t itself can overflow while every weighted
# information integrand remains finite. No domain tail is omitted.
kappa_weighted_scores <- function(logp, logq, k, h, logroot) {
  if (h == 0) {
    w <- log(-logp)
    logB <- w
  } else {
    log_abs_expm1 <- function(v) {
      answer <- numeric(length(v))
      large <- v > 36
      answer[large] <- v[large] + log1p(-exp(-v[large]))
      answer[!large] <- log(abs(expm1(v[!large])))
      answer
    }
    w <- log_abs_expm1(h * logp) - log(abs(h))
    logB <- log_abs_expm1(-h * logp) - log(abs(h))
  }
  upper <- logq < -36
  w[upper] <- logq[upper]
  logB[upper] <- logq[upper]
  root <- exp(logroot)
  Broot <- exp(logB + logroot)
  Aroot <- (1 - k) * root - (1 - h) * Broot
  ARroot <- (1 - k) * exp(-k * w + logroot) - (1 - h) * exp(logB - k * w + logroot)
  if (k == 0) {
    scale_score <- -root - w * Aroot
    shape_score <- -w * root - w * w * Aroot / 2
  } else {
    scale_score <- -root + (ARroot - Aroot) / k
    shape_score <- -w * root + ((1 - k * w) * Aroot - ARroot) / k^2
    small <- abs(k * w) < .05
    scale_score[small] <- -root[small] + Aroot[small] * expm1(-k * w[small]) / k
    shape_score[small] <- -w[small] * root[small] - Aroot[small] * expm1_minus_x(-k * w[small]) / k^2
  }
  if (h == 0) {
    hscore <- (-logp - logp * logp / 2) * root
  } else {
    hscore <- (-logp * root - (1 - h) * Broot) / h
    small <- abs(h * logp) < .05
    hscore[small] <- Broot[small] - expm1_minus_x(-h * logp[small]) / h^2 * root[small]
  }
  cbind(ARroot, scale_score, shape_score, hscore)
}

fixed_log_density <- function(x, parameters, family) {
  xi <- parameters[1]; alpha <- parameters[2]; k <- parameters[3]
  y <- (x - xi) / alpha
  if (alpha <= 0 || (k != 0 && 1 - k * y <= 0)) return(-Inf)
  w <- if (k == 0) -y else log1p(-k * y) / k
  if (family == "GNO") {
    z <- -w
    return(-log(alpha) + k * z + dnorm(z, log = TRUE))
  }
  h <- if (family == "GLO") -1 else if (family == "GEV") 0 else parameters[4]
  t <- exp(w)
  if (h != 0 && 1 - h * t <= 0) return(-Inf)
  logF <- if (h == 0) -t else log1p(-h * t) / h
  -log(alpha) + (1 - k) * w + (1 - h) * logF
}

quantile <- function(p, parameters, family) {
  xi <- parameters[1]; alpha <- parameters[2]; k <- parameters[3]
  if (family == "GNO") z <- qnorm(p) else {
    h <- if (family == "GLO") -1 else if (family == "GEV") 0 else parameters[4]
    t <- if (h == 0) -log(p) else -expm1(h * log(p)) / h
    z <- -log(t)
  }
  xi + alpha * if (k == 0) z else -expm1(-k * z) / k
}

# Five-point fixed-x differences, at two step sizes; observations are never
# recomputed inside the differentiation stencil.
derivative_check <- function(family, k, h) {
  parameters <- c(0.7, 1.3, k, if (family == "K4") h)
  worst <- 0
  agreement <- 0
  for (p in c(.1, .5, .9)) {
    x <- quantile(p, parameters, family)
    expected <- if (family == "GNO") normal_scores(qnorm(p), k) else {
      kappa_scores(log(p), log1p(-p), k, h)[, seq_along(parameters), drop = FALSE]
    }
    expected[1:2] <- expected[1:2] / parameters[2]
    for (j in seq_along(parameters)) {
      difference <- function(step) {
        values <- vapply(c(-2, -1, 1, 2), function(multiple) {
          candidate <- parameters
          candidate[j] <- candidate[j] + multiple * step
          fixed_log_density(x, candidate, family)
        }, 0.0)
        (values[1] - 8 * values[2] + 8 * values[3] - values[4]) / (12 * step)
      }
      step <- 1e-4 * (1 + abs(parameters[j]))
      coarse <- difference(step)
      fine <- difference(step / 2)
      worst <- max(worst, abs(fine - expected[j]) / max(1, abs(expected[j])))
      agreement <- max(agreement, abs(fine - coarse) / max(1, abs(fine)))
    }
  }
  stopifnot(is.finite(worst), worst < 2e-7)
  c(derivative_relative_error = worst, derivative_step_agreement = agreement)
}

integral <- function(f, lower, upper) {
  answer <- integrate(f, lower, upper, subdivisions = 1500L,
                      rel.tol = relative_tolerance, abs.tol = absolute_tolerance,
                      stop.on.error = TRUE)
  stopifnot(is.finite(answer$value), is.finite(answer$abs.error))
  c(value = answer$value, error = answer$abs.error)
}

# Alternative full-domain parametrizations: p=u^power/2 and q=u^power/2.
# The production algorithm is neither loaded nor called. Different powers give
# independent error-sensitive coordinates; GNO additionally uses direct z-space.
probability_integral <- function(family, k, h, i, j, power) {
  halves <- lapply(c(FALSE, TRUE), function(upper_tail) {
    integral(function(u) {
      logtail <- power * log(u) - log(2)
      logother <- log1p(-exp(logtail))
      lp <- if (upper_tail) logother else logtail
      lq <- if (upper_tail) logtail else logother
      logjac <- log(power / 2) + (power - 1) * log(u)
      weighted <- if (family == "GNO") {
        z <- if (upper_tail) qnorm(logtail, lower.tail = FALSE, log.p = TRUE) else qnorm(logtail, log.p = TRUE)
        normal_scores(z, k) * exp(logjac / 2)
      } else kappa_weighted_scores(lp, lq, k, h, logjac / 2)
      if (j == 0) weighted[, i] * exp(logjac / 2) else weighted[, i] * weighted[, j]
    }, 0, 1)
  })
  Reduce(`+`, halves)
}

normal_integral <- function(k, i, j) {
  integral(function(z) {
    scores <- normal_weighted_scores(z, k)
    if (j == 0) scores[, i] * exp(-z * z / 4) / (2 * pi)^.25 else scores[, i] * scores[, j]
  }, -Inf, Inf)
}

cases <- list()
add <- function(family, k, h) cases[[length(cases) + 1L]] <<- list(family = family, k = k, h = h)
for (k in c(0, -1e-6, 1e-6, -.2, .2, -1, 1, -2, 2)) add("GNO", k, NA_real_)
for (k in c(0, -1e-6, 1e-6, -.2, .2, -.4, .4, -.45, .45, -.49, .49)) add("GLO", k, -1)
for (k in c(0, -1e-6, 1e-6, -.2, .2, -1, .45, .49)) add("GEV", k, 0)
for (pair in list(c(0,0), c(1e-6,0), c(-1e-6,0), c(0,1e-6), c(0,-1e-6), c(0,-1), c(.2,-1), c(-.2,-1),
                  c(.2,.2), c(-.2,-.2), c(.4,-1), c(-.4,-1), c(-1,-.2),
                  c(-.2,-2), c(.49,.2), c(.2,.49), c(-.49,-1), c(-.2,-2.45))) add("K4", pair[1], pair[2])

rows <- list(); summaries <- list()
append_row <- function(family, k, h, quantity, i, j, value, error, agreement) {
  rows[[length(rows) + 1L]] <<- data.frame(family = family, kappa = k, hondo = h,
    location = 0, scale = 1, sample_size = 1, quantity = quantity, row = i,
    column = j, value = value, absolute_error_estimate = error,
    alternate_coordinate_difference = agreement)
}

for (case in cases) {
  family <- case$family; k <- case$k; h <- case$h
  dimension <- if (family == "K4") 4L else 3L
  check <- derivative_check(family, k, h)
  information <- matrix(0, dimension, dimension)
  errors <- information; alternatives <- information
  score_mean <- score_errors <- numeric(dimension)
  # The high-power map makes approach-to-boundary cases well conditioned at u=0.
  near_boundary <- family != "GNO" && max(k, h, k*h) > .4
  first_power <- if (near_boundary) 128 else 16
  second_power <- if (near_boundary) 192 else 24
  for (i in seq_len(dimension)) {
    mean <- if (family == "GNO") normal_integral(k, i, 0) else probability_integral(family, k, h, i, 0, first_power)
    mean2 <- probability_integral(family, k, h, i, 0, second_power)
    score_mean[i] <- mean["value"]; score_errors[i] <- mean["error"]
    append_row(family, k, h, "score_mean", i, 0, mean["value"], mean["error"], abs(mean["value"] - mean2["value"]))
    for (j in i:dimension) {
      a <- if (family == "GNO") normal_integral(k, i, j) else probability_integral(family, k, h, i, j, first_power)
      b <- probability_integral(family, k, h, i, j, second_power)
      information[i,j] <- information[j,i] <- a["value"]
      errors[i,j] <- errors[j,i] <- a["error"]
      alternatives[i,j] <- alternatives[j,i] <- b["value"]
    }
  }
  eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
  covariance <- solve(information)
  alternate_covariance <- solve(alternatives)
  delta <- norm(information - alternatives, "2")
  error_norm <- norm(errors, "2")
  covariance_error <- norm(covariance, "2")^2 * error_norm / (1 - norm(covariance, "2") * error_norm)
  residual <- max(abs(information %*% covariance - diag(dimension)))
  stopifnot(min(eigenvalues) > 0, is.finite(covariance_error), covariance_error > 0,
    norm(covariance, "2") * error_norm < 1, max(abs(score_mean)) < 5e-8,
    delta < 2e-7 * max(1, norm(information, "2")), residual < 1e-8)
  for (i in seq_len(dimension)) for (j in seq_len(dimension)) {
    append_row(family, k, h, "information", i, j, information[i,j], errors[i,j], abs(information[i,j] - alternatives[i,j]))
    append_row(family, k, h, "covariance", i, j, covariance[i,j], covariance_error, abs(covariance[i,j] - alternate_covariance[i,j]))
  }
  if (family == "GNO" && k == 0) {
    exact <- matrix(c(7/6,0,1/3, 0,1/2,0, 1/3,0,2/3), 3, 3)
    stopifnot(max(abs(covariance - exact)) < 1e-10)
    for (i in 1:3) for (j in 1:3) append_row(family, k, h, "exact_covariance", i, j, exact[i,j], 0, abs(covariance[i,j] - exact[i,j]))
  }
  summaries[[length(summaries) + 1L]] <- data.frame(family = family, k = k, h = h,
    min_eigenvalue = min(eigenvalues), condition = max(eigenvalues) / min(eigenvalues),
    max_mean = max(abs(score_mean)), error_norm = error_norm, coordinate_difference = delta,
    covariance_error = covariance_error, inverse_residual = residual,
    derivative_error = check[1], derivative_agreement = check[2])
  cat(family, "k=", k, "h=", h, "minimum eigenvalue=", min(eigenvalues), "coordinate delta=", delta, "\n")
}

# Additional analytical moment references. GNO is a shifted/reflected lognormal.
# Near zero, GLO uses power-series algebra performed on coefficient arrays before
# evaluation, rather than subtracting nearly equal floating-point Gamma products.
normal_moments <- function(k) {
  if (k == 0) return(c(0, 1, 0, 3))
  t <- k*k; v <- expm1(t); relative <- expm1(t)/t
  c(-expm1(t/2)/k, exp(t/2)*sqrt(relative),
    -k*(v+3)*sqrt(relative), 3+v*(16+v*(15+v*(6+v))))
}

logistic_moments <- function(k) {
  if (k == 0) return(c(0, pi/sqrt(3), 0, 21/5))
  if (abs(k) > .01) {
    r <- 1:4; B <- pi*r*k/sin(pi*r*k)
    V <- B[2]-B[1]^2
    return(c((1-B[1])/k, sqrt(V)/abs(k),
      sign(k)*(-B[3]+3*B[1]*B[2]-2*B[1]^3)/V^1.5,
      (B[4]-4*B[1]*B[3]+6*B[1]^2*B[2]-3*B[1]^4)/V^2))
  }
  degree <- 14L
  reciprocal_sinc <- numeric(degree+1L); reciprocal_sinc[1] <- 1
  sinc <- (-1)^(0:degree)*pi^(2*(0:degree))/factorial(2*(0:degree)+1)
  for (n in 1:degree) reciprocal_sinc[n+1] <- -sum(sinc[2:(n+1)]*rev(reciprocal_sinc[1:n]))
  multiply <- function(a,b) {
    answer <- numeric(degree+1L)
    for (n in 0:degree) answer[n+1] <- sum(a[1:(n+1)]*rev(b[1:(n+1)]))
    answer
  }
  B <- lapply(1:4, function(r) reciprocal_sinc*r^(2*(0:degree)))
  B11 <- multiply(B[[1]],B[[1]])
  V <- B[[2]]-B11
  N3 <- -B[[3]]+3*multiply(B[[1]],B[[2]])-2*multiply(B11,B[[1]])
  N4 <- B[[4]]-4*multiply(B[[1]],B[[3]])+6*multiply(B11,B[[2]])-3*multiply(B11,B11)
  horner <- function(a,t) Reduce(function(value,coefficient) value*t+coefficient, rev(a), init=0)
  t <- k*k
  # Exact orders of vanishing are removed symbolically, not inferred by tolerance.
  variance <- horner(V[-1],t)
  c(-k*horner(reciprocal_sinc[-1],t), sqrt(variance),
    k*horner(N3[-c(1,2)],t)/variance^1.5,
    horner(N4[-c(1,2)],t)/variance^2)
}

stopifnot(abs(logistic_moments(.00010001)[4]-4.2000018682670925) < 3e-13,
          abs(logistic_moments(.0002)[3]+.0017412478719479963) < 3e-15)
for (family in c("GNO","GLO")) for (k in c(0,-1e-6,1e-6,-.00010001,.00010001,-.0002,.0002,-.01,.01,-.2,.2)) {
  values <- if (family == "GNO") normal_moments(k) else logistic_moments(k)
  for (i in 1:4) append_row(family,k,if (family=="GLO") -1 else NA_real_,"moment",i,0,
    values[i],5e-13*max(1,abs(values[i])),NA_real_)
}

result <- do.call(rbind, rows)
formatted <- result
for (column in names(formatted)) if (is.numeric(formatted[[column]])) {
  formatted[[column]] <- ifelse(is.na(formatted[[column]]),"NA",sprintf("%.17g",formatted[[column]]))
}
write.table(formatted, target, sep = ",", quote = FALSE, row.names = FALSE, na = "NA", eol = "\n")
summary <- do.call(rbind, summaries)
header <- c("# Independent generalized-family Fisher information oracle", "",
  paste("Generated with", R.version.string, "on", R.version$platform, "."), "",
  "The generator uses base R only. It does not load Numerics or any production numerical helper.",
  "All CSV matrices use location=0, scale=1, sample_size=1 and parameter order xi, alpha, kappa[, hondo].",
  "GNO denotes Hosking's transformed normal, not the symmetric exponential-power distribution.",
  "GLO fixes Kappa Four hondo=-1; GEV fixes hondo=0. A fixed-shape family inverts its own information block.", "",
  "## Method and acceptance checks", "",
  sprintf("R integrate (QUADPACK): relative tolerance %.3g; absolute tolerance %.3g; at most 1500 subdivisions.", relative_tolerance, absolute_tolerance),
  "GNO uses direct integration over the entire normal latent line, with an independent probability-coordinate cross-check.",
  "GLO, GEV and K4 pair p=u^a/2 and 1-p=u^a/2 over the complete u interval [0,1]. Powers 16/24 are cross-checked; boundary-approach cases use 128/192.",
  "Log probabilities and log survival probabilities are retained separately. Below exp(-36), log(t)=log(survival) is used to binary64 relative precision; no probability tail is discarded.",
  "Analytical scores are checked against five-point fixed-observation derivatives of an independently expressed log density at probabilities .1, .5 and .9, using two step sizes.",
  "Finite matrices, near-zero mean scores, positive eigenvalues, inverse residuals and agreement between integration coordinates are required before writing the fixture.",
  "The covariance error column is a conservative matrix-norm perturbation estimate from the numerical integration error estimates; it is not a rigorous interval bound.",
  "The alternate-coordinate difference is independent numerical corroboration, not a rigorous bound. Values are binary64 R references, not arbitrary-precision goldens.",
  "The exact GNO kappa=0 covariance is [[7/6,0,1/3],[0,1/2,0],[1/3,0,2/3]] with all three parameters estimated.", "",
  "The additional moment rows use row=1 mean, row=2 standard deviation, row=3 skewness, row=4 non-excess kurtosis. GNO uses analytical lognormal moments; GLO uses trigonometric moments, with factored coefficient-array series for abs(kappa)<=.01.",
  "Moment absolute_error_estimate is an arithmetic comparison allowance of 5e-13*max(1,abs(value)), not a quadrature error estimate. The GLO near-zero series is cross-checked against independent 70-digit Decimal values at kappa=.00010001 and .0002.", "",
  "## Mathematical domain", "",
  "GNO information is finite for every finite kappa. GLO uses abs(kappa)<1/2. K4 uses kappa<1/2, hondo<1/2 and kappa*hondo<1/2. GEV uses kappa<1/2.",
  "These are local asymptotic-information conditions, not guarantees of global-MLE existence, estimator convergence, or finite-sample coverage.",
  "MLE covariance is not L-moment or product-moment covariance. The fixture does not authorize substitution between estimators.", "",
  "## Primary-source mapping and references", "",
  "- [Hosking lmom R defining quantiles](https://raw.githubusercontent.com/cran/lmom/master/R/lmom.r): quagno(p,c(xi,alpha,k)), quaglo(p,c(xi,alpha,k)), quakap(p,c(xi,alpha,k,h)), quagev(p,c(xi,alpha,k)).",
  "- [Park and Kim (2007), Fisher information matrix for a four-parameter kappa distribution](https://doi.org/10.1016/j.spl.2007.03.002). The present scores/domain are independently derived; inaccessible full-paper formulae were not treated as verified values.",
  "- [Wang and Flournoy (2015), local likelihood estimation for the three-parameter lognormal](https://doi.org/10.1016/j.spl.2015.05.021): finite Fisher information does not imply bounded global likelihood.", "",
  "## Numerical summary", "",
  "| Family | kappa | hondo | Minimum eigenvalue | Condition number | Max mean score | Information error norm | Coordinate delta norm | Covariance error estimate | Fixed-x derivative relative error |",
  "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
table <- apply(summary, 1, function(r) sprintf("| %s | %.6g | %.6g | %.6g | %.6g | %.3g | %.3g | %.3g | %.3g | %.3g |",
    r[1], as.numeric(r[2]), as.numeric(r[3]), as.numeric(r[4]), as.numeric(r[5]), as.numeric(r[6]), as.numeric(r[7]), as.numeric(r[8]), as.numeric(r[9]), as.numeric(r[11])))
writeLines(c(header, table, "", paste(nrow(result), "CSV rows from", nrow(summary), "parameter cases.")), notes, useBytes = TRUE)
cat("Wrote", target, "and", notes, "\n")
