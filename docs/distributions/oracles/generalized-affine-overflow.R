# Independent four-case Hosking GNO/GLO tail reference.
# Base R only; no Numerics assembly or production helper is loaded.
# Run: Rscript --vanilla generalized-affine-overflow.R
# Inputs are the binary64 constants used by the C# tests. The standardized
# ratio x/alpha exceeds binary64 range, so the defining support factor is
# assembled from its logarithm before evaluating the normal or logistic law.

cat("# ", R.version.string, "\n", sep = "")
cat("# platform=", R.version$platform, "; locale=", Sys.getlocale(), "\n", sep = "")
cat("family,x,alpha,kappa,latent,log_tail,log_pdf,pdf\n")
alpha <- 1e-200
for (family in c("GNO", "GLO")) {
  for (direction in c(-1, 1)) {
    x <- direction * 1e200
    kappa <- -direction * 20
    log_product <- log(abs(kappa)) + log(abs(x)) - log(alpha)
    log_support <- log_product + log1p(exp(-log_product))
    latent <- -log_support / kappa
    if (family == "GNO") {
      log_tail <- pnorm(abs(latent), lower.tail = FALSE, log.p = TRUE)
      log_pdf <- dnorm(latent, log = TRUE) + kappa * latent - log(alpha)
    } else {
      # Defining logistic tail/density, independently cross-checked by stats.
      log_tail <- -abs(latent) - log1p(exp(-abs(latent)))
      log_pdf <- -abs(latent) - 2 * log1p(exp(-abs(latent))) + kappa * latent - log(alpha)
      stopifnot(abs(log_tail - plogis(abs(latent), lower.tail = FALSE, log.p = TRUE)) < 1e-13)
      stopifnot(abs(log_pdf - (dlogis(latent, log = TRUE) + kappa * latent - log(alpha))) < 1e-12)
    }
    stopifnot(is.finite(latent), is.finite(log_tail), is.finite(log_pdf))
    cat(paste(c(family, sprintf("%.17g", c(x, alpha, kappa, latent, log_tail, log_pdf, exp(log_pdf)))), collapse = ","), "\n", sep = "")
  }
}
