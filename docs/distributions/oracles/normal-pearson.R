# Independent R stats oracle for the approved Normal/Pearson family repairs.
# Run with R 4.4.3. No Numerics assembly is loaded; CSV is frozen for runtime-independent tests.
options(digits = 17)
arguments <- commandArgs(trailingOnly = TRUE)
directory <- if (length(arguments)) arguments[[1]] else "."
qpearson <- function(p, theta) {
  mu <- theta[1]; sigma <- theta[2]; skew <- theta[3]
  if (skew == 0) return(qnorm(p, mu, sigma))
  shape <- 4/skew^2
  mu + sigma * skew/2 * (qgamma(p, shape, lower.tail = skew > 0) - shape)
}
qfamily <- function(p, family, theta, base) {
  if (family == "PearsonTypeIII") return(qpearson(p, theta))
  if (family == "LogPearsonTypeIII") return(base^qpearson(p, theta))
  if (family == "LogNormal") return(qlnorm(p, theta[1] * log(base), theta[2] * log(base)))
  variance <- log1p((theta[2]/theta[1])^2)
  qlnorm(p, log(theta[1]) - variance/2, sqrt(variance))
}
gradient <- function(p, family, theta, base) {
  vapply(seq_along(theta), function(i) {
    step <- 1e-4 * max(1, abs(theta[i]))
    at <- function(offset) { candidate <- theta; candidate[i] <- candidate[i] + step * offset; qfamily(p, family, candidate, base) }
    (at(-2) - 8*at(-1) + 8*at(1) - at(2))/(12*step)
  }, 0.0)
}
covariance <- function(family, theta, base, n = 100) {
  if (family == "LogNormal") return(diag(c(theta[2]^2, theta[2]^2/2))/n)
  if (family == "LnNormal") {
    variance <- log1p((theta[2]/theta[1])^2); s <- sqrt(variance)
    J <- matrix(c(theta[1], theta[1]*s, theta[2], theta[2]*s + theta[1]^2*exp(variance)*s/theta[2]), 2, byrow = TRUE)
    return(J %*% diag(c(variance, variance/2)) %*% t(J)/n)
  }
  sigma <- theta[2]; skew <- theta[3]
  if (skew == 0) return(diag(c(sigma^2, sigma^2/2, 6))/n)
  shape <- 4/skew^2; beta <- sigma*skew/2
  information <- matrix(c(1/(beta^2*(shape-2)), 1/beta^2, 1/(beta*(shape-1)),
                          1/beta^2, shape/beta^2, 1/beta,
                          1/(beta*(shape-1)), 1/beta, trigamma(shape)), 3, byrow = TRUE)
  J <- matrix(c(1, shape, beta, 0, sign(beta)*sqrt(shape), abs(beta)/(2*sqrt(shape)), 0, 0, -sign(beta)/shape^1.5), 3, byrow = TRUE)
  J %*% solve(information) %*% t(J)/n
}
rows <- list()
add <- function(family, theta, base, p) {
  grad <- gradient(p, family, theta, base)
  cov <- covariance(family, theta, base)
  count <- length(theta)
  rows[[length(rows)+1]] <<- data.frame(family=family, mu=theta[1], sigma=theta[2], gamma=if(count==3) theta[3] else 0,
    base=base, probability=p, quantile=qfamily(p,family,theta,base), gradient_mu=grad[1], gradient_sigma=grad[2],
    gradient_gamma=if(count==3) grad[3] else 0, cov11=cov[1,1], cov12=cov[1,2], cov13=if(count==3) cov[1,3] else 0,
    cov22=cov[2,2], cov23=if(count==3) cov[2,3] else 0, cov33=if(count==3) cov[3,3] else 0,
    variance_mle=as.numeric(t(grad)%*%cov%*%grad))
}
for (skew in c(-1.2,-0.2,0,0.2,1.2)) for (p in c(0.01,0.83,0.99)) {
  add("PearsonTypeIII", c(2,3,skew), exp(1), p)
  add("LogPearsonTypeIII", c(0.3,0.2,skew), exp(1), p)
}
for (p in c(0.01,0.5,0.83,0.99)) {
  add("LnNormal", c(10,2), exp(1), p)
  add("LnNormal", c(1,1), exp(1), p)
  for (base in c(2,exp(1),10)) add("LogNormal", c(0.3,0.2), base, p)
}
write.csv(do.call(rbind,rows), file.path(directory,"normal-pearson.csv"), row.names=FALSE, quote=FALSE)
sink(file.path(directory,"normal-pearson-evidence.txt"))
cat(R.version.string, "\n")
cat("Normal logCDF(-40):", pnorm(-40,log.p=TRUE), "\n")
cat("Normal logPDF(0), sigma=1e308:", dnorm(0,sd=1e308,log=TRUE), "\n")
cat("Logistic logPDF(-1000):", dlogis(-1000,log=TRUE), "\n")
cat("LogNormal moment sigma, mean10 sd2 base10:", sqrt(log1p(0.2^2))/log(10), "\n")
cat("PIII(0,1,-2), quantile 1e-20:", qpearson(1e-20,c(0,1,-2)), "\n")
cat("LnNormal indirect-MoM median variance(1,1), n100:", log(2)/200, "\n")
cat("LnNormal(10,2) exact physical mean derivative at median = 135/(26*sqrt(26)):",135/(26*sqrt(26)),"\n")
cat("LnNormal direct physical-moment estimator is a different estimator; its median variance would be 0.00875.\n")
cat("LP3(0,.8,1), base e, finite mean and SD:")
raw <- function(r) if(1-0.4*r <= 0) Inf else exp(-1.6*r)*(1-0.4*r)^-4
print(c(mean=raw(1),sd=sqrt(raw(2)-raw(1)^2),third=raw(3),fourth=raw(4)))
cat("Legacy LP3 default and (1,1,1) corrected modes:\n")
lp3mode <- function(mu,sigma,skew,base=10) {
  t <- sigma*log(base); beta <- sigma*skew/2*log(base)
  if(skew==0) return(exp(mu*log(base)-t*t))
  exp(mu*log(base)-(t*t+beta)/(1+beta))
}
print(c(default=lp3mode(3,.5,0), custom=lp3mode(1,1,1)))
cat("Asymptotic MLE covariance uses independent inversion of location/scale/shape Fisher information, transformed to public coordinates.\n")
cat("Quantile gradients use five-point finite differences of R stats qgamma/qlnorm, with public parameter perturbations.\n")
cat("Independent review overflow regressions (defining lognormal identities):\n")
cat("LP3 zero-skew Sigma16 skew:",(exp(256)+2)*sqrt(expm1(256)),"\n")
cat("LP3 zero-skew Sigma12 kurtosis:",3+expm1(576)+2*expm1(432)+3*expm1(288),"\n")
cat("LogNormal Mu370 Sigma1e-10 base e median variance n100:",exp(740+log(1e-22)),"\n")
cat("Nonzero LP3 large-moment checks from its gamma moment generating function, log normalized before exponentiation:\n")
for (skew in c(-0.001,0.001)) for (s in c(12,16)) {
  shape <- 4/skew^2; beta <- s*skew/2
  d <- function(r) shape*(r*log1p(-beta)-log1p(-r*beta))
  logvar <- log(expm1(d(2)))
  third <- exp(d(3)-1.5*logvar)*(1-3*exp(d(2)-d(3))+2*exp(-d(3)))
  fourth <- exp(d(4)-2*logvar)*(1-4*exp(d(3)-d(4))+6*exp(d(2)-d(4))-3*exp(-d(4)))
  cat("Sigma",s,"Gamma",skew,"Skewness",third,"Kurtosis",fourth,"\n")
}
sink()
evidence_path <- file.path(directory,"normal-pearson-evidence.txt")
writeLines(sub("[[:blank:]]+$", "", readLines(evidence_path)), evidence_path)
