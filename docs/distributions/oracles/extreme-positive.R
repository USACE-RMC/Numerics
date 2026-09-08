# Offline, independent oracle generator. Requires only R 4.4.3 base/stats.
# Run from the repository root; never calls Numerics or changes production/tests.
options(digits=17, warn=1)
stopifnot(as.character(getRversion()) == "4.4.3")
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[1] else "docs/distributions/oracles/extreme-positive.csv"
rows <- list()
fmt <- function(x) if (is.na(x)) "NaN" else sprintf("%.17g", x)
add <- function(family, case, scale=1, shape=0, xi=0, x=NA_real_, p=NA_real_,
                n=NA_real_, quantity, value, oracle, rel=2e-11, abs=0,
                error=0, status="finite") {
  if (is.na(value) && status == "finite") status <- "undefined"
  if (is.infinite(value)) status <- "infinite"
  rows[[length(rows)+1L]] <<- data.frame(family, case, xi=fmt(xi), scale=fmt(scale),
    shape=fmt(shape), x=fmt(x), p=fmt(p), n=fmt(n), quantity, value=fmt(value),
    absolute_tolerance=fmt(abs), relative_tolerance=fmt(rel),
    estimated_absolute_error=fmt(error), oracle, status, stringsAsFactors=FALSE)
}
exprel <- function(t) {
  ans <- expm1(t)/t
  small <- abs(t)<1e-3
  z <- t[small]
  ans[small] <- 1+z*(1/2+z*(1/6+z*(1/24+z*(1/120+z*(1/720+z/5040)))))
  ans
}
exprel2 <- function(t) {
  ans <- (expm1(t)-t)/(t*t)
  small <- abs(t)<1e-3
  z <- t[small]
  ans[small] <- 1/2+z*(1/6+z*(1/24+z*(1/120+z*(1/720+z*(1/5040+z/40320)))))
  ans
}
log1mexp <- function(x) ifelse(x < -log(2), log1p(-exp(x)), log(-expm1(x)))
standard_q <- function(family, logp, k) {
  z <- switch(family, GEV=-log(-logp), GPA=-log1mexp(logp), Gumbel=-log(-logp))
  if (family=="Gumbel") z else z*exprel(-k*z)
}
gradient <- function(family,p,a,k) {
  if (family=="Exponential") return(c(1,-log1p(-p)))
  if (family=="Gumbel") return(c(1,-log(-log(p))))
  if (family=="Weibull") {
    lt <- log(-log1p(-p)); q <- a*exp(lt/k)
    return(c(q/a,-q*lt/k^2))
  }
  z <- if(family=="GEV") -log(-log(p)) else -log1p(-p)
  # Independent derivative of z*(exp(-k*z)-1)/(-k*z).
  c(1,z*exprel(-k*z),a*z*z*(exprel2(-k*z)-exprel(-k*z)))
}
# Sixth-order Richardson estimates using positive log-parameter perturbations.
# Error estimates are empirical convergence indicators, not rigorous bounds.
differentiate <- function(f) {
  hs <- .04/2^(0:11)
  d <- vapply(hs,function(h)(f(-2*h)-8*f(-h)+8*f(h)-f(2*h))/(12*h),0.0)
  r <- (16*d[-1]-d[-length(d)])/15
  err <- abs(diff(r))
  good <- which(is.finite(err) & is.finite(r[-1]))
  if (!length(good)) return(c(value=NA_real_,error=Inf))
  j <- good[which.min(err[good])]+1
  c(value=r[j],error=max(err[j-1],32*.Machine$double.eps*abs(r[j])))
}
gamma_quantile <- function(p,k) if(p<=.5) qgamma(log(p),shape=k,log.p=TRUE) else qgamma(log1p(-p),shape=k,lower.tail=FALSE,log.p=TRUE)
gamma_shape_derivative <- function(p,k) {
  q <- gamma_quantile(p,k)
  if (!is.finite(q) || q==0) return(c(value=NA,error=NA,implicit=NA))
  d <- differentiate(function(u) log(gamma_quantile(p,k*exp(u))))
  value <- (q/k)*d["value"]
  err <- abs(q/k)*d["error"]
  lower <- p<=.5
  dp <- differentiate(function(u) pgamma(q,shape=k*exp(u),lower.tail=lower,log.p=TRUE))
  lp <- if(lower) log(p) else log1p(-p)
  implicit <- (if(lower)-1 else 1)*exp(lp-dgamma(q,shape=k,log=TRUE))/k*dp["value"]
  err <- max(err,abs(value-implicit))
  c(value=unname(value),error=unname(err),implicit=unname(implicit))
}
# Native R stats probability and log-probability oracles in standardized units.
for (family in c("Exponential","Gamma","Weibull")) {
  shapes <- switch(family,Exponential=0,Gamma=c(.001,.01,.1,.5,1,2,10,100,1e4),Weibull=c(.1,.5,1,2,10,100))
  for(k in shapes) for(a in c(1,1e-200,1e200)) {
    zs <- if(family=="Gamma") c(0,1e-10,.1,1,max(1,k),1000) else c(0,1e-20,.1,1,10,1000)
    for(z in unique(zs)) {
      lp <- switch(family,Exponential=dexp(z,log=TRUE),Gamma=dgamma(z,shape=k,log=TRUE),Weibull=dweibull(z,shape=k,log=TRUE))-log(a)
      lc <- switch(family,Exponential=pexp(z,log.p=TRUE),Gamma=pgamma(z,shape=k,log.p=TRUE),Weibull=pweibull(z,shape=k,log.p=TRUE))
      ls <- switch(family,Exponential=pexp(z,lower.tail=FALSE,log.p=TRUE),Gamma=pgamma(z,shape=k,lower.tail=FALSE,log.p=TRUE),Weibull=pweibull(z,shape=k,lower.tail=FALSE,log.p=TRUE))
      log_density_repaired <- FALSE; log_cdf_repaired <- FALSE
      if(family=="Weibull" && z>0 && is.finite(z)) {
        log_power <- k*log(z)
        defining_lp <- log(k)+(k-1)*log(z)-exp(log_power)-log(a)
        if(is.infinite(lp) && lp<0 && is.finite(defining_lp)) {lp<-defining_lp;log_density_repaired<-TRUE}
        # When z^k underflows, log(1-exp(-z^k)) equals k*log(z) to binary64 precision.
        if(is.infinite(lc) && lc<0 && is.finite(log_power) && log_power<0) {lc<-log_power;log_cdf_repaired<-TRUE}
      }
      for(nm in c("LogPDF","LogCDF","LogCCDF")) add(family,"R-standardized-evaluation",a,k,x=a*z,quantity=nm,value=switch(nm,LogPDF=lp,LogCDF=lc,LogCCDF=ls),oracle=if((nm=="LogPDF" && log_density_repaired)||(nm=="LogCDF" && log_cdf_repaired))"Weibull-defining-log-formula-R-underflow" else "R-stats",abs=if(nm=="LogPDF")2e-11 else 0)
    }
  }
}
for(family in c("Exponential","Gamma","Weibull")) {
  shapes <- switch(family,Exponential=0,Gamma=c(.001,.01,.1,.5,1,2,10,100,1e4),Weibull=c(.1,.5,1,2,10,100))
  for(k in shapes) for(p in c(1e-12,1e-6,.01,.5,.99,1-1e-12)) {
    q <- switch(family,Exponential=qexp(p),Gamma=gamma_quantile(p,k),Weibull=qweibull(p,shape=k))
    add(family,"R-quantile",shape=k,p=p,quantity="InverseCDF",value=q,oracle="R-stats",status=if(q==0)"underflowed" else "finite")
    if(family=="Gamma") {
      d <- gamma_shape_derivative(p,k)
      add(family,"exact-quantile-gradient",shape=k,p=p,quantity="GradientScale",value=q,oracle="R-qgamma-scale-identity",status=if(q==0)"underflowed" else "finite")
      # Do not freeze a zero shape derivative when R qgamma has underflowed.
      add(family,"exact-quantile-gradient",shape=k,p=p,quantity="GradientShape",value=d["value"],oracle="R-qgamma-log-shape-Richardson-crosschecked-pgamma",rel=5e-8,error=d["error"],status=if(q==0)"unresolved-quantile-underflow" else "finite")
      if(q>0 && !(is.finite(d["value"]) && d["value"]>0 && d["error"]<=5e-8*abs(d["value"]))) stop(sprintf("Gamma derivative unresolved k=%g p=%.17g value=%.17g error=%.17g implicit=%.17g",k,p,d["value"],d["error"],d["implicit"]))
    } else {
      g <- gradient(family,p,1,k)
      for(j in seq_along(g)) add(family,"analytic-quantile-gradient",shape=k,p=p,quantity=paste0("Gradient",j),value=g[j],oracle="defining-quantile-derivative")
    }
  }
}
# GEV/GPA use Hosking's kappa convention; SciPy mappings: genextreme(c=k), genpareto(c=-k).
for(family in c("GEV","GPA","Gumbel")) {
  ks <- if(family=="Gumbel")0 else c(-1,-.5,-.2,-1e-6,0,1e-6,.2,.5,1,2)
  for(k in ks) {
    for(p in c(1e-12,.01,.5,.99,1-1e-12)) {
      q <- standard_q(family,log(p),k)
      add(family,"defining-quantile",shape=k,p=p,quantity="InverseCDF",value=q,oracle="Hosking-quantile-expm1")
      g <- gradient(family,p,1,k)
      for(j in seq_along(g)) add(family,"analytic-quantile-gradient",shape=k,p=p,quantity=paste0("Gradient",j),value=g[j],oracle="defining-quantile-derivative")
    }
    lo <- if(family=="GPA")0 else if(k<0)1/k else -Inf
    hi <- if(k>0)1/k else Inf
    for(nm in c("Minimum","Maximum")) add(family,"support",shape=k,quantity=nm,value=if(nm=="Minimum")lo else hi,oracle="mathematical-support")
    xs <- unique(c(lo,hi,if(family=="Gumbel")c(-10,0,10) else c(0,.1)))
    for(x in xs) {
      if(x<lo || x>hi) {lp <- -Inf; lc <- if(x<lo)-Inf else 0; ls <- if(x<lo)0 else -Inf}
      else if(x==hi && is.finite(hi)) {lp<-if(k<1)-Inf else if(k==1)0 else Inf;lc<-0;ls<--Inf}
      else if(x==lo && family!="GPA") {lp<--Inf;lc<--Inf;ls<-0}
      else {
        z <- if(k==0)x else -log1p(-k*x)/k
        if(family=="GPA") {lp<--(1-k)*z;ls<--z;lc<-log1mexp(ls)}
        else {lp<--(1-k)*z-exp(-z);lc<--exp(-z);ls<-log1mexp(lc)}
        if(is.nan(lp) && is.infinite(x))lp<--Inf
      }
      for(nm in c("LogPDF","LogCDF","LogCCDF")) add(family,"boundary-and-evaluation",shape=k,x=x,quantity=nm,value=switch(nm,LogPDF=lp,LogCDF=lc,LogCCDF=ls),oracle="defining-density-and-one-sided-limits")
    }
  }
}
# Subnormal scale reproduction without forming an overflowing reciprocal scale.
for(family in c("Exponential","Weibull"))add(family,"subnormal-scale-log-density",scale=1e-320,shape=if(family=="Weibull")1 else 0,x=1e-320,quantity="LogPDF",value=-1-log(1e-320),oracle="unit-exponential-density-plus-log-scale",abs=2e-11)
# Defining Hosking tail formulas with log(1-kappa*(x-xi)/scale) evaluated
# from the physical coordinates, never from an overflowing standardized value.
tail_cases <- list(c(xi=0,a=1e-200,k=-20,x=1e200),
                   c(xi=0,a=1,k=-20,x=1e308),
                   c(xi=-1e308,a=1e-200,k=-20,x=1e308),
                   c(xi=0,a=1e-200,k=-1e-200,x=1e200),
                   c(xi=0,a=1e-200,k=0,x=1e200))
for(family in c("GEV","GPA")) {
  cases <- tail_cases
  if(family=="GEV") cases <- c(cases,list(c(xi=0,a=1e-200,k=20,x=-1e200),c(xi=1e308,a=1e-200,k=20,x=-1e308)))
  for(test in cases) {
    xi<-test["xi"];a<-test["a"];k<-test["k"];x<-test["x"]
    magnitude<-max(abs(x),abs(xi))
    log_difference<-log(magnitude)+log(abs(x/magnitude-xi/magnitude))
    if(k==0)y<-Inf else {
      log_product<-log(abs(k))+log_difference-log(a)
      log_support<-max(0,log_product)+log1p(exp(-abs(log_product)))
      y<--log_support/k
    }
    if(family=="GPA") {lp<--(1-k)*y-log(a);ls<--y;lc<-log1mexp(ls)}
    else {lp<--(1-k)*y-exp(-y)-log(a);lc<--exp(-y);ls<-if(y>36)-y else log1mexp(lc)}
    for(nm in c("LogPDF","LogCDF","LogCCDF"))add(family,"overflowing-tail-transform",a,k,xi=xi,x=x,quantity=nm,value=switch(nm,LogPDF=lp,LogCDF=lc,LogCCDF=ls),oracle="R-defining-log-Hosking-physical-coordinates")
  }
}
# Existing rounded Weibull MLE coefficients, contracted analytically before
# the coordinate magnitudes (lambda=kappa=1e100) can hide finite contributions.
log_t<-log(log(2))
add("Weibull","MLE-mismatched-range-variance",scale=1e100,shape=1e100,p=.5,n=100,
    quantity="QuantileVariance",value=(1.108665-2*.257022*log_t+.607927*log_t^2)/100,
    oracle="R-analytical-Weibull-variance-retained-coefficients")
# At kappa >= 1e16, integrate the scaled Fisher-residual kernel in u=kappa*t.
# The first omitted integrated Bernoulli term has magnitude 1/(30*kappa^3).
for(test in list(c(a=1,k=1e16,n=100),c(a=1e150,k=1e16,n=100),
                 c(a=1e-150,k=1e16,n=100),c(a=1,k=1e155,n=2e9))) {
  a<-test["a"];k<-test["k"];n<-test["n"]
  residual<-integrate(function(u)exp(-u)*(u/2+u*u/(12*k)),0,Inf,rel.tol=1e-13,abs.tol=1e-14)$value
  vals<-c(exp(2*log(a)-log(n)+log(1/residual+1/k)),
          -exp(log(a)+log(k)-log(n)-log(residual)),
          exp(2*log(k)-log(n)-log(residual)))
  for(j in 1:3)add("Gamma","MLE-large-shape-Fisher-residual",a,k,n=n,quantity=c("Covariance11","Covariance12","Covariance22")[j],value=vals[j],oracle="R-integrate-scaled-Fisher-residual-kernel",rel=2e-10)
}
# At the median q=k-1/3+O(1/k), dq/dk=1+O(1/k^2). The defining
# positive Fisher/MoM factor therefore gives k/n with relative correction O(1/k).
for(method in c("MLE","MoM"))add("Gamma",paste0(method,"-large-shape-median-variance"),shape=1e16,p=.5,n=100,
    quantity="QuantileVariance",value=1e16/100,oracle="R-defining-concentrated-Gamma-median-variance-limit",rel=2e-10)
add("Weibull","scale-rescued-median",scale=1e200,shape=.0004,p=.5,
    quantity="Median",value=exp(log(1e200)+log(log(2))/.0004),oracle="R-defining-log-Weibull-median")
for(k in c(1e200,1e308)) {
  r<-1/k
  vals<-c((2*(r-1)/(r+3))*sqrt(k)*sqrt(r+2),
          k*(3*(r+2)*(3*r*r-r+2)/((r+3)*(r+4))))
  for(j in 1:2)add("GPA","large-positive-shape-scaled-higher-moments",shape=k,
      quantity=c("Skewness","Kurtosis")[j],value=vals[j],oracle="R-Hosking-rational-moments-in-reciprocal-shape")
}
for(family in c("GEV","GPA"))for(k in c(-1e-12,1e-12,-1e-8,1e-8,-1e-4,1e-4)) {
  add(family,"nonzero-small-shape-support",shape=k,quantity="Minimum",value=if(family=="GPA")0 else if(k<0)1/k else -Inf,oracle="exact-nonzero-shape-support")
  add(family,"nonzero-small-shape-support",shape=k,quantity="Maximum",value=if(k>0)1/k else Inf,oracle="exact-nonzero-shape-support")
  for(p in c(.01,.5,.99)) {
    add(family,"nonzero-small-shape-quantile",shape=k,p=p,quantity="InverseCDF",value=standard_q(family,log(p),k),oracle="Hosking-quantile-expm1")
    g<-gradient(family,p,1,k)
    for(j in 1:3)add(family,"nonzero-small-shape-gradient",shape=k,p=p,quantity=paste0("Gradient",j),value=g[j],oracle="defining-quantile-derivative")
  }
}
# Full probability integration for existence-aware GEV/GPA centered moments.
integrate_pair <- function(fun) integrate(function(u){lt=8*log(u)-log(2);jac=4*u^7;(fun(lt)+fun(log1p(-exp(lt))))*jac},0,1,rel.tol=1e-10,abs.tol=1e-12,subdivisions=1000)
for(family in c("GEV","GPA")) for(k in c(-1,-.5,-.2,-1e-6,0,1e-6,.2,.5,1,2)) {
  vals <- rep(NA_real_,4); errs <- rep(0,4)
  if(k> -1) {
    q <- function(lp)standard_q(family,lp,k)
    r <- integrate_pair(q);vals[1]<-r$value;errs[1]<-r$abs.error
    if(k> -.5) {
      r<-integrate_pair(function(lp)(q(lp)-vals[1])^2);vals[2]<-sqrt(r$value);errs[2]<-r$abs.error/(2*vals[2])
      if(k> -1/3){r<-integrate_pair(function(lp)((q(lp)-vals[1])/vals[2])^3);vals[3]<-r$value;errs[3]<-r$abs.error}
      if(k> -.25){r<-integrate_pair(function(lp)((q(lp)-vals[1])/vals[2])^4);vals[4]<-r$value;errs[4]<-r$abs.error}
    }
  }
  for(j in 1:4)add(family,"full-probability-central-moments",shape=k,quantity=c("Mean","StandardDeviation","Skewness","Kurtosis")[j],value=vals[j],oracle="R-integrate-full-probability-Hosking-quantile",rel=2e-8,abs=2e-9,error=errs[j])
}
# Large positive GEV shapes have finite standardized ratios even when raw gamma moments overflow.
for(k in c(50,100,200)) for(a in c(1,1e-200,1e200)) {
  l<-lgamma(1+(1:4)*k)
  lv<-l[2]+log1p(-exp(2*l[1]-l[2]))
  vals<-c(-exp(log(a)+l[1]+log1p(-exp(-l[1]))-log(k)),
          exp(log(a)+lv/2-log(k)),
          -exp(l[3]+log1p(-3*exp(l[1]+l[2]-l[3])+2*exp(3*l[1]-l[3]))-1.5*lv),
          exp(l[4]+log1p(-4*exp(l[1]+l[3]-l[4])+6*exp(2*l[1]+l[2]-l[4])-3*exp(4*l[1]-l[4]))-2*lv))
  for(j in 1:4)add("GEV","large-positive-shape-lgamma-moments",a,k,quantity=c("Mean","StandardDeviation","Skewness","Kurtosis")[j],value=vals[j],oracle="R-lgamma-normalized-analytical-moments",rel=2e-10)
}
# Weibull X=scale*T^(1/kappa), T unit exponential. Small kappa needs lgamma
# moments, while large kappa is integrated after centering and dividing by c=1/kappa.
for(k in c(.005,.01,.02)) for(a in c(1,1e-200,1e200)) {
  c<-1/k;l<-lgamma(1+(1:4)*c);lv<-l[2]+log1p(-exp(2*l[1]-l[2]))
  vals<-c(exp(log(a)+l[1]),exp(log(a)+lv/2),
          exp(l[3]+log1p(-3*exp(l[1]+l[2]-l[3])+2*exp(3*l[1]-l[3]))-1.5*lv),
          exp(l[4]+log1p(-4*exp(l[1]+l[3]-l[4])+6*exp(2*l[1]+l[2]-l[4])-3*exp(4*l[1]-l[4]))-2*lv))
  for(j in 1:4)add("Weibull","small-shape-lgamma-moments",a,k,quantity=c("Mean","StandardDeviation","Skewness","Kurtosis")[j],value=vals[j],oracle="R-lgamma-normalized-analytical-power-moments",rel=2e-10)
}
for(k in c(1e3,1e6,1e12,1e200)) {
  c<-1/k;q<-function(lp)expm1(c*log(-lp))/c
  m<-integrate_pair(q)$value;v<-integrate_pair(function(lp)(q(lp)-m)^2)$value
  skew<-integrate_pair(function(lp)((q(lp)-m)/sqrt(v))^3)$value
  kurt<-integrate_pair(function(lp)((q(lp)-m)/sqrt(v))^4)$value
  for(a in c(1,1e-200,1e200)) {
    vals<-c(exp(log(a)+lgamma(1+c)),exp(log(a)+log(c)+log(v)/2),skew,kurt)
    for(j in 1:4)add("Weibull","large-shape-full-probability-moments",a,k,quantity=c("Mean","StandardDeviation","Skewness","Kurtosis")[j],value=vals[j],oracle="R-integrate-centered-scaled-exponential-power",rel=2e-8)
  }
}
for(family in c("Exponential","Gamma","Gumbel","Weibull")) for(a in c(1,1e-200,1e200)) {
  ks<-switch(family,Exponential=0,Gamma=c(.1,1,2,100),Gumbel=0,Weibull=c(.5,1,2,10))
  for(k in ks) {
    vals <- switch(family,Exponential=c(a,a,2,9),Gamma=c(a*k,a*sqrt(k),2/sqrt(k),3+6/k),Gumbel=c(-digamma(1)*a,a*pi/sqrt(6),1.1395470994046487,5.4),Weibull={g<-gamma(1+(1:4)/k);v<-g[2]-g[1]^2;c(a*g[1],a*sqrt(v),(g[3]-3*g[1]*g[2]+2*g[1]^3)/v^1.5,(g[4]-4*g[1]*g[3]+6*g[1]^2*g[2]-3*g[1]^4)/v^2)})
    for(j in 1:4)add(family,"analytic-scale-separated-moments",a,k,quantity=c("Mean","StandardDeviation","Skewness","Kurtosis")[j],value=vals[j],oracle="R-gamma-and-analytical-moments",rel=if(family=="Weibull")2e-8 else 2e-11)
  }
}
covrows <- function(family,k,C,method,a=1,n=100)for(i in 1:nrow(C))for(j in i:ncol(C))add(family,paste0(method,"-covariance"),a,k,n=n,quantity=paste0("Covariance",i,j),value=C[i,j],oracle=method,rel=2e-9,abs=1e-12)
for(n in c(2,10,100,50000))covrows("Exponential",0,diag(c(1/n^2,(n-1)/n^2)),"actual-MLE-order-statistics",n=n)
for(k in c(.001,.1,1,2,100,1e4)) {
  t<-trigamma(k); d<-k*t-1
  covrows("Gamma",k,matrix(c(t,-1,-1,k),2)/(100*d),"MLE-R-trigamma")
  J<-matrix(c(-1/k,2,1/k,-1),2);S<-matrix(c(k,2*k,2*k,2*k*(k+3)),2)
  covrows("Gamma",k,J%*%S%*%t(J)/100,"MoM-independent-moment-Jacobian")
}
euler <- -digamma(1); B<-6/pi^2; A<-1+B*(1-euler)^2; C<-B*(1-euler)
covrows("Gumbel",0,matrix(c(A,C,C,B),2)/100,"MLE-exact-Euler-pi")
for(k in c(.5,1,2,10))covrows("Weibull",k,matrix(c(A/k^2,C,C,B*k^2),2)/100,"MLE-exact-Euler-pi")
for(k in c(-.2,0,.2,.49)) {
  covrows("GPA",k,matrix(c(1/((100+2*k)*(100+k)^2)*100,0,0,0,2*(1-k)/100,(1-k)/100,0,(1-k)/100,(1-k)^2/100),3),"MLE-Hosking-Wallis")
}
for(k in c(-.2,0,.2,1,2)) {
  den<-(1+2*k)*(1+3*k)*(1+4*k)
  vscale<-2*(1+k)^2*(1+6*k+12*k^2)/(100*den)
  vshape<-(1+k)^2*(1+2*k)^2*(1+k+6*k^2)/(100*den)
  cross<-(1+k)^2*(1+2*k)*(1+4*k+12*k^2)/(100*den)
  C<-matrix(c(100/((100+2*k)*(100+k)^2),0,0,0,vscale,cross,0,cross,vshape),3)
  covrows("GPA",k,C,"MoM-Hosking-Wallis")
}
gev_scores <- function(lp,k) {
  z<--log(-lp);r<-1+lp
  A<-z*exprel(k*z);B<-z*z*exprel2(k*z)
  cbind((r-k)*exp(k*z),-1+(r-k)*A,z+(k-r)*B)
}
for(k in c(-.2,-1e-6,0,1e-6,.2,.49)) {
  I<-matrix(0,3,3)
  for(i in 1:3)for(j in i:3){z<-integrate_pair(function(lp){s<-gev_scores(lp,k);s[,i]*s[,j]});I[i,j]<-I[j,i]<-z$value}
  stopifnot(min(eigen(I,symmetric=TRUE,only.values=TRUE)$values)>0)
  covrows("GEV",k,solve(I)/100,"MLE-independent-score-outer-product-integral")
}
for(family in c("GEV","GPA"))for(k in c(.5,1,2))add(family,"MLE-nonregular-domain",shape=k,n=100,quantity="CovarianceDefined",value=0,oracle="regular-Fisher-domain-kappa-less-than-half",status="outside-uncertainty-domain")
for(k in c(-1,-.3,-.25))add("GPA","MoM-missing-fourth-moment",shape=k,n=100,quantity="CovarianceDefined",value=0,oracle="Hosking-Wallis-fourth-moment-domain",status="outside-uncertainty-domain")
result <- do.call(rbind,rows)
write.table(result,out,sep=",",row.names=FALSE,col.names=TRUE,quote=TRUE,na="NaN",eol="\n")
cat(R.version.string,"\nRows:",nrow(result),"\n")
cat("Gamma finite derivative rows:",sum(result$family=="Gamma" & result$quantity=="GradientShape" & result$status=="finite"),"\n")
cat("Gamma unresolved underflow rows:",sum(result$status=="unresolved-quantile-underflow"),"\n")
