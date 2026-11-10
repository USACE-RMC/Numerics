# gstat oracle values for Kriging.
#
# Writes kriging-gstat.json with:
#   * the sp::meuse data (x, y, log(zinc)) and the first 10 points of sp::meuse.grid;
#   * the spherical variogram of log(zinc) fitted by fit.variogram, starting from
#     vgm(0.6, "Sph", 900, 0.05) (the fitted parameters are inputs of the tests);
#   * krige() predictions and variances at the grid points for simple (beta = mean of log(zinc)),
#     ordinary (log(zinc) ~ 1) and universal (log(zinc) ~ x + y) kriging, all with the global
#     neighborhood, and krige.cv() leave-one-out predictions and variances for ordinary kriging;
#   * the same for exponential and Gaussian variogram models with fixed parameters
#     (ordinary kriging only);
#   * the ordinary kriging example of Isaaks and Srivastava (1989, chapter 12): seven samples, the
#     covariance C(h) = 10 exp(-0.3 h), estimate at (65, 137).
#
# Usage: Rscript kriging_gstat_oracle.R kriging-gstat.json
suppressMessages({ library(sp); library(gstat) })
args <- commandArgs(trailingOnly = TRUE)
data(meuse); data(meuse.grid)
meuse$lz <- log(meuse$zinc)
coordinates(meuse) <- ~x + y
grid <- meuse.grid[1:10, c("x", "y")]
coordinates(grid) <- ~x + y

v <- fit.variogram(variogram(lz ~ 1, meuse), vgm(0.6, "Sph", 900, 0.05))
sph <- list(nugget = v$psill[1], psill = v$psill[2], range = v$range[2])
run <- function(model, formula, beta = NULL) {
  k <- if (is.null(beta)) krige(formula, meuse, grid, model = model, debug.level = 0)
       else krige(formula, meuse, grid, model = model, beta = beta, debug.level = 0)
  list(mean = k$var1.pred, variance = k$var1.var)
}
cv <- function(model) {
  k <- krige.cv(lz ~ 1, meuse, model = model, verbose = FALSE)
  list(mean = k$var1.pred, variance = k$var1.var)
}
fixed <- function(type, psill, range, nugget) {
  m <- vgm(psill, type, range, nugget)
  list(model = type, nugget = nugget, psill = psill, range = range,
       ordinary = run(m, lz ~ 1), crossValidation = cv(m))
}

pts <- data.frame(x = c(61, 63, 64, 68, 71, 73, 75), y = c(139, 140, 129, 128, 140, 141, 128),
                  v = c(477, 696, 227, 646, 606, 791, 783))
coordinates(pts) <- ~x + y
target <- data.frame(x = 65, y = 137); coordinates(target) <- ~x + y
isk <- krige(v ~ 1, pts, target, model = vgm(10, "Exp", 1 / 0.3, 0), debug.level = 0)

out <- list(
  R = as.character(getRversion()), gstat = as.character(packageVersion("gstat")),
  sp = as.character(packageVersion("sp")),
  meuse = list(x = coordinates(meuse)[, 1], y = coordinates(meuse)[, 2], logZinc = meuse$lz),
  grid = list(x = coordinates(grid)[, 1], y = coordinates(grid)[, 2]),
  spherical = list(nugget = sph$nugget, psill = sph$psill, range = sph$range,
                   simpleMean = mean(meuse$lz),
                   simple = run(v, lz ~ 1, beta = mean(meuse$lz)),
                   ordinary = run(v, lz ~ 1),
                   universal = run(v, lz ~ x + y),
                   crossValidation = cv(v)),
  exponential = fixed("Exp", 0.6, 300, 0.05),
  gaussian = fixed("Gau", 0.6, 500, 0.05),
  isaaksSrivastava = list(mean = isk$var1.pred, variance = isk$var1.var))
writeLines(jsonlite::toJSON(out, digits = NA, auto_unbox = TRUE, pretty = TRUE), args[1], useBytes = TRUE)
