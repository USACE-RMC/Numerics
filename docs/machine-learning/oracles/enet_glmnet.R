# Reference values for ElasticNetRegression from R glmnet.
# Run from the repository root:  Rscript docs/machine-learning/oracles/enet_glmnet.R
# Writes docs/machine-learning/oracles/enet-glmnet.json.
#
# Data: the prostate cancer training set (67 cases, 8 predictors, response lpsa) from
# docs/example-data/prostate-cancer.csv (Hastie, Tibshirani and Friedman, 2009, Section 3.2.1).
# glmnet stops when the largest weighted squared coefficient change is below thresh, so thresh = 1e-24
# (changes below 1e-12) gives reference solutions converged far beyond the test tolerances.
# cv.glmnet uses a fixed foldid so the folds are reproducible outside R.
suppressPackageStartupMessages({ library(glmnet); library(jsonlite) })
d <- read.csv("docs/example-data/prostate-cancer.csv", comment.char = "#")
train <- d[d$train == 1, ]
X <- as.matrix(train[, 1:8])
y <- train$lpsa
foldid <- rep(1:5, length.out = nrow(X))

path <- function(fit) {
  list(lambda = fit$lambda, a0 = unname(fit$a0), beta = unname(as.matrix(fit$beta)),
       dev_ratio = fit$dev.ratio, df = fit$df)
}
out <- list(
  r = R.version.string,
  glmnet = as.character(packageVersion("glmnet")),
  foldid = foldid,
  lasso = path(glmnet(X, y, alpha = 1, thresh = 1e-24)),
  enet = path(glmnet(X, y, alpha = 0.5, thresh = 1e-24)),
  ridge = path(glmnet(X, y, alpha = 0, thresh = 1e-24)),
  enet_no_intercept = path(glmnet(X, y, alpha = 0.5, intercept = FALSE, thresh = 1e-24)),
  enet_unstandardized = path(glmnet(X, y, alpha = 0.5, standardize = FALSE, thresh = 1e-24)),
  enet_penalty_factor = path(glmnet(X, y, alpha = 0.5, penalty.factor = c(0, 1, 1, 1, 1, 1, 1, 2), thresh = 1e-24))
)
cv <- cv.glmnet(X, y, alpha = 0.5, foldid = foldid, thresh = 1e-24)
out$cv <- list(lambda = cv$lambda, cvm = cv$cvm, cvsd = cv$cvsd,
               lambda_min = cv$lambda.min, lambda_1se = cv$lambda.1se,
               index_min = cv$index[1], index_1se = cv$index[2])

# cv.glmnet fits each fold on its own penalty sequence and interpolates the predictions to the
# full-data sequence. ElasticNetRegression fits each fold exactly at the full-data penalties, so the
# exact reference refits the folds with lambda = the full-data sequence and applies the cv.glmnet
# fold-weighted mean and standard error.
lam <- cv$lambda
raw <- matrix(0, 5, length(lam))
nf <- tabulate(foldid)
for (f in 1:5) {
  fit <- glmnet(X[foldid != f, ], y[foldid != f], alpha = 0.5, lambda = lam, thresh = 1e-24)
  pred <- predict(fit, X[foldid == f, , drop = FALSE])
  raw[f, ] <- colMeans((y[foldid == f] - pred)^2)
}
cvm <- apply(raw, 2, weighted.mean, w = nf)
cvsd <- sqrt(apply(sweep(raw, 2, cvm)^2, 2, weighted.mean, w = nf) / (5 - 1))
imin <- which(cvm <= min(cvm))[1]
i1se <- which(cvm <= cvm[imin] + cvsd[imin])[1]
out$cv_exact <- list(lambda = lam, cvm = cvm, cvsd = cvsd, index_min = imin, index_1se = i1se)

# Wide case (p > n): the first 6 cases, all 8 predictors.
out$wide <- path(glmnet(X[1:6, ], y[1:6], alpha = 0.5, thresh = 1e-24))

writeLines(toJSON(out, digits = I(17), auto_unbox = TRUE), "docs/machine-learning/oracles/enet-glmnet.json")
cat("lasso path length", length(out$lasso$lambda), "\n")
cat("cv lambda.min", cv$lambda.min, "lambda.1se", cv$lambda.1se, "\n")
