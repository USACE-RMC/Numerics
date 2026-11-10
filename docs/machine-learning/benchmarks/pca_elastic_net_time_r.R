# Times R prcomp and glmnet on the benchmark data written by pca_elastic_net_generate_data.py.
# Usage:  Rscript pca_elastic_net_time_r.R <data-folder> [repeats]
# Writes pca_elastic_net_timings_r.json (median seconds per case) with versions and the BLAS and LAPACK libraries.
# Cases, matched to the C# benchmarks:
#   pca_fit     prcomp(X, retx = FALSE)
#   enet_path   glmnet(X, y, alpha = 0.5) with defaults (100 penalties, thresh = 1e-7)
#   enet_cv     cv.glmnet(X, y, alpha = 0.5, nfolds = 10), sequential
suppressPackageStartupMessages(library(glmnet))
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
repeats <- if (length(args) > 1) as.integer(args[2]) else 7L
txt <- paste(readLines(file.path(root, "sizes.json"), warn = FALSE), collapse = "")
grab <- function(key) {
  body <- sub(paste0('.*"', key, '": \\[(\\[[^"]*?\\])\\].*'), "\\1", txt, perl = TRUE)
  v <- as.integer(unlist(regmatches(body, gregexpr("[0-9]+", body))))
  matrix(v, ncol = 2, byrow = TRUE)
}
read_bin <- function(name, n) readBin(file.path(root, name), "double", n = n, size = 8, endian = "little")
# The elapsed-time clock has a resolution of about 10 ms on Windows, so fast cases are repeated until one
# sample takes at least 0.2 s and the per-call time is reported.
median_time <- function(f, reps) {
  f()
  inner <- 1L
  repeat {
    t <- system.time(for (i in seq_len(inner)) f())[["elapsed"]]
    if (t >= 0.2 || inner >= 4096L) break
    inner <- inner * 4L
  }
  t <- numeric(reps)
  for (i in seq_len(reps)) t[i] <- system.time(for (j in seq_len(inner)) f())[["elapsed"]] / inner
  median(t)
}
cases <- list()
pca <- grab("pca")
for (k in seq_len(nrow(pca))) {
  n <- pca[k, 1]; p <- pca[k, 2]
  X <- matrix(read_bin(sprintf("pca_%dx%d.bin", n, p), n * p), nrow = n, byrow = TRUE)
  cases[[sprintf("pca_fit_%dx%d", n, p)]] <- median_time(function() prcomp(X, retx = FALSE), repeats)
  cat("pca", n, p, cases[[sprintf("pca_fit_%dx%d", n, p)]], "\n")
}
enet <- grab("enet")
for (k in seq_len(nrow(enet))) {
  n <- enet[k, 1]; p <- enet[k, 2]
  X <- matrix(read_bin(sprintf("enet_%dx%d_X.bin", n, p), n * p), nrow = n, byrow = TRUE)
  y <- read_bin(sprintf("enet_%dx%d_y.bin", n, p), n)
  reps <- if (n * p < 1e7) repeats else max(3L, repeats %/% 2L)
  cases[[sprintf("enet_path_%dx%d", n, p)]] <- median_time(function() glmnet(X, y, alpha = 0.5), reps)
  cases[[sprintf("enet_cv_%dx%d", n, p)]] <- median_time(function() cv.glmnet(X, y, alpha = 0.5, nfolds = 10), reps)
  cat("enet", n, p, cases[[sprintf("enet_path_%dx%d", n, p)]], cases[[sprintf("enet_cv_%dx%d", n, p)]], "\n")
}
esc <- function(s) gsub("\\\\", "/", s)
# R on Windows ships its own reference BLAS and LAPACK (Rblas.dll, Rlapack.dll) and reports them as "".
blas <- extSoftVersion()[["BLAS"]]
if (!nzchar(blas)) blas <- "R reference BLAS (Rblas)"
lapack <- La_library()
if (!nzchar(lapack)) lapack <- "R reference LAPACK (Rlapack)"
json <- paste0('{\n "r": "', R.version.string, '",\n "glmnet": "', as.character(packageVersion("glmnet")),
  '",\n "blas": "', esc(blas), '",\n "lapack": "', esc(lapack),
  '",\n "cases": {\n', paste(sprintf('  "%s": %.6f', names(cases), unlist(cases)), collapse = ",\n"), '\n }\n}\n')
writeLines(json, file.path(root, "pca_elastic_net_timings_r.json"))
cat(json)
