# Reference values for PrincipalComponentAnalysis from R prcomp (stats).
# Run from the repository root:  Rscript docs/machine-learning/oracles/pca_prcomp.R
# Reads the random matrix written by pca_oracle.py (run it first) and prints sdev^2 and rotation.
# prcomp does not fix the sign of each axis, so rotation columns and scores agree with
# PrincipalComponentAnalysis up to a sign per component.
iris4 <- read.csv("docs/example-data/iris-dataset.csv", comment.char = "#")[, 1:4]
A <- as.matrix(read.csv("docs/machine-learning/oracles/pca-random-matrix.csv", header = FALSE))
show <- function(label, p) {
  cat(label, "\n")
  cat("sdev^2:", formatC(p$sdev^2, digits = 17, format = "g"), "\n")
  cat("rotation (column-major):", formatC(as.vector(p$rotation), digits = 17, format = "g"), "\n")
  cat("x[1:3,] (column-major):", formatC(as.vector(p$x[1:3, ]), digits = 17, format = "g"), "\n\n")
}
cat("R", R.version.string, "\n\n")
show("iris covariance", prcomp(iris4))
show("iris correlation", prcomp(iris4, scale. = TRUE))
show("random covariance", prcomp(A))
show("random correlation", prcomp(A, scale. = TRUE))
