# Offline oracle generator for Hartigan's dip statistic and dip test.
# Requires R 4.4.3 with the CRAN package diptest 0.77-2 (GPL; called as a black box only,
# no code or tables are copied into Numerics). Run from the repository root; never calls
# Numerics or changes production or test files.
options(digits = 17, warn = 1)
stopifnot(as.character(getRversion()) == "4.4.3")
stopifnot(as.character(packageVersion("diptest")) == "0.77.2")
suppressPackageStartupMessages(library(diptest))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "docs/statistics/oracles/dip-test/dip-test.csv"

fmt <- function(x) sprintf("%.17g", x)
rows <- list()
add <- function(case, x, quantity, value, oracle) {
  rows[[length(rows) + 1L]] <<- data.frame(case = case, n = length(x),
    values = paste(fmt(x), collapse = " "), quantity = quantity, value = fmt(value),
    oracle = oracle, stringsAsFactors = FALSE)
}
dipcase <- function(case, x) add(case, x, "dip", dip(x), "diptest::dip")

# Statistic: small, tied, constant and empty samples (Test_HartiganDip_MatchesR)
dipcase("empty", numeric(0))
dipcase("single", 5)
dipcase("two", c(1, 2))
dipcase("four", c(4, 1, 3, 2))
dipcase("constant5", rep(3, 5))
dipcase("ties10", c(1, 1, 1, 2, 2, 5, 5, 5, 5, 9))
dipcase("ties15", c(0.5, 1.5, 1.5, 2, 2, 2, 2.5, 3, 7, 7.5, 8, 8, 8.5, 9, 9))
dipcase("halfsteps20", c(1, 0.5, -1, 0.5, -0.5, -1, 1, 1, -1, 0.5, 2.5, 1, -0.5, 0.5, -1.5, 0.5, -2, 2.5, -1, 1))

# A randomized two-component sample of 20 values (seed 20261105 of the validation set)
mixture20 <- c(-0.2600038160874425, -0.16285100167438976, -0.90087458650539542, -0.23135597425869309,
  0.86195506130218913, 0.3763467681724586, -0.17715607751019694, -0.085177164177269976,
  -0.43612827351821809, 1.4237823458859504, 3.416678901465473, 3.8604758569957012,
  3.6775664316148138, 3.0153586596756199, 2.0927864607817059, 1.2240458252098512,
  3.6884804177102644, 2.4817555475005166, 1.5391226578204371, 3.8898272439175168)
dipcase("mixture20", mixture20)

# Bimodal fixture: statistic, table p-value and simulated p-value (Test_HartiganDipTest_MatchesR)
bimodal30 <- c(-2.1, -1.7, -1.5, -1.2, -1.1, -0.9, -0.6, -0.4, -0.3, -0.1, 0.2, 0.4, 0.6, 0.9, 1.3,
  4.1, 4.4, 4.6, 4.8, 5.0, 5.1, 5.3, 5.6, 5.9, 6.2, 6.5, 6.8, 7.0, 7.4, 8.0)
dipcase("bimodal30", bimodal30)
add("bimodal30", bimodal30, "p_table", dip.test(bimodal30)$p.value, "diptest::dip.test")
set.seed(3)
add("bimodal30", bimodal30, "p_simulated_B100000_seed3",
    dip.test(bimodal30, simulate.p.value = TRUE, B = 100000)$p.value, "diptest::dip.test(simulate)")

# Normal quantiles: minimum dip 1/(2n) and p = 1
quantiles40 <- qnorm(ppoints(40))
dipcase("quantiles40", quantiles40)
add("quantiles40", quantiles40, "p_table", dip.test(quantiles40)$p.value, "diptest::dip.test")

# Three values: p = 1
add("three", c(1, 2, 5), "p_table", dip.test(c(1, 2, 5))$p.value, "diptest::dip.test")

table <- do.call(rbind, rows)
con <- file(out, open = "wb")
write.csv(table, con, row.names = FALSE, quote = TRUE, eol = "\n")
close(con)
cat("wrote", nrow(table), "rows to", out, "\n")

# Optional randomized validation set (not committed):
# Rscript --vanilla dip-test.R <csv> <random-set-path> writes 400 samples (n = 1..2000; normal,
# two-component, rounded, exponential, uniform and rounded mixtures) with their dips.
if (length(args) >= 2) {
  set.seed(20261105)
  con <- file(args[2], "w")
  for (k in 1:400) {
    n <- sample(c(1:12, 15, 20, 30, 50, 100, 200, 500, 2000), 1); type <- k %% 6
    x <- switch(type + 1, rnorm(n), c(rnorm(n %/% 2), rnorm(n - n %/% 2, 3)), round(rnorm(n) * 2) / 2,
                rexp(n), runif(n), round(c(rnorm(n %/% 3, -2), rnorm(n - n %/% 3, 2)), 1))
    writeLines(paste(fmt(dip(x)), paste(fmt(x), collapse = " ")), con)
  }
  close(con)
}
