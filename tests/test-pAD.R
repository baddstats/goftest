# pAD warns when input cannot be parsed as numeric
options(warn = 2) # turn warning into error
observed_warning <- NULL
observed_warning <- tryCatch(
  cvm.test(letters, "pnorm", mean = 0, sd = 1),
  error = function(e) {
    observed_warning <- 1
  }
)
if (is.null(observed_warning)) stop("pAD input validation failed")
options(warn = 0)

# pAD returns a similar value compared to an empirical distribution
# This dataset contains values retrieved from original table
# P. A. Lewis’s 1961 paper, "Distribution of the Anderson-Darling Statistic" (Annals of Mathematical Statistics, 32, 1118–1124)
empirical_cdf <- read.csv("anderson_darling_empirical_cdf.csv")

result_n8_z0.375 <- pAD(0.375, n = 8)
empirical_n8_z0.375 <- empirical_cdf[empirical_cdf[["z"]] == 0.375, "n8"]
if (round(result_n8_z0.375, digits = 3) != empirical_n8_z0.375) {
  stop("reference table and empirical distribution have different value")
}

result_n7_z0.425 <- pAD(0.425, n = 7)
empirical_n7_z0.425 <- empirical_cdf[empirical_cdf[["z"]] == 0.425, "n7"]
if (round(result_n7_z0.425, digits = 3) != empirical_n7_z0.425) {
  stop("reference table and empirical distribution have different value")
}

# pAD returns 1 for q = Inf
if (pAD(Inf) != 1) {
  stop("probability was different than expected when using infinite quantile")
}

# pAD works with default values
result <- pAD(1.5)
if (result < 0 || result > 1) stop("pAD returning wrong probability value")

# pAD works with n = Inf and fast = FALSE (exact algorithm)
result_fast   <- pAD(1.5, fast = TRUE)
result_slow   <- pAD(1.5, fast = FALSE)
if (result_slow < 0 || result_slow > 1) stop("pAD returning wrong probability value")

expected_tolerance <- 1e-4
observed_tolerance <- abs(result_fast - result_slow) / result_slow
if (observed_tolerance > expected_tolerance) stop ("difference between result slow and fast is larger than expected")

# pAD lower.tail = FALSE returns the survival probability
if (pAD(1.5, lower.tail = FALSE) != (1 - pAD(1.5))) stop ("lower.tail now working as expected")


# qAD returns a similar quantile compared to an empirical distribution
# This dataset contains values retrieved from original table
# P. A. Lewis’s 1961 paper, "Distribution of the Anderson-Darling Statistic" (Annals of Mathematical Statistics, 32, 1118–1124)
empirical_quantiles <- read.csv("anderson_darling_empirical_quantile.csv")

result_n2_z0.90 <- qAD(0.90, 2)
empirical_n2_z0.90 <- empirical_quantiles[empirical_quantiles[["z"]] == 0.90, "n2"]
if (round(result_n2_z0.90, digits = 2) != round(empirical_n2_z0.90, digits = 2)) {
  stop("quantile distribution different than the expected values")
}

result_n5_z0.95 <- qAD(0.95, 5)
empirical_n5_z0.95 <- empirical_quantiles[empirical_quantiles[["z"]] == 0.95, "n5"]
if (round(result_n5_z0.95, digits = 2) != round(empirical_n5_z0.95, digits = 2)) {
stop("quantile distribution different than the expected values")
}

# qAD lower tail works
if (qAD(0.1, lower.tail = FALSE) != qAD(0.9)) {
  stop ("qAD lower.tail returning non expected values")
}
