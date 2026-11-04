# tests for CvM provability and quantile functions
# returns a valid qCvM quantile"
q <- qCvM(0.9)
if (isFALSE(is.numeric(q) && q > 0)) stop("not a valid quantile")
expected_tolerance <- 1e-3
observed_tolerance <- abs(pCvM(q) - 0.9) / pCvM(q)
if (observed_tolerance > expected_tolerance) stop("quantile different than expected value")

# returns matching qCvM quantiles for lower.tail inversion"
if (qCvM(0.1, lower.tail = FALSE) != qCvM(0.9)) stop("lower.tail qCvM not working as expected")

# The tests below are ported from sciPy repo:
# https://github.com/scipy/scipy/blob/1e9a17b001c24e833eed24622741db8af92dc198/scipy/stats/tests/test_hypotests.py#L108
# the expected values of the cdfs are taken from Table 1 in
# Csorgo / Faraway: The Exact and Asymptotic Distribution of
# Cramér-von Mises Statistics, 1996.
# TestCvm
dat <- list(
  list(n = 4, x = c(0.02983, 0.04111, 0.12331, 0.94251), ref = c(0.01, 0.05, 0.5, 0.999)),
  list(n = 10, x = c(0.02657, 0.03830, 0.12068, 0.56643), ref = c(0.01, 0.05, 0.5, 0.975)),
  list(n = 1000, x = c(0.02481, 0.03658, 0.11889, 1.16120), ref = c(0.01, 0.05, 0.5, 0.999)),
  list(n = Inf, x = c(0.02480, 0.03656, 0.11888, 1.16204), ref = c(0.01, 0.05, 0.5, 0.999))
)

expected_tolerance <- 1e-3
lapply(dat, function(d) {
  observed_tolerance <- abs(pCvM(d$x, n = d$n) - d$ref) / pCvM(d$x, n = d$n)
  if (any(observed_tolerance > expected_tolerance)) stop("probability different than expected value")
})

# test_cdf_support: CDF has support on [1/(12*n), n/3]
# returns 0 and 1 at the boundaries of the CDF support
n <- 533
if (pCvM(1 / (12 * n), n = n) != 0) stop("pCvM not within boundaries")
if (pCvM(n / 3, n = n) != 1) stop("pCvM not within expected boundaries")

n <- 27
if (pCvM(1 / (12 * (n + 1)), n = n) != 0) stop("pCvM not within expected boundaries")
if (pCvM((n + 1) / 3, n = n) != 1) stop("pCvM not within expected boundaries")

# test_cdf_large_n: asymptotic CDF and CDF for large samples should be close
# asymptotic CDF closely matches CDF for large n
x <- c(0.02480, 0.03656, 0.11888, 1.16204, 100)
expected_tolerance <- 1e-3
observed_tolerance <- abs(pCvM(x, n = 10000) - pCvM(x, n = Inf)) / pCvM(x, n = 10000)
if (any(observed_tolerance > expected_tolerance)) stop("pCvM not matching large samples")

# test_large_x: for large x and n, the series converges slowly;
# cdf should be in (0.99999, 1] for x = 333.3 and n = 1000
# returns a value just below 1 for large x near the upper support boundary
x <- 333.3
p_finite <- pCvM(x, n = 1000)
p_asymp  <- pCvM(x, n = Inf)
if (p_finite < 0.99999) stop ("value not closer to 1")
if (p_asymp < 0.99990) stop ("value not closer to 1")
