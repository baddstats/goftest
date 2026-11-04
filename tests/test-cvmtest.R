library(goftest)
# Argument validation for cvm.test
# It throws error if 'x' is not the expected type
observed_error <- NULL
observed_error <- tryCatch(
  cvm.test(letters, "pnorm", mean = 0, sd = 1),
  error = function(e) {
    observed_error <- 1
  }
)
if (is.null(observed_error)) {
  stop("cvm.test argument validation test on `x` failed")
}

# It throws error if 'estimated' is not the expected type
observed_error <- NULL
observed_error <- tryCatch(
  cvm.test(rnorm(10), "pnorm", mean = 0, sd = 1, estimated = "random"),
  error = function(e) {
    observed_error <- 1
  }
)
if (is.null(observed_error)) {
  stop("cvm.test argument validation test on `estimated` failed")
}

# Throws error when CDF returns values outside [0, 1]
bad_cdf <- function(x, ...) rep(2, length(x))
observed_error <- NULL
observed_error <- tryCatch(
  cvm.test(1:10, bad_cdf),
  error = function(e) {
    observed_error <- 1
  }
)
if (is.null(observed_error)) {
  stop("cvm.test argument validation test on CDF values failed")
}

# The cvm.test should not reject the null hypothesis
set.seed(1)
p_value_null_hypothesis <- 0.05

# Does not reject with default parameters
values <- rnorm(1000)
result <- cvm.test(values, "pnorm")
if (p_value_null_hypothesis > result$p.value) {
  stop("Null hypothesis was rejected")
}

# Does not reject with defined parameters in normal distribution
values <- rnorm(1000, mean = 5, sd = 1)
result <- cvm.test(values, "pnorm", mean = 5, sd = 1)
if (p_value_null_hypothesis > result$p.value) {
  stop("Null hypothesis was rejected with defined parameters")
}

# Does not reject with estimated parameters in normal distribution
values <- rnorm(2000, mean = 5, sd = 1)
result <- cvm.test(
  values, "pnorm", mean = mean(values), sd = sd(values), estimated = TRUE
)
if (p_value_null_hypothesis > result$p.value) {
  stop("Null hypothesis was rejected with estimated parameters")
}

# Does not reject with log normal distribution and default parameters
values <- rlnorm(1000)
result <- cvm.test(values, "plnorm")
if (p_value_null_hypothesis > result$p.value) {
  stop("Null hypothesis was rejected with lognormal distribution")
}

# The cvm.test should reject the null hypothesis
set.seed(1)
p_value_null_hypothesis <- 0.05

# Does reject with default parameters with different distributions
values <- rnorm(1000)
result <- cvm.test(values, "plnorm")
if (p_value_null_hypothesis < result$p.value) {
  stop("Null hypothesis was not rejected when testing different distributions")
}

# Does reject with different distributions and defined parameters
values <- rnorm(1000, mean = 5, sd = 10)
result <- cvm.test(values, "plnorm", mean = 5, sd = 10)
if (p_value_null_hypothesis < result$p.value) {
  stop(
    paste(
      "Null hypothesis was not rejected when testing different distributions",
      "with defined parameters"
    )
  )
}
