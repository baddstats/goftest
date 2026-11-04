library(goftest)
# recogniseCdf returns NULL for an unknown distribution name
if (!is.null(recogniseCdf("random"))) stop("recogniseCdf does not return NULL when expected")

# recogniseCdf returns correct distribution function for a known distribution name
if (recogniseCdf("gamma") != "Gamma distribution") stop("incorrect returned distribution")
if (recogniseCdf("binom") != "binomial distribution") stop("incorrect returned distribution")
if (recogniseCdf("cauchy") != "Cauchy distribution") stop("incorrect returned distribution")

# recogniseCdf removes 'p' on the provided string
if (recogniseCdf("pgamma") != "Gamma distribution") stop("incorrect returned distribution")
if (!is.null(recogniseCdf("poison"))) stop ("not processing correctly distribution name")
