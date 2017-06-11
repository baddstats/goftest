\name{goftest-package}
\alias{goftest-package}
\alias{goftest}
\docType{package}
\title{
  Classical Goodness-of-Fit Tests
}
\description{
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises
  and Anderson-Darling tests of goodness-of-fit
  for continuous univariate distributions, using modern
  algorithms to compute the null distributions.
}
\details{
  The \pkg{goftest} package contains implementations of the
  classical \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises
  and Anderson-Darling tests of goodness-of-fit
  for continuous univariate distributions.

  The \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test
  is performed by \code{\link{cvm.test}}. The cumulative distribution
  function of the null distribution of the test statistic
  is computed by \code{\link{pCvM}}
  using the algorithm of \ifelse{latex}{\out{Cs\"org\H{o}}}{Csorgo}
  and Faraway (1996). The quantiles are computed by \code{\link{qCvM}}
  by root-finding.

  The Anderson-Darling test is performed by 
  \code{\link{ad.test}}. The cumulative distribution
  function of the null distribution of the test statistic
  is computed by \code{\link{pAD}}
  using the algorithm of Marsaglia and Marsaglia (2004).
  The quantiles are computed by \code{\link{qAD}} by root-finding.
}
\author{
  Adrian Baddeley, Julian Faraway, John Marsaglia, George Marsaglia.

  Maintainer: Adrian Baddeley <adrian.baddeley@uwa.edu.au>
}
\references{
  \ifelse{latex}{\out{Cs\"org\H{o}}}{Csorgo}, S. and Faraway, J.J. (1996)
  The exact and asymptotic distributions of
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises statistics.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{58}, 221--234.

  Marsaglia, G. and Marsaglia, J. (2004)
  Evaluating the Anderson-Darling Distribution.
  \emph{Journal of Statistical Software} \bold{9} (2), 1--5.
  February 2004. 
  \url{http://www.jstatsoft.org/v09/i02}
}
\keyword{package}
\keyword{htest}
\seealso{
  \code{\link[stats]{ks.test}}
}
\examples{
  x <- rnorm(10, mean=2, sd=1)
  cvm.test(x, "pnorm", mean=2, sd=1)
  ad.test(x, "pnorm", mean=2, sd=1)
}
