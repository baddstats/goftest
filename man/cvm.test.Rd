\name{cvm.test}
\alias{cvm.test}
\title{
  Cramer-Von Mises Test of Goodness-of-Fit
}
\description{
  Performs the
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test
  of goodness-of-fit to a specified continuous univariate
  probability distribution.
}
\usage{
cvm.test(x, null = "punif", ..., nullname)
}
\arguments{
  \item{x}{
    Numeric vector of data values.
  }
  \item{null}{
    A function, or a character string giving the name of a function,
    to compute the cumulative distribution function for the
    null distribution.
  }
  \item{\dots}{
    Additional arguments for the cumulative distribution function.
  }
  \item{nullname}{
    Optional character string describing the null distribution.
    The default is \code{"uniform distribution"}.
  }
}
\details{
  This command performs the
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test
  of goodness-of-fit to the distribution specified by the argument
  \code{null}. It is assumed that the values in \code{x} are
  independent and identically distributed random values, with some
  cumulative distribution function \eqn{F}.
  The null hypothesis is that \eqn{F} is the function
  specified by the argument \code{null}, while the alternative
  hypothesis is that \eqn{F} is some other function.
}
\value{
  An object of class \code{"htest"} representing the result of
  the hypothesis test.
}
\references{
  \ifelse{latex}{\out{Cs\"org\H{o}}}{Csorgo}, S. and Faraway, J.J. (1996)
  The exact and asymptotic distributions of
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises statistics.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{58}, 221--234.
}
\author{
  Adrian Baddeley.
}
\seealso{
  \code{\link{pCvM}} for the null distribution of the test statistic.
}
\examples{
x <- rnorm(10, mean=2, sd=1)
cvm.test(x, "pnorm", mean=2, sd=1)
}
\keyword{htest}

