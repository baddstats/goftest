\name{recogniseCdf}
\alias{recogniseCdf}
\title{
  Explanatory Name of Distribution Function
}
\description{
  Recognises many standard cumulative distribution functions and
  returns a string describing the distribution.
}
\usage{
recogniseCdf(s="punif")
}
\arguments{
  \item{s}{
    A single character string giving the \emph{name}
    of an \R function that calculates cumulative
    probabilities.
  }
}
\details{
  The list of recognised distribution functions
  includes all those available in the \pkg{stats} package
  and in \pkg{goftest}.
}
\value{
  Character string, or \code{NULL} if the name is not recognised.
}
\author{
  Adrian Baddeley.
}
\seealso{
  \code{\link{pAD}}
}
\examples{
   recogniseCdf("punif")
   recogniseCdf("pt")
}
\keyword{distribution}
