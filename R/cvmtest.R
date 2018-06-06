##
## cvmtest.R
##
## Cramer-von Mises test
##
## $Revision: 1.7 $ $Date: 2018/06/06 08:25:46 $
##

cvm.test <- function(x, null="punif", ..., estimated=FALSE, nullname) {
  xname <- deparse(substitute(x))
  nulltext <- deparse(substitute(null))
  if(is.character(null)) nulltext <- null
  if(missing(nullname) || is.null(nullname)) {
    reco <- recogniseCdf(nulltext)
    nullname <- if(!is.null(reco)) reco else 
                paste("distribution", sQuote(nulltext))
  }
  stopifnot(is.numeric(x))
  x <- as.vector(x)
  F0 <- getCdf(null)
  U <- F0(x, ...)
  n <- length(U)
  if(any(U < 0 | U > 1))
    stop("null distribution function returned values outside [0,1]")
  if(!estimated || n <= 4) {
    #' simple null hypothesis
    z <- do.goftest.CvM(U)
    PVAL <- z$pvalue
    STATISTIC <- z$omega2
    names(STATISTIC) <- "omega2"
  } else {
    #' composite - use Braun (1980)
    first <- sample(n, ceiling(n/2), replace=TRUE)
    z1 <- do.goftest.CvM(U[first])
    z2 <- do.goftest.CvM(U[-first])
    PVAL <- 1 - (1 - min(z1$pvalue, z2$pvalue))^2
    STATISTIC <- max(z1$omega2, z2$omega2)
    names(STATISTIC) <- "omega2max"
  }
  METHOD <- c("Cramer-von Mises test of goodness-of-fit",
              if(estimated) "(with Braun's adjustment)" else NULL,
              paste("Null hypothesis:", nullname))

  extras <- list(...)
  parnames <- intersect(names(extras), names(formals(F0)))
  if(length(parnames) > 0) {
    pars <- extras[parnames]
    pard <- character(0)
    for(i in seq_along(parnames))
      pard[i] <- paste(parnames[i], "=", paste(pars[[i]], collapse=" "))
    pard <- paste("with",
                  ngettext(length(pard), "parameter", "parameters"),
                  "  ", 
                  paste(pard, collapse=", "))
    METHOD <- c(METHOD, pard)
  }

  coda <- paste("Parameters assumed to",
                if(estimated) "have been estimated from data" else "be fixed")
  METHOD <- c(METHOD, coda)

  out <- list(statistic = STATISTIC,
               p.value = PVAL,
               method = METHOD,
               data.name = xname)
  class(out) <- "htest"
  return(out)
}

#' not exported

do.goftest.CvM <- function(U) {
  U <- sort(U)
  n <- length(U)
  k <- seq_len(n)
  omega2 <- 1/(12 * n) + sum((U - (2*k - 1)/(2*n))^2)
  pvalue <- pCvM(omega2, n=n, lower.tail=FALSE)
  return(list(omega2=omega2, pvalue=pvalue))
}
