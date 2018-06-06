##
## andarl.R
##
##  Anderson-Darling test and null distribution
##
## $Revision: 1.10 $ $Date: 2018/06/06 08:25:51 $
##

ad.test <- function(x, null="punif", ..., estimated=FALSE, nullname) {
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
  n <- length(x)
  F0 <- getCdf(null)
  U <- F0(x, ...)
  if(any(U < 0 | U > 1))
    stop("null distribution function returned values outside [0,1]")
  if(!estimated || n <= 4) {
    #' simple null hypothesis
    z <- do.goftest.AD(U)
    PVAL <- z$pvalue
    STATISTIC <- z$adstat
    names(STATISTIC) <- "An"
  } else {
    #' composite - use Braun (1980)
    first <- sample(n, ceiling(n/2), replace=TRUE)
    z1 <- do.goftest.AD(U[first])
    z2 <- do.goftest.AD(U[-first])
    STATISTIC <- max(z1$adstat, z2$adstat)
    names(STATISTIC) <- "AnMax"
    PVAL <- 1 - (1 - min(z1$pvalue, z2$pvalue))^2
  }
  METHOD <- c("Anderson-Darling test of goodness-of-fit",
              if(estimated) "(with Braun's adjustment)" else NULL,
              paste("Null hypothesis:", nullname))

  extras <- list(...)
  parnames <- intersect(names(extras), names(formals(F0)))
  if(length(parnames) > 0) {
    pars <- extras[parnames]
    pard <- character(length(parnames))
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


do.goftest.AD <- function(U) {
  ## Internal: call Marsaglia C code
  U <- sort(U)
  n <- length(U)
  z <- .C(CgofADtestR,
          x = as.double(U),
          n = as.integer(n),
          adstat = as.double(numeric(1)),
          pvalue = as.double(numeric(1)),
	  PACKAGE="goftest"
          )
  return(list(adstat=z$adstat, pvalue=z$pvalue))
}

pAD <- function(q, n=Inf, lower.tail=TRUE, fast=TRUE) {
  q <- as.numeric(q)
  p <- rep(NA_real_, length(q))
  if(any(ones <- is.infinite(q) & (q == Inf)))
    p[ones] <- 1
  if(any(zeroes <- (is.finite(q) & q <= 0) | (is.infinite(q) & (q == -Inf))))
    p[zeroes] <- 0
  ok <- is.finite(q) & (q > 0)
  nok <- sum(ok)
  if(nok > 0) {
    if(is.finite(n)) {
      z <- .C(CgofADprobN,
              a       = as.double(q[ok]),
              na      = as.integer(nok),
              nsample = as.integer(n),
              prob    = as.double(numeric(nok)),
	      PACKAGE="goftest")
      p[ok] <- z$prob
    } else if(fast) {
      ## fast version adinf()
      z <- .C(CgofADprobApproxInf,
              a    = as.double(q[ok]),
              na   = as.integer(nok),
              prob = as.double(numeric(nok)),
	      PACKAGE="goftest")
      p[ok] <- z$prob
    } else {
      ## slow, accurate version ADinf()
      z <- .C(CgofADprobExactInf,
              a    = as.double(q[ok]),
              na   = as.integer(nok),
              prob = as.double(numeric(nok)),
	      PACKAGE="goftest")
      p[ok] <- z$prob
    }
      
  }
  if(!lower.tail)
    p <- 1 - p
  return(p)
}

qAD <- local({

  f <- function(x, N, P, Fast) {
    pAD(x, N, fast=Fast) - P
  }
    
  qAD <- function(p, n=Inf, lower.tail=TRUE, fast=TRUE) {
    ## quantiles of null distribution of Anderson-Darling test statistic
    stopifnot(all(p >= 0))
    stopifnot(all(p <= 1))
    if(!lower.tail) p <- 1-p
    ans <- rep(NA_real_, length(p))
    for(i in which(p >= 0 & p < 1)) 
      ans[i] <- uniroot(f, c(0, 1), N=n, P=p[i], Fast=fast, extendInt="up")$root
    return(ans)
  }

  qAD
})


  

