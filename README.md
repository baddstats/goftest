# goftest

This repository contains the sources for the
contributed `R` package `goftest`, which performs
goodness-of-fit tests.

The package implements the classical Cramer-von Mises
and Anderson-Darling tests of goodness-of-fit
for continuous univariate distributions, using modern
algorithms to compute the null distributions.

The package was written by Adrian Baddeley, Julian Faraway,
John Marsaglia and George Marsaglia.

This package is free open source under the GNU Public Licence GPL 2|3.

## Installation

The current official release of `goftest` is available
on [CRAN](http://cran.r-project.org/web/packages/goftest)
and can be downloaded and installed automatically
using the R command `install.packages`. 

The code in this repository is the development version,
which may be newer than the official release.
The easiest way to install the development version of `goftest` 
from github is through the `remotes` package:

```R
require(remotes)
install_github('baddstats/goftest')
```

If you don't have `remotes` installed you should first run

```R
install.packages('remotes')
```

## Bug reports 

Users of `goftest` are encouraged to report bugs here 
(go to *issues* in the menu above, 
and press *new issue* to start a new bug report
or feature request).



