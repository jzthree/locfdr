\name{locfdr}

\alias{locfdr}

\title{Computation of Local False Discovery Rates}

\description{
  Compute local false discovery rates, following the definitions and
  description in Efron (2004) JASA, Volume 99, pages 96--104.
}

\usage{
locfdr(zz, bre = 120, df = 7, pct = 1/1000, pct0 = 2/3, nulltype = 1, type = 0)
}

\arguments{
  \item{zz}{A vector of summary statistics, one for each case under
    simultaneous consideration. In a microarray experiment there would be
    one component of \eqn{zz} for each gene, perhaps a \eqn{t}-statistic
    comparing gene expression levels under two different conditions. The
    calculations assume a large number of cases, say at least length(zz)
    \eqn{>} 100.  
  }
  
  \item{bre}{Number of breaks in the discretization of the \eqn{z}-score
    axis, set to 120 by default. This can also be a vector of
    breakpoints fully describing the discretization. 
  }
  
  \item{df}{Number degrees of freedom for fitting the estimated density
    \eqn{f(z)}, set to 7 by default.
  }
  
  \item{pct}{Excluded tail proportions of zz's when fitting \eqn{f(z)};
    pct=\eqn{1/1000} by default; pct=0 includes full range of zz's; pct
    can also be a 2-vector, describing the fitting range.
  }
  
  \item{pct0}{Included proportion of range used in fitting null density
    \eqn{f0(z)}; default value is \eqn{2/3}.
  }
  
  \item{nulltype}{Type of null hypothesis assumed in estimating
    \eqn{f0(z)}; 0 is theoretical null \eqn{N(0,1)} [which assumes that
    the original zz scores have been scaled to have a \eqn{N(0,1)}
    distribution under the null hypothesis]; 1 is the empirical null
    [which assumes a \eqn{N(a,b)} null hypothesis, with a=zmax and
    b=sig2 estimated from the central part of the \eqn{f(z)} fit]; 2 is
    a "split normal" version of 1 in which the \eqn{f0(z)} is allowed to
    have different scales on the two sides of the maximum.] The default
    is nulltype=1.
  }
  
  \item{type}{Type of fitting used for \eqn{f(z)}; 0 is a natural
    spline, 1 is a polynomial, in either case with degrees of freedom df
    [so total degrees of freedom df+1 including the intercept.] The
    default is type=0. 
  }
}

\value{
  A list with three components.
  
  \item{fdr}{the estimated local false discovery rates for each case} 
  \item{f0.p0}{The estimated mean and variance(s) for \eqn{f0(z)}
    assuming nulltype 1 or 2, and also the estimated proportion \eqn{p0}
    of null cases}
  
  \item{mat}{A matrix summarizing the estimates of \eqn{f(z)},
    \eqn{f0(z)}, \eqn{fdr(z)}, etc. at the midpoints "\eqn{z}" of the
    break discretization. These are convenient for comparisons and
    plotting; mat includes fdr from nulltype 1 or 2 as specified,
    estimates of the usual tail-area False Discovery Rates, Fdrleft and
    Fdrright, and also fdrtheo and f0theo if nulltype=0 is specified.}
}

\references{Efron, B. (2004) \emph{Large-scale simultaneous hypothesis
    testing: the choice of a null hypothesis}, JASA, Vol. 99, pp 96-104
}
\author{Bradley Efron}

\examples{
## HIV data example
data(hivdata)
w <- locfdr(hivdata)
print(w)

## Second Simulation Example

}

\keyword{htest}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
