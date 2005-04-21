\name{locfdr}

\alias{locfdr}

\title{Computation of Local False Discovery Rates}

\description{
  Compute local false discovery rates, following the definitions and
  description in Efron (2004) JASA, Volume 99, pages 96--104 and
  Efron, B (2005) "Local false discovery rates"
  \url{http://www-stat.stanford.edu/~brad/papers/}. 
}

\usage{
locfdr(zz, bre = 120, df = 7, pct = 1/1000, pct0 = 1/3, nulltype = 1,
type = 0, plot = 1)
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
    \eqn{f(z)}; df=7 by default. Larger values of df may be required if
    \eqn{f(z)} has sharp bends or other irregularities. A warning is
    issued if the fitted curve does not adequately match the histogram
    counts. It is a good idea to use the plot option to view the
    histogram and fitted curve.
  }
  
  \item{pct}{Excluded tail proportions of zz's when fitting \eqn{f(z)};
    pct=\eqn{1/1000} by default; pct=0 includes full range of zz's; pct
    can also be a 2-vector, describing the fitting range.
  }
  
  \item{pct0}{Included proportion of zz distribution  used in fitting
    null density \eqn{f0(z)} is range [pct0, 1-pct0]; default pct0=1/3;
    pct0 can be a 2-vector, eg pct0=c(.25,.60).
  }
  
  \item{nulltype}{Type of null hypothesis assumed in estimating
    \eqn{f0(z)}; 0 is theoretical null \eqn{N(0,1)} [which assumes that the
    original zz scores have been scaled to have a \eqn{N(0,1)} distribution
    under the null hypothesis]; 1 is the empirical null [which assumes a
    \eqn{N(a,b)} null hypothesis, with \eqn{a=zmax} and \eqn{b=sig^2}
    estimated from the central part of the \eqn{f(z)} fit]; 2 is a "split
    normal" version of 1 in which the \eqn{f0(z)} is allowed to have
    different scales on the two sides of the maximum.] The default is
    nulltype=1. Note that the output includes most of the results for both
    nulltype=0 and nulltype=1 or 2 no matter which choice is made here, the
    exception being the indiviual results "fdr" below.
  }

  \item{type}{Type of fitting used for \eqn{f(z)}; 0 is a natural
    spline, 1 is a polynomial, in either case with degrees of freedom df
    [so total degrees of freedom df+1 including the intercept.] The
    default is type=0. 
  }
  
  \item{plot}{Number of plots desired; plot=1 gives single plot showing
    histogram of zz and fitted densities \eqn{f(z)} and \eqn{f0(z)};
    colored
    histogram bars indicate non-null "thinned counts", see Section 5,
    2nd reference above; square dots on the x-axis indicate threshold
    z-values for fdr <= .2. plot=2 also gives plot of fdr, and the right
    and left tail area Fdr curves; plot=3 gives instead the f1 cdf of
    the estimated fdr curve, as in figure 6 of the second reference above.
  }

}
\value{
  A list with three components.
  \item{fdr} the estimated local false discovery rates for each case,
  using the selected options for type and nulltype.

  \item{fp0} the estimated mean and variance(s) for \eqn{f0(z)} assuming 
  nulltype 1 or 2, the estimated proportion p0 of
  null cases, and also p0theo, the p0 estimate for nulltype=0,
  labelled "p0theo".( The "theo" suffix indicates nulltype=0
  results below also.)

  \item{Efdr} the expected false discovery rate for the Non-Null cases,
  a measure of the experiment's power as described in Section 5
  of the second reference above. Large values of Efdr, say Efdr>.4,
  indicate low power. Overall Efdr and right and left values are
  given, both for nulltype =1 or 2 and for nulltype=0.

  \item{cdf1} a 2x99 matrix giving the estimated cdf of fdr under the
  non-null distribution f1. Large values of the cdf for small fdr
  values indicate good power; see Section 5 of the second reference
  above.
  
  \item{mat} A matrix summarizing the estimates of \eqn{f(z)},
  \eqn{f0(z)}, \eqn{fdr(z)}, etc. at the midpoints "z." of the break
  discretization. These are convenient for comparisons and plotting; mat
  includes fdr from nulltype 1 or 2 as specified, estimates of the usual
  tail-area False Discovery Rates, Fdrleft and Fdrright, and also
  fdrtheo and f0theo. Notice that fp0 and mat contain the information
  for nulltype 0 and either nulltype 1 or 2, no matter which nulltype
  has been selected. The choice of nulltype does affect the 10th column
  of mat, "lfdrse", an estimate of standard error for the curve
  \eqn{log(fdr)}. The 11th column of mat is an estimate "f1" of density
  for the non-null z-scores. Column "counts" gives the histogram counts
  for zz. 

}

\details{
  The standard error estimate lfdrse assumes independence of the
  zz values, and should usually be considered as a lower bound on
  the true standard errors.
  
  The density estimates f, f0 , f0theo are scaled to add up to
  approximately the number of zz's. The non-null density \eqn{f1} is scaled
  to add up to approximately (1-p0) times the number of zz's.
}

\references{Efron, B. (2004) \emph{Large-scale simultaneous hypothesis
    testing: the choice of a null hypothesis}, JASA, Vol. 99, pp 96-104

  Efron, B. (2005). \emph{Local False Discovery Rates},
  \url{http://www-stat.stanford.edu/~brad/papers/} 
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
 