\name{locfdr}

\alias{locfdr}

\title{Local False Discovery Rate Calculation}

\description{
  Compute local false discovery rates, following the definitions and
  description in Efron (2004) JASA, Volume 99, pages 96--104 and
  Efron, B (2005) "Local false discovery rates" and Efron, B. (2005)
  "Correlation and large-scale simultaneous significance testing" 
  \url{http://www-stat.stanford.edu/~brad/papers/}. 
}

\usage{
locfdr(zz, bre=120, df=7, pct=0, pct0=1/4, nulltype=1, type=0, plot=1,
	mult, main=" ", sw=0)
}

\arguments{
  \item{zz}{A vector of summary statistics, one for each case under
  simultaneous consideration. In a microarray experiment there would be
  one component of \eqn{zz} for each gene, perhaps a
  \eqn{t}-statistic comparing gene expression levels under two
  different conditions. Results may be improved by transforming
  \eqn{zz} so that its components are theoretically distributed
  as \eqn{N(0,1)} under the null hypothesis, for example via
  \code{z[i] = qnorm(pt(t[i],df))} when using
  \eqn{t}-statistics. This is especially important when the
  theoretical null option is invoked (see below). Recentering
  and rescaling \eqn{zz} may be necessary if its central
  histogram looks very far removed from mean 0 and variance 1.
  The calculations  assume a large number of cases, say at least
  \code{length(zz)} exceeding 200.}

\item{bre} {Number of breaks in the discretization of the \eqn{z}-score axis,
        or a vector of breakpoints fully describing the
        discretization.  If \code{length(zz)} is small, such as when the
        number of cases is less than about 1000, set bre to a number
        lower than the default of 120.}

\item{df} {Degrees of freedom for fitting the estimated
        density \eqn{f(z)}. Larger values of df may be required if
        \eqn{f(z)} has sharp bends or other irregularities. A warning
        is issued if the fitted curve does not adequately match the
        histogram counts. It is a good idea to use the plot option to
        view the histogram and fitted curve.}

\item{pct} {Excluded tail proportions of \eqn{zz}'s when fitting
	   \eqn{f(z)}. \code{pct=0} includes full range of \eqn{zz}'s. pct can
	   also be a 2-vector, describing the fitting range.}

\item{pct0} {Proportion of the \eqn{zz} distribution used in fitting the
	    null density \eqn{f0(z)} by central matching.  If a 2-vector,
	    e.g. \code{pct0=c(0.25,0.60)}, the range [pct0[1],
	    pct0[2]] is used.  If a scalar, [pct0, 1-pct0] is used.}

\item{nulltype} {Type of null hypothesis assumed in estimating \eqn{f0(z)},
	for use in the fdr calculations; 0 is 
	theoretical null \eqn{N(0,1)} [which assumes that the original \eqn{zz}
	scores have been scaled to have a \eqn{N(0,1)} distribution under
	the null hypothesis]; 1 (the default) is the empirical null
	with parameters estimated by maximum likelihood; 2 is the
        empirical null
	with parameters estimated by central matching (see
	second reference); 3 is a "split normal" version of 2, in
	which \eqn{f0(z)} is allowed to have different
	scales on the two sides of the maximum.  Unless \eqn{sw==2} or
	\eqn{3}, the theoretical, maximum
	likelihood, and central matching estimates all will be output in
	the matrix fp0, and both the theoretical and the specified
	nulltype will be used in the calculations output in mat, but
	only the specified nulltype is used in the calculation of the
	output fdr (local fdr estimates for every case).}

\item{type} {Type of fitting used for \eqn{f(z)}; 0 is a natural spline, 1 is
	a polynomial, in either case with degrees of freedom df [so
	total degrees of freedom including the intercept is \code{df+1}.]}

\item{plot} {Plots desired. plot=0 gives no plots. plot=1 gives single
	    plot showing the histogram of \eqn{zz} and fitted
	    densities \eqn{f(z)} and \eqn{f0(z)}; colored histogram
	    bars indicate estimated non-null counts; yellow triangles
	    on the x-axis indicate threshold z-values for \code{fdr <=
	    0.2}. plot=2 also gives plot of fdr, and the right and
	    left tail area Fdr curves; plot=3 gives instead the f1 cdf
	    of the estimated fdr curve, as in figure 4 of the second
	    reference; plot=4 gives all three plots.}

\item{mult} {Optional scalar multiple (or vector of multiples) of the
	    sample size for calculation of the corresponding
	    hypothetical Efdr value(s).}

\item{main} {Main heading for the histogram plot when \code{plot>0}.}

\item{sw} {Determines the type of output desired.  \eqn{sw=2} gives a list
	  consisting of the last 5 values listed below. \eqn{sw=3} gives
	  the square matrix of dimension bre-1 representing the influence
	  function of log(fdr), i.e. the derivative of log(fdr) (for
	  each bin) with respect to the bin counts.  Any other value
	  of sw returns a list consisting of the first 5 (6 if mult is
	  supplied) values listed below.}
}

\details{
	The standard error estimate lfdrse assumes independence of the
	\eqn{zz} values and should usually be considered as a lower bound on
	the true standard errors. See the third reference.

	The density estimates f, f0, f0theo are scaled to add up to
	approximately the number of \eqn{zz}'s. The non-null density
	f1 is scaled to add up to approximately (1-p0) times the number of
	\eqn{zz}'s. i.e. the estimated number of non-null \eqn{zz}'s.	
}

\value{
\item{fdr} {the estimated local false discovery rate for each case,
	using the selected options for type and nulltype.}

\item{fp0} {the estimated parameters delta (mean of f0), sigma
	(standard deviation of f0), and p0, along with
	their standard errors.  If \code{nulltype<3}, fp0 is a 5 by 3 matrix,
	with columns representing delta, sigma, and p0 and rows
	representing nulltypes and estimate vs. standard error.  If
	\code{nulltype==3}, a fourth column represents the sigma estimate for
	the right side of f0.}

\item{Efdr} {the expected false discovery rate for the non-null cases,
	a measure of the experiment's power as described in Section 3
	of the second reference. Large values of Efdr, say \code{Efdr>0.4},
	indicate low power. Overall Efdr and right and left values are
	given, both for the specified nulltype and for nulltype 0.  If
	\code{nulltype==0}, values are given for nulltypes 1 and 0.}

\item{cdf1} {a 99x2 matrix giving the estimated cdf of fdr under the
	non-null distribution f1. Large values of the cdf for small fdr
	values indicate good power; see Section 3 of the second
	reference.  Set plot to 3 or 4 to see the cdf plot.}

\item{mat} {A matrix summarizing the estimates of \eqn{f(z)},
	\eqn{f0(z)}, \eqn{fdr(z)}, etc. at the \eqn{bre-1} midpoints "x"
	of the break discretization. These are convenient for
	comparisons and plotting; mat includes fdr from nulltype 1, 2,
	or 3 as specified, estimates of the usual tail-area False
	Discovery Rates, Fdrleft and Fdrright, and also fdrtheo and
	f0theo, the fdr and f0 estimates assuming the theoretical null
	density \eqn{N(0,1)}. If \code{nulltype==0}, the fdr and f0
	columns of mat are calculated using nulltype 1. The 10th column
	of mat, "lfdrse", is an estimate of standard error for the curve
	log(fdr) and is calculated based on the specified nulltype. The
	11th column of mat is an estimate p1f1 of the subdensity for
	the non-null \eqn{z}-scores. Column "counts" gives the histogram
	counts for \eqn{zz}.}

\item{mult} {If the argument mult was supplied, vector of the ratios of
	    hypothetical Efdr for the supplied multiples of the sample
	    size to Efdr for the actual sample size.}

\item{pds} {The estimates of p0, delta, and sigma.}

\item{x} {The bin midpoints.}

\item{f} {The values of \eqn{f(z)} at the bin midpoints.}

\item{pds.} {The derivative of the estimates of p0 (when
	    \code{nulltype==1}) or \code{log(p0)} (when
	    \code{nulltype==0} or \code{2}), delta, and sigma with
	    respect to the bin counts.}

\item{stdev} {The delta-method estimates of the standard deviations of
	     the p0, delta, and sigma estimates.}
}

\references{
  Efron, B. (2004) "Large-scale simultaneous hypothesis testing: the choice of
  a null hypothesis", Jour Amer Stat Assoc, \bold{99}, pp. 96--104
  
  Efron, B. (2006) "Size, Power, and False Discovery Rates"
  
  Efron, B. (2006) "Correlation and Large-Scale Simultaneous Significance
  Testing"
  
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
	     