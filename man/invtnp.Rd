\name{invtnp}
\alias{invtnp}
\alias{invtnp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{invtnp}
\description{
Performs the inverse lifting transform on a detail and scaling coefficient vector with grid \env{X} and corresponding coefficients \env{coeff}. There is a unique correspondence between the grid values and the function values.
}
\usage{
invtnp(X, coeff, lengths, lengthsremove, pointsin, removelist, 
neighbrs, schemehist, interhist, nadd = length(X) - 2, 
intercept = TRUE, neighbours = 1, closest = FALSE, LocalPred = LinearPred)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{data vector of the grid used in the transform.}
  \item{coeff}{vector of detail and scaling coefficients in the wavelet decomposition of the signal.}
  \item{lengths}{vector of interval lengths to be used in the update step of the transform. This is of length \env{pointsin}.}
  \item{lengthsremove}{vector of interval lengths corresponding to the points removed during the forward transform.}
  \item{pointsin}{indices into \env{X} of the scaling coefficients in the wavelet decomposition.}
  \item{removelist}{a vector of indices into \env{X} of the lifted coefficients during the transform (in the order of removal). }
  \item{neighbrs}{a list of indices into \env{X}. Each list entry gives the indices of the neighbours of the removed point used at that particular step of the forward transform.}
  \item{schemehist}{a vector of character strings indicating the type of regression used at each step of the forward transform. This is NULL apart from when \env{AdaptNeigh} is to be used in the transform.}
  \item{interhist}{a boolean vector indicating whether or not an intercept was used in the regression curve at each step of the forward transform. This is NULL apart from when \env{AdaptNeigh} is to be used in the transform.}
  \item{nadd}{The number of steps to perform of the inverse transform. This corresponds to (\code{length(X)-nkeep}) in the forward transform.}
  \item{intercept}{Boolean value for whether or not an intercept is used in the prediction step of the transform.}
  \item{neighbours}{the number of neighbours in the computation of the predicted value.}
  \item{closest}{a boolean value showing whether or not the neighbours were symmetrical (FALSE) about the removed point during the transform.}
  \item{LocalPred}{The type of regression to be performed. Possible options are \code{LinearPred}, \code{QuadPred}, \code{CubicPred}, \code{AdaptPred} and \code{AdaptNeigh}.}
}
\details{
 This algorithm reconstructs an estimate of a function/signal from information about detail and scaling coefficients in its wavelet decomposition. 
Step One. Extract information about the first point to be added in the transform from the \bold{last} entries in \env{removelist}, \env{lengthsremove} and \env{neighbrs}. Use this information to discover the correct placement of this point in relation to the indices in \env{pointsin}. 

Step Two. Using the information about the prediction scheme used in the "forward" transform, use the corresponding version of \env{LocalPred} to obtain prediction weights and value for the lifted point. 

Step Three. "Undo" the update step of the transform, and then the prediction step of the transform. The vector of scaling and detail coefficients, as well as the interval lengths are modified accordingly. 

Step Four. Remove the added point from \env{removelist}. Update \env{pointsin} and \env{lengths} to contain the added point.

Step Five. Return to step 1 but in the identification of the next point to add, the second to last entries in (the original) \env{removelist}, \env{lengthsremove} and \env{neighbrs} are used to indentify the point and then place it in \env{pointsin}. 

The algorithm continues like this until as many points as desired are added. 


}
\value{
\item{X}{data vector of the grid used in the transform.}
\item{coeff}{vector of signal function values after inversion.}
\item{lengths}{vector of interval lengths at the start of the transform. This should be the same as \code{intervals(X)}.}
\item{lengthsremove}{vector of interval lengths corresponding to the points added during the transform.}
\item{pointsin}{indices into \env{X} of the scaling coefficients in the wavelet decomposition. These are the indices of the \env{X} values which remain after all points in \env{removelist} have been added.  For a straight forward-inverse transform implementation, this should be \code{order(X)}.}
\item{removelist}{a vector of indices into \env{X} of the lifted coefficients during the transform (in the reverse order of how they were added).}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}

\seealso{\code{\link{AdaptNeigh}}, \code{\link{AdaptPred}}, \code{\link{CubicPred}}, \code{\link{fwtnp}}, \code{\link{invtnpmp}}, \code{\link{LinearPred}}, \code{\link{QuadPred}}, \code{\link{UndoPointsUpdate}}     }

\examples{
#
# Generate some one-dimensional data: 500 observations.
x2<-runif(500)
f2<-make.signal2("bumps",x=x2)
#
# perform the forward transform...
out<-fwtnp(x2,f2,LocalPred=AdaptPred)
#
# and now invert using the information from out...
#
fhat<-invtnp(x2,out$coeff,out$lengths,out$lengthsremove,out$pointsin,out$removelist,
 out$neighbrs,out$schemehist,out$interhist,LocalPred=AdaptPred)
#
# Now compare the original signal with the reconstruction.
sum(abs(f2-fhat$coeff))
# 
}
\keyword{methods}% at least one, from doc/KEYWORDS



