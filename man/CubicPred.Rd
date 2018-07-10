\name{CubicPred}
\alias{CubicPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CubicPred}
\description{
This function performs the prediction lifting step using a cubic regression curve given a configuration of neighbours.   
}
\usage{
CubicPred(pointsin, X, coeff, nbrs, remove, intercept, 
neighbours)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pointsin}{The indices of gridpoints still to be removed.}
  \item{X}{the vector of grid values.}
  \item{coeff}{the vector of detail and scaling coefficients at that step of the transform.}
  \item{nbrs}{the indices (into \env{X}) of the neighbours to be used in the prediction step.}
  \item{remove}{the index (into \env{X}) of the point to be removed.}
  \item{intercept}{Boolean value for whether or not an intercept is used in the prediction step of the transform. }
  \item{neighbours}{the number of neighbours in the computation of the predicted value. This is not actually used specifically in \command{CubicPred}, since this is known already from \env{nbrs}.}
}
\details{
  The procedure performs cubic regression using the given neighbours using an intercept if chosen. The regression coefficients (\env{weights}) are used to predict the new function value at the removed point. If there are not enough neighbours to generate a cubic regression curve, the order of prediction is decreased until it is possible (i.e. to \command{QuadPred}, then \command{LinearPred}).  
}
\value{
\item{Xneigh}{matrix of \env{X} values corresponding to the neighbours of the removed point. The matrix consists of columns \eqn{X[nbrs],X[nbrs]^2,X[nbrs]^3} augmented with a column of ones if an intercept is used. Refer to any reference on linear regression for more details.}
\item{mm}{the matrix from which the prediction is made. In terms of \env{Xneigh}, it is \cr \eqn{(Xneigh^T Xneigh)^{-1} Xneigh^T} .}
\item{bhat}{The regression coefficients used in prediction.}
\item{weights}{the prediction weights for the neighbours.}
\item{pred}{the predicted function value obtained from the regression.}
\item{coeff}{vector of (modified) detail and scaling coefficients to be used in the update step of the transform.}

}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }
%\note{The \pkg{Matrix} package is needed for this function. }

\seealso{\code{\link{fwtnp}}, \code{\link{LinearPred}}, \code{\link{QuadPred}}     }

\examples{
#
# Generate some doppler data: 500 observations.
#
tx <- runif(500)
ty<-make.signal2("doppler",x=tx)
#
# Compute the neighbours of point 173 (2 neighbours on each side)
#
out<-getnbrs(tx,173,order(tx),2,FALSE)

#
# Perform cubic prediction based on the neighbours (without intercept) 
#
cp<-CubicPred(order(tx),tx,ty,out$nbrs,173,FALSE,2)
#
cp$bhat

#
#the coefficients which define the cubic regression curve
#
cp$pred

#
#the predicted value from the regression curve
#
}
\keyword{regression}% at least one, from doc/KEYWORDS

