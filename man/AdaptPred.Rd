\name{AdaptPred}
\alias{AdaptPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AdaptPred}
\description{This function performs the prediction lifting step over intercept and regression order.
}
\usage{
AdaptPred(pointsin, X, coeff, nbrs, remove, intercept, 
neighbours)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pointsin}{The indices of gridpoints still to be removed.}
  \item{X}{the vector of grid values}
  \item{coeff}{the vector of detail and scaling coefficients at that step of the transform.}
  \item{nbrs}{the indices (into \env{X}) of the neighbours to be used in the prediction step. Note that the value to this input is not important, since the procedure checks the neighbourhoods structure in the minimisation step anyway, but is for standardisation of arguments to the non-adaptive prediction schemes. }
  \item{remove}{the index (into \env{X}) of the point to be removed.}
  \item{intercept}{Boolean value for whether or not an intercept is used in the prediction step of the transform. (Note that this is actually a dummy argument, since it is not necessary for the computation of the detail coefficient in \command{AdaptPred}(the intercept is part of the adaptiveness), though is used for standardising its arguments with other prediction schemes for use in the \command{fwtnp} function).}
  \item{neighbours}{the number of neighbours to be considered in the computation of predicted values and detail coefficients.}
}
\details{
     The procedure performs adaptive regression (through AdaptPred) over the three types of regression and also over intercept.  The combination (type of regression, intercept) is chosen which gives the smallest detail coefficient (in absolute value).   
}
\value{
  \strong{results.} This is a ten item list giving the regression information chosen from the detail coefficient minimisation:
  \item{Xneigh}{matrix of \env{X} values corresponding to the neighbours of the removed point. The matrix consists of columns \eqn{1,X[nbrs],X[nbrs]^2,...} depending on the order of the prediction used and whether or not an intercept is used. Refer to any reference on linear regression for more details. }
  \item{mm}{the matrix from which the prediction is made. In terms of \env{Xneigh}, it is \cr \eqn{(Xneigh^T Xneigh)^{-1} Xneigh^T} . }
  \item{bhat}{The regression coefficients used in prediction. }
  \item{weights}{the prediction weights for the neighbours. }
  \item{pred}{the predicted function value obtained from the regression. }
  \item{coeff}{vector of (modified) detail and scaling coefficients to be used in the update step of the transform. }
  \item{int}{if TRUE, an intercept was used in the regression. }
  \item{scheme}{a character vector denoting the type of regression used in the prediction ("Linear", "Quad" or "Cubic"). }
  \item{details}{a vector of the detail coefficients from which \command{AdaptPred} selects the minimum value. There are six entries. The first three entries represent the detail coefficients from regression with no intercept in increasing order of prediction. The second three details are values for regression with intercept. }
  \item{minindex}{the index into details (\env{results[[9]]}) which produces the minimum value.}
}

%\note{The \pkg{Matrix} package is needed for this function.}
\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}


\seealso{\code{\link{AdaptNeigh}}, \code{\link{CubicPred}}, \code{\link{fwtnp}}, \code{\link{LinearPred}}, \code{\link{QuadPred}}    }

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
# Perform the adaptive lifting step 
#
ap<-AdaptPred(order(tx),tx,ty,out$nbrs,173,FALSE,2)
#
#the detail coefficient:
ap[[3]]

#and let's check the scheme used:
ap[[4]]

ap[[5]]



}
\keyword{regression}% at least one, from doc/KEYWORDS









