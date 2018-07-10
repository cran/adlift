\name{PointsUpdate}
\alias{PointsUpdate}
\alias{PointsUpdate2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PointsUpdate}
\description{
This function performs the update lifting step using a given configuration of neighbours and boundary handling.   
}
\usage{
PointsUpdate(X, coeff, nbrs, index, remove, pointsin, weights, lengths)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the vector of grid values.}
  \item{coeff}{the vector of detail and scaling coefficients at that step of the transform.}
  \item{nbrs}{the indices (into \env{X}) of the neighbours to be used in the lifting step.}
  \item{index}{the indices into \env{pointsin} of \env{nbrs}, the neighbours of \env{remove}.}
  \item{remove}{the index (into \env{X}) of the point to be removed.}
  \item{pointsin}{The indices of gridpoints still to be removed.}
  \item{weights}{the prediction weights obtained from the regression in the prediction step of the transform.}
  \item{lengths}{the vector of interval lengths at the present step of the transform (to be updated).}
}
\details{
The procedure performs a minimum norm update lifting step. Firstly the interval lengths are updated using the coefficients obtained. Secondly, the scaling and detail coefficient vector is modified using the new interval lengths.
}
\value{
  \item{coeff}{vector of (modified) detail and scaling coefficients to be used in the next step of the transform.}
  \item{lengths}{the vector of interval lengths after the update step of the transform.}
  \item{r}{the index into \env{pointsin} of \env{remove}.}
  \item{N}{length(\env{pointsin}).}
  \item{weights}{The regression coefficients used in prediction.}
  \item{alpha}{the update weights used to update \env{lengths} and \env{coeff}.}

}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}

\seealso{\code{\link{AdaptNeigh}}, \code{\link{AdaptPred}}, \code{\link{CubicPred}}, \code{\link{fwtnp}}, \code{\link{LinearPred}}, \code{\link{QuadPred}}, \code{\link{UndoPointsUpdate}}      }
\examples{
#
# Generate some blocks data: 100 observations.
#
x <- runif(100)
y <-make.signal2("blocks",x=x)
#
#find initial interval lengths...
#
I<-intervals(x,"reflect")
lengths<-lengthintervals(x,I,neighbours=2,closest=FALSE)
#
#perform prediction step...
p<-AdaptNeigh(order(x),x,y,32,5,TRUE,2)
#
#
u<-PointsUpdate(x,p$results[[6]],p$newinfo[[2]],p$newinfo[[3]],5,order(x),p$results[[4]],lengths)
#
#and here are the updated coefficients...
u$coeff
# 


}
\keyword{arith}% at least one, from doc/KEYWORDS

