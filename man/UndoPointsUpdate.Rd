\name{UndoPointsUpdate}
\alias{UndoPointsUpdate}
\alias{UndoPointsUpdate2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{UndoPointsUpdate}
\description{
This function undoes the update lifting step in the inverse transform.
}
\usage{
UndoPointsUpdate(X, coeff, nbrs, index, remove, r, N, pointsin, gamweights,
 lengths, lengthrem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the vector of grid values.}
  \item{coeff}{the vector of detail and scaling coefficients at that step of the transform.}
  \item{nbrs}{the indices (into \env{X}) of the neighbours to be used in the lifting step.}
  \item{index}{the indices into \env{pointsin} of \env{nbrs}, the neighbours of \env{remove}, the point to be added.}
  \item{remove}{the index (into \env{X}) of the point to be added.}
  \item{r}{the index into \env{pointsin} of the added point, \env{remove}.}
  \item{N}{length(\env{pointsin}).}
  \item{pointsin}{The indices of gridpoints still to be added.}
  \item{gamweights}{the prediction weights obtained from the regression in the prediction step of the transform.}
  \item{lengths}{the vector of interval lengths at the present step of the transform.}
  \item{lengthrem}{the interval length associated to the point to be added.}
}
\details{
 This procedure uses minimum norm update coefficients to invert the update step of the transform. The prediction weights are used to change the interval lengthsm before the update weights are used to modify \env{coeff}.   
}
\value{
  \item{coeff}{vector of (modified) detail and scaling coefficients to be used later in the transform.}
  \item{lengths}{vector of interval lengths after inverting the update step of the transform.}
  \item{alpha}{the weights used to modify \env{lengths} and \env{coeff}.}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }
\seealso{\code{\link{AdaptNeigh}}, \code{\link{AdaptPred}}, \code{\link{CubicPred}}, \code{\link{invtnp}}, \code{\link{LinearPred}}, \code{\link{PointsUpdate}},
 \code{\link{QuadPred}}     }
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
u<-PointsUpdate(x,p$results[[6]],p$newinfo[[2]],p$newinfo[[3]],5,order(x),p$results[[4]],
lengths)
#
p2<-setdiff(order(x),5)
a<-which(order(x)==5)
l2<-lengths[setdiff(1:100, a)]
#
#remove the lifted coefficient
#
#now undo the update step...
#
undo<-UndoPointsUpdate(x,u$coeff,p$newinfo[[2]],p$newinfo[[3]],5,a,99,p2,
p$results[[4]],l2,lengths[a])
#

}
\keyword{arith}% at least one, from doc/KEYWORDS


