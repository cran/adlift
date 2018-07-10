\name{intervals}
\alias{intervals}
\alias{intervals2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{intervals}
\description{
This function constructs the intervals around the grid values to be used as scaling integrals during the transform
}
\usage{
intervals(X, initboundhandl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The vector of gridpoints.}
  \item{initboundhandl}{the interval construction at the boundary. Takes the value \code{"reflect"} for intervals symmetric about the endpoints or \code{"stop"} if the endpoint intervals are limited to the edges of the dataset, i.e. the intervals end at the first and last gridpoints respectively.}
}
\details{
The function constructs the intervals by sorting the observed gridpoints. The endpoints of the intervals are found as the midpoints between consecutive (sorted) gridpoints. In this way the intervals are not necessarily centered around the gridpoints. The first and last intervals are then modified according to \env{initboundhandl} (see above). These intervals represent the support of the initial scaling functions associated to each gridpoint.
}
\value{
\item{intervals}{a vector of length (\code{length(X)+1}) with the \env{X} values of the endpoints of the intervals (including the edges).}
\item{order}{\code{order(X)} (the sorted observation order).}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\seealso{\code{\link{lengthintervals}}     }
\examples{
x2<-runif(50)
x2
#
intervals(x2,"reflect")
#
#check that the gridpoints are between the interval vector entries...
#
}
\keyword{arith}% at least one, from doc/KEYWORDS

