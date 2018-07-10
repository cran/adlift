\name{lengthintervals}
\alias{lengthintervals}
%\alias{lengthintervals2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{lengthintervals}
\description{
This function constructs the vector of interval lengths from a vector of interval endpoints.   
}
\usage{
lengthintervals(X, I, type = "midpoints", neighbours, closest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The vector of gridpoints.}
  \item{I}{a vector of interval endpoints. This is of length \code{length(X)+1}. }
  \item{type}{a character string, either \code{"midpoints"} or \code{"average"}, denoting the way of computing the interval lengths, if \env{closest=TRUE}. If \code{"average"}, then the average neighbour distance is associated as the interval lengths to the gridpoints; otherwise the lengths are associated from the interval vector, \env{I} in the obvious way : right endpoint - left endpoint.}
  \item{neighbours}{the number of neighbours to be used in the prediction step of the transform. This is only used if \env{closest=TRUE}, since it specifies how many distances to average over when \env{type="average"}. }
  \item{closest}{indicates whether the neighbourhood structure to be used in the transform is symmetrical or not. When combined with \env{type="average"}, enables the option of average closest neighbour distance as the associated interval lengths to the gridpoints.}
}
\details{
The function computes the interval lengths by finding the differences between the consecutive entries of the supplied interval vector \env{I}. In the case of the choice of average closest neighbour distance interval association, the method uses the function \command{getnbrs} to find the initial neighbours of each gridpoint to compute the average distances. 
}
\value{
  \item{lengths}{a vector of \code{length(X)} with the intervals lengths associated to the gridpoints.}
  \item{initialnbrs}{a matrix with columns \code{order(X)}, possibly together with the neighbour indices into \env{X} of each gridpoint, if \env{type="average"}.}
  \item{initialindex}{If \env{closest=TRUE} and \env{type="average"}, a matrix of dimension \code{length(X)} x \code{neighbours}, showing the indices into \code{order(X)} of the neighbours of each gridpoint. Otherwise is NULL.
}
}
\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\seealso{\code{\link{getnbrs}}, \code{\link{intervals}}     }
\examples{
input<-runif(10)
#gridpoint vector
#
I<-intervals(input,"reflect")
#create the interval endpoint vector using the input
#
lengthintervals(input,I,"average",3,TRUE)
#
#computes 'intervals' based on 3 closest neighbour distance averages
#
}
\keyword{arith}% at least one, from doc/KEYWORDS










