\name{getnbrs}
\alias{getnbrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getnbrs}
\description{
This function uses the user's neighbourhood configuration input to find the neighbours of the lifted datapoint to be used in the prediction step of the transform.
}
\usage{
getnbrs(X, remove, pointsin, neighbours, closest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The vector of gridpoints.}
  \item{remove}{the index (into \env{X}) of the point to be removed.}
  \item{pointsin}{The indices of gridpoints still to be removed.}
  \item{neighbours}{the number of neighbours to find for use in prediction.}
  \item{closest}{Boolean argument: If FALSE, the neighbours selected are the ones on both sides of the removed point.}
}
\details{
  The function uses the value of \env{neighbours} and \env{closest} to choose the neighbours to return. If \env{closest} is FALSE, \env{pointsin} is used to find \env{neighbours} indices on both sides of the index of the removed point (\env{remove}). If \env{closest} is TRUE, then the function uses the gridpoint vector (\env{X}) to calculate distances from the removed point to \env{neighbours} neighbours on each side of the removed point (if they exist) and then uses this information to choose the closest \env{neighbours} ones, recording where they lie in relation to the removed point, and accordingly their index can be obtained. If the removed point is on the boundary, then by choice, we take only one neighbour.
}
\value{
\item{nbrs}{the indices of the neighbours corresponding to the specified configuration.}
%\item{r}{the index into \env{pointsin} of the removed point.}
%\item{leftneigh}{the number of neighbours on the left of the removed point.}
%\item{rightneigh}{the number of neighbours on the right of the removed point.}
\item{index}{the indices into \env{pointsin} of the neighbours}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\seealso{\code{\link{fwtnp}}     }

\examples{
x1<-runif(20)
#
x1
#
y1<-make.signal2("bumps",x=x1)
#
y1
#
order(x1)
#
# shows where the points lie in relation to each other.
#
neigh<-getnbrs(x1,3,order(x1),4,TRUE)
#
neigh$nbrs
#
# these are the indices of the 4 closest neighbours to point 3.
#
}
\keyword{arith}% at least one, from doc/KEYWORDS

