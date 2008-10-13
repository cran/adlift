\name{pts}
\alias{pts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pts}
\description{
This function constructs the gridpoints (\env{X} values) from a vector of lengths, with optional starting point.   
}
\usage{
pts(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{a vector containing the \env{X} starting point in the first entry, and the distances between observed values. \code{input[1]} can be NA.}
}
\details{
     The function constructs the grid values using the first entry in input as the first gridpoint, consecutively adding the given lengths in the rest of input to find the other gridpoints. If \code{input[1]} is NA, then the starting point of the gridpoint vector \env{X} is taken to be \code{sum(lengths)/(length(lengths)+1)}.   
}
\value{
  \item{lengths}{a vector of length \code{(length(input)-1)} with the distances between the constructed gridpoints(this is \code{input[2:length(input)]}).}
  \item{X}{the constructed gridpoint vector.}
}

\author{Matt Nunes (\email{matt.nunes@bristol.ac.uk}), Marina Popa (\email{Marina.Popa@bristol.ac.uk})}
\seealso{\code{\link{fwtnp}}, \code{\link{fwtnpmp}}     }
\examples{
y<-runif(20)
y
#
y1<-c(0,y)
y2<-c(NA,y)
#
#
pts(y1)$X
#the gridpoints, with specified startpoint of zero.
#
#
pts(y2)$X
#the grid vector with unspecified startpoint.
# 

}
\keyword{arith}% at least one, from doc/KEYWORDS

