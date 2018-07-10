\name{as.row}
\alias{as.row}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{as.row}
\description{
This function returns a given vector as a row (with dimension). }
\usage{
as.row(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{any vector or array}
}
\details{
\env{x} can either be a vector with no dimension attributes (a list of values), a vector with dimensions, or a matrix/array. If \env{x} is a matrix/array, the function gives \env{x} if \command{ncol(x)} is greater than or equal to \command{nrow(x)}, or its transpose if \command{ncol(x)} is less than or equal to \command{nrow(x)}. For any input, the input is given non-null dimensions.  
}
\value{
  \item{y}{a vector identical to \env{x}, but given as a row.}
 }

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}

\seealso{\code{\link{as.column}}     }

\examples{
X<-0:5
#
X
#
as.row(X)
#
#puts input into row (matrix)
#
Y<-matrix(0:5,6,1)
#
Y
#
as.row(Y)
#
#input forced into a row.
#
}
\keyword{manip}% at least one, from doc/KEYWORDS

