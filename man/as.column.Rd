\name{as.column}
\alias{as.column}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{as.column}
\description{This function returns a given vector as a column (with dimension).
}
\usage{
as.column(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{any vector or array.}
}
\details{
\env{x} can either be a vector with no dimension attributes (a list of values), a vector with dimensions, or a matrix/array. If \env{x} is a matrix/array, the function gives \env{x} if \command{ncol(x)} is less than or equal to \command{nrow(x)}, or its transpose if \command{ncol(x)} is greater than or equal to \command{nrow(x)}. For any input, the input is given non-null dimensions. 
}
\value{
\item{y}{a vector identical to \env{x}, but given as a column.} 
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}

\seealso{\code{\link{as.row}}    }
\examples{
vector<-1:8
#
vector          
#
#...vector has no dimension attributes
# 
as.column(vector)        
#
#...gives output dimension of (8,1)
#
A<-matrix(c(6,2,2,10,6,17),3,2)
#
#
as.column(A)

#
#the function has no effect on F
#
F<-t(A)
F
#now has dimension (2,3)...
#
as.column(A)
#
#the output is made to have more rows than columns

}
\keyword{manip}% at least one, from doc/KEYWORDS

