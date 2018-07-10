\name{Rmatsolve}
\alias{Rmatsolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rmatsolve}
\description{
This function calculates matrix inverses for symmetric matrices.
}
\usage{
Rmatsolve(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{a (symmetric) matrix.}
}

\details{
This function uses the eigenvalue decomposition of a matrix \env{m} to work out its inverse.  The function is used here since standard matrix inverse algorithms do not cope well with matrices which are near singular (this often happens in the regression stages of the forward transforms.
}
\value{
\item{inv}{the matrix inverse of \env{m}.}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}

\examples{
#
#create a 4x4 matrix
m<-matrix(runif(16),4,4)

temp<-crossprod(m)

#i.e. temp is t(m)%*%m

inv<-Rmatsolve(temp)

}
\keyword{array}% at least one, from doc/KEYWORDS

