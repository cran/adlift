\name{condno}
\alias{condno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{condno}
\description{
This function uses a specified norm to compute the condition number of a matrix representation of a wavelet transform.    
}
\usage{
condno(W, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{a matrix which represents a wavelet transform.}
  \item{type}{a character string denoting which norm to use when computing the condition number. Possible values are \env{"l1"}, or one of the standard norm types, \env{"F"} (Frobenius norm), \env{"i"} (infinity norm), \env{"m"} (max modulus of a matrix) or \env{"1"} (1-norm). 
}
}
\details{
     The function computes the condition number as \env{condno} = \eqn{||W||*||W^{-1}||}. 
}
\value{
  \item{condno}{the condition number of the matrix \env{W}. }
}

\author{Matt Nunes (\email{matt.nunes@bristol.ac.uk}), Marina Popa (\email{Marina.Popa@bristol.ac.uk})}
\note{The matrix \env{W} must be invertible.}

\examples{
#create test signal data
#
x<-runif(100)
y<-make.signal2("blocks",x=x)
#
a<-transmatdual(x,y,Pred=AdaptNeigh,neigh=2)
#
#computes the transition matrix for the specified options
#
W<-a$Wnew
#
condno(W,"F")
#
condno(W,"l1")
#
condno(W,"1")
# 
}
\keyword{array}% at least one, from doc/KEYWORDS
\keyword{algebra}%   ONLY ONE   keyword per line