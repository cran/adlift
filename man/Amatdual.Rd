\name{Amatdual}
\alias{Amatdual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Amatdual}
\description{Combines filter matrices to produce a refinement matrix A for a wavelet transform. 
}
\usage{
Amatdual(steps, pointsin, removelist, nbrs, weights, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{steps}{a value indicating which refinement matrix to construct. It refers to the number of points already removed during the transform.}
  \item{pointsin}{The indices of gridpoints still to be removed. }
  \item{removelist}{a vector of indices into env{X} of the lifted coefficients during the transform (in the order of removal). }
  \item{nbrs}{indices of the neighbours used in the last step of the decomposition.}
  \item{weights}{the prediction weights obtained from the regression in the prediction step of the transform.}
  \item{alpha}{the update weights used to update \env{lengths} and \env{coeff}.}
}
\details{The function uses the prediction and update weights to construct the filter matrices \env{Hdual} and \env{Gdual}. Combining these two matrices results in the refinement matrix \env{Adual}.
}
\value{
\item{Adual}{the refinement matrix for the particular step of the transform. }
\item{Hdual}{the high-pass filter matrix for the current step of the transform.}
\item{Gdual}{the low-pass filter matrix for the current step of the transform.}
\item{o}{the indices of \env{nbrs} into the vector of \env{pointsin} and the \env{steps} removed points of the transform.}
\item{alpha}{the update weights used to update \env{lengths} and \env{coeff}.}
\item{weights}{the prediction weights obtained from the regression in the prediction step of the transform.}
 }

\author{Matt Nunes (\email{matt.nunes@bristol.ac.uk}), Marina Popa (\email{Marina.Popa@bristol.ac.uk})}
\seealso{\code{\link{transmatdual}}     }
\examples{
#
x<-runif(256)
y<-make.signal2("doppler",x=x)
a<-fwtnp(x,y,LocalPred=AdaptNeigh,neighbours=2)
#
A<-Amatdual(90,a$pointsin,a$removelist,a$neighbrs[[90]],a$gamlist[[90]],a$alphalist[[90]])
#
A$Adual
#
#the 90th refinement matrix for the transform above.
#
}
\keyword{array}% at least one, from doc/KEYWORDS
