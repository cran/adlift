\name{matcond}
\alias{matcond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{matcond}
\description{
Works out two alternative condition numbers for the transform associated to the prediction scheme given in the arguments to the function. 
}
\usage{
matcond(x, f, Pred, neigh, int, clo, keep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of grid values. Can be of any length, not necessarily equally spaced.}
  \item{f}{A vector of function values corresponding to \env{x}. Must be of the same length as \env{x}.}
  \item{Pred}{The type of regression to be performed. Possible options are \code{\link{LinearPred}}, \code{\link{QuadPred}}, \code{\link{CubicPred}}, \code{\link{AdaptPred}} and \code{\link{AdaptNeigh}}.}
  \item{neigh}{The number of neighbours over which the regression is performed at each step. If \env{clo} is false, then this in fact denotes the number of neighbours on each side of the removed point.}
  \item{int}{Indicates whether or not the regression curve includes an intercept.}
  \item{clo}{Refers to the configuration of the chosen neighbours. If \env{clo} is false, the neighbours will be chosen symmetrically around the removed point. Otherwise, the closest neighbours will be chosen. }
  \item{keep}{The number of scaling coefficients to be kept in the final representation of the initial signal. This must be at least two.}
}
\details{
   The function uses the transform matrices to work out their norms and singular value decompositions.  Condition numbers are calculated by \eqn{||T_j ||*||T_j^{-1} ||} and \code{svd$d[1]/svd$d[nrow(T_j)]} respectively.
}
\value{
  \item{cno}{the condition numbers for the augmented transform matrices, calculated using the Frobenius norm (see condno).}
  \item{v}{the condition numbers for the augmented transform matrices, calculated using the ratio between the largest to the smallest singular values in the singular value decomposition of the augmented matrices.}
  \item{a}{the transform matrix information for the transform (output from fwtnp).} 
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\seealso{\code{\link{condno}}, \code{\link{fwtnp}}     }
\examples{
x1<-runif(256)
y1<-make.signal2("doppler",x=x1)
#
m<-matcond(x1,y1,AdaptNeigh,2,TRUE,TRUE,2)
#
m$cno
#
m$v
# shows the two different condition number measures for the matrix associated
# to the transform performed.
#
}
\keyword{array}% at least one, from doc/KEYWORDS
\keyword{algebra}%   ONLY ONE   keyword per line





