\name{basisfns}
\alias{basisfns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{basisfns}
\description{
This function plots all mother and father wavelets associated with a given wavelet transform.   
}
\usage{
basisfns(x, f, pred, neigh, int, clo, keep, plot.f = FALSE, 
plot.bas = FALSE, separate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a gridpoint vector.}
  \item{f}{the vector of associated function values.}
  \item{pred}{The type of regression to be performed. Possible options are \link{LinearPred}, \link{QuadPred}, \link{CubicPred}, \link{AdaptPred} and \link{AdaptNeigh}.}
  \item{neigh}{The number of neighbours over which the regression is performed at each step. If \env{closest} is false, then this in fact denotes the number of neighbours on each side of the removed point.}
  \item{int}{Indicates whether or not the regression curve includes an intercept.}
  \item{clo}{Refers to the configuration of the chosen neighbours. If \env{closest} is false, the neighbours will be chosen symmetrically around the removed point. Otherwise, the closest neighbours will be chosen.}
  \item{keep}{The number of scaling coefficients to be kept in the final representation of the initial signal. This must be at least two.}
  \item{plot.f}{a boolean value indicating whether to plot the original function or not. If so, the signal is plotted with vertical coloured lines, showing which prediction method was used on the different parts of the signal. The plot also shows which gridpoints correspond to scaling functions.}
  \item{plot.bas}{ subset of \code{1:length(f)}, denoting which basis functions to plot. Each basis function is colour-coded according to which prediction scheme was used in the lifting of the corresponding gridpoint.}
  \item{separate}{a boolean argument indicating if the basis functions should be plotted on a single graphsheet.}
}
\details{
 The procedure constructs \env{W}, the matrix representation of the forward transform specified in the arguments to the function, and then uses the inverse matrix to calculate the vectors of basis function values: to work out the basis function values, one inverts the transform with a delta vector, with a one in the position corresponding to the basis function required. Since this is equivalent to pre-multiplying the delta vector by the matrix representation for the inverse transform \eqn{(W^{-1})}, the basis function values are precisely the columns of \eqn{W^{-1}}. The procedure then plots the basis functions (each on a separate graphsheet, if chosen), colour coded according to the prediction scheme used or whether it is a scaling function. }
\value{
\item{out}{the output from the forward transform which is specified in the arguments to this function}
%\item{maxv}{a vector of maximum values attained by the wavelet basis functions used in the specified wavelet transform. }
%\item{removelist}{the vector of point indices removed during the transform.}
$\item{pointsin}{the vector of indices of points still to be removed.}
\item{schhist}{a character string vector of the prediction scheme used for the prediction of each gridpoint (in the order of \env{x}).}
\item{inthist}{vector of boolean values indicating whether an intercept was used in the prediction steps during the transform (in the order of \env{x}). }
\item{basmat}{a matrix of wavelet basis function values. The row \env{i} represents the function values corresponding to the grid for the basis function associated to the gridpoint \eqn{i}.}
%\item{schemehist}{a character string vector of the prediction scheme used for the prediction of each gridpoint (in the order of removal). }
%\item{interhist}{vector of boolean values indicating whether an intercept was used in the prediction steps during the transform (in the order of removal). }
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}
\note{If \env{plot.bas=T}, since the function produces one graph for each gridpoint, R or Splus is likely to exceed the total number of open devices for large datasets.}

\seealso{\code{\link{fwtnp}},\code{\link{fwtnp}}     }
\examples{
#create test signal data
#
x<-runif(100)
y<-make.signal2("blocks",x=x)
#
#perform procedure...
#
a<-basisfns(x,y,AdaptNeigh,2,TRUE,TRUE,2,FALSE,c(1,14,15),FALSE)
#
#this produces plots of three basis functions all on one graph.
}
\keyword{graphs}% at least one, from doc/KEYWORDS

