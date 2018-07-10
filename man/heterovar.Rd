\name{heterovar}
\alias{heterovar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{heterovar}
\description{
Estimates individual wavelet coefficient variances using a sliding window approach.
}
\usage{
heterovar(y, detail, al)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a vector of the gridpoints of \env{removelist} after executing the forward transform, in the order of the gridpoint vector.}
  \item{detail}{the vector of detail coefficients after the forward transform has been performed, in the order of the gridpoint vector.}
  \item{al}{The list of indices into \env{removelist} divided into artificial levels.}
}
\details{
The function works out the interval endpoints for each gridpoint in \env{removelist}, based on an initial window length of one fifth of the range of \env{y}, and then adjusts them so that they lie within the range of \env{y}. The indices of the \env{removelist} points inside these intervals are then compared against the indices of the first artificial level for the data. These new indices are then used to compute the individual coefficient variances, based on the detail values of the new indices. If any of the window indices list entries contains less than four values, then the initial window length is increased by 5\% and the process redone, until each window contains at least four coefficients. 
}
\value{
\item{ep1}{a two-column matrix with the (true) endpoints of the windows from which to calculate the coefficient variances (according to the specified window length). }
\item{ep2}{a two-column matrix with the endpoints of the windows from which to calculate the coefficient variances (adjusted to be of the window length and within the range of \env{y}). }
\item{idlist}{a list of indices into \env{y} showing the points each interval contains.}
\item{newidlist}{a list of indices into \env{y} showing the points each interval contains, which are also in the first artificial level.}
\item{dlist}{a list of detail coefficients which correspond to the indices in \env{newidlist}.}
\item{varvec}{a vector of median absolute deviation values (from the median) for the coefficients in \env{dlist}.}
\item{varvec1}{a vector of median absolute deviation values (from the median), centered at zero, for the coefficients in \env{dlist}.}
\item{varvec2}{a vector of median absolute deviation values (from the median), centered at zero, for the coefficients in \env{dlist}.}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\seealso{\code{\link{denoisehetero}}     }

\examples{
x1<-runif(256)
#
y1<-make.signal2("doppler",x=x1)
#
fwd<-fwtnp(x1,y1,LocalPred=AdaptNeigh,neighbours=2)
#
y<-fwd$lengthsremove
rem<-fwd$removelist
al<-artlev(y,rem)
#
yrem<-x1[sort(rem)]
detail<-fwd$coeff[sort(rem)]
#
h<-heterovar(yrem,detail,al)
#
h$varvec[1:10]
#
#the first ten coefficient variances to be used in the normalisation of the detail 
#coefficients
}
\keyword{arith}% at least one, from doc/KEYWORDS

