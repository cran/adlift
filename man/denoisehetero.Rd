\name{denoisehetero}
\alias{denoisehetero}
\alias{denoisehetero2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{denoisehetero}
\description{
Denoises the inputted signal using artificial levels noise variance estimation and bayesian thresholding, using heteroscedastic (estimated) noise variances.  
}
\usage{
denoisehetero(x, f, pred, neigh, int, clo, keep, rule = "median",
returnall=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of grid values. Can be of any length, not necessarily equally spaced. }
  \item{f}{A vector of function values corresponding to \env{x}. Must be of the same length as \env{x}.}
  \item{pred}{The type of regression to be performed. Possible options are \command{\link{LinearPred}}, \command{\link{QuadPred}}, \command{\link{CubicPred}}, \command{\link{AdaptPred}} and \command{\link{AdaptNeigh}}.}
  \item{neigh}{The number of neighbours over which the regression is performed at each step. If clo is false, then this in fact denotes the number of neighbours on each side of the removed point.}
  \item{int}{Indicates whether or not the regression curve includes an intercept.}
  \item{clo}{Refers to the configuration of the chosen neighbours. If \env{clo} is false, the neighbours will be chosen symmetrically around the removed point. Otherwise, the closest neighbours will be chosen.}
  \item{keep}{The number of scaling coefficients to be kept in the final representation of the initial signal. This must be at least two.}
  \item{rule}{The type of bayesian thresholding used in the procedure. Possible values are \env{"mean"}, \env{"median"} (posterior mean or median thresholding) or "hard or \cr "soft" (hard or soft thresholding). }
  \item{returnall}{Indicates whether the function returns useful variables or just the denoised datapoints.}
}
\details{
The function uses the transform matrix to normalise the detail coefficients produced from the forward transform, so that they can be used in the bayesian thresholding procedure EbayesThresh. The coefficients are divided into artificial levels, and the first (largest)level is used to estimate the noise variances of the coefficients, based on the MAD of those coefficients falling in a sliding window around each gridpoint.  EbayesThresh is then used to threshold the coefficients. The resulting new coefficients are then unnormalised and the transform inverted to obtain an estimate of the true (unnoisy) signal.    
}
\value{
  If \env{returnall=FALSE}, the estimate of the function after denoising.  If \env{returnall=TRUE}, a list with components:\cr

  \item{fhat}{the estimate of the function after denoising.}
  \item{fhat1}{the estimate of the function after denoising, using the alternate variance estimate of MAD, centered at zero.}
  \item{fhat2}{the estimate of the function after denoising, using the alternate variance estimate of the median of the absolute values of the detail coefficients.}
  \item{w}{the matrix associated to the wavelet transform.}
  \item{indsd}{the individual coefficient variances introduced by the transform.}
  \item{al}{the artificial levels used to estimate the noise variance.}
  \item{sd}{the standard deviation of the noise.}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}

\seealso{\code{\link{denoise}},  \code{\link{heterovar}}     }
\examples{
x1<-runif(256)
y1<-make.signal2("doppler",x=x1)
n1<-rnorm(256,0,.1)
z1<-y1+n1
#
est1<-denoisehetero(x1,z1,AdaptNeigh,1,TRUE,TRUE,2)
traceback()
sum(abs(y1-est1))
#
#the error between the true signal and the denoised version. 


}
\keyword{smooth}% at least one, from doc/KEYWORDS
\keyword{regression}%   ONLY ONE   keyword per line


