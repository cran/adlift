\name{make.signal2}
\alias{make.signal2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{make.signal2}
\description{
This function computes signal function values based on a grid input.
}
\usage{
make.signal2(name, x, snr = Inf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{a character string of the test signal to create.}
  \item{x}{a vector of gridpoints.}
  \item{snr}{optional argument to scale the function values according 
  to a signal-to-noise ratio.}
  \item{\dots}{any additional arguments.}
}
\details{
    This function is based on the \command{make.signal} function included in the \bold{S-Plus wavelets} module, except that the \env{x} vector can be irregular. As well as the signals included for the original version (e.g. the Donoho/Johnstone test signals), a piecewise polynomial can be sampled.
}
\value{
  \item{z}{the signal function values.}
}
\note{The test signals have domain [0,1], so the grid vector \env{x} must have values within this interval.}
\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\examples{
#create grid vector
#
xgrid<-rnorm(50)
xgrid
#
pp<-make.signal2("ppoly",x=xgrid)
#
#piecewise polynomial data vector 
#
plot(sort(xgrid),pp[order(xgrid)],type="l")
# 
}
\keyword{datagen}% at least one, from doc/KEYWORDS

