\name{modjitter}
\alias{modjitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{modjitter}
\description{
   This function jitters grid values by a proportion of the regular distance between consecutive gridpoints and then alters it to lie in [0,1].    
}
\usage{
modjitter(x, amount)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector to be jittered (e.g. a gridpoint vector). }
  \item{amount}{a value of how much to jitter the vector (expressed as a proportion of the regular gridpoint distance, d). }
}
\details{
    The function uses dojitter to jitter the gridpoint vector by (\code{amount*d}) . The endpoints are fixed to be zero and one, and the corresponding \env{jx} values to \code{x[2]} and \code{x[length(x)-1]} are randomised again in the intervals \code{[0,x[2]+amount*d]} and \code{[x[length(x)-1]-amount*d,1]} respectively.  
}
\value{
  \item{jx}{the jittered version of \env{x}}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\seealso{\code{\link{dojitter}}, \code{\link{make.signal2}}    }
\examples{
#create grid vector
#
xgrid<-seq(0,1,length=51)
#
xgrid
#
#a regularly-spaced grid on [0,1]
#
modjitter(xgrid,1)
#
#jitters xgrid with a maximum change of .02, keeping endpoints of zero and one 
}
\keyword{arith}% at least one, from doc/KEYWORDS
\keyword{manip}%   ONLY ONE   keyword per line
