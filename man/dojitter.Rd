\name{dojitter}
\alias{dojitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{dojitter}
\description{This function adds a random uniform vector of the same length as the input to modify the input.
}
\usage{
dojitter(x, amount = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector to be jittered (e.g. a gridpoint vector). }
  \item{amount}{a value of how much to jitter the vector \env{x}.}
}
\details{
The function creates \code{length(x)} samples from a uniform[\env{-amount},\env{amount}], and adds these to the original vector \env{x}. If \env{amount=0}, the new vector \env{jx} is the same as the original vector.
}
\value{
  \item{jx}{the jittered version of \env{x}}

}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}


\seealso{\code{\link{modjitter}}     }
\examples{
#create grid vector
#
xgrid<-seq(0,1,length=51)
#
xgrid
#
#a regularly-spaced grid
#
dojitter(xgrid,.01)
#
#a jittered grid.
# 
}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{arith}%   ONLY ONE   keyword per line
