\name{mcdata}
\alias{mcdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mcdata}
\description{
Loads motorcycle dataset.   
}
\usage{
mcdata()
}
%- maybe also 'usage' for other objects documented here.
\details{The function reads the motorcycle data (Hardle) in from the text file included with this package.
}

\value{
  \item{time}{time after a simulated impact with motorcycles.}
  \item{accel}{head acceleration of a PTMO (post mortem human test object). }
}

\author{Matt Nunes (\email{matt.nunes@bristol.ac.uk}), Marina Popa (\email{Marina.Popa@bristol.ac.uk})}
\references{Hardle, W.(1990). \emph{Applied Nonparametric Regression}. Cambridge University Press.}

\seealso{\code{\link{motorcycledata}}}

\examples{
mcdata()
#
time
accel
#
}
\keyword{datasets}% at least one, from doc/KEYWORDS

