\name{artlev}
\alias{artlev}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{artlev}
\description{
     This function splits the coefficients into levels according to increasing quantiles of the removed interval lengths.    
}
\usage{
artlev(y, rem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a vector of the removed interval lengths (in the order of \env{removelist}). }
  \item{rem}{vector of indices of the removed points (from the output of the forward transform). }
}
\details{
     The function finds the median of the removed interval lengths, and takes all \env{pointsin} indices with removed interval lengths at most this value as the first artificial level. These indices are now not considered in later groups. The cut-off value, q, is now increased to the 75th percentile, and the indices at most this value are grouped into the second level. The procedure is continued with successive percentiles (1+q)/2 until all indices are grouped. At each stage, the level size is checked to ensure it has at least 10 elements, and if not, the level is taken together with the next level (i.e. the present percentile is ignored, and increased to the q value).   
}
\value{
    \item{p}{a list of the grouped indices of \env{removelist} (in decreasing group size) indicating thresholding groups.}
  }

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}

\seealso{\code{\link{denoise}},     }
\examples{
#create test signal data
#
x<-runif(100)
y<-make.signal2("blocks",x=x)
#
#perform forward transform...
#
out<-fwtnp(x,y,LocalPred=AdaptNeigh,neighbours=2)
#
al<-artlev(out$lengthsremove,out$removelist)
#
#
# the indices of removelist split into levels:
al
#
}
\keyword{arith}% at least one, from doc/KEYWORDS
%\keyword{lists}%   ONLY ONE   keyword per line









