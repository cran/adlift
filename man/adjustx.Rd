\name{adjustx}
\alias{adjustx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{adjustx}
\description{
This function produces new grid values to cope with data with repeated grid values according to the method chosen to deal with it. }
\usage{
adjustx(x, f, type = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of the original (repeated) gridpoints.}
  \item{f}{the vector of function values associated to the grid vector \env{X}.}
  \item{type}{The method used to cope with the multiple points. \env{"mean"} averages all function values with the same grid value. The \env{"jitter"} option adds a small amount to all but one of each repeated grid value, and associates the function values to these new gridpoints. In this way, the each gridpoint value corresponds uniquely to the function values.}
}
\details{
The function compares \env{x} to \command{unique(x)} to find the occurences of repeated grid values, and stores the information in groups. In the \env{"jitter"} case, this is then used to modify the original gridpoints by adding an epsilon to the repeated values. In the case of \env{type="mean"}, the new gridpoints are, in fact \command{unique(x)}, and the information is used to average the groups of original function values to construct \env{sepf}. 

}
\value{
\item{sepx}{the vector of new gridpoints.}
\item{sepf}{the function values associated to \env{sepx}. }
\item{groups}{a list of indices into \env{x} showing where the original repeated grid values occured.}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }


\seealso{ \code{\link{fwtnp}},     }

\examples{
#read in the motorcycle crash data 
#
data(motorcycledata)

#
dim(motorcycledata)

#check data.
#
times<-motorcycledata$time
accel<-motorcycledata$accel

a<-adjustx(times,accel,"mean")
#
#note the repeated values in the original grid data
#
#display new data vectors
a$sepx
#
a$sepf
# 
#and now the new adjusted data has length 94.
#
}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{arith}%   ONLY ONE   keyword per line
