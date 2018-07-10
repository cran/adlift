\name{findadds}
\alias{findadds}
\title{Finds minimum number of inversion steps}
\description{This function finds the minimum inversion steps to perform 
to fully reconstruct (a subset of) data
}
\usage{
findadds(rem, neighbrs, po, index = 1:(length(rem) + length(po)))
}
\arguments{
  \item{rem}{the \code{removelist} variable from a lifting 
  decomposition.  See \code{\link{fwtnp}}.}
  \item{neighbrs}{A list of neighbour indices corresponding to lifting 
  steps in a decomposition.  See \code{\link{fwtnp}}.}
  \item{po}{The \code{pointsin} variable in a lifting decomposition 
  (the index into the unlifted datapoints). See \code{\link{fwtnp}}.}
  \item{index}{a vector of indices into the original data, indicating 
  which points should be fully reconstructed during inversion.}
}
\details{
This function gives a computational shortcut to get datapoint 
information in certain inversion cases.  In some circumstances,when 
inverting, you might only be interested in the inverted coefficients for 
a subset of timepoints.  In this case, it is not necessary to do a full 
inversion to look at the desired coefficients;  the function uses the 
neighbourhood and removal order of the forward transform information and 
notes: 1) when the desired points were lifted (if at all) and 2) when 
the desired points were used as neighbours in prediction (if applicable).  
The number of inversion steps needed for each index individually is then 
taken as the maximum for these two conditions to be met.  Inverting the 
transform with this number will yield the correct inverted coefficient.  
Note that to get the correct coefficients for all \code{index}, the number of 
inversion steps is \code{max(adds)}.
}
\value{\code{adds}: a vector corresponding to \code{index}, each element 
of which is 
the number of inversion steps needed for that datapoint to be fully 
reconstructed.}
\author{
Matt Nunes (\email{nunesrpackages@gmail.com})
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{invtnp}}
}
\examples{
#create data:
x<-runif(256)

f<-make.signal2("bumps",x=x)

#do forward transform:
fwd<-fwtnp(x,f)

#I want to invert enough so that points 1:3 are reconstructed.

adds<-findadds(fwd$removelist,fwd$neighbrs,fwd$pointsin,1:3)
adds

#now reconstruct...
fhat<-invtnp(fwd,f,nadd=max(adds))

#...and check that the desired points are the same:
f[1:3]
fhat[1:3]

}
\keyword{methods}
