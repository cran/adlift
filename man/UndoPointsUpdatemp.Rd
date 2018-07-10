\name{UndoPointsUpdatemp}
\alias{UndoPointsUpdatemp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{UndoPointsUpdatemp}
\description{
This function undoes the update lifting step in the multiple observation inverse transform.
}
\usage{
UndoPointsUpdatemp(X, coeff, nbrs, newnbrs, index, remove, r, N, pointsin,
 gamweights, lengths, lengthrem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the vector of grid values.}
  \item{coeff}{the vector of detail and scaling coefficients at that step of the transform.}
  \item{nbrs}{the indices (into \env{X}) of the neighbours to be used in the lifting step.}
  \item{newnbrs}{as nbrs, but repeated according to the multiple point structure of the grid.}
  \item{index}{the indices into \env{pointsin} of \env{nbrs}, the neighbours of \env{remove}, the point to be added.}
  \item{remove}{the index (into \env{X}) of the point to be added.}
  \item{r}{the index into \env{pointsin} of the added point, \env{remove}.}
  \item{N}{length(\env{pointsin}).}
  \item{pointsin}{The indices of gridpoints still to be added.}
  \item{gamweights}{the prediction weights obtained from the regression in the prediction step of the transform.}
  \item{lengths}{the vector of interval lengths at the present step of the transform.}
  \item{lengthrem}{the interval length associated to the point to be added.}
}
\details{
 This procedure uses minimum norm update coefficients to invert the update step of the transform. The prediction weights are used to change the interval lengthsm before the update weights are used to modify \env{coefflist}.   
}
\value{
  \item{coeff}{vector of (modified) detail and scaling coefficients to be used later in the transform.}
  \item{lengths}{vector of interval lengths after inverting the update step of the transform.}
  \item{alpha}{the weights used to modify \env{lengths} and \env{coeff}.}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }
\seealso{\code{\link{AdaptNeighmp}}, \code{\link{AdaptPredmp}}, \code{\link{CubicPredmp}}, \code{\link{invtnpmp}}, \code{\link{LinearPredmp}}, \code{\link{PointsUpdatemp}}, \code{\link{QuadPredmp}}     }
\examples{
#read in data with multiple values...

data(motorcycledata)
times<-motorcycledata$time
accel<-motorcycledata$accel
short<-adjustx(times,accel,"mean")
X<-short$sepx
coeff<-short$sepx
g<-short$g

coefflist<-list()
for (i in 1:length(g)){
coefflist[[i]]<-accel[g[[i]]]
}

I<-intervals(X,"reflect")
lengths<-lengthintervals(X,I,neighbours=2,closest=TRUE)

#work out neighbours of point to be removed (31)

out<-getnbrs(X,31,order(X),2,TRUE)
nbrs<-out$n

nbrs

newnbrs<-NULL
for (i in 1:length(nbrs)){
newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
}

#work out repeated neighbours using g...
newnbrs

p<-AdaptNeighmp(order(X),X,coefflist,coeff,nbrs,newnbrs,31,TRUE,2,"ave",g)

nbrs<-p$newinfo[[3]]
newnbrs<-NULL
for (i in 1:length(nbrs)){
newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
}
coefflist[[31]]<-p$results[[6]][31]

u<-PointsUpdatemp(X,coefflist,p$newinfo[[2]],newnbrs,p$newinfo[[3]],31,order(X),p$results[[4]],
lengths)

p2<-setdiff(order(X),31)
a<-which(order(X)==31)
l2<-lengths[setdiff(1:length(X), a)]
#
#remove the lifted coefficient
#
#now undo the update step...
#
undo<-UndoPointsUpdatemp(X,coeff,newnbrs,p$newinfo[[2]],p$newinfo[[3]],31,
a,length(X)-1,p2,p$results[[4]],l2,lengths[a])
#

}
\keyword{arith}% at least one, from doc/KEYWORDS








