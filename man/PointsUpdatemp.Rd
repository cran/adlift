\name{PointsUpdatemp}
\alias{PointsUpdatemp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PointsUpdatemp}
\description{
This function performs the update lifting step using a given configuration of neighbours and boundary handling.   
}
\usage{
PointsUpdatemp(X, coeff, nbrs, newnbrs, index, remove, pointsin, 
weights, lengths)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the vector of grid values.}
  \item{coeff}{the vector of detail and scaling coefficients at that step of the transform.}
  \item{nbrs}{the indices (into \env{X}) of the neighbours to be used in the lifting step.}
  \item{newnbrs}{as nbrs, but repeated according to the multiple point structure of the grid.}
  \item{index}{the indices into \env{pointsin} of \env{nbrs}, the neighbours of \env{remove}.}
  \item{remove}{the index (into \env{X}) of the point to be removed.}
  \item{pointsin}{The indices of gridpoints still to be removed.}
  \item{weights}{the prediction weights obtained from the regression in the prediction step of the transform.}
  \item{lengths}{the vector of interval lengths at the present step of the transform (to be updated).}
}
\details{
The procedure performs a minimum norm update lifting step. Firstly the interval lengths are updated using the coefficients obtained. Secondly, the scaling and detail coefficient list is modified using the new interval lengths.
}
\value{
  \item{coeff}{vector of (modified) detail and scaling coefficients to be used in the next step of the transform.}
  \item{lengths}{the vector of interval lengths after the update step of the transform.}
  \item{r}{the index into \env{pointsin} of \env{remove}.}
  \item{N}{length(\env{pointsin}).}
  \item{weights}{The regression coefficients used in prediction.}
  \item{alpha}{the update weights used to update \env{lengths} and \env{coeff}.}
}

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\seealso{\code{\link{AdaptNeighmp}}, \code{\link{AdaptPredmp}}, \code{\link{CubicPredmp}}, \code{\link{fwtnpmp}}, \code{\link{LinearPredmp}},
 \code{\link{QuadPredmp}}, \code{\link{UndoPointsUpdatemp}}      }

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
nbrs
newnbrs<-NULL
for (i in 1:length(nbrs)){
newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
}
newnbrs

coefflist[[31]]<-p$results[[6]][31]

u<-PointsUpdatemp(X,coefflist,p$newinfo[[2]],newnbrs,p$newinfo[[3]],31,
order(X),p$results[[4]],lengths)
#
#and here is the updated coefficient list...
u$coeff
}
\keyword{arith}% at least one, from doc/KEYWORDS








