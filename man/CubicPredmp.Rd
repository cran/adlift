\name{CubicPredmp}
\alias{CubicPredmp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CubicPredmp}
\description{
This function performs the prediction lifting step using a cubic regression curve given a configuration of neighbours, for multiple point data.
}
\usage{
CubicPredmp(pointsin, X, coefflist, coeff, nbrs, newnbrs, remove, 
intercept, neighbours, mpdet, g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pointsin}{The indices of gridpoints still to be removed.}
  \item{X}{the vector of grid values.}
  \item{coeff}{the vector of detail and scaling coefficients at that step of the transform.}
  \item{coefflist}{the list of detail and multiple scaling coefficients at that step of the transform.}
  \item{nbrs}{the indices (into \env{X}) of the neighbours to be used in the prediction step.}
  \item{newnbrs}{as nbrs, but repeated according to the multiple point structure of the grid.}
  \item{remove}{the index (into \env{X}) of the point to be removed.}
  \item{intercept}{Boolean value for whether or not an intercept is used in the prediction step of the transform.}
  \item{neighbours}{the number of neighbours in the computation of the predicted value. This is not actually used specifically in \command{CubicPredmp}, since this is known already from \env{nbrs}.}
  \item{mpdet}{how the mutiple point detail coefficients are computed.  Possible values are "ave", in which the multiple detail coefficients produced when performing the multiple predictions are averaged, or "min", where the overall minimum detail coefficient is taken.  Note that this is taken to standardise the input when \env{LocalPredmp} is called.}
  \item{g}{the group structure of the multiple point data.  Note that this is taken to standardise the input when \env{LocalPredmp} is called.}  
}
\details{
The procedure performs cubic regression using the given neighbours using an intercept if chosen. The regression coefficients (\env{weights}) are used to predict the new function value at the removed point.   
}
\value{
\item{Xneigh}{matrix of X values corresponding to the neighbours of the removed point. The matrix consists of the column X[newnbrs] augmented with a column of ones if an intercept is used. Refer to any reference on linear regression for more details.}
\item{mm}{the matrix from which the prediction is made. In terms of \env{Xneigh}, it is \cr \eqn{(Xneigh^T Xneigh)^{-1} Xneigh^T}.}
\item{bhat}{The regression coefficients used in prediction.}
\item{weights}{the prediction weights for the neighbours.}
\item{pred}{the predicted function value obtained from the regression.}
\item{coeff}{vector of (modified) detail and scaling coefficients to be used in the update step of the transform.}
}

%\note{The \pkg{Matrix} is needed for this function.}
\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight}

\seealso{\code{\link{fwtnpmp}}, \code{\link{LinearPredmp}}, \code{\link{QuadPredmp}}     }

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

CubicPredmp(order(X),X,coefflist,coeff,nbrs,newnbrs,31,TRUE,2,"ave",g)

}
\keyword{regression}% at least one, from doc/KEYWORDS

