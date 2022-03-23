\name{AdaptNeighmp}
\alias{AdaptNeighmp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AdaptNeighmp}
\description{
This function performs the prediction lifting step over neighbourhoods and interpolation schemes, for multiple point data.
}
\usage{
AdaptNeighmp(pointsin, X, coefflist, coeff, nbrs, newnbrs, 
remove, intercept, neighbours, mpdet, g)
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
  \item{neighbours}{the number of neighbours in the computation of the predicted value. This is not actually used specifically in \command{AdaptNeighmp}, since this is known already from \env{nbrs}.}
  \item{mpdet}{how the mutiple point detail coefficients are computed.  Possible values are "ave", in which the multiple detail coefficients produced when performing the multiple predictions are averaged, or "min", where the overall minimum detail coefficient is taken.  }
  \item{g}{the group structure of the multiple point data.  }  
}
\details{
The procedure performs adaptive regression (through \command{AdaptPred}) over the three types of regression and also over the 3*\env{neighbours} configurations of neighbours. The combination (type of regression, configuration of neighbours) is chosen which gives the smallest detail coefficient (in absolute value). 
}
\value{
  \strong{results.}  This is a ten item list giving the regression information chosen from the detail coefficient minimisation (i.e, the information supplied to \command{AdaptNeigh} by \command{AdaptPred}): 
  \item{Xneigh}{matrix of \env{X} values corresponding to the neighbours of the removed point. The matrix consists of columns \eqn{1,X[nbrs],X[nbrs]^2,...} depending on the order of the prediction used and whether or not an intercept is used. Refer to any reference on linear regression for more details. }
  \item{mm}{the matrix from which the prediction is made. In terms of \env{Xneigh}, it is \cr \eqn{({Xneigh}^T Xneigh)^{-1} {Xneigh}^T} . }
  \item{bhat}{The regression coefficients used in prediction. }
  \item{weights}{the prediction weights for the neighbours. }
  \item{pred}{the predicted function value obtained from the regression. }
  \item{coeff}{vector of (modified) detail and scaling coefficients to be used in the update step of the transform. }
  \item{int}{if TRUE, an intercept was used in the regression. }
  \item{scheme}{a character vector denoting the type of regression used in the prediction ("Linear", "Quad" or "Cubic"). }
  \item{details}{a vector of the detail coefficients from which \command{AdaptPred} selects the minimum value. There are six entries. The first three entries represent the detail coefficients from regression with no intercept in increasing order of prediction. The second three details are values for regression with intercept. }
\item{minindex}{the index into details (\env{results[[9]]}) which produces the minimum value.}

\strong{newinfo.} A six item list containing extra information to be used in the main transform procedure (\command{fwtnp}) obtained from the minimisation in \command{AdaptNeigh}: 
\item{clo}{boolean value telling the configuration of the neighbours which produce the overall minimum detail coefficient. }
\item{totalminindex}{the index into \env{mindetails} (below) indicating the overall minimum detail coefficient produced by the procedure. }
\item{nbrs}{the indices into \env{X} of the neighbours used in the best prediction scheme.}
\item{index}{the indices into \env{pointsin} of the neighbours used in the best prediction. }
\item{mindetails}{a vector of 3*\env{neighbours} entries giving the minimum details produced by each call of \command{AdaptPred} in \command{AdaptNeigh} (for the different number and configuration of neighbours). }
\item{minindices}{vector of 3*\env{neighbours} entries giving the index (out of 6) of the schemes which produce the best predictions by each call of \cr \command{AdaptPred} in \command{AdaptNeigh}.}
}

%\note{The \pkg{Matrix} package is needed for this function.} 

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina Knight }

\seealso{\code{\link{AdaptPredmp}}, \code{\link{fwtnpmp}}     }

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

AdaptNeighmp(order(X),X,coefflist,coeff,nbrs,newnbrs,31,TRUE,2,"ave",g)

}
\keyword{regression}% at least one, from doc/KEYWORDS

