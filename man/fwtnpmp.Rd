\name{fwtnpmp}
\alias{fwtnpmp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fwtnpmp}
\description{
Performs the lifting transform on a signal with grid \env{input} and corresponding function values \env{f}, where f has multiple points, that is, more than one function value for (some of) the grid values. 
}
\usage{
fwtnpmp(input, f,  nkeep = 2, intercept = TRUE,
 initboundhandl = "reflect", neighbours = 1,
 closest = FALSE, LocalPred = LinearPredmp, mpdet="ave")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{A vector of grid values. Can be of any length, not necessarily equally spaced.}
  \item{f}{A vector of function values corresponding to \env{input}. Must be of the same length as \env{input}.}
  \item{nkeep}{The number of scaling coefficients to be kept in the final representation of the initial signal. This must be at least two. }
  \item{intercept}{Indicates whether or not the regression curve includes an intercept.}
  \item{initboundhandl}{variable specifying how to handle the boundary at the start of the transform.  Possible values are \code{"reflect"} - the intervals corresponding to the first and last datapoints are taken to have the respective grid values as midpoints; and \code{"stop"} - the first and last intervals have the first and last grid values (respectively) as outer endpoints.}
  \item{neighbours}{The number of neighbours over which the regression is performed at each step. If closest is false, then this in fact denotes the number of neighbours on each side of the removed point.}
  \item{closest}{Refers to the configuration of the chosen neighbours. If \env{closest} is false, the neighbours will be chosen symmetrically around the removed point. Otherwise, the closest neighbours will be chosen.}
  \item{LocalPred}{The type of regression to be performed. Possible options are \code{LinearPredmp}, \code{QuadPredmp}, \code{CubicPredmp}, \code{AdaptPredmp} and \code{AdaptNeighmp}.}
  \item{mpdet}{how the mutiple point detail coefficients are computed.  Possible values are "ave", in which the multiple detail coefficients produced when performing the multiple predictions are averaged, or "min", where the overall minimum detail coefficient is taken.}
}
\details{
Given \eqn{n} points on a line, \env{input}, with multiple \env{f} values, this algorithm computes a lifting transform of the (\env{input},\env{f}) data.

Step One. Order the grid values so that corresponding intervals can be constructed, using the average function value at multiple points. 

Step Two. Compute "integrals" for each point. For each point its integral is the length of the interval associated to the gridpoint. 

Step Three. Identify the point to remove as that with the smallest integral. Generally, we remove points in order of smallest to largest integral. The integrals of neighbours of removed points change at each step. 

Step Four(a). The neighbours of the removed point are identified using the specified neighbour configuration. The values of \env{f} at the removed point are predicted using the specified regression curve over the neighbours, unless an adaptive procedure is chosen. In this case, the algorithm adjusts itself.  If the removed point has multiple point neighbours, the extra points are used in the regression.  The difference between the removed point(s) \env{f} value and the prediction is computed: these are the wavelet coefficient for the removed point.  When the removed point is itself a multiple point, this will produce multiple detail coefficients at that point.  \env{mpdet} says how the final detail coefficient for that point is recorded (either averaged or the minimum).  The detail replaces the function value in the vector coeff at the removed point's location. In this way wavelet coefficients gradually overwrite (scaling) function values in \env{coeff}. 

Step Four(b). The integrals and the scaling function values (other \env{coeff} and \env{coefflist} values) of neighbours of the removed point are updated. The values of the rest of the scaling coefficients are unaffected. 

Step Five. Return to step 3 but in the identification of a point to remove the updated integrals are used. 

The algorithm continues until as many points as desired are removed.    
}
\value{
\item{x}{data vector of the grid used in the transform.}
\item{coeff}{vector of detail and scaling coefficients in the wavelet decomposition of the signal.}
\item{coefflist}{list of detail and scaling coefficients.  Should be the same as \env{coeff}, apart from possible multiple points at the scaling function values.}
\item{origlengths}{vector of initial interval lengths corresponding to the gridpoints.}
\item{lengths}{vector of (updated) interval lengths at the end of the transform. This is of length \env{nkeep}.}
\item{lengthsremove}{vector of interval lengths corresponding to the points removed during the transform (in \env{removelist}).}
%\item{weights}{the prediction (regression) weights used in the last step of the transform.}
\item{pointsin}{indices into \env{X} of the scaling coefficients in the wavelet decomposition. These are the indices of the \env{X} values which remain after all points in \env{removelist} have been predicted and removed. This has length \env{nkeep}.}
\item{removelist}{a vector of indices into \env{X} of the lifted coefficients during the transform (in the order of removal).}
\item{neighbrs}{a list of indices into \env{X}. Each list entry gives the indices of the neighbours of the removed point used at that particular step of the transform.}
%\item{nbrs}{indices of the neighbours used in the last step of the decomposition.}
%\item{r}{the index of the removed point as it appears in env{pointsin}.}
%\item{N}{\code{length(pointsin)}}
%\item{alpha}{the update weights for the last step of the transform.}
%\item{pred}{the predicted function value of the last removed point. }
%\item{bhat}{the regression coefficients in the last prediction step.}
\item{neighbours}{the user-specified number of neighbours used in the prediction step of the transform.}
\item{gamlist}{a list of all the prediction weights used at each step of the transform.}
%\item{predlist}{a list of the predicted functions values used at each step of the transform.}
\item{alphalist}{a list of the update coefficients used in the update step of the decomposition.}
%\item{lengthlist}{a list showing how the interval lengths change over the course of the decomposition.}
\item{schemehist}{a vector of character strings indicating the type of regression used at each step of the transform.}
\item{interhist}{a boolean vector indicating whether or not an intercept was used in the regression curve at each step.}
%\item{detailslist}{a vector of the detail coefficients produced from the prediction steps of the transform.}
\item{clolist}{a boolean vector showing whether or not the neighbours were symmetrical \cr (FALSE) about the removed point during the transform. This is NULL except when \code{LocalPred=AdaptNeigh}.}
%\item{minindexlist}{a vector indicating which configuration of neighbours was used at each step. Each component is a value between 1 and 3*\env{neighbours}. If the value is between 1 and \env{neighbours}, then value symmetrical neighbours were used. Otherwise, (value - \env{neighbours}) closest neighbours were used. This is NULL except when \code{LocalPred=AdaptNeigh}.}
%\item{indiceslist}{a list indicating the prediction scheme (Linear, Quad, Cubic and intercept TRUE or FALSE chosen by \command{AdaptPred}) corresponding to each of the details in \env{mindetailslist}. If the value is between 1 and 3, then no intercept was used in the prediction. Larger values indicate that an intercept was used in the prediction. Values of 1 and 4 represent Linear prediction, 2 and 5 represent Quadratic prediction, 3 and 6 represent Cubic prediction. This is NULL except when \code{LocalPredmp=AdaptNeighmp}.}
%\item{history}{a list of pairs of indices from \env{indiceslist} and \env{minindexlist}. The first number is the type of regression used, and the second indicates the neighbour configuration which results in the minimum detail coefficient.}
\item{g}{a list desscribing the group structure (indices) of the initial function values.}
\item{mp}{a boolean vector of which of the groups are actually multiple points.}
} 

\author{Matt Nunes (\email{nunesrpackages@gmail.com}), Marina.Knight}

\seealso{\code{\link{AdaptNeighmp}}, \code{\link{AdaptPredmp}}, \code{\link{CubicPredmp}}, \code{\link{fwtnp}}, \code{\link{invtnpmp}}, \code{\link{LinearPredmp}}, \code{\link{QuadPredmp}}     }

\examples{
#read in multiple point data...

data(motorcycledata)
times<-motorcycledata$time
accel<-motorcycledata$accel

out<-fwtnpmp(times,accel,LocalPred=AdaptPredmp,neighbours=2)
out$coeff

#these are the detail coefficients of the transform.

}
\keyword{methods}% at least one, from doc/KEYWORDS

