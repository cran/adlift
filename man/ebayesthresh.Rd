\name{ebayesthresh}
\alias{beta.cauchy}
\alias{beta.laplace}
\alias{cauchy.medzero}
\alias{cauchy.threshzero}
\alias{ebayesthresh}
\alias{ebayesthresh.wavelet}
\alias{isotone}
\alias{laplace.threshzero}
\alias{negloglik.laplace}
\alias{postmean.cauchy}
\alias{postmean.cauchy2}
\alias{postmean.laplace}
\alias{postmean}
\alias{postmed.cauchy}
\alias{postmed.laplace}
\alias{postmed}
\alias{tfromw}
\alias{tfromx}
\alias{threshld}
\alias{vecbinsolv}
\alias{wandafromx}
\alias{wfromt}
\alias{wfromx}
\alias{wmonfromx}
\alias{qmonfromx}
\alias{zetafromx}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ebayesthresh}
\description{
Bayesian thresholding based on heavy-tailed prior.
}
\details{
These functions perform Bayesian thresholding.  For more information, see help by Silverman (see references below).
}
\references{
Johnstone, I.M. and Silverman, B.W. (2002) EbayesThresh:  R and S-PLUS software for Empirical Bayes thresholding (Submitted for publication).


Johnstone, I.M. and Silverman, B.W. (2004) Needles and hay in haystacks: Empirical Bayes estimates of possibly sparse sequences.  \emph{Ann. Statist.}, \bold{32}, 1594--1649.
    
Bernard Silverman's website: \url{http://www.stats.ox.ac.uk/~silverma/ebayesthresh/}.
}
%\author{   who you are   }
%\note{The new function \code{postmean.cauchy2} approximates \code{postmean} at small \env{x} values, since the original function \code{postmean.cauchy} does not cope well with this.  To make use of this function, \code{postmean} has been changed from the original function (supplied with ebayesthresh) to call \code{postmean.cauchy2} instead of \code{postmean.cauchy}.
}

\seealso{\code{\link{denoise}},\code{\link{denoisehetero}}, \code{\link{denoiseheteromp}}, \code{\link{denoiseheteroprop}}     }
}
\keyword{methods}% at least one, from doc/KEYWORDS

