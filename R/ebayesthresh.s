"ebayesthresh"<-
function(x, prior = "laplace", a = 0.5, bayesfac = FALSE, sdev = 1, verbose = FALSE, 
	threshrule = "median")
{
#  Given a vector of data x, find the marginal maximum likelihood estimator
#   of the mixing weight w, and apply an appropriate thresholding rule using
#   this weight.
#  If the prior is laplace and a=NA, then the scale factor is also found by MML.
#  The thresholding rules allowed are "median", "mean", "hard", "soft" and "none";
#   if "none" is used, then only the parameters are worked out.
#  If hard or soft thresholding is used, the argument "bayesfac" specifies
#   whether to use the bayes factor threshold or the posterior median threshold.
#  If verbose=T then the routine returns a list with several arguments, including
#   muhat which is the result of the thresholding.
#  If verbose=F then only muhat is returned.
#  It is assumed that the standard deviation of the data is sdev; if sdev=NA, then
#   it is estimated using the function mad(x).
#
#  find the standard deviation if necessary and estimate the parameters
	if(is.na(sdev)) sdev <- mad(x, center = 0)
	x <- x/sdev
	pr <- substring(prior, 1, 1)
	if((pr == "l") & is.na(a)) {
		pp <- wandafromx(x)
		w <- pp$w
		a <- pp$a
	}
	else w <- wfromx(x, prior = prior, a = a)	#
	if(pr != "m" | verbose) {
		tt <- tfromw(w, prior = prior, bayesfac = bayesfac, a = a)
		tcor <- sdev * tt
	}
	if(threshrule == "median")
		muhat <- postmed(x, w, prior = prior, a = a)
	if(threshrule == "mean")
		muhat <- postmean(x, w, prior = prior, a = a)
	if(threshrule == "hard")
		muhat <- threshld(x, tt)
	if(threshrule == "soft")
		muhat <- threshld(x, tt, hard = FALSE)
	if(threshrule == "none") muhat <- NA	#
# Now return desired output
#
	muhat <- sdev * muhat
	if(!verbose)
		return(muhat)
	retlist <- list(muhat = muhat, x = x, threshold.sdevscale = tt, 
		threshold.origscale = tcor, prior = prior, w = w, a = a, 
		bayesfac = bayesfac, sdev = sdev, threshrule = threshrule)
	if(pr == "c")
		retlist <- retlist[-7]
	if(threshrule == "none")
		retlist <- retlist[-1]
	return(retlist)
}
"ebayesthresh.wavelet" <- 
function(x.dwt, vscale = "independent", smooth.levels = Inf, prior = "laplace",
	a = 0.5, bayesfac = FALSE, threshrule = "median", package = "spluswavelets"
	)
{
#
#  Carry out the empirical Bayes smoothing approach for the detail
#  coefficients of the wavelet transform  xdwt, produced using the 
#  S+Wavelets module or WaveThresh package
#
#  If smooth.levels is set to a smaller number than n.levels then 
#  only smooth.levels levels are processed.
#
#  The parameter vscale controls the estimation of the variances, 
#  as follows:
#     if vscale is a scalar, it is used as the standard deviation 
#     	of the wavelet coefficients at every level.
#     if vscale = "independent", the standard deviation is estimated 
#    	using the coefficients at the highest level, and the same 
#	value is used at every level
#     if vscale = "level" then level-dependent standard deviations 
#	are estimated
#
#  The other arguments (prior, a, bayesfac, threshrule) have the 
#  same meaning as for ebayesthresh.  If prior="laplace" and a=0, 
#  then a is estimated separately at each level
#
#  To translate this routine to work with other wavelet software, 
#  bear in mind the following features of the S+WAVELETS wavelet 
#  transform:
#  	
#  nlevs is the number of levels of detail computed within the 
#  transform.  The detail is then held in vectors which are 
#  accessible as x.dwt[[2]], x.dwt[[3]], ... , x.dwt[[nlevs+1]], 
#  with x.dwt[[2]] being the coarsest level and x.dwt[[nlevs+1]] 
#  the finest.
#
#  The routine reconstruct inverts the transform, automatically 
#  doing the right thing if the nondecimated transform is being used
#
#  In the present version, set package="wavethresh" to use the wavethresh
#   option
#
#  Grateful acknowledgements to Stuart Barber for the wavethresh implementation
#
 pkg <- casefold(substring(package, 1, 1))
 if(pkg == "s") {
	nlevs <- attributes(x.dwt)$n.levels
	slevs <- min(nlevs, smooth.levels)
	if(is.character(vscale)) {
		vs <- substring(vscale, 1, 1)
		if(vs == "i")
			vscale <- mad(x.dwt[[nlevs + 1]])
		if(vs == "l")
			vscale <- NA
	}
	for(j in ((nlevs - slevs + 2):(nlevs + 1)))
		x.dwt[[j]] <- ebayesthresh(as.vector(x.dwt[[j]]), prior,
			a, bayesfac, vscale, FALSE, threshrule)
	x <- reconstruct(x.dwt)
 }
 else {
	print("Using WaveThresh")
	nlevs <- x.dwt$nlevels
	slevs <- min(nlevs - 1, smooth.levels)
	if(is.character(vscale)) {
		vs <- substring(vscale, 1, 1)
		if(vs == "i")
			vscale <- mad(accessD(x.dwt, level = nlevs -
				1))
		if(vs == "l")
			vscale <- NA
	}
	for(j in (nlevs - slevs):(nlevs - 1)) {
		x.dwt <- putD(x.dwt, level = j, v = ebayesthresh(
			accessD(x.dwt, level = j), prior, a, bayesfac,
			vscale, FALSE, threshrule))
		print(paste("Done level", j))
	}
	if(x.dwt$type == "station")
		x <- AvBasis(convert(x.dwt))
	else x <- wr(x.dwt)
}
return(x)
}
"tfromw"<-
function(w, prior = "laplace", bayesfac = FALSE, a = 0.5)
{
#  given the vector of weights w, find the threshold or vector of
#   thresholds corresponding to these weights, under the specified prior.
#
#  if bayesfac=T the Bayes factor thresholds are found, otherwise the posterior median
#   thresholds are found.
#
#  if the Laplace prior is used, a gives the value of the scale factor
#
	pr <- substring(prior, 1, 1)
	if(bayesfac) {
		z <- 1/w - 2
		if(pr == "l")
			tt <- vecbinsolv(z, beta.laplace, 0, 10, a = a)
		if(pr == "c")
			tt <- vecbinsolv(z, beta.cauchy, 0, 10)
	}
	else {
		zz <- rep(0, length(w))
		if(pr == "l")
			tt <- vecbinsolv(zz, laplace.threshzero, 0, 10, w = w, 
				a = a)
		if(pr == "c")
			tt <- vecbinsolv(zz, cauchy.threshzero, 0, 10, w = w)
	}
	return(tt)
}
"tfromx"<-
function(x, prior = "laplace", bayesfac = FALSE, a = 0.5)
{
#  given the data x, the prior, and any other parameters, 
#   find the threshold
#   corresponding to the marginal maximum likelihood estimator
#   of the mixing weight.
#
	if ( prior=="laplace" && is.na (a) ) { wa <- wandafromx( x)
	w <- wa$w
	a <- wa$a 
	}  	else	w <- wfromx(x, prior, a = a)
	t <- tfromw(w, prior, a = a, bayesfac = bayesfac)
	return(t)
}
"wfromt"<-
function(tt, prior = "laplace", a = 0.5)
{
# find the weight that has posterior median threshold tt, 
#
	pr <- substring(prior, 1, 1)
	if(pr == "l")
		wi <- (a * pnorm(tt - a))/dnorm(tt - a) - beta.laplace(tt, a)
	if(pr == "c") {
		dnz <- dnorm(tt)
		wi <- 1 + (pnorm(tt) - tt * dnz - 1/2)/(sqrt(pi/2) * dnz * tt^2
			)
		wi[!is.finite(wi)] <- 1
	}
	return(1/wi)
}
"wfromx"<-
function(x, prior = "laplace", a = 0.5)
{
#  given the vector of data x and the function betaf
#   which finds beta(x), 
#  find the value of w that zeroes S(w) in the
#  range 
#
#  additional arguments to betaf can be passed through ...
#
#  works by successive bisection, carrying out nits harmonic bisections
#   of the original interval between wlo and 1.  
#  
#
	pr <- substring(prior, 1, 1)
	tuniv <- sqrt(2 * log(length(x)))
	wlo <- wfromt(tuniv, prior, a)
	if(pr == "l") {
		beta <- beta.laplace(x, a)
	}
	if(pr == "c") {
		beta <- beta.cauchy(x)
	}
	whi <- 1
	beta <- pmin(beta, 1e20) 
	shi <- sum(beta/(1 + beta))
	if(shi >= 0)
		return(w = 1)
	slo <- sum(beta/(1 + wlo * beta))
	if(slo <= 0)
		return(w = wlo)
	for(j in (1:30)) {
		wmid <- sqrt(wlo * whi)
		smid <- sum(beta/(1 + wmid * beta))
		if(smid == 0)
			return(w = wmid)
		if(smid > 0) {
			wlo <- wmid
		}
		else {
			whi <- wmid
		}
	}
	return(w = sqrt(wlo * whi))
}
"wandafromx"<-
function(x)
{
#  finds the marginal max lik estimators of w and a, using a bivariate optimization
#
#  The threshold is constrained to lie between 0 and sqrt ( 2 log (n))
#
#  If running R, the routine optim is used; in S-PLUS the routine is nlminb
#
	thi <- sqrt(2 * log(length(x)))
	lo <- c(0,0.04)
	hi <- c(thi,3)
	startpar <- c(1,0.5)
	if (exists("optim")) {
 		uu <- optim(startpar, negloglik.laplace, method="L-BFGS-B",
			lower = lo, upper = hi, xx = x)
               	uu <- uu$par
		}
	else {uu <- nlminb(startpar, negloglik.laplace, lower = lo, upper = hi, xx = x)
	uu <- uu$parameters}
	a <- uu[2]
	w <- wfromt(uu[1], a = a)
	return(w, a)
}
"postmean"<-
function(x, w, prior = "laplace", a = 0.5)
{
#
# find the posterior mean for the appropriate prior for 
#   given x and w and a, assuming the error variance
#   is 1.  
#
	pr <- substring(prior, 1, 1)
	if(pr == "l")
		mutilde <- postmean.laplace(x, w, a = a)
	if(pr == "c")
		mutilde <- postmean.cauchy2(x, w)
	return(mutilde)
}
"postmean.cauchy"<-
function(x, w)
{
#
#  Find the posterior mean for the quasi-Cauchy prior with mixing weight w
#   given data x, which may be a scalar or a vector.
#
	muhat <- x
	ind <- (x == 0)
	x <- x[!ind]
	ex <- exp( - x^2/2)
	z <- w * (x - (2 * (1 - ex))/x)
	z <- z/(w * (1 - ex) + (1 - w) * ex * x^2)
	muhat[!ind] <- z
	return(muhat)
}
"postmean.laplace"<-
function(x, w, a = 0.5)
{
#
# find the posterior mean for the double exponential prior for 
#   given x, w and a, assuming the error variance
#   is 1.
#
#  only allow a < 20. 
	a <- min(a, 20)	#
#  First find the odds of zero and the shrinkage factor
#
	wpost <- 1 - (1 - w)/(1 + w * beta.laplace(x, a))	
	#  now find the posterior mean conditional on being non-zero
#
	sx <- sign(x)
	x <- abs(x)
	cp1 <- pnorm(x - a)
	dp1 <- dnorm(x - a)
	cp2 <- pnorm( - x - a)
	dp2 <- dnorm(x + a)
	ef <- exp(pmin(2 * a * x, 100))
	postmeancond <- ((x - a) * cp1 + dp1 + ef * ((x + a) * cp2 - dp2))/(cp1 +
		ef * cp2)	#
#  calculate posterior mean and return
#
	mutilde <- sx * wpost * postmeancond
	return(mutilde)
}
"postmed"<-
function(x, w, prior = "laplace", a = 0.5)
{
#
# find the posterior median for the appropriate prior for 
#   given x and w and a, assuming the error variance
#   is 1.  
#
	pr <- substring(prior, 1, 1)
	if(pr == "l")
		muhat <- postmed.laplace(x, w, a = a)
	if(pr == "c")
		muhat <- postmed.cauchy(x, w)
	return(muhat)
}
"postmed.cauchy"<-
function(x, w)
{
#
# find the posterior median of the Cauchy prior with
#   mixing weight w, pointwise for each of the data points x
#
	nx <- length(x)
	zest <- rep(NA, length(x))
	w <- rep(w, length.out = nx)
	ax <- abs(x)
	j <- (ax < 20)
	zest[!j] <- ax[!j] - 2/ax[!j]
	if(sum(j) > 0) {
		zest[j] <- vecbinsolv(zf = rep(0, sum(j)), fun = cauchy.medzero,
			tlo = 0, thi = max(ax[j]), z = ax[j], w = w[j])
	}
	zest[zest < 1e-007] <- 0
	zest <- sign(x) * zest
	return(zest)
}
"postmed.laplace"<-
function(x, w, a = 0.5)
{
#
# find the posterior median for the Laplace prior for 
#   given x (possibly vector), w and a, assuming the error variance
#   is 1.
#
#  only allow a < 20. 
	a <- min(a, 20)	#
#  Work with the absolute value of x, and for x > 25 use the approximation
#    to dnorm(x-a)*beta.laplace(x, a)
#
	sx <- sign(x)
	x <- abs(x)
	xma <- x - a
	zz <- (dnorm(x - a) * (1/w + beta.laplace(x, a)))/a
	zz[x > 25] <- 0.5
	mucor <- qnorm(pmin(zz, 1))
	muhat <- sx * pmax(0, x - a - mucor)
	return(muhat)
}
"threshld"<-
function(x, t, hard = TRUE)
{
#
#  threshold the data x using threshold t
#  if hard=T use hard thresholding
#  if hard=F use soft thresholding
	if(hard) z <- x * (abs(x) >= t) else {
		z <- sign(x) * pmax(0, abs(x) - t)
	}
	return(z)
}
"vecbinsolv"<-
function(zf, fun, tlo, thi, nits = 30, ...)
{
#  given a monotone function fun, and a vector of values
#  zf find a vector of numbers t such that f(t) = zf.
#  The solution is constrained to lie on the interval (tlo, thi)
#
#  The function fun may be a vector of increasing functions 
#
#  Present version is inefficient because separate calculations
#  are done for each element of z, and because bisections are done even
#  if the solution is outside the range supplied
#    
#  It is important that fun should work for vector arguments.
#   Additional arguments to fun can be passed through ...
#
#  Works by successive bisection, carrying out nits harmonic bisections
#   of the  interval between tlo and thi
#
	nz <- length(zf)
	tlo <- rep(tlo, nz)
	thi <- rep(thi, nz)	#
#  carry out nits bisections
#
	for(jj in (1:nits)) {
		tmid <- (tlo + thi)/2
		fmid <- fun(tmid, ...)
		indt <- (fmid <= zf)
		tlo[indt] <- tmid[indt]
		thi[!indt] <- tmid[!indt]
	}
	tsol <- (tlo + thi)/2
	return(tsol)
}
"negloglik.laplace"<-
function(xpar, xx)
{
#
#  marginal negative log likelihood function for laplace prior
#  
#  xx data
#  xpar vector of two parameters:
#      xpar[1] :  threshold
#      xpar[2] :  scale factor a
#
	a <- xpar[2]
	w <- wfromt(xpar[1], a = a)
	loglik <- sum(log(1 + w * beta.laplace(xx, a)))
	return( - loglik)
}
"beta.cauchy"<-
function(x)
{
#
#   Find the function beta
#    for the mixed normal prior with Cauchy tails.  It is assumed that the 
#    noise variance is equal to one.  
#
	phix <- dnorm(x)
	j <- (x != 0)
	beta <- x
	beta[!j] <- -1/2
	beta[j] <- (dnorm(0)/phix[j] - 1)/x[j]^2 - 1
	return(beta)
}
"beta.laplace"<-
function(x, a = 0.5)
{
#
#  The function beta for the Laplace prior with parameter a
#
	x <- abs(x)
	a <- min(a, 35)
	xpa <- x + a
	xma <- x - a
	rat1 <- 1/xpa
	rat1[xpa < 35] <- pnorm( - xpa[xpa < 35])/dnorm(xpa[xpa < 35])
	rat2 <- pnorm(xma)/dnorm(xma)
	beta <- (a/2) * (rat1 + rat2) - 1
	return(beta)
}
"laplace.threshzero"<-
function(x, w, a = 0.5)
{
#
# the function that has to be zeroed to find the threshold with the Laplace
#    prior.  
#  only allow a < 20. 
	a <- min(a, 20)	#
	z <- pnorm(x - a) - (dnorm(x - a) * (1/w + beta.laplace(x, a)))/a
	return(z)
}
"cauchy.medzero"<-
function(x, z, w)
{
# the objective function that has to be zeroed, component by component, 
# to find the 
#  posterior median when the quasi-Cauchy prior is used. 
#   x is the parameter vector, z is the data vector, w
# is the weight
#   x and z may be scalars
#
	hh <- z - x
	dnhh <- dnorm(hh)
	yleft <- pnorm(hh) - z * dnhh + ((z * x - 1) * dnhh * pnorm( - x))/
		dnorm(x)
	yright2 <- 1 + exp( - z^2/2) * (z^2 * (1/w - 1) - 1)
	return(yright2/2 - yleft)
}
"cauchy.threshzero"<-
function(z, w)
{
# the objective function that has to be zeroed
# to find the Cauchy
#  threshold.  z is the putative threshold vector, w
# is the weight
#   w can be a vector
#
	y <- pnorm(z) - z * dnorm(z) - 1/2 - (z^2 * exp( - z^2/2) * (1/w - 1))/
		2
	return(y)
}
"wmonfromx"<-
function(xd, prior = "laplace", a = 0.5, tol = 1e-008, maxits = 20)
{
#
#   Find the monotone marginal maximum likelihood estimate of the mixing weights
#    for the Laplace prior with parameter a.  It is assumed that the 
#    noise variance is equal to one.
#
#  Find the beta values and the minimum weight
#
	pr <- substring(prior, 1, 1)
	nx <- length(xd)
	wmin <- wfromt(sqrt(2 * log(length(xd))), prior, a)
	winit <- 1
	if(pr == "l")
		beta <- beta.laplace(xd, a)
	if(pr == "c") beta <- beta.cauchy(xd)	#
#   now conduct iterated weighted least squares isotone regression
#    
	w <- rep(winit, length(beta))
	for(j in (1:maxits)) {
		aa <- w + 1/beta
		ps <- w + aa
		ww <- 1/aa^2
		wnew <- isotone(ps, ww, increasing = FALSE)
		wnew <- pmax(wmin, wnew)
		wnew <- pmin(1, wnew)
		zinc <- max(abs(range(wnew - w)))
		w <- wnew
		if(zinc < tol)
			return(w)
	}
	cat("Warning: more iterations required to achieve convergence \n")
	return(w)
}
"isotone"<-
function(x, wt = rep(1, length(x)), increasing = FALSE)
{
#
#   find the weighted least squares isotone fit to the 
#   sequence x, the weights given by the sequence wt
#
#   if increasing == T the curve is set to be increasing, 
#   otherwise to be decreasing
#
#   the vector ip contains the indices on the original scale of the
#   breaks in the regression at each stage
#
	nn <- length(x)
	if(nn == 1)
		return(x)
	if(!increasing)
		x <-  - x
	ip <- (1:nn)
	dx <- diff(x)
	nx <- length(x)
	while((nx > 1) && (min(dx) < 0)) {
#
#  do single pool-adjacent-violators step
#
#  find all local minima and maxima
#
		jmax <- (1:nx)[c(dx <= 0, FALSE) & c(TRUE, dx > 0)]
		jmin <- (1:nx)[c(dx > 0, TRUE) & c(FALSE, dx <= 0)]	
#	
#  do pav step for each pair of maxima and minima
#
#  add up weights within subsequence that is pooled
#  set first element of subsequence to the weighted average
#  the first weight to the sum of the weights within the subsequence
#    and remainder of the subsequence to NA
#
		for(jb in (1:length(jmax))) {
			ind <- (jmax[jb]:jmin[jb])
			wtn <- sum(wt[ind])
			x[jmax[jb]] <- sum(wt[ind] * x[ind])/wtn
			wt[jmax[jb]] <- wtn
			x[(jmax[jb] + 1):jmin[jb]] <- NA
		}
#
#  clean up within iteration, eliminating the parts of sequences that were set
#   to NA
#
		ind <- !is.na(x)
		x <- x[ind]
		wt <- wt[ind]
		ip <- ip[ind]
		dx <- diff(x)
		nx <- length(x)
	}
# 
#  final cleanup:  reconstruct z at all points by repeating the pooled values
#    the appropriate number of times
#
	jj <- rep(0, nn)
	jj[ip] <- 1
	z <- x[cumsum(jj)]
	if(!increasing)
		z <-  - z
	return(z)
}
"zetafromx"<-
function(xd, cs, pilo=NA, prior = "laplace", a = 0.5)
{
#
#  given a sequence xd, a vector of scale factors cs and
#  a lower limit pilo, find the marginal maximum likelihood
#  estimate of the parameter zeta such that the prior prob
#  is of the form median( pilo, zeta*cs, 1)
#
#  if pilo=NA then it is calculated according to the sample size
#  to corrrespond to the universal threshold
#  
#
#  Find the beta values and the minimum weight if necessary
#
	pr <- substring(prior, 1, 1)
	nx <- length(xd)
	if (is.na(pilo)) pilo <- wfromt(sqrt(2 * log(length(xd))), prior, a)
	if(pr == "l")
		beta <- beta.laplace(xd, a)
	if(pr == "c") beta <- beta.cauchy(xd)
#
#  Find jump points zj in derivative of log likelihood as function
#    of z, and other preliminary calculations
#
	zs1 <- pilo/cs
	zs2 <- 1/cs
	zj <- sort(unique(c(zs1, zs2)))
	cb <- cs * beta
	mz <- length(zj)
	zlmax <- NULL	#
#  Find left and right derivatives at each zj
#   and check which are local minima
#  Check internal zj first
#
	lmin <- rep(F, mz)
	for(j in (2:(mz - 1))) {
		ze <- zj[j]
		cbil <- cb[(ze > zs1) & (ze <= zs2)]
		ld <- sum(cbil/(1 + ze * cbil))
		if(ld <= 0) {
			cbir <- cb[(ze >= zs1) & (ze < zs2)]
			rd <- sum(cbir/(1 + ze * cbir))
			lmin[j] <- (rd >= 0)
		}
	}
#
#  Deal with the two end points in turn, finding right deriv
#   at lower end and left deriv at upper
#
#  In each case the corresponding end point is either a local min
#   or a local max depending on the sign of the relevant deriv
#
	cbir <- cb[zj[1] == zs1]
	rd <- sum(cbir/(1 + zj[1] * cbir))
	if(rd > 0) lmin[1] <- TRUE	else zlmax <- zj[1]
	cbil <- cb[zj[mz] == zs2]
	ld <- sum(cbil/(1 + zj[mz] * cbil))
	if(ld < 0) lmin[mz] <- TRUE else zlmax <- zj[mz]	#
#  Flag all local minima and do a binary search between them to find the local maxima
#
	zlmin <- zj[lmin]
	nlmin <- length(zlmin)
	if(nlmin >= 2) for(j in (2:nlmin)) {
			zlo <- zlmin[j - 1]
			zhi <- zlmin[j]
			ze <- (zlo + zhi)/2
			zstep <- (zhi - zlo)/2
			for(nit in (1:10)) {
				cbi <- cb[(ze >= zs1) & (ze <= zs2)]
				likd <- sum(cbi/(1 + ze * cbi))
				zstep <- zstep/2
				ze <- ze + zstep * sign(likd)
			}
			zlmax <- c(zlmax, ze)
		}
#
#  Evaluate all local maxima and find global max
#
	nlmax <- length(zlmax)
	zm <- rep(NA, nlmax)
	for(j in (1:nlmax)) {
		pz <- pmax(zs1, pmin(zlmax[j], zs2))
		zm[j] <- sum(log(1 + cb * pz))
	}
	zmax <- zlmax[zm == max(zm)]
	return(zmax)
}



