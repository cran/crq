coef.crq <- function(object, taus = 1:4/5, ...)
	{
	# Extract coefficients from the crq solution array 
	
	if(min(taus) < 0 || max(taus) > 1) stop("taus out of range [0,1]")
	taus <- sort(taus)
	S <- object$sol
	r <- S[1, ]
	r <- c(r[1],r)
	r <- (r[-1]+r[-length(r)])/2
	B <- S[-1,]
	J <- length(r)
	ts <- taus[taus < max(r)]
	bin <- cut(ts,r,label=FALSE)
	wgt <- (ts - r[bin])/(r[bin + 1] - r[bin])
	coef <- wgt * B[,bin] + (1-wgt) * B[,bin - 1]
	nna <- length(taus) - length(ts)
	if(nna > 0)
		coef <- cbind(coef, matrix(NA,nrow(coef),nna))
	taulabs <- paste("tau=", format(round(taus, 3)))
	dimnames(coef)[[2]] <- taulabs
	coef
	}
