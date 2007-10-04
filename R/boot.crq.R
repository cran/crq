boot.crq <- function(x, y, c, taus, R=100,  mboot,  bmethod = "xy-pair", ...) 
	{
	n <- length(y)
	p <- ncol(x)
	mboot <- if(missing(mboot)) mboot <- n
	A <- array(0,dim=c(p+1,length(taus),R))
	
	for (i in 1:R){
		if(bmethod == "xy-pair"){
	    		w <- table(sample(1:n,mboot,replace=TRUE))
    			s <- as.numeric(names(w))
    			w <- as.numeric(w)
			yb <- y[s]
			xb <- x[s,]
			cb <- c[s]
			}
		else if(bmethod == "Bose"){
			w <- rexp(n)
			yb <- y
			xb <- x
			cb <- c
			}
		else
			stop("invalid method for boot.crq") 
	    	a <- crq.fit(xb,yb,cb, weights = w, method = "grid", ... )
		if((i %% floor(R/10)) == 0 ) 
			cat(paste("bootstrap roughly ",100*(i/R)," percent complete\n"))
		A[,,i] <- coef.crq(a,taus)
	        }
	A <- A[-nrow(A),,] # Delete Qbar entries
	list(A = A, n = length(y), mboot = mboot)
	}
