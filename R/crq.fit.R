crq.fit <- function(x, y, cen, weights, method, mw = 20, ginit = 0.001, gstep = NULL)
{
      p <- ncol(x)
      n <- length(y) 
      cen <- 1 - cen #!!! Fortran routine wants censoring indicator flipped.
      mp <- n+5+sum(cen)
      gstep <- ifelse(is.null(gstep),min(0.01,1/(2*n^{.7})),gstep)
      if(length(weights)){
		if (any(weights < 0)) 
		        stop("negative weights not allowed")
    		contr <- attr(x, "contrasts")
    		x <- x * weights
    		y <- y * weights
		}
      if(method == "pivot"){
		nsol <- 3*n
		if(mw < 1) stop("mw must be positive for pivot method")
		}
      else if(method == "grid"){
		#nsol <- length(seq(ginit,1,by=gstep)) #too small for some reason!
		nsol <- 3*n
		mw = -1
		}
      else 
		stop("Unrecognized crq estimation method")
	z <- .Fortran("crq",
		as.integer(n),
		as.integer(p),
                as.integer(mp), 
		as.integer(p+2),
		as.double(x),
		as.double(y),
		as.integer(cen),
		as.double(ginit),
                as.integer(mw),
                as.double(gstep),
		ift = integer(1),
		h = integer(p),
                xh = double(p*p),
		wa = double(mp*p),
		wb = double(mp),
		wc = double(mp*(p+2)),
		wd = double(mp),
		we = double(mp),
		wf = double(p),
		iflag = integer(mp),
                as.integer(nsol),
		sol = double(nsol*(p+2)),
                lsol = integer(1),
		icen = integer(n),
		tcen = double(n),
		lcen = integer(1),
		PACKAGE = "crq")
	nw <- z$h[1]
	flag <- z$ift
	msg <- switch(flag,
		paste("Error in input dimensions, n,p,mw "),
		paste("Error in input dimensions, n,p,mw "),
		paste("Error in input dimensions, n,p,mw "),
		paste("Less than p=",p,"observations above tau = 0 solution"),
		paste("Possible degeneracy at",nw,"tau values.",
          		"$tau.degen: first mp =", n + 5 + sum(cen)," such tau values"), 
		paste("Number of pivots to be saved in sol > nsol.",
          		"Redefine nsol: use nsol < n to save for tau = i/(nsol-1)"),
		paste("Error with partial return: possible degeneracies",
         		"Max number of rq calls exceeded: dither x or increase mw"),
		paste("Premature stop: degeneracy or infinite basis element in rq"))
	if(flag > 0)
		ifelse(flag <= 3,stop(msg),warning(msg))
	J <- z$lsol
	B <- matrix(z$sol, nrow=p+2, ncol=nsol, byrow=FALSE)[,1:J]
	dimnames(B) <- list(c("tau",dimnames(x)[[2]],"Qbar"),NULL)
	ic <- z$icen
	sp <- (1:n)[ic == 1]
	tsp <- z$tcen[sp]
        t1 <- z$wd[1:nw]
	a<-list(sol=B, Isplit=sp, tausplit=tsp,status=ic)
	return(a) 
	}
