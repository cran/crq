plot.summary.crqs <-
function (x, nrow = 3, ncol = 3, CoxPHit = NULL, ...) {
    taus <- function(x) x$tau
    xx <- unlist(lapply(x, taus))
    coef <- lapply(x, coefficients)
    p <- nrow(coef[[1]])
    k <- ncol(coef[[1]])
    if(k != 6) stop("summary.crqs object has wrong column dimension")
    m <- length(xx)
    blab <- dimnames(coef[[1]])[[1]]
    a <- array(unlist(coef), c(p, k, m))
    if(length(CoxPHit))
	CoxQTE <- QTECox(CoxPHit)
    par(mfrow = c(nrow, ncol))
    for (i in 2:p) {
            b  <- a[i, 1, ]
            bl <- a[i, 2, ]
            bu <- a[i, 3, ]
        plot(rep(xx, 2), c(bl, bu), xlab = "", ylab = "", type = "n")
        title(paste(blab[i]), cex = 0.75)
        polygon(c(xx, rev(xx)), c(bl, rev(bu)), col = "LightSkyBlue")
        points(xx, b, cex = 0.5, pch = "o", col = "blue")
        lines(xx, b, col = "blue")
        abline(h = 0)
        if(length(CoxPHit)) {
	    lines(CoxQTE$taus,CoxQTE$QTE[,i-1],col="red")
            }
    }
}
