summary.crq <-
function (object, taus = 1:4/5, alpha = .1, se = "boot", ...) 
{
    mt <- terms(object)
    m <- model.frame(object)
    Y <- model.response(m)
    y <- Y[,1]
    cen  <- Y[,2]
    x <- model.matrix(mt, m, contrasts = object$contrasts)
    coef <- coef(object,taus)
    coef <- coef[-nrow(coef),] #delete Qbar row
    B <- boot.crq(x, y, cen, taus, ...)
    sqmn <- sqrt(B$mboot/B$n)
    fact <-   qnorm(1 - alpha/2)/qnorm(.75)
    B <- apply(B$A, 1:2, quantile, probs = 1:3/4, na.rm = TRUE)
    L <- B[2,,] - fact * (B[3,,] - B[2,,]) * sqmn
    U <- B[2,,] - fact * (B[1,,] - B[2,,]) * sqmn
    D <- (U - L)/(2*fact*qnorm(.75))
    T <- coef/D
    P <- 2 * (1 - pnorm(abs(T)))
    G <- list()
    cnames <- c("Value","Lower Bd","Upper Bd","Std Error","T Value","Pr(>|t|)")
    for(i in 1:length(taus)){
	tab <- cbind(coef[,i],L[,i],U[,i],D[,i],T[,i],P[,i])
	dimnames(tab)[[2]] <- cnames
	G[[i]] <- list(tau = taus[i], coefficients = tab)
	}
    class(G) <- "summary.crqs"
    G
   }
print.crq <- function(x, ...)
    print(coef(x, ...), ...)
print.summary.crqs <- function(x, ...)
    lapply(x,print.summary.crq)
print.summary.crq <- function (x, digits = max(5, .Options$digits - 2), ...) {
    coef <- x$coefficients
    tau <- x$tau
    cat("\ntau: ")
    print(format(round(tau, digits = digits)), quote = FALSE, ...)
    cat("\nCoefficients:\n")
    print(format(round(coef, digits = digits)), quote = FALSE, ...)
    invisible(x)
}
