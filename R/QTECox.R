
QTECox <- function(x, smooth = TRUE){
# compute quantile treatment effect for a Cox PH fit
g <- survfit(x)
if(smooth)
	g <- supsmu(1-g$surv,g$time)
taus <- (g$x[-1] + g$x[-length(g$x)])/2
Qhat <- (g$y[-1] + g$y[-length(g$y)])/2

dQ <- diff(Qhat)/diff(taus)
taus <- taus[-1]; Qhat <- Qhat[-1]
QTE <- outer(dQ * (1-taus) * log(1-taus) / Qhat, coef(x))
list(QTE = QTE , taus = taus)
}
