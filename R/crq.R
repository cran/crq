crq <- function (formula, data, subset, weights, na.action,
    method = "grid", contrasts = NULL, ...)
{
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
        names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    if (method == "model.frame")
        return(mf)
    mt <- attr(mf, "terms")
    Y <- model.extract(mf, "response")
    if(!inherits(Y,"Surv")) 
	stop("Response must be a survival object")
    if(attr(Y,"type") != "right")
	stop("Only right censoring Surv objects are allowed")
    y <- Y[,1]
    cen <-  Y[,2]
    x <- model.matrix(mt, mf, contrasts)
    weights <- model.weights(mf)
    fit <- crq.fit(x, y, cen, weights, method, ...)  
 

class(fit) <- "crq"
fit$terms <- mt
fit$call <- call
fit$formula <-  formula(mt) 
attr(fit, "na.message") <- attr(m, "na.message")
fit
}
