
splineGenericEval <- function(FUN, x, knots, pars){
    splFun <- match.fun(paste("spline",FUN, sep="_"))
    anyAD <- inherits(x,"advector") || inherits(pars,"advector")
    ans <- splFun(advector(x), knots, advector(pars))
    if(anyAD)
        return(ans)
    Value(ans)
}

##' @export
bcspline <- function(x, knots, pars){
    splineGenericEval("bcspline",x=x, knots=knots, pars=pars)
}

##' @export
ibcspline <- function(x, knots, pars){
    splineGenericEval("ibcspline",x=x, knots=knots, pars=pars)
}

##' @export
ibcdspline <- function(x, knots, pars){
    splineGenericEval("ibcdspline",x=x, knots=knots, pars=pars)
}


##' @export
ibcispline <- function(x, knots, pars){
    splineGenericEval("ibcispline",x=x, knots=knots, pars=pars)
}


##' @export
iibcspline <- function(x, knots, pars){
    splineGenericEval("iibcspline",x=x, knots=knots, pars=pars)
}

##' @export
iibcdspline <- function(x, knots, pars){
    splineGenericEval("iibcdspline",x=x, knots=knots, pars=pars)
}


##' @export
iibcispline <- function(x, knots, pars){
    splineGenericEval("iibcispline",x=x, knots=knots, pars=pars)
}
