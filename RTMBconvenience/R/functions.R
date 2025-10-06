

##' @export
simplify <- function(x){
    if(!is.list(x)) return(x)
    if(length(x) == 1) return(x[[1]])
    r <- AD(numeric(length(x)))
    for(i in seq_along(x))
        r[i] <- AD(x[[i]])
    r
}

##' @export
## Value <- function(x){
##     if(is(x,"advector"))
##         return(Im(unclass(x)))
##     x
## }
Value <- function(x){
    if(is(x,"advector"))
        return(RTMB:::getValues(x))
    x
}

##' @importFrom RTMB cbind.advector
##' @export
safe_cbind <- function(...){
    args <- list(...)
    isAD <- sapply(args,function(x) is(x,"advector"))
    if(any(isAD))
        return(do.call(RTMB::cbind.advector,args))
    return(do.call(cbind,args))
}

##' @importFrom RTMB rbind.advector
##' @export
safe_rbind <- function(...){
    args <- list(...)
    isAD <- sapply(args,function(x) is(x,"advector"))
    if(any(isAD))
        return(do.call(RTMB::rbind.advector,args))
    return(do.call(rbind,args))
}

##' @export
safe_apply <- function(x, MARGIN, FUN, ...){
    if(is.data.frame(x)){
        if(MARGIN == 1){
            xl <- split(x, seq_len(nrow(x)))
            return(lapply(xl, FUN, ...))
        }else if(MARGIN == 2){
            xl <- as.list(x)
            return(lapply(xl, FUN, ...))
        }else{
            stop("Wrong margin for a data frame")
        }
    }
    indx <- slice.index(x, MARGIN)
    xl <- split(x, indx)
    lapply(xl, FUN, ...)
}

##' @export
safe_aggregate <- function(x, by, FUN){
    r0 <- lapply(split(x,do.call("paste",c(by,list(sep=":")))),FUN)
    cl <- lapply(by,class)
    ii <- as.list(as.data.frame(do.call("rbind",strsplit(names(r0),":"))))
    ii <- lapply(seq_along(ii), function(qq) as(ii[[qq]],cl[[qq]]))
    if(!is.null(names(by)))
        names(ii) <- names(by)
    oo <- do.call(order,rev(ii))
    c(list(x = simplify(r0[oo])),lapply(ii,function(y)y[oo]))
}


##' @export
squeeze <- function(u, eps = .Machine$double.eps){
    (1.0 - eps) * (u - 0.5) + 0.5
}

##' @export
lgamma <- function(x) UseMethod("lgamma")

##' @export
lgamma.default <- function(x){
    (.Primitive("lgamma"))(x)    
}
##' @export
lgamma.advector <- function(x){
    RTMB:::Math1(x,"lgamma")
}


##' @export
toRowLogPropMatrix <- function(x){
    y <- cbind(x,0)
    ys <- simplify(safe_apply(y,1,logspace_sum))
    y - ys[row(y)]
}

##' @export
logspace_1m <- function(logx){
    log1p(-exp(logx))
}

##' @export
undim <- function(x){
    dim(x) <- NULL
    x
}
