##' @useDynLib RTMBconvenience, .registration=TRUE
NULL

##' @export
pnorm5 <- function(x, mu, sigma, lower_tail = TRUE, log_p = FALSE) UseMethod("pnorm5")
##' @export
pnorm5.default <- function(x, mu, sigma, lower_tail = TRUE, log_p = FALSE){
    Value(pnorm5_ad(RTMB::advector(x), RTMB::advector(mu), RTMB::advector(sigma), lower_tail, log_p))
}
##' @export
pnorm5.advector <- function(x, mu, sigma, lower_tail = TRUE, log_p = FALSE){
    pnorm5_ad(x, mu, sigma, lower_tail, log_p)
}

##' @export
log_ipnorm <- function(x, y, mu, sigma) UseMethod("log_ipnorm")
##' @export
log_ipnorm.default <- function(x, y, mu, sigma){
    Value(log_ipnorm_ad(RTMB::advector(x),RTMB::advector(y),RTMB::advector(mu),RTMB::advector(sigma)))
}
##' @export
log_ipnorm.advector <- function(x, y, mu, sigma){
    log_ipnorm_ad(x,y,mu,sigma)
}


##' @export
pt2 <- function(x, df, lower_tail = TRUE, log_p = FALSE) UseMethod("pt2")
pt2.default <- function(x, df, lower_tail = TRUE, log_p = FALSE){
    Value(pt_ad(RTMB::advector(x), RTMB::advector(df), lower_tail, log_p))
}
##' @export
pt2.advector <- function(x, df, lower_tail = TRUE, log_p = FALSE){
    pt_ad(x, df, lower_tail, log_p)
}
##' @export 
pde_scheme <- function(x) UseMethod("pde_scheme")
##' @export 
pde_scheme.default <- function(x){
    Value(pde_scheme_ad(RTMB::advector(x)))
}
##' @export 
pde_scheme.advector <- function(x){
    pde_scheme_ad(x)
}

##' @export
logspace_add <- function(x, y) UseMethod("logspace_add")
##' @export
logspace_add.default <- function(x, y){
    Value(logspace_add_ad(RTMB::advector(x),RTMB::advector(y)))
}
##' @export
logspace_add.advector <- function(x, y){
    logspace_add_ad(x,y)
}

        

##' @export
logspace_sum <- function(x) UseMethod("logspace_sum")
##' @export
logspace_sum.default <- function(x){
    Value(logspace_sum_ad(RTMB::advector(x)))
}
##' @export
logspace_sum.advector <- function(x, y){
    logspace_sum_ad(x,y)
}



##' @export
logspace_sub <- function(x, y) UseMethod("logspace_sub")
##' @export
logspace_sub.default <- function(x, y){
    Value(logspace_sub_ad(RTMB::advector(x),RTMB::advector(y)))
}
##' @export
logspace_sub.advector <- function(x, y){
    logspace_sub_ad(x,y)
}




##' @export
quantreg <- function(x, tau) UseMethod("quantreg")
##' @export
quantreg.default <- function(x, tau){
    Value(quantreg_ad(RTMB::advector(x),RTMB::advector(tau)))
}
##' @export
quantreg.advector <- function(x, tau){
    quantreg_ad(x,tau)
}



##' @export
login_log_besselI <- function(logx, nu) UseMethod("login_log_besselI")
##' @export
login_log_besselI.default <- function(logx, nu){
    Value(login_log_besselI_ad(RTMB::advector(logx),RTMB::advector(nu)))
}
##' @export
login_log_besselI.advector <- function(x, y){
    login_log_besselI_ad(x,y)
}



##' @export
log_besselI <- function(x, nu) UseMethod("log_besselI")
##' @export
log_besselI.default <- function(x, nu){
    Value(log_besselI_ad(RTMB::advector(x),RTMB::advector(nu)))
}
##' @export
log_besselI.advector <- function(x, nu){
    log_besselI_ad(x,nu)
}



##' @export
log_MarcumQ <- function(a, b, nu) UseMethod("log_MarcumQ")
##' @export
log_MarcumQ.default <- function(a, b, nu){
    Value(log_MarcumQ_ad(RTMB::advector(a), RTMB::advector(b), RTMB::advector(nu)))
}
##' @export
log_MarcumQ.advector <- function(a, b, nu){
    log_MarcumQ_ad(a,b, nu)
}


##' @export
login_log_MarcumQ <- function(loga, logb, nu) UseMethod("login_log_MarcumQ")
##' @export
login_log_MarcumQ.default <- function(loga, logb, nu){
    Value(login_log_MarcumQ_ad(RTMB::advector(loga), RTMB::advector(logb), RTMB::advector(nu)))
}
##' @export
login_log_MarcumQ.advector <- function(loga, logb, nu){
    login_log_MarcumQ_ad(loga,logb, nu)
}


##' @export
log_Marcum1mQ <- function(a, b, nu) UseMethod("log_Marcum1mQ")
##' @export
log_Marcum1mQ.default <- function(a, b, nu){
    Value(log_Marcum1mQ_ad(RTMB::advector(a), RTMB::advector(b), RTMB::advector(nu)))
}
##' @export
log_Marcum1mQ.advector <- function(a, b, nu){
    log_Marcum1mQ_ad(a,b, nu)
}


##' @export
login_log_Marcum1mQ <- function(loga, logb, nu) UseMethod("login_log_Marcum1mQ")
##' @export
login_log_Marcum1mQ.default <- function(loga, logb, nu){
    Value(login_log_Marcum1mQ_ad(RTMB::advector(loga), RTMB::advector(logb), RTMB::advector(nu)))
}
##' @export
login_log_Marcum1mQ.advector <- function(loga, logb, nu){
    login_log_Marcum1mQ_ad(loga,logb, nu)
}

##' @export
logdrice <- function(logx, lognu, logsigma) UseMethod("logdrice")
##' @export
logdrice.default <- function(logx, lognu, logsigma){
    Value(logdrice_ad(RTMB::advector(logx), RTMB::advector(lognu), RTMB::advector(logsigma)))
}
##' @export
logdrice.advector <- function(logx, lognu, logsigma){
    logdrice_ad(logx, lognu, logsigma)
}


##' @export
logprice <- function(logx, lognu, logsigma, lower_tail = TRUE) UseMethod("logprice")
##' @export
logprice.default <- function(logx, lognu, logsigma, lower_tail = TRUE){
    Value(logprice_ad(RTMB::advector(logx), RTMB::advector(lognu), RTMB::advector(logsigma), lower_tail))
}
##' @export
logprice.advector <- function(logx, lognu, logsigma, lower_tail = TRUE){
    logprice_ad(logx, lognu, logsigma, lower_tail)
}

