##' @useDynLib RTMBconvenience, .registration=TRUE
NULL

## ##' @export
## pnorm5 <- function(x, mu, sigma, lower_tail = TRUE, log_p = FALSE) UseMethod("pnorm5")
## ##' @export
## pnorm5.default <- function(x, mu, sigma, lower_tail = TRUE, log_p = FALSE){
##     Value(pnorm5_ad(RTMB::advector(x), RTMB::advector(mu), RTMB::advector(sigma), lower_tail, log_p))
## }
## ##' @export
## pnorm5.advector <- function(x, mu, sigma, lower_tail = TRUE, log_p = FALSE){
##     pnorm5_ad(x, mu, sigma, lower_tail, log_p)
## }

##' @export
setGeneric("pnorm5", function(x,mu,sigma,lower_tail=TRUE,log_p=FALSE){
              (pnorm5_ad(advector(x),advector(mu),advector(sigma), lower_tail, log_p))
          })
##' @export
setMethod("pnorm5",
          signature(x = "num", mu = "num", sigma="num",lower_tail="ANY",log_p="ANY"),
          function(x,mu,sigma,lower_tail=TRUE,log_p=FALSE){
              Value(pnorm5_ad(advector(x),advector(mu),advector(sigma), lower_tail, log_p))
          })
##' @export
setMethod("pnorm5",
          signature(x = "ad", mu="ad", sigma="ad",lower_tail="ANY",log_p="ANY"),
          function(x,mu,sigma,lower_tail=TRUE,log_p=FALSE){
              pnorm5_ad(advector(x),advector(mu),advector(sigma), lower_tail, log_p)
          })
##' @export
setMethod("pnorm5",
          signature(x = "advector",mu="advector",sigma="advector",lower_tail="ANY",log_p="ANY"),
          function(x,mu,sigma,lower_tail=TRUE,log_p=FALSE){
              pnorm5_ad(x,mu,sigma, lower_tail, log_p)
          })


##' @export
setGeneric("qnorm5", function(x,mu,sigma,lower_tail=TRUE,log_p=FALSE){
              (qnorm5_ad(advector(x),advector(mu),advector(sigma), lower_tail, log_p))
          })
##' @export
setMethod("qnorm5",
          signature(x = "num", mu = "num", sigma="num",lower_tail="ANY",log_p="ANY"),
          function(x,mu,sigma,lower_tail=TRUE,log_p=FALSE){
              Value(qnorm5_ad(advector(x),advector(mu),advector(sigma), lower_tail, log_p))
          })
##' @export
setMethod("qnorm5",
          signature(x = "ad", mu="ad", sigma="ad",lower_tail="ANY",log_p="ANY"),
          function(x,mu,sigma,lower_tail=TRUE,log_p=FALSE){
              qnorm5_ad(advector(x),advector(mu),advector(sigma), lower_tail, log_p)
          })
##' @export
setMethod("qnorm5",
          signature(x = "advector",mu="advector",sigma="advector",lower_tail="ANY",log_p="ANY"),
          function(x,mu,sigma,lower_tail=TRUE,log_p=FALSE){
              qnorm5_ad(x,mu,sigma, lower_tail, log_p)
          })


## ##' @export
## log_ipnorm <- function(x, y, mu, sigma) UseMethod("log_ipnorm")
## ##' @export
## log_ipnorm.default <- function(x, y, mu, sigma){
##     Value(log_ipnorm_ad(RTMB::advector(x),RTMB::advector(y),RTMB::advector(mu),RTMB::advector(sigma)))
## }
## ##' @export
## log_ipnorm.advector <- function(x, y, mu, sigma){
##     log_ipnorm_ad(x,y,mu,sigma)
## }

##' @export
setGeneric("log_ipnorm", function(x,y,mu,sigma){
              (log_ipnorm_ad(advector(x),advector(y),advector(mu),advector(sigma)))
          })
##' @export
setMethod("log_ipnorm",
          signature(x = "num", y = "num", mu = "num", sigma="num"),
          function(x, y, mu, sigma){
              Value(log_ipnorm_ad(advector(x),advector(y),advector(mu),advector(sigma)))
          })
##' @export
setMethod("log_ipnorm",
          signature(x = "ad", y = "ad", mu="ad", sigma="ad"),
          function(x,y,mu,sigma){
              log_ipnorm_ad(advector(x),advector(y),advector(mu),advector(sigma))
          })
##' @export
setMethod("log_ipnorm",
          signature(x = "advector",y = "advector",mu="advector",sigma="advector"),
          function(x,y,mu,sigma){
              log_ipnorm_ad(x,y,mu,sigma)
          })

        


## ##' @export
## pt2 <- function(x, df, lower_tail = TRUE, log_p = FALSE) UseMethod("pt2")
## ##' @export
## pt2.default <- function(x, df, lower_tail = TRUE, log_p = FALSE){
##     Value(pt_ad(RTMB::advector(x), RTMB::advector(df), lower_tail, log_p))
## }
## ##' @export
## pt2.advector <- function(x, df, lower_tail = TRUE, log_p = FALSE){
##     pt_ad(x, df, lower_tail, log_p)
## }

##' @export
setGeneric("pt2", function(x,df,lower_tail=TRUE,log_p=FALSE){
              (pt_ad(advector(x),advector(df), lower_tail, log_p))
          })
##' @export
setMethod("pt2",
          signature(x = "num", df = "num", lower_tail="ANY",log_p="ANY"),
          function(x,df,lower_tail=TRUE,log_p=FALSE){
              Value(pt_ad(advector(x),advector(df), lower_tail, log_p))
          })
##' @export
setMethod("pt2",
          signature(x = "ad", df = "ad", lower_tail="ANY",log_p="ANY"),
          function(x,df,lower_tail=TRUE,log_p=FALSE){
              pt_ad(advector(x),advector(df), lower_tail, log_p)
          })
##' @export
setMethod("pt2",
          signature(x = "advector",df = "advector", lower_tail="ANY",log_p="ANY"),
          function(x,df,lower_tail=TRUE,log_p=FALSE){
              pt_ad(x,df, lower_tail, log_p)
          })


##' @export
setGeneric("qt2", function(p,df,lower_tail=TRUE,log_p=FALSE){
              (qt_ad(advector(p),advector(df), lower_tail, log_p))
          })
##' @export
setMethod("qt2",
          signature(p = "num", df = "num", lower_tail="ANY",log_p="ANY"),
          function(p,df,lower_tail=TRUE,log_p=FALSE){
              Value(qt_ad(advector(p),advector(df), lower_tail, log_p))
          })
##' @export
setMethod("qt2",
          signature(p = "ad", df = "ad", lower_tail="ANY",log_p="ANY"),
          function(p,df,lower_tail=TRUE,log_p=FALSE){
              qt_ad(advector(p),advector(df), lower_tail, log_p)
          })
##' @export
setMethod("qt2",
          signature(p = "advector",df = "advector", lower_tail="ANY",log_p="ANY"),
          function(p,df,lower_tail=TRUE,log_p=FALSE){
              qt_ad(p,df, lower_tail, log_p)
          })



## ##' @export 
## pde_scheme <- function(x) UseMethod("pde_scheme")
## ##' @export 
## pde_scheme.default <- function(x){
##     Value(pde_scheme_ad(RTMB::advector(x)))
## }
## ##' @export 
## pde_scheme.advector <- function(x){
##     pde_scheme_ad(x)
## }

##' @export
setGeneric("logspace_add", function(x,y){
              (logspace_add_ad(advector(x),advector(y)))
          })
##' @export
setMethod("logspace_add",
          signature(x = "num", y = "num"),
          function(x, y){
              Value(logspace_add_ad(advector(x),advector(y)))
          })
##' @export
setMethod("logspace_add",
          signature(x = "ad", y = "ad"),
          function(x,y){
              logspace_add_ad(advector(x),advector(y))
          })
##' @export
setMethod("logspace_add",
          signature(x = "advector",y = "advector"),
          function(x,y){
              logspace_add_ad(x,y)
          })

        

##' @export
setGeneric("logspace_sum", function(x){
              (logspace_sum_ad(advector(x)))
          })
##' @export
setMethod("logspace_sum",
          signature(x = "num"),
          function(x){
              Value(logspace_sum_ad(advector(x)))
          })
##' @export
setMethod("logspace_sum",
          signature(x = "ad"),
          function(x){
              logspace_sum_ad(advector(x))
          })
##' @export
setMethod("logspace_sum",
          signature(x = "advector"),
          function(x){
              logspace_sum_ad(x)
          })


## ##' @export
## logspace_sub <- function(x, y) UseMethod("logspace_sub")
## ##' @export
## logspace_sub.default <- function(x, y){
##     Value(logspace_sub_ad(RTMB::advector(x),RTMB::advector(y)))
## }
## ##' @export
## logspace_sub.advector <- function(x, y){
##     logspace_sub_ad(x,y)
## }
##' @export
setGeneric("logspace_sub", function(x,y){
              (logspace_sub_ad(advector(x),advector(y)))
          })
##' @export
setMethod("logspace_sub",
          signature(x = "num", y = "num"),
          function(x, y){
              Value(logspace_sub_ad(advector(x),advector(y)))
          })
##' @export
setMethod("logspace_sub",
          signature(x = "ad", y = "ad"),
          function(x,y){
              logspace_sub_ad(advector(x),advector(y))
          })
##' @export
setMethod("logspace_sub",
          signature(x = "advector",y = "advector"),
          function(x,y){
              logspace_sub_ad(x,y)
          })




##' @export
setGeneric("quantreg_loss", function(x,tau){
              (quantreg_ad(advector(x),advector(tau)))
          })
##' @export
setMethod("quantreg_loss",
          signature(x = "num", tau = "num"),
          function(x, tau){
              Value(quantreg_ad(advector(x),advector(tau)))
          })
##' @export
setMethod("quantreg_loss",
          signature(x = "ad", tau = "ad"),
          function(x,tau){
              quantreg_ad(advector(x),advector(tau))
          })
##' @export
setMethod("quantreg_loss",
          signature(x = "advector",tau = "advector"),
          function(x,tau){
              quantreg_ad(x,tau)
          })



## ##' @export
## login_log_besselI <- function(logx, nu) UseMethod("login_log_besselI")
## ##' @export
## login_log_besselI.default <- function(logx, nu){
##     Value(login_log_besselI_ad(RTMB::advector(logx),RTMB::advector(nu)))
## }
## ##' @export
## login_log_besselI.advector <- function(x, y){
##     login_log_besselI_ad(x,y)
## }

##' @export
setGeneric("login_log_besselI", function(logx,nu){
              (login_log_besselI_ad(advector(logx),advector(nu)))
          })
##' @export
setMethod("login_log_besselI",
          signature(logx = "num", nu = "num"),
          function(logx, nu){
              Value(login_log_besselI_ad(advector(logx),advector(nu)))
          })
##' @export
setMethod("login_log_besselI",
          signature(logx = "ad", nu = "ad"),
          function(logx,nu){
              login_log_besselI_ad(advector(logx),advector(nu))
          })
##' @export
setMethod("login_log_besselI",
          signature(logx = "advector",nu = "advector"),
          function(logx,nu){
              login_log_besselI_ad(logx,nu)
          })



## ##' @export
## log_besselI <- function(x, nu) UseMethod("log_besselI")
## ##' @export
## log_besselI.default <- function(x, nu){
##     Value(log_besselI_ad(RTMB::advector(x),RTMB::advector(nu)))
## }
## ##' @export
## log_besselI.advector <- function(x, nu){
##     log_besselI_ad(x,nu)
## }
##' @export
setGeneric("log_besselI", function(x,nu){
              (log_besselI_ad(advector(x),advector(nu)))
          })
##' @export
setMethod("log_besselI",
          signature(x = "num", nu = "num"),
          function(x, nu){
              Value(log_besselI_ad(advector(x),advector(nu)))
          })
##' @export
setMethod("log_besselI",
          signature(x = "ad", nu = "ad"),
          function(x,nu){
              log_besselI_ad(advector(x),advector(nu))
          })
##' @export
setMethod("log_besselI",
          signature(x = "advector",nu = "advector"),
          function(x,nu){
              log_besselI_ad(x,nu)
          })



## ##' @export
## log_MarcumQ <- function(a, b, nu) UseMethod("log_MarcumQ")
## ##' @export
## log_MarcumQ.default <- function(a, b, nu){
##     Value(log_MarcumQ_ad(RTMB::advector(a), RTMB::advector(b), RTMB::advector(nu)))
## }
## ##' @export
## log_MarcumQ.advector <- function(a, b, nu){
##     log_MarcumQ_ad(a,b, nu)
## }
##' @export
setGeneric("log_MarcumQ", function(a,b,nu){
              (log_MarcumQ_ad(advector(a),advector(b),advector(nu)))
          })
##' @export
setMethod("log_MarcumQ",
          signature(a = "num",b="num", nu = "num"),
          function(a,b, nu){
              Value(log_MarcumQ_ad(advector(a),advector(b),advector(nu)))
          })
##' @export
setMethod("log_MarcumQ",
          signature(a = "ad", b="ad", nu = "ad"),
          function(a,b,nu){
              log_MarcumQ_ad(advector(a),advector(b),advector(nu))
          })
##' @export
setMethod("log_MarcumQ",
          signature(a = "advector", b = "advector", nu = "advector"),
          function(a,b,nu){
              log_MarcumQ_ad(a,b,nu)
          })


##' @export
## login_log_MarcumQ <- function(loga, logb, nu) UseMethod("login_log_MarcumQ")
## ##' @export
## login_log_MarcumQ.default <- function(loga, logb, nu){
##     Value(login_log_MarcumQ_ad(RTMB::advector(loga), RTMB::advector(logb), RTMB::advector(nu)))
## }
## ##' @export
## login_log_MarcumQ.advector <- function(loga, logb, nu){
##     login_log_MarcumQ_ad(loga,logb, nu)
## }
##' @export
setGeneric("login_log_MarcumQ", function(loga,logb,nu){
              (login_log_MarcumQ_ad(advector(loga),advector(logb),advector(nu)))
          })
##' @export
setMethod("login_log_MarcumQ",
          signature(loga = "num",logb="num", nu = "num"),
          function(loga,logb, nu){
              Value(login_log_MarcumQ_ad(advector(loga),advector(logb),advector(nu)))
          })
##' @export
setMethod("login_log_MarcumQ",
          signature(loga = "ad", logb="ad", nu = "ad"),
          function(loga,logb,nu){
              login_log_MarcumQ_ad(advector(loga),advector(logb),advector(nu))
          })
##' @export
setMethod("login_log_MarcumQ",
          signature(loga = "advector", logb = "advector", nu = "advector"),
          function(loga,logb,nu){
              login_log_MarcumQ_ad(loga,logb,nu)
          })


## ##' @export
## log_Marcum1mQ <- function(a, b, nu) UseMethod("log_Marcum1mQ")
## ##' @export
## log_Marcum1mQ.default <- function(a, b, nu){
##     Value(log_Marcum1mQ_ad(RTMB::advector(a), RTMB::advector(b), RTMB::advector(nu)))
## }
## ##' @export
## log_Marcum1mQ.advector <- function(a, b, nu){
##     log_Marcum1mQ_ad(a,b, nu)
## }
##' @export
setGeneric("log_Marcum1mQ", function(a,b,nu){
              (log_Marcum1mQ_ad(advector(a),advector(b),advector(nu)))
          })
##' @export
setMethod("log_Marcum1mQ",
          signature(a = "num",b="num", nu = "num"),
          function(a,b, nu){
              Value(log_Marcum1mQ_ad(advector(a),advector(b),advector(nu)))
          })
##' @export
setMethod("log_Marcum1mQ",
          signature(a = "ad", b="ad", nu = "ad"),
          function(a,b,nu){
              log_Marcum1mQ_ad(advector(a),advector(b),advector(nu))
          })
##' @export
setMethod("log_Marcum1mQ",
          signature(a = "advector", b = "advector", nu = "advector"),
          function(a,b,nu){
              log_Marcum1mQ_ad(a,b,nu)
          })


## ##' @export
## login_log_Marcum1mQ <- function(loga, logb, nu) UseMethod("login_log_Marcum1mQ")
## ##' @export
## login_log_Marcum1mQ.default <- function(loga, logb, nu){
##     Value(login_log_Marcum1mQ_ad(RTMB::advector(loga), RTMB::advector(logb), RTMB::advector(nu)))
## }
## ##' @export
## login_log_Marcum1mQ.advector <- function(loga, logb, nu){
##     login_log_Marcum1mQ_ad(loga,logb, nu)
## }
##' @export
setGeneric("login_log_Marcum1mQ", function(loga,logb,nu){
              (login_log_Marcum1mQ_ad(advector(loga),advector(logb),advector(nu)))
          })
##' @export
setMethod("login_log_Marcum1mQ",
          signature(loga = "num",logb="num", nu = "num"),
          function(loga,logb, nu){
              Value(login_log_Marcum1mQ_ad(advector(loga),advector(logb),advector(nu)))
          })
##' @export
setMethod("login_log_Marcum1mQ",
          signature(loga = "ad", logb="ad", nu = "ad"),
          function(loga,logb,nu){
              login_log_Marcum1mQ_ad(advector(loga),advector(logb),advector(nu))
          })
##' @export
setMethod("login_log_Marcum1mQ",
          signature(loga = "advector", logb = "advector", nu = "advector"),
          function(loga,logb,nu){
              login_log_Marcum1mQ_ad(loga,logb,nu)
          })

## ##' @export
## logdrice <- function(logx, lognu, logsigma) UseMethod("logdrice")
## ##' @export
## logdrice.default <- function(logx, lognu, logsigma){
##     Value(logdrice_ad(RTMB::advector(logx), RTMB::advector(lognu), RTMB::advector(logsigma)))
## }
## ##' @export
## logdrice.advector <- function(logx, lognu, logsigma){
##     logdrice_ad(logx, lognu, logsigma)
## }
setGeneric("logdrice", function(logx,lognu,logsigma){
              (logdrice_ad(advector(logx),advector(lognu),advector(logsigma)))
          })
##' @export
setMethod("logdrice",
          signature(logx = "num",lognu="num", logsigma = "num"),
          function(logx,lognu, logsigma){
              Value(logdrice_ad(advector(logx),advector(lognu),advector(logsigma)))
          })
##' @export
setMethod("logdrice",
          signature(logx = "ad", lognu="ad", logsigma = "ad"),
          function(logx,lognu,logsigma){
              logdrice_ad(advector(logx),advector(lognu),advector(logsigma))
          })
##' @export
setMethod("logdrice",
          signature(logx = "advector", lognu = "advector", logsigma = "advector"),
          function(logx,lognu,logsigma){
              logdrice_ad(logx,lognu,logsigma)
          })


## ##' @export
## logprice <- function(logx, lognu, logsigma, lower_tail = TRUE) UseMethod("logprice")
## ##' @export
## logprice.default <- function(logx, lognu, logsigma, lower_tail = TRUE){
##     Value(logprice_ad(RTMB::advector(logx), RTMB::advector(lognu), RTMB::advector(logsigma), lower_tail))
## }
## ##' @export
## logprice.advector <- function(logx, lognu, logsigma, lower_tail = TRUE){
##     logprice_ad(logx, lognu, logsigma, lower_tail)
## }
##' @export
setGeneric("logprice", function(logx,lognu,logsigma, lower_tail = TRUE){
              (logprice_ad(advector(logx),advector(lognu),advector(logsigma), lower_tail))
          })
##' @export
setMethod("logprice",
          signature(logx = "num",lognu="num", logsigma = "num", lower_tail = "ANY"),
          function(logx,lognu, logsigma, lower_tail = TRUE){
              Value(logprice_ad(advector(logx),advector(lognu),advector(logsigma), lower_tail))
          })
##' @export
setMethod("logprice",
          signature(logx = "ad", lognu="ad", logsigma = "ad", lower_tail = "ANY"),
          function(logx,lognu,logsigma, lower_tail = TRUE){
              logprice_ad(advector(logx),advector(lognu),advector(logsigma), lower_tail)
          })
##' @export
setMethod("logprice",
          signature(logx = "advector", lognu = "advector", logsigma = "advector", lower_tail = "ANY"),
          function(logx,lognu,logsigma, lower_tail = TRUE){
              logprice_ad(logx,lognu,logsigma, lower_tail)
          })



##' @export
setGeneric("pnchisq", function(x,df,ncp, lower_tail = TRUE, log_p=FALSE){
              (pnchisq_ad(advector(x),advector(df),advector(ncp), lower_tail, log_p))
          })
##' @export
setMethod("pnchisq",
          signature(x = "num",df="num", ncp = "num", lower_tail = "ANY", log_p="ANY"),
          function(x,df,ncp, lower_tail = TRUE, log_p=FALSE){
              Value(pnchisq_ad(advector(x),advector(df),advector(ncp), lower_tail, log_p))
          })
##' @export
setMethod("pnchisq",
          signature(x = "ad",df="ad", ncp = "ad", lower_tail = "ANY", log_p="ANY"),
          function(x,df,ncp, lower_tail = TRUE, log_p=FALSE){
              pnchisq_ad(advector(x),advector(df),advector(ncp), lower_tail, log_p)
          })
##' @export
setMethod("pnchisq",
          signature(x = "advector",df="advector", ncp = "advector", lower_tail = "ANY", log_p="ANY"),
          function(x,df,ncp, lower_tail = TRUE, log_p=FALSE){
              pnchisq_ad(advector(x),advector(df),advector(ncp), lower_tail, log_p)
          })





##' @export
setGeneric("pde_scheme", function(x){
              pde_scheme_ad(advector(x))
          })
##' @export
setMethod("pde_scheme",
          signature(x = "num"),
          function(x){
              Value(pde_scheme_ad(advector(x)))
          })
##' @export
setMethod("pde_scheme",
          signature(x = "ad"),
          function(x){
              pde_scheme_ad(advector(x))
          })
##' @export
setMethod("pde_scheme",
          signature(x = "advector"),
          function(x){
              pde_scheme_ad(x)
          })

##' @export
sparse_solve <- function(x,b){
    if(is(x,"adsparse")){
        return(sparse_solve_ad(x,advector(b)))
    }else{
        return(sparse_solve_NotAD(as(x,"sparseMatrix"),b))
    }
}




##' @export
setGeneric("logspace_1m", function(x){
              logspace_1m_ad(advector(x))
          })
##' @export
setMethod("logspace_1m",
          signature(x = "num"),
          function(x){
              Value(logspace_1m_ad(advector(x)))
          })
##' @export
setMethod("logspace_1m",
          signature(x = "ad"),
          function(x){
              logspace_1m_ad(advector(x))
          })
##' @export
setMethod("logspace_1m",
          signature(x = "advector"),
          function(x){
              logspace_1m_ad(x)
          })


##' @export
setGeneric("logspace_1p", function(x){
              logspace_1p_ad(advector(x))
          })
##' @export
setMethod("logspace_1p",
          signature(x = "num"),
          function(x){
              Value(logspace_1p_ad(advector(x)))
          })
##' @export
setMethod("logspace_1p",
          signature(x = "ad"),
          function(x){
              logspace_1p_ad(advector(x))
          })
##' @export
setMethod("logspace_1p",
          signature(x = "advector"),
          function(x){
              logspace_1p_ad(x)
          })




##' @export
setGeneric("continuationRatioLogitToLogProbability", function(x){
              continuationRatioLogitToLogProbability_ad(advector(x))
          })
##' @export
setMethod("continuationRatioLogitToLogProbability",
          signature(x = "num"),
          function(x){
              Value(continuationRatioLogitToLogProbability_ad(advector(x)))
          })
##' @export
setMethod("continuationRatioLogitToLogProbability",
          signature(x = "ad"),
          function(x){
              continuationRatioLogitToLogProbability_ad(advector(x))
          })
##' @export
setMethod("continuationRatioLogitToLogProbability",
          signature(x = "advector"),
          function(x){
              continuationRatioLogitToLogProbability_ad(x)
          })
