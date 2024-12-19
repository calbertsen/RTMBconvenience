#' @export
makeAdvDiff <- function(x0,f,g){
    D <- function(X,Time){
        v <- g(X,Time)
        0.5 * v %*% t(v)
    }
    J_tp <- MakeTape(function(p) D(tail(p,-1),head(p,1)),c(0,x0))$jacfun()$atomic()
    J <- function(X,Time){
        RTMB::matrix(J_tp(AD(c(Time,X))),ncol=length(X)+length(Time))[,-1]
    }
    mu <- function(X,Time){
        vJ <- array(J(AD(X),AD(Time)),dim = c(length(X),length(X),length(X)))
        gradD <- simplify(lapply(seq_along(X), function(i) Reduce("+",lapply(seq_along(X), function(j) vJ[i,j,j]))))
        f(X,Time) + gradD
    }
    list(D=D,mu=mu,J=J)
}

#' @export
buildPDESystem <- function(mu,D,rho,S, grd, phi0, usePointer = FALSE){
    ## Check input

    ## Defaults for rho, S0, S1
    if(missing(rho))
        rho <- function(...) AD(1)
    if(missing(S))
        S <- function(...) AD(0)
    numCells <- ncol(grd@Gcentroids)

    S_atomic <- MakeTape(function(p){
        phi0i <- p[1]
        t0 <- p[2]
        i <- p[3]
        X <- RTMB::DataEval(function(i) grd@Gcentroids[,i], i)
        S(X = X, Time = t0, Cell = i, Phi = phi0i)
    }, c(.1,0,1))
    S_atomic$simplify()
    S_atomic <- S_atomic$atomic()
    
    buildEquations <- function(t0,t1, phi0){
        dt <- t1 - t0
        rho_P <- simplify(lapply(seq_len(numCells), function(i) rho(X = grd@Gcentroids[,i], Time = t, Cell = i)))
        Sstar <- simplify(lapply(seq_len(numCells), function(i)  S_atomic(c(phi0[i],t0,i)) ))
        dS <- simplify(lapply(seq_len(numCells), function(i) S_atomic$jacobian(c(phi0[i],t0,i))[1,1]))
        S0 <- Sstar - dS * phi0
        S1 <- dS
        S0 <- AD(rep(0, numCells))
        S1 <- AD(rep(0,numCells))
        innerProd <- function(x, y) (x %*% y)[1,1]
        norm <- function(x) sqrt(innerProd(x,x))
        makeOneCell <- function(i){                    
            ga <- grd@Gareas[i]                
            T1 <- rho_P[i] - S1[i] * dt ## Diagonal
            makeOneNeighbour <- function(j){
                nj <- grd@neighbours[[i]][j]
                snv <- grd@Snvec[[i]][,j]
                ## Advection
                muval <- mu(grd@Scentroids[[i]][,j], t0)
                mf <- innerProd(muval,snv) * grd@Sarea[[i]][j]
                ## Diffusion
                Dval <- D(grd@Scentroids[[i]][,j], t0)
                dcf <- grd@Gcentroids[,nj] - grd@Gcentroids[,i]
                Sf <- Dval %*% snv
                SfNorm <- norm(Sf)
                dcfNorm <- norm(dcf)
                ## Convection
                Ff <- mf
                Df <- SfNorm / dcfNorm
                Pf <- mf / SfNorm * dcfNorm
                af <- pde_scheme(Pf)
                aA <- Df - (1 - af) * Ff
                ## rbind(c(i,nj,-aA * dt / ga),
                ##       c(i,i,(aA + Ff) * dt / ga))
                c(-aA * dt / ga, (aA + Ff) * dt / ga)
            }
            do.call("c",c(list(T1),lapply(seq_along(grd@neighbours[[i]]), makeOneNeighbour)))
        }
        X <- do.call("c",lapply(seq_along(grd@neighbours), makeOneCell))
        I <- rep(seq_along(grd@Sarea), times = 2 * sapply(grd@Sarea,length) + 1)
        J <- do.call("c",lapply(seq_along(grd@Sarea), function(i) c(i,sapply(grd@neighbours[[i]],function(j) c(j,i)))))
        X2 <- safe_aggregate(X,list(I=I,J=J), sum)    
        p <- c(0,cumsum(sapply(split(X2$J,cumsum(c(0,diff(X2$J)))),length)))
        dp <- diff(p)
        pp <- rep(seq_along(dp),dp)
        if(RTMB:::ad_context()){
            A <- new("adsparse",
                     x = X2$x,
                     i = as.integer(X2$I)-1L,
                     p = as.integer(p),
                     Dim = c(length(grd@Sarea),length(grd@Sarea)))
        }else{
            A <- new("dgCMatrix",
                     x=X2$x,
                     i=as.integer(X2$I)-1L,
                     p=as.integer(p),
                     Dim = c(length(grd@Sarea),length(grd@Sarea))
                     )
        }
        b <- S0 * dt + rho_P * phi0
        return(list(A=A,b=b))
    }
    Eq0 <- buildEquations(0,0.1,phi0)
    if(usePointer){
        sparse_solver_ptr <- RTMBconvenience:::sparse_solve_ptr_NotAD(Eq0$A)
        sparse_solver_update <- function(A) RTMBconvenience:::sparse_solve_ptr_update_NotAD(sparse_solver_ptr,A)
        sparse_solver_eval <- function(ph) as.vector(RTMBconvenience:::sparse_solve_ptr_eval_NotAD(sparse_solver_ptr, matrix(ph,ncol=1)))
        step <- function(t0, t1, phi0){
            Eq <- buildEquations(t0,t1,phi0)
            sparse_solver_update(Eq$A)
            sparse_solver_eval(Eq$b)
        }
    }else{ ## This should be used with AD in RTMB
        sparse_solver_ptr <- NULL
        sparse_solver_update <- NULL
        sparse_solver_eval <- NULL
        step <- function(t0, t1, phi0){
            Eq <- buildEquations(t0,t1,phi0)
            RTMB::solve(Eq$A,RTMB::matrix(Eq$b,ncol=1))[,1]
        }   
    }
    solve <- function(Times, phi0){
        res <- RTMB::matrix(0,length(phi0), length(Times))
        res[,1] = phi0
        for(i in tail(seq_along(Times),-1)){
            res[,i] <- step(Times[i-1],Times[i],res[,i-1])
        }
        res
    }  
    r <- list(buildEquations = buildEquations,
              sparse_solver_ptr = sparse_solver_ptr,
              sparse_solver_update = sparse_solver_update,
              sparse_solver_eval = sparse_solver_eval,
              step = step,
              solve = solve)
    class(r) <- "PDE_Solver"
    r
}

## buildA <- function(t, dt, rho, mu, D, S0, S1, grd){
##     innerProd <- function(x, y) (x %*% y)[1,1]
##     norm <- function(x) sqrt(innerProd(x,x))
##     makeOneCell <- function(i){
##         T1 <- 1 ## Diagonal
##         makeOneNeighbour <- function(j){
##             nj <- grd@neighbours[[i]][j]
##             snv <- grd@Snvec[[i]][,j]
##             ga <- grd@Gareas[i]
##             ## Advection
##             muval <- mu(grd@Scentroids[[i]][,j], t)
##             mf <- innerProd(muval,snv) * grd@Sarea[[i]][j]
##             ## Diffusion
##             Dval <- D(grd@Scentroids[[i]][,j], t)
##             dcf <- grd@Gcentroids[,nj] - grd@Gcentroids[,i]
##             Sf <- Dval %*% snv
##             SfNorm <- norm(Sf)
##             dcfNorm <- norm(dcf)
##             ## Convection
##             Ff <- mf
##             Df <- SfNorm / dcfNorm
##             Pf <- mf / SfNorm * dcfNorm
##             af <- ps(Pf)
##             aA <- Df - (1 - af) * Ff
##             ## rbind(c(i,nj,-aA * dt / ga),
##             ##       c(i,i,(aA + Ff) * dt / ga))
##             c(-aA * dt / ga, (aA + Ff) * dt / ga)
##         }
##         do.call("rbind",c(list(T1),lapply(seq_along(grd@neighbours[[i]]), makeOneNeighbour)))
##     }
##     X <- do.call("c",lapply(seq_along(grd@neighbours), makeOneCell))
##     I <- rep(seq_along(grd@Sarea), times = 2 * sapply(grd@Sarea,length) + 1)
##     J <- do.call("c",lapply(seq_along(grd@Sarea), function(i) c(i,sapply(grd@neighbours[[i]],function(j) c(j,i)))))
##     X2 <- safe_aggregate(X,list(I=I,J=J), sum)    
##     p <- c(0,cumsum(sapply(split(X2$J,cumsum(c(0,diff(X2$J)))),length)))
##     dp <- diff(p)
##     pp <- rep(seq_along(dp),dp)
##     if(RTMB:::ad_context()){
##         A <- new("adsparse",
##                  x = X2$x,
##                  i = as.integer(X2$I)-1L,
##                  p = as.integer(p),
##                  Dim = c(length(grd@Sarea),length(grd@Sarea)))
##     }else{
##         A <- new("dgCMatrix",
##                  x=X2$x,
##                  i=as.integer(X2$I)-1L,
##                  p=as.integer(p),
##                  Dim = c(length(grd@Sarea),length(grd@Sarea))
##                  )
##     }
##     return(A)
## }


## buildG <- function(t, dt, mu, D, grd){
##     innerProd <- function(x, y) (x %*% y)[1,1]
##     norm <- function(x) sqrt(innerProd(x,x))
##     makeOneCell <- function(i){
##         ##T1 <- c(i,i,0)
##         T1 <- 0
##         makeOneNeighbour <- function(j){
##             nj <- grd@neighbours[[i]][j]
##             snv <- grd@Snvec[[i]][,j]
##             ga <- grd@Gareas[i]
##             ## Advection
##             muval <- mu(grd@Scentroids[[i]][,j], t)
##             mf <- innerProd(muval,snv) * grd@Sarea[[i]][j]
##             ## Diffusion
##             Dval <- D(grd@Scentroids[[i]][,j], t)
##             dcf <- grd@Gcentroids[,nj] - grd@Gcentroids[,i]
##             Sf <- Dval %*% snv
##             SfNorm <- norm(Sf)
##             dcfNorm <- norm(dcf)
##             ## Convection
##             Ff <- mf
##             Df <- SfNorm / dcfNorm
##             Pf <- mf / SfNorm * dcfNorm
##             af <- ps(Pf)
##             aA <- Df - (1 - af) * Ff
##             ## rbind(c(i,nj,aA / ga),
##             ##       c(i,i,-(aA + Ff) / ga))
##             c(aA / ga,-(aA + Ff) / ga)
##         }
##         do.call("c",c(list(T1),lapply(seq_along(grd@neighbours[[i]]), makeOneNeighbour)))
##     }
##     X <- do.call("c",lapply(seq_along(grd@neighbours), makeOneCell))
##     I <- rep(seq_along(grd@Sarea), times = 2 * sapply(grd@Sarea,length) + 1)
##     J <- do.call("c",lapply(seq_along(grd@Sarea), function(i) c(i,sapply(grd@neighbours[[i]],function(j) c(j,i)))))
##     X2 <- safe_aggregate(X,list(I=I,J=J), sum)    
##     p <- c(0,cumsum(sapply(split(X2$J,cumsum(c(0,diff(X2$J)))),length)))
##     dp <- diff(p)
##     pp <- rep(seq_along(dp),dp)
##     if(RTMB:::ad_context()){
##         G <- new("adsparse",
##                  x = X2$x,
##                  i = as.integer(X2$I)-1L,
##                  p = as.integer(p),
##                  Dim = c(length(grd@Sarea),length(grd@Sarea)))
##     }else{
##         G <- new("dgCMatrix",
##                  x=X2$x,
##                  i=as.integer(X2$I)-1L,
##                  p=as.integer(p),
##                  Dim = c(length(grd@Sarea),length(grd@Sarea))
##                  )
##     }
##     return(G)
## }

## SolveSDE <- function(grd, times, F, G, x0, ...){
##     if(is.numeric(grd)){
##         ## Vector of 1D grid cell endpoints
##         grd <- Grid1D(grd)
##     }else if(is(grd,"sfc")){
##         ## sf map
##         ndim <- st_dimension(grd)[1]
##         grd <- Grid2D(grd,...)
##     }
##     #if(!is(grd, "Grid"
## }
