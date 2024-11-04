#' @export
UKF <- function(x0, P0, f, h, Q, R, Y,W0 = 0.5){
    Nx <- length(x0)    
    W <- c(W0,rep((1-W0)/(2*Nx),2*Nx))
    Xk <- matrix(0,Nx,nrow(Y)+1)
    Xk[,1] <- x0
    Pk <- array(0,dim=c(Nx,Nx,nrow(Y)+1))
    Pk[,,1] <- P0
    nll <- 0
    ## Loop over time
    for(i in 1:nrow(Y)){
        ## Make sigma points
        cP <- chol(((Nx / (1-W0)) * Pk[,,i]))
        S <- do.call(safe_cbind,list(Xk[,i,drop=FALSE], Xk[,i] +  cP, Xk[,i] - cP))
        ## Predict state
        Xp <- safe_apply(S, 2, f)
        Xkp1 <- Reduce("+",lapply(seq_along(W), function(j) Xp[[j]] * W[j]))
        Pkp1 <- Reduce("+",lapply(seq_along(W), function(j) ((Xp[[j]] - Xkp1) %*% t((Xp[[j]] - Xkp1))) * W[j])) + Q
        ## Predict obs
        Yp <- lapply(Xp, h)
        Zk <- Reduce("+",lapply(seq_along(W), function(j) Yp[[j]] * W[j]))
        Pzz <- Reduce("+",lapply(seq_along(W), function(j) ((Yp[[j]] - Zk) %*% t((Yp[[j]] - Zk))) * W[j])) + R
        Pxz <- Reduce("+",lapply(seq_along(W), function(j) ((Xp[[j]] - Xkp1) %*% t((Yp[[j]] - Zk))) * W[j]))
        Pzzi <- solve(Pzz)
        K <- Pxz %*% Pzzi
        ## Update
        Xk[,i+1] <- Xkp1 + K %*% (Y[i,] - Zk)
        Pk[,,i+1] <- Pkp1 - K %*% Pzz %*% t(K)
        ## log-Likelihood
        nll <- nll + 0.5 * (log(2 * pi * det(Pzz)) + ((Y[i,] - Zk) %*% Pzzi %*% (Y[i,] - Zk))[1,1])
    }
    list(nll = nll,
         states = t(Xk[,-1,drop=FALSE]),
         var = Pk[,,-1,drop=FALSE])
}
