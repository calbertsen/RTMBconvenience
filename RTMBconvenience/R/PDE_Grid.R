setOldClass("sfc")


validGridObject <- function(object) TRUE

#' @export
setClass("Grid",
         slots = c(Gcentroids = "matrix",
                   Gareas = "numeric",
                   neighbours = "list",
                   Scentroids = "list",
                   Sarea = "list",
                   Snvec = "list",
                   gI = "list",
                   vI = "list",
                   time = "numeric",
                   Gindex = "matrix",
                   geometryIndex = "integer",
                   geometry = "list"),
         sealed = TRUE,
         validity = validGridObject)

get_centroid_in <- function(x, forceCentroidInCell = TRUE){
    cent <- st_centroid(x)
    if(forceCentroidInCell){
        centInCell <- diag(st_intersects(cent,x, FALSE))    
        cent[!centInCell] <- st_point_on_surface(x[!centInCell])
    }
    cent
}

get_coordvec <- function(x) as.numeric(st_coordinates(x))

#' @export
GridCombine <- function(g1,g2){
    dim1 <- nrow(g1@Gcentroids)
    dim2 <- nrow(g2@Gcentroids)
    ix1 <- seq_len(ncol(g1@Gcentroids))
    ix2 <- seq_len(ncol(g2@Gcentroids))
    ixN <- expand.grid(ix1,ix2)

    Gcentroids <- rbind(g1@Gcentroids[,ixN[,1]],
                        g2@Gcentroids[,ixN[,2]])
    Gareas <- g1@Gareas[ixN[,1]] * g2@Gareas[ixN[,2]]

    neigh <- lapply(seq_len(nrow(ixN)), function(i){
        ## if(diagonalNeighbours){
        ##     setdiff(intersect(which(ixN[,1] %in% (g1@neighbours[[ixN[i,1]]]+1)),
        ##                       which(ixN[,2] %in% (g2@neighbours[[ixN[i,2]]]+1))),
        ##             i)
        ## }else{
        setdiff(sort(c(intersect(which(ixN[,1] %in% (g1@neighbours[[ixN[i,1]]]+1)),
                                 which(ixN[,2] == ixN[i,2])),
                       intersect(which(ixN[,2] %in% (g2@neighbours[[ixN[i,2]]]+1)),
                                 which(ixN[,1] == ixN[i,1])))),
                i)    
        ## }
    })

    Sarea <- sapply(seq_along(neigh), function(i){
        sapply(seq_along(neigh[[i]]), function(j){
            ni <- as.numeric(ixN[i,])
            nn <- as.numeric(ixN[neigh[[i]][j],])
            k1 <- which(g1@neighbours[[ni[1]]]+1 == nn[1])
            k2 <- which(g2@neighbours[[ni[2]]]+1 == nn[2])
            if(length(k1) == 0){
                A1 <- g1@Gareas[ni[1]]
                A2 <- g2@Sarea[[ni[2]]][k2]
            }else if(length(k2) == 0){
                A1 <- g1@Sarea[[ni[1]]][k1]
                A2 <- g2@Gareas[ni[2]]
            }else{
                stop("Diagonal neighbours not implemented")
            }
            A1 * A2
        })
    })
    Scentroids <- sapply(seq_along(neigh), function(i){
        sapply(seq_along(neigh[[i]]), function(j){
            ni <- as.numeric(ixN[i,])
            nn <- as.numeric(ixN[neigh[[i]][j],])
            k1 <- which(g1@neighbours[[ni[1]]]+1 == nn[1])
            k2 <- which(g2@neighbours[[ni[2]]]+1 == nn[2])
            if(length(k1) == 0){
                ## Neighbour in g2
                return(c(g1@Gcentroids[,ni[1]], g2@Scentroids[[ni[2]]][,k2]))
            }else if(length(k2) == 0){
                ## Neighbour in g1
                return(c(g1@Scentroids[[ni[1]]][,k1],g2@Gcentroids[,ni[2]]))
            }else{
                stop("Diagonal neighbours not implemented")
            }              
        })
    })
    Snvec <- sapply(seq_along(neigh), function(i){
        sapply(seq_along(neigh[[i]]), function(j){
            ni <- as.numeric(ixN[i,])
            nn <- as.numeric(ixN[neigh[[i]][j],])
            k1 <- which(g1@neighbours[[ni[1]]]+1 == nn[1])
            k2 <- which(g2@neighbours[[ni[2]]]+1 == nn[2])
            if(length(k1) == 0){
                ## Neighbour in g2
                return(c(rep(0,dim1), g2@Snvec[[ni[2]]][,k2]))
            }else if(length(k2) == 0){
                ## Neighbour in g1
                return(c(g1@Snvec[[ni[1]]][,k1],rep(0,dim2)))
            }else{
                stop("Diagonal neighbours not implemented")
            }              
        })
    })

    ## gI <- lapply(seq_len(ncol(Gcentroids)), function(i){
    ##     sapply(seq_along(neigh[[i]]), function(j){           
    ##         ## norm(Scentroids[[i]][,j] - Gcentroids[,neigh[[i]][j]],"2") / norm(Gcentroids[,i] - Gcentroids[,neigh[[i]][j]],"2")
    ##         ef <- Scentroids[[i]][,j] / norm(Scentroids[[i]][,j],"2")
    ##         dCf <- Scentroids[[i]][,j] - Gcentroids[,i]
    ##         dfF <- Gcentroids[,neigh[[i]][j]] - Scentroids[[i]][,j]
    ##         a <- sum(dCf * ef)
    ##         b <- sum(dfF * ef)
    ##         a / (a+b)
    ##     })
    ## })
    gI <- lapply(seq_len(ncol(Gcentroids)), function(i){
        sapply(seq_along(neigh[[i]]), function(j){           
            ## norm(Scentroids[[i]][,j] - Gcentroids[,neigh[[i]][j]],"2") / norm(Gcentroids[,i] - Gcentroids[,neigh[[i]][j]],"2")
            dCF <- Gcentroids[,neigh[[i]][j]] - Gcentroids[,i]
            ef <- dCF / norm(dCF,"2")
            dCf <- Scentroids[[i]][,j] - Gcentroids[,i]
            dfF <- Gcentroids[,neigh[[i]][j]] - Scentroids[[i]][,j]
            a <- norm(dCf,"2")
            b <- norm(dfF,"2")
            ## a <- sum(dCf * ef)
            ## b <- sum(dfF * ef)
            ## if(a+b == 0)
            ##     return(0.5)
            a / (a+b)
        })
    })
    vI <- lapply(seq_len(ncol(Gcentroids)), function(i){
        sapply(seq_along(neigh[[i]]), function(j){
            Gareas[i] / (Gareas[i] + Gareas[neigh[[i]][j]])            
        })
    })

    new("Grid",
        Gcentroids = Gcentroids,
        Gareas = as.numeric(Gareas),
        neighbours = lapply(neigh,function(x)x),
        Scentroids = Scentroids,
        Sarea = Sarea,
        Snvec = Snvec,
        gI = gI,
        vI = vI,
        time = 0,
        Gindex = rbind(g1@Gindex[,ixN[,1]],
                       g2@Gindex[,ixN[,2]]),
        geometryIndex = c(g1@geometryIndex,max(g1@geometryIndex) + g2@geometryIndex),
        geometry = c(g1@geometry,g2@geometry))
}

#' @export
Grid1D <- function(map, cellSize, time = NA, clip = FALSE){
    if(is.numeric(map)){
        ## If range
        if(length(map) == 1)
            stop("wrong class, map must be numeric with length >= 2 or sfc")
        if(length(map) == 2)
            map <- unique(c(map[1],seq(map[1], map[2], by = cellSize),map[2]))
        ## If seq, don't do anything
        ## Convert to sfc
        grid <- st_sfc(sapply(tail(seq_along(map),-1), function(i){
            xx <- cbind(map[c(i-1,i,i,i-1,i-1)], c(-0.5,0.5)[c(1,1,2,2,1)])
            st_polygon(list(xx))
        }, simplify = FALSE))
    }else if(is(map,"sfc")){
        ## If sfc
        if(missing(cellSize)){
            grid <- map
        }else{
            gridraw <- st_make_grid(map, cellsize = c(cellSize,diff(st_bbox(map)[c(2, 4)])), square = TRUE)
            fullCellArea <- as.numeric(st_area(gridraw[1]))
            if(clip){
                grid <- st_intersection(gridraw, map)
            }else{
                grid <- gridraw[sapply(st_intersects(gridraw,map),length)>0]
            }
        }
    }else{
        stop("wrong class, x must be numeric with length >= 2 or sfc")
    }

    grid <- grid[as.numeric(st_area(grid)) > 0]
    Gcentroids <- get_centroid_in(grid)
    Gareas <- as.numeric(st_area(grid))

    neigh <- st_intersects(grid,grid, sparse = TRUE)
    for(i in 1:length(neigh)) neigh[[i]] <- neigh[[i]][neigh[[i]] != i]
    Sarea <- sapply(1:length(neigh), function(i){
        sapply(neigh[[i]], function(j)as.numeric(st_length(st_intersection(grid[i],grid[j]))))
    })
    for(i in 1:length(neigh)) neigh[[i]] <- neigh[[i]][Sarea[[i]] > 0]
    for(i in 1:length(Sarea)) Sarea[[i]] <- Sarea[[i]][Sarea[[i]] > 0]

    Scentroids <- sapply(1:length(neigh), function(i){
        sapply(neigh[[i]], function(j)st_centroid(st_intersection(grid[i],grid[j])))
    })
    Snvec <- sapply(1:length(neigh), function(i){
        lapply(neigh[[i]], function(j){
            lsec <- st_intersection(grid[i],grid[j])
            l <- st_coordinates(lsec)[,c("X","Y")]
            l <- l[c(1,nrow(l)),]
            if(diff(l[,2])==0){
                z1 <- 0
                z2 <- 1
            }else{
                z1 <- 1 # ifelse(diff(l[,1])==0,0,1)#-diff(l[,2])/diff(l[,1]))
                z2 <- -diff(l[,1])/diff(l[,2])
            }
            ss <- sign(st_coordinates(Gcentroids[j] - Gcentroids[i]))
            ns <- 1
            if((ss[1] != 0 && ss[1] != sign(z1)) || (ss[1] == 0 && ss[2] != sign(z2))){
                ns <- -1
            }
            st_point(ns * c(z1,z2) / sqrt(z1^2 + z2^2))
        })
    })

    ## gI <- lapply(seq_along(grid), function(i){
    ##     sapply(seq_along(neigh[[i]]), function(j){
    ##         norm(Scentroids[[i]][[j]] - Gcentroids[[neigh[[i]][j]]],"2") / norm(Gcentroids[[i]] - Gcentroids[[neigh[[i]][j]]],"2")
    ##     })
    ## })
    gI <- lapply(seq_along(Gcentroids), function(i){
        sapply(seq_along(neigh[[i]]), function(j){           
            ## norm(Scentroids[[i]][,j] - Gcentroids[,neigh[[i]][j]],"2") / norm(Gcentroids[,i] - Gcentroids[,neigh[[i]][j]],"2")
            dCF <- get_coordvec(Gcentroids[neigh[[i]][j]]) - get_coordvec(Gcentroids[[i]])
            ef <- dCF / norm(dCF,"2")
            dCf <- get_coordvec(Scentroids[[i]][[j]]) - get_coordvec(Gcentroids[[i]])
            dfF <- get_coordvec(Gcentroids[neigh[[i]][j]]) - get_coordvec(Scentroids[[i]][[j]])
            a <- norm(dCf,"2")
            b <- norm(dfF,"2")
            ## a <- sum(dCf * ef)
            ## b <- sum(dfF * ef)
            ## if(a+b == 0)
            ##     return(0.5)
            a / (a+b)
        })
    })
    vI <- lapply(seq_along(grid), function(i){
        sapply(seq_along(neigh[[i]]), function(j){
            Gareas[i] / (Gareas[i] + Gareas[neigh[[i]][j]])
        })
    })

    tim <- time
    if(is.na(time)){
      #warning("Time not supplied. Using a single time point,")
        tim <- 0
    }

    new("Grid",
        Gcentroids = t(st_coordinates(Gcentroids))[1,,drop = FALSE],
        Gareas = as.numeric(Gareas),
        neighbours = lapply(neigh,function(x)x),
        Scentroids = lapply(Scentroids,function(x)sapply(x,st_coordinates)[1,,drop = FALSE]),
        Sarea = Sarea,
        Snvec = lapply(Snvec,function(x)sapply(x,st_coordinates)[1,,drop = FALSE]),
        gI = gI,
        vI = vI,
        time = tim,
        Gindex = matrix(seq_len(length(Gareas)),nrow=1),
        geometryIndex = 1L,
        geometry = list(grid))

}

#' @export
Grid2D <- function(map, cellSize, square = FALSE, time = NA, clip = FALSE, minArea = 0.1, centroidInCell = TRUE){
    if(!is(map,"sfc"))
        stop("wrong class, map must be sfc")

    if(missing(cellSize)){
        grid <- map
    }else{
        gridraw <- st_make_grid(map, cellsize = cellSize, square = square)
        fullCellArea <- as.numeric(st_area(gridraw[1]))
        if(clip){
            grid <- st_intersection(gridraw, map)
        }else{
            grid <- gridraw[sapply(st_intersects(st_centroid(gridraw),map),length)>0]
        }
    }
    stag <- as.numeric(st_area(grid))    
    grid <- grid[stag / max(stag) > minArea]
    if(length(grid)==0) stop("minArea too large. No cells left.")
    
    ## Check for neighbour (linestring intersection)
    hasN <- sapply(seq_along(grid),function(i)max(as.numeric(st_length(st_intersection(grid[i],grid[-i])))) > 0)
    grid <- grid[hasN]
    if(length(grid)==0) stop("minArea too large. No cells left.")
    ## Remove neighbours with only points
    

    Gcentroids <- get_centroid_in(grid, centroidInCell)

    Gareas <- as.numeric(st_area(grid))

    neigh <- st_intersects(grid,grid, sparse = TRUE)
    for(i in 1:length(neigh)) neigh[[i]] <- neigh[[i]][neigh[[i]] != i]
    Sarea <- sapply(1:length(neigh), function(i){
        sapply(neigh[[i]], function(j)as.numeric(st_length(st_intersection(grid[i],grid[j]))))
    })
    for(i in 1:length(neigh)) neigh[[i]] <- neigh[[i]][Sarea[[i]] > 0]
    for(i in 1:length(Sarea)) Sarea[[i]] <- Sarea[[i]][Sarea[[i]] > 0]

    Scentroids <- sapply(1:length(neigh), function(i){
        sapply(neigh[[i]], function(j)st_centroid(st_intersection(grid[i],grid[j])))
    })
    Snvec <- sapply(1:length(neigh), function(i){
        lapply(neigh[[i]], function(j){
            lsec <- st_intersection(grid[i],grid[j])
            l <- st_coordinates(lsec)[,c("X","Y")]
            l <- l[c(1,nrow(l)),]
            if(diff(l[,2])==0){
                z1 <- 0
                z2 <- 1
            }else{
                z1 <- 1 # ifelse(diff(l[,1])==0,0,1)#-diff(l[,2])/diff(l[,1]))
                z2 <- -diff(l[,1])/diff(l[,2])
            }
            ss <- sign(st_coordinates(Gcentroids[j] - Gcentroids[i]))
            ns <- 1
            if((ss[1] != 0 && ss[1] != sign(z1)) || (ss[1] == 0 && ss[2] != sign(z2))){
                ns <- -1
            }
            st_point(ns * c(z1,z2) / sqrt(z1^2 + z2^2))
        })
    })

    ## gI <- lapply(seq_along(grid), function(i){
    ##     sapply(seq_along(neigh[[i]]), function(j){
    ##         norm(Scentroids[[i]][[j]] - Gcentroids[[neigh[[i]][j]]],"2") / norm(Gcentroids[[i]] - Gcentroids[[neigh[[i]][j]]],"2")
    ##     })
    ## })
   ## gI <- lapply(seq_along(grid), function(i){
   ##      sapply(seq_along(neigh[[i]]), function(j){           
   ##          ## norm(Scentroids[[i]][,j] - Gcentroids[,neigh[[i]][j]],"2") / norm(Gcentroids[,i] - Gcentroids[,neigh[[i]][j]],"2")
   ##          ef <- get_coordvec(Scentroids[[i]][[j]]) / norm(get_coordvec(Scentroids[[i]][[j]]),"2")
   ##          dCf <- get_coordvec(Scentroids[[i]][[j]]) - get_coordvec(Gcentroids[[i]])
   ##          dfF <- get_coordvec(Gcentroids[neigh[[i]][j]]) - get_coordvec(Scentroids[[i]][[j]])
   ##          a <- sum(dCf * ef)
   ##          b <- sum(dfF * ef)
   ##          a / (a+b)
   ##      })
    ##  })
       gI <- lapply(seq_along(Gcentroids), function(i){
        sapply(seq_along(neigh[[i]]), function(j){           
            ## norm(Scentroids[[i]][,j] - Gcentroids[,neigh[[i]][j]],"2") / norm(Gcentroids[,i] - Gcentroids[,neigh[[i]][j]],"2")
            dCF <- get_coordvec(Gcentroids[neigh[[i]][j]]) - get_coordvec(Gcentroids[[i]])
            ef <- dCF / norm(dCF,"2")
            dCf <- get_coordvec(Scentroids[[i]][[j]]) - get_coordvec(Gcentroids[[i]])
            dfF <- get_coordvec(Gcentroids[neigh[[i]][j]]) - get_coordvec(Scentroids[[i]][[j]])
            a <- norm(dCf,"2")
            b <- norm(dfF,"2")
            ## a <- sum(dCf * ef)
            ## b <- sum(dfF * ef)
            ## if(a+b == 0)
            ##     return(0.5)
            a / (a+b)
        })
    })
    vI <- lapply(seq_along(grid), function(i){
        sapply(seq_along(neigh[[i]]), function(j){
            Gareas[i] / (Gareas[i] + Gareas[neigh[[i]][j]])
        })
    })

    tim <- time
    if(is.na(time)){
      #warning("Time not supplied. Using a single time point,")
        tim <- 0
    }

    new("Grid",
        Gcentroids = t(st_coordinates(Gcentroids)),
        Gareas = as.numeric(Gareas),
        neighbours = lapply(neigh,function(x)x),
        Scentroids = lapply(Scentroids,function(x)sapply(x,st_coordinates)),
        Sarea = Sarea,
        Snvec = lapply(Snvec,function(x)sapply(x,st_coordinates)),
        gI = gI,
        vI = vI,
        time = tim,
        Gindex = matrix(seq_len(length(Gareas)),nrow=2,ncol=length(Gareas),byrow=TRUE),
        geometryIndex = c(1L,1L),
        geometry = list(grid))
}

#' @export
setMethod("show", "Grid", function(object){
    show(object@geometry)
})


#' @export
setMethod("plot", "Grid", function(x, y,..., set_mfrow = TRUE){
    if(set_mfrow)
        par(mfrow=n2mfrow(length(x@geometry)))
    for(i in seq_along(x@geometry))
    plot(x@geometry[[i]], ...)
})
