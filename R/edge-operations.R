#' @title Edge matrix operations
#' @description Edge matrix operations needed for ips algorithms
#' @name emat-operations
#' @note These functions may well be removed from the package in future relases
#' @details An emat with p edges is represented by a 2 x p matrix.
#'
#' @param emat,emat1,emat2 Edge matrix (a 2 x p matrix)
#' @param dim A vector with dimensions
#' @param index A vector of integers
#' @param prob Probability of any edge being present.
#' @param type Output type.
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @examples
#' emat1 <- emat_saturated_model(3:4)
#' emat2 <- emat_saturated_model(1:4)
#' emat_complement(emat1, emat2)
#' emat3 <- emat_saturated_model(2:4)
#' emat_compare(emat1, emat3)

#' @rdname emat-operations
#' @export
emat_compare <- function(emat1, emat2){
    ## emat1, emat2 : 2 x p matrices
    ## common: Edges in emat1 AND emat2
    ## emat1only: Edges in emat1 that are not in emat2
    ## emat2only: Edges in emat2 that are not in emat1    
    emat1 <- colmat2list(emat1)
    emat2 <- colmat2list(emat2)

    common=emat1[(sapply(emat1, is_inset, emat2))]
    emat1only=emat1[!sapply(emat1, is_inset, emat2)]
    emat2only=emat2[!sapply(emat2, is_inset, emat1)]

    common=do.call(cbind, common)
    emat1only=do.call(cbind, emat1only)
    emat2only=do.call(cbind, emat2only)
    out  <- list(common=common, emat1only=emat1only, emat2only=emat2only)
    out
}

#' @rdname emat-operations
#' @export
emat_complement <- function(emat1, emat2){
    ## emat1, emat2 : 2 x p matrices
    ## Find complement of emat1 in emat2
    gl1 <- colmat2list(emat1)
    gl2 <- colmat2list(emat2)
    emat2[, !sapply(gl2, is_inset, gl1)]
}

#' @rdname emat-operations
#' @export
emat_sort <- function(emat1){
    id <- emat1[1,] > emat1[2,]
    emat1[1:2, id] <- emat1[2:1, id]
    o <- order(emat1[1,])
    emat1[, o]
}

#' @rdname emat-operations
#' @export
order_rows <- function(emat){
    b <- emat[1,] > emat[2,]
    emat[1:2, b] <- emat[2:1, b]
    emat
}

## #################################################################


#' @rdname emat-operations
#' @export
emat_saturated_model <- function(index){
    t(names2pairs(index, result="matrix"))
}

#' @rdname emat-operations
#' @export
model_saturated <- function(index, type="emat"){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_saturated_model(index)
    switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )        
}




## #################################################################

#' @rdname emat-operations
#' @export
emat_random_tree <- function(index, prob=0){

    emat_random_tree0 <- function(index){
        ## Based on https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence
        N <- max(index)
        emat <- matrix(0, nrow=2, ncol=N-1)
        degree <- rep(1, N)
        S <- sample(N, N-2, replace=TRUE)
        for (i in S){
            degree[i] <- degree[i] + 1
        }
        
        k <- 1
        for (i in S){
            for (j in seq_len(N)){
                ## cat(sprintf("i %d j %d\n", i, j))
                if (degree[j] == 1){
                    edge <- c(j, i)
                    ## cat("edge: \n"); print(edge)
                    emat[, k] <- edge
                    k <- k + 1
                    degree[i] <- degree[i] - 1
                    degree[j] <- degree[j] - 1
                    break
                }
            }    
        }
        emat[, k] <- which(degree == 1)    
        emat
    }

    if ((prob < 0) || (prob > 1)) stop("prob must be in [0, 1]")
    if (prob < 1e-6)
        return(emat_random_tree0(index))

    emat1 <- emat_random_tree(index)
    emat2 <- emat_random_model(index, prob)
    ee   <- cbind(order_rows(emat1), order_rows(emat2))
    emat <- t.default(unique(t.default(ee)))
    emat
}

#' @rdname emat-operations
#' @export
model_random_tree <- function(index, prob=0, type="emat"){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_random_tree(index, prob)
    switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )        
}


## #################################################################

#' @rdname emat-operations
#' @export
emat_rectangular_grid <- function(dim){

    if (length(dim)==1)
        dim <- rep(dim, 2)

    nr.nc <- dim
    nr <- nr.nc[1]
    nc <- nr.nc[2]
    
    ## first row of edges
    rr <- rbind(
        1:(nc-1),
        2:nc
    )
    
    ## first column of edges
    bb <- 0:(nr-1)*nc + 1
    cc <- rbind(
        bb[1:(nr-1)],
        bb[2:nr]
    )
        
    ## all rows of edges
    roe <- lapply((0:(nr-1))*nc, function(r) rr + r)  |>  do.call(cbind, args=_)
    ## all columns of edges
    coe <- lapply((0:(nc-1)), function(c) cc + c)  |>  do.call(cbind, args=_)
    
    out <- cbind(roe, coe)
    out
}

#' @rdname emat-operations
#' @export
model_rectangular_grid <- function(dim, type="emat"){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_rectangular_grid(dim)
    switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, prod(dim))},
           "amat" ={as_emat2amat(em, prod(dim))}
           )        

}

## #################################################################

#' @rdname emat-operations
#' @export
emat_line_model <- function(index){
  unname(rbind(index[-length(index)], index[-1]))
}

#' @rdname emat-operations
#' @export
model_line <- function(index, type="emat"){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_line_model(index)
    switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )        
}


## #################################################################

#' @rdname emat-operations
#' @export
emat_star_model <- function(index){
    unname(rbind(index[1], index[-1]))
}

#' @rdname emat-operations
#' @export
model_star <- function(index, type="emat"){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_star_model(index)
    switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )        
}


## #################################################################

#' @rdname emat-operations
#' @export
emat_loop_model <- function(index, prob=0){

    emat_loop_model0 <- function(index){
        if (length(index) == 2) matrix(index)
        else unname(rbind(index, c(index[-1], index[1])))
    }
    
    if ((prob < 0) || (prob > 1)) stop("prob must be in [0, 1]")
    if (prob < 1e-6)
        return(emat_loop_model0(index))
    
    emat1 <- emat_loop_model0(index)
    emat2 <- emat_random_model(index, prob)
    ee   <- cbind(order_rows(emat1), order_rows(emat2))
    emat <- t.default(unique(t.default(ee)))
    emat
}

#' @rdname emat-operations
#' @export
model_loop <- function(index, prob=0, type="emat"){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_loop_model(index, prob)
    switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )
}


## #################################################################

#' @rdname emat-operations
#' @export
emat_random_model <- function(index, prob=0.1){
    N <- max(index)
    g     <- igraph::erdos.renyi.game(N, prob)
    emat  <- t(igraph::get.edgelist(g))
    emat    
}

#' @rdname emat-operations
#' @export
model_random <- function(index, prob=0.1, type="emat"){
    type <- match.arg(type, c("emat", "glist", "cliq", "amat"))
    em <- emat_random_model(index, prob)
    switch(type,
           "emat" ={em},
           "glist"={as_emat2glist(em)},
           "cliq" ={as_emat2cq(em, max(index))},
           "amat" ={as_emat2amat(em, max(index))}
           )    
}




## #' @rdname emat-operations
## #' @export
## emat_grid_model <- function(nr.nc){
##     if (length(nr.nc==1))
##         nr.nc <- rep(nr.nc, 2)
    
##     nr <- nr.nc[1]
##     nc <- nr.nc[2]
    
##     ## first row of edges
##     rr <- rbind(
##         1:(nc-1),
##         2:nc
##     )
    
##     ## first column of edges
##     bb <- 0:(nr-1)*nc + 1
##     cc <- rbind(
##         bb[1:(nr-1)],
##         bb[2:nr]
##     )

##     . <- NULL
##     ## all rows of edges
##     roe <- lapply((0:(nr-1))*nc, function(r) rr + r) %>% do.call(cbind, .)
##     ## all columns of edges
##     coe <- lapply((0:(nc-1)), function(c) cc + c) %>% do.call(cbind, .)
    
##     grid <- cbind(roe, coe)
##     grid
## }

## ' @rdname emat-operations
## ' @export
## emat_cq_random_model <- function(index, prob=0.1){
    ## g     <- igraph::erdos.renyi.game(max(index), prob)
    ## emat  <- t(igraph::get.edgelist(g))
    ## cq <- as_emat2cq(emat, max(index))
    ## list(emat=emat, cq=cq)
## }
