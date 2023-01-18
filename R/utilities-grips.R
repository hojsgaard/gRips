
#' Genrate matrix of N(0, 1) variables
#'
#' @param n.obs Number of rows
#' @param nvar Number of columns
#' @param seed Seed for random number generator
#' 
#' @export
generate_n01 <- function(n.obs, nvar, seed=2022){
    set.seed(seed)
    d <-matrix(rnorm(n.obs * nvar), nrow=n.obs)
    d
}


## FIXME : Need a solve_fun 

solve_qr <- function(A, b=NULL){
    if (is.null(b))
        solve.qr(qr(A))
    else
        solve.qr(qr(A)) %*% b
}

solve_ginv <- function(A, b=NULL){
    if (is.null(b))
        MASS::ginv(A)
    else
        MASS::ginv(A) %*% b
}

solveM <- function(A, b=NULL){
    if (is.null(b))
        Matrix::solve(A)
    else
        Matrix::solve(A, b)
}

ips_logL <- function(S, K, nobs){
    nobs*(log(det(K)) - sum(S*K))/2
}


## UTILITIES

.get.time <- function(){
    Sys.time()
}

.get.diff.time <- function(t0, units="secs"){
  unit_ <- match.arg(units, c("secs", "millisecs"))
  d <- as.numeric(difftime(Sys.time(), t0, units="secs"))
  if (identical(units, "millisecs"))
      d <- d * 1000
  d
}


get_zero_edges <- function(emat, dim){
  if (is.null(emat)) emat <- matrix(NA, 2, 0)
  emat <- cbind(emat, emat[2:1,])
  amat <- matrix(0L, nrow=dim, ncol=dim)
  amat[t.default(emat)] <- 1L
  amat <- upper.tri(amat, diag=TRUE) + amat
  out <- which(amat==0L, arr.ind = TRUE)

  t.default(out)
}


#' @title Utes
#' @description Utes
#' @name utes
#' @param emat Edge matrix (2 x p)
#' @param Sigma,S Symmetrix positive definite matrices
#' @param eps Threshold
#' 

## # do backquote
## dobq <- function(fnlist){  ## To doBy
##    lapply(fnlist, function(g) bquote(.(g)()))
## }


#' @export
#' @rdname utes
does_fit <- function(Sigma, S, emat, eps=1e-4){
  v <- max_abs_diff_on_Emat(Sigma, S, emat)
  cat("max deviation : ", v, "\n")
  v < eps
}


#' Parameter equality
#'
#' @param p2a,p2b Parameter lists.
#' @param eps A small number.
#' 
#' @export
#' @rdname utilities
parmEq <- function(p2a, p2b, eps=1e-4){
  max(abs(p2a$K - p2b$K)) / max(abs(p2a$K)) < eps && max(abs(p2a$Sigma - p2b$Sigma)) / max(abs(p2a$Sigma)) < eps
}


#' Impose zeros in matrix entries which do not correspond to an edge.
#'
#' @param emat Edge matrix (2 x p matrix)
#' @param K Matrix; typically a concentration matrix.
#'
#' @export
impose_zero <- function(emat, K){
    emc <- as_emat_complement(emat, nrow(K))
    if (ncol(emc) == 0){
        return(K)  
    } 
    emc <- t.default(emc)
    K[rbind(emc, emc[2:1,])]  <- 0
    return(K)
}


#' @title Utilities
#'
#' @description Utilities
#'
#' @name utilities
#'
#' @param emat Edge matrix (2 x p)
#' @param amat Adjacency matrix
#' @param glist Generator list
#' @param elist Edge list (list of pairs)
#' @param nvar Number of variables
#' @param K Concentration matrix
#' @param x An object (to print)
#' @param d Number of columns in output.
#'

# Convert edges to cliques
#' @export
#' @rdname utilities
as_emat2cq <- function(emat, nvar=NULL){
    ## if (ncol(emat) == 0) return(as.list(1:nvar))
    if (ncol(emat) == 0) return(list())
    if (is.null(nvar)){
        nvar <- max(emat)
    }
    am <- matrix(0, nrow=nvar, ncol=nvar)
    
    am[t(emat)] <- 1  
    am <- am + t(am)

    ## print(am)
    ## am <<- am
    cq <- as(am, "graphNEL")  %>% as("dgCMatrix")  %>% gRbase::maxCliqueMAT() 
    cq <- lapply(cq$maxCliques, "as.numeric")
    cq
}


#' @export
#' @rdname utilities
as_emat_complement <- function(emat, nvar){ # Those edges NOT in emat
  am <- matrix(0, nrow=nvar, ncol=nvar)
  am[t(emat)] <- 1
  am[t(emat[2:1, ,drop=FALSE])] <- 1
  am <- (am == 0) * lower.tri(am)
  t.default(which(am==1, arr.ind=TRUE))
}


#' @export
#' @rdname utilities
as_emat2amat <-function(emat, d=-1){
    vn <- sort(unique(c(emat)))
    if (d == -1) d=length(vn)
    ##print(vn)
    M <- matrix(0, nrow=d, ncol=d)
    for (j in 1:ncol(emat)){
        e <- emat[,j]
        M[e[1], e[2]] <- M[e[2], e[1]] <- 1
    }
    M
}

#' @export
#' @rdname utilities
as_emat2elist <- function(emat){
    if (ncol(emat) > 0)
        split(t.default(emat), 1:(ncol(emat)))
    else
        list()
}

#' @export
#' @rdname utilities
as_elist2emat <- function(elist){

    idx <- sapply(elist, length) > 1
    elist <- elist[idx]

    if (length(elist) > 0){
        . <- NULL 
        ed <- lapply(elist, combn_prim, 2)  %>% do.call(cbind, .)  %>% t
        ## ed <- lapply(elist, combn_prim, 2)  %>% {function(zz) {do.call(cbind, zz)}}  %>% t
        b <- ed[, 1] > ed[, 2]
        ed[b,] <- ed[b, 2:1]
        ed <- ed[!duplicated.matrix(ed),,drop=FALSE]
        Emat <- t.default(ed)
    } else
        Emat <- matrix(NA, ncol=0, nrow=2)
    Emat
}

#' @export
#' @rdname utilities
as_glist2emat <- function(glist){
  if (!inherits(glist, "list")) stop("'glist' must be a list\n")
  if (length(glist) == 0)
    return(matrix(NA, nrow=2, ncol=0))
  glist <- glist[sapply(glist, function(x) length(x)>1)]
  if (length(glist) == 0)
    return(matrix(NA, nrow=2, ncol=0))

  g2 <- lapply(glist, function(x) names2pairs(x))
  g2 <- unlist(g2, recursive=FALSE)
  g2 <- do.call(rbind, g2)
  t.default(g2[!duplicated(g2),])
}

#' @export
#' @rdname utilities
as_emat2graph  <- function(emat){
    as(as_emat2amat(emat), "graphNEL")
}

#' @export
#' @rdname utilities
as_amat2edges <- function(amat, eps=1e-4){
    amat[lower.tri(amat, diag=TRUE)] <- 0
    e <- which(abs(amat) > eps, arr.ind = T)
    t.default(e)
}


#' @export
#' @rdname utilities
as_amat2emat <- function(amat, eps=1e-4){
    amat[lower.tri(amat, diag=TRUE)] <- 0
    e <- which(abs(amat) > eps, arr.ind = T)
    rownames(e) <- NULL
    colnames(e) <- NULL
    t.default(e)
}

#' @export
#' @rdname utilities
as_emat2glist <- function(emat){
    gRbase::colmat2list(emat)
}

#' @export
#' @rdname utilities
as_glist2out_edges <- function(glist){
    u <- gRbase::ug(glist)
    n <- gRbase::nonEdgeList(u)
    do.call(rbind, lapply(n, as.numeric))
}

#' @export
#' @rdname utilities
as_K2amat <- function(K, eps=1e-4){
    MM <- zapsmall(K)
    MM <- 1 * (abs(MM) > eps)
    diag(MM) <- 0
    MM
}

#' @export
#' @rdname utilities
as_K2edges <- as_amat2edges

#' @export
#' @rdname utilities
as_K2graph <- function(K){
    as(as_K2amat(zapsmall(K)), "graphNEL")
}

#' @export
#' @rdname utilities
as_sparse <- function(K){
    as(zapsmall(K), "dgCMatrix")}

.form2emat <- function(form){
    if (is.matrix(form) && nrow(form) == 2) return(form)  ## FIXME
    else
        if (length(form) > 0){
            em <- do.call(rbind, lapply(form, function(g)
                if (length(g) > 1) t.default(combn_prim(g, 2))))
            em <- t.default(unique(em))
        } else
            em  <- matrix(nrow=2, ncol=0)
    em
}



































### Internal .functions

.sol <- function(A){ ## FIXME : silly name
    .solve2x2 <- function(A){
        ##d <- A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1] ## FIXME: Handle d=0 case
        d <- A[1] * A[4] - A[2] * A[3] ## FIXME: Handle d=0 case
        #B <- A
        B <- c(0, 0, 0, 0)
        B[1] <- A[4]
        B[2] <- -A[2]
        B[3] <- -A[3]
        B[4] <- A[1]
        B <- B * (1 / d)
        dim(B) = c(2L, 2L)
        B
    }
    
    if (dim(A)[1] == 2) .solve2x2(A)
    else solve.default(A) ## FIXME: Could do chol2inv(chol(A))
}

.checkS <- function(S){
    if (!is.matrix(S)) 
        stop("Second argument is not a matrix!")
    
    if (dim(S)[1] != dim(S)[2]) 
        stop("Second argument is not a square matrix!")
}

.penalizeS <- function(S, lambda=0){
    S + diag(lambda, nrow(S))
}


## Just once
.make_rc <- function(glist, S){
    if (!is.matrix(glist))
        glist <- do.call(cbind, glist)
    vn    <- unique.default(sort(glist))
    amat  <- matrix(0L, nrow=length(vn), ncol=length(vn))
    amat[t.default(glist)] <- 1L
    amat <- amat + t.default(amat)
    dimnames(amat) <- dimnames(S)
    lt   <- lower.tri(amat)
    rrr  <- which.arr.index(lt * amat != 0)
    colnames(rrr) <- c("row", "col")
    rrr
}


.getDiag <- function(S){
    r <- nrow(S)
    S[1 + (r+1) * (0:(r-1))]
}

#' @export
#' @rdname utilities
#' @param object Model object.
#' @param ... Additional arguments; currently not used.
#' @param k Penalty parameter for calculating AIC; only k=2 gives genuine AIC.
logLik.gips_fit_class <- function(object, ...){

  nobs <- nobs(object)# unname(object$details["nobs"])
  if (is.na(nobs))
    stop("'nobs' not given; can not compute log L\n")
  
  trKS  <- unname(object$details$trKS)
  nparm <- unname(object$details$dim)
  
  ##str(list(nobs, trKS))
  out  <- c(nobs * (log(det(object$K)) - trKS) / 2)
  
  ## FIXME: Perhaps need number of parameters in specified model.
  attr(out, "nobs")  = nobs
  attr(out, "nparm") = nparm
  attr(out, "df")    = nparm
  
  class(out) <- "logLik"
  out    
}


#' @method AIC gips_fit_class
#' @export
#' @rdname utilities
AIC.gips_fit_class <- function(object, ..., k=2){
    ll <- logLik(object) 
    -2 * as.numeric(ll) + k * attr(ll, "nparm")
}

#' @method BIC gips_fit_class
#' @export
#' @rdname utilities
BIC.gips_fit_class <- function(object, ...){
    ll <- logLik(object) 
    -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "nparm")    
}

#' @export
#' @rdname utilities
print.gips_fit_class <- function(x, ...){
  ## cat("Method: ", x$method, "\n")
    ## xx <- c(method=x$method, eng=x$engine, x$details)
  ## print(as.data.frame(x$details))
    xx <- x$details
    print(as.data.frame(xx))
  invisible(x)
}

#' @export
#' @rdname utilities
summary.gips_fit_class <- function(object, ...){
  as.data.frame(object$details)
}

#' @export
#' @rdname utilities
glance.gips_fit_class <- function(x, ...){
    as.data.frame(x$details)
}






## logLik.gips_fit_class <- function(object, ...){
##     n <- unname(object$details["n"])
##     if (is.null(n))
##         stop("'n' not given; can not compute log L\n")
##     out <- n * (log(det(object$K)) - sum(object$K * object$S)) / 2
##     di <- sum((abs(object$K) > 1e-12)[upper.tri(object$K, diag=T)])

##     attr(out, "nobs") = n
##     attr(out, "nparm") = di
    
##     class(out) <- "logLik"
##     out
## }



