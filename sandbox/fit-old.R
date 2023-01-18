#' @title Fit graphical Gaussian models
#' 
#' @description Functions for fitting graphical Gaussian models
#'     (ggms); deliberately with a simple interface. These are not
#'     intended to be called directly by the user but by more proper
#'     "modelling functions".
#'
#' @name old-fit-ggm
#'
#' @param S Sample covariance matrix.
#' @param edges Generators of model; a list of integer vectors or a 2
#'     x p matrix of integers.
#' @param nobs Number of observations.
#' @param zero FIXME
#' @param iter Maximum number of iterations.
#' @param eps Convergence criterion.
#' @param engine Either \code{"R"} or \code{"cpp"}.
#' @param check Convergence check; either 1 or 2. FIXME
#' @param aux List specifying various restrictions on models.
NULL


#' @rdname old-fit-ggm
#' @export
gips18 <- function (S, edges=NULL, nobs=NULL, zero=NULL,
                    iter=1000, eps=1e-6, engine="cpp", check=1, aux=list())
{ ## FIXME: engine: WHERE 
    engine <- match.arg(tolower(engine), c("cpp", "r"))

    ## Only difference between 17 and 18
    fitfun <- if (identical(engine, "cpp")) .gips18_c else .gips18_r
    edges  <- .form2emat(edges) ## FIXME Needed???; 
    ## End Only...
    
    t0 <- proc.time()
    aux <- .finalize_aux(aux, p=nrow(S))

    ## FIXME: Synes ikke at virke
    if (identical(aux$method, "mtp2")) edges  <- edges[, S[t(edges)] > 0]

    out <- fitfun(S, nobs=nobs, edges=edges, iter=iter, eps=eps, check=check, aux=aux)
    .finalize_fit(out, S=S, t0=t0, method="gips", fitfun=fitfun)    
}


## ## 'unified ips' (c++); requires edges
## #' @rdname old-fit-ggm
## #' @export
## gips17 <- function (S, edges=NULL, nobs=NULL, zero=NULL,
##                     iter=1000, eps=1e-6, engine="cpp", aux=list())
## { ## FIXME: engine WHERE and No check argumuent

##     engine <- match.arg(tolower(engine), c("cpp", "r"))

##     ## Only difference between 17 and 18
##     fitfun <- if (identical(engine, "cpp")) .gips17_c else .gips17_r
##     edges  <- .form2glist(edges) ## FIXME Needed???; Only difference between 17 and 18
##     ## End Only...
    
##     t0 <- proc.time()
##     aux <- .finalize_aux(aux, p=nrow(S))

##     ## FIXME: Synes ikke at virke
##     if (identical(aux$method, "mtp2")) edges  <- edges[, S[t(edges)] > 0]

##     out <- fitfun(S, nobs=nobs, edges=edges, iter=iter, eps=eps, aux=aux)
##     .finalize_fit(out, S=S, t0=t0, method="gips", fitfun=fitfun)    
## } 

## -------------------------------------------------------------------
## Interfaces
## -------------------------------------------------------------------

.gips18_c <- function (S, nobs=NA, edges=NULL, zero=NULL,
                       iter=1000L, eps=1e-6, check=1, aux=list()){
    K     = .initK(S)        
    Sigma = .initSigma(S)    
    .c_gips18_fit_(S=S, edges=edges, K=K, Sigma=Sigma, 
                   iter=iter, eps=eps, check=check,
                   lambda=aux$lambda, method=aux$method, pos=aux$pos)
    list(K=K, Sigma=Sigma, iter=iter, eps=eps, check=check, nobs=nobs)
}

.gips18_r <- function(S, nobs=NA, edges=NULL, zero=NULL, 
                      iter=1000L, eps=1e-6, check=1, aux=list()){
    .r_gips18_fit_(S, edges=edges, iter=iter, eps=eps,
                  lambda=aux$lambda, method=aux$method, pos=aux$pos)    
}


## .gips17_c <- function(S, nobs=NA, edges=NULL, zero=NULL,
##                       iter=1000L, eps=1e-6, aux=list()){
##     edges <- .form2glist(edges) ## Necessary???    
##     .c_gips17_fit_(S=S, edges=edges, iter=iter, eps=eps,
##                    lambda=aux$lambda, method=aux$method, pos=aux$pos)
## }

## ## edges: 2xp matrix
## .gips17_r <- function (S, nobs=NA, edges=NULL, zero=NULL,
##                        iter=1000L, eps=1e-6, aux=list()){
##     .r_gips17_fit_(S=S, edges=edges, iter=iter, eps=eps, 
##                    lambda=aux$lambda, method=aux$method, pos=aux$pos)
## }


## -------------------------------------------------------------------
## Workers
## -------------------------------------------------------------------

.r_gips18_fit_ <- function(S, edges,
                          iter=1000, eps=1e-6, lambda=0, method="ips", pos=FALSE){  ## FIXME: Strange order...
    it <- 0
    t0 <- proc.time()
    parm <- .parmInit(S)  

    model.fits <- .does_model_fit_to_data(edges, parm, S, eps=eps, method=method)
    
    while((!model.fits) && (it < iter)){
        it <- it + 1
        parm <- .inner_loop18(S, edges, parm, eps=eps,
                              method=method, lambda=lambda, pos=pos)
        model.fits <- .does_model_fit_to_data(edges, parm, S, eps=eps, method=method)
    }
    
    parm$iter <- it
    parm$time <- (proc.time() - t0)[3]
    parm$eps  <- eps
    parm 
}



## .r_gips17_fit_ <- function (S, edges,
##                             iter=1000L, eps=1e-6, lambda=0, method="ips", pos=FALSE) {    
##     it <- 0
##     t0 <- proc.time()
##     parm <- .parmInit(S)  
##     repeat {
##         Sigma.old <- parm$Sigma ## FIXME Careful here
##         it <- it + 1L
##         parm <- .inner_loop_17(S, edges, parm, eps=eps, 
##                                method=method, lambda=lambda, pos=pos)
##         diff <- .maxKdiff(Sigma.old, parm$Sigma)
##         if (diff < eps || it == iter) break
##     }
    
##     parm$iter <- it
##     parm$time <- (proc.time() - t0)[3]    
##     parm$eps  <- eps
##     parm
## }

## ---------------------------------------------------------------

.inner_loop18 <- function(S, edges, parm,
                          eps=1e-6, method, lambda, pos){
    for (j in 1:ncol(edges)){
        cc <- edges[, j]
        edge.fits <- .does_edge_fit_to_data(cc, parm, S, eps=eps, method=method)
        if (!edge.fits) {
            parm <- .update_cc_parm(cc, parm, S, eps=eps,
                                method=method, lambda=lambda, pos=pos) 
        }        
    }
    parm
}


## .inner_loop_17 <- function(S, edges, parm,
##                            eps=1e-6, method, lambda, pos){    
##     for (j in 1:ncol(edges)){
##         cc <- edges[, j]
##         parm <- .update_cc_parm(cc, parm, S, eps=eps,
##                             method=method, lambda=lambda, pos=pos)
##     }
##     parm
## }


















