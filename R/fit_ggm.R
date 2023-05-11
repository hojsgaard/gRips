#' @title Fit Gaussian graphical models
#' @description Fit Gaussian graphical models using various algorithms.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name fit_ggm
#'
#' @param S Sample covariance matrix.
#' @param edges Generators of model; a list of integer vectors or a 2
#'     x p matrix of integers.
#' @param nobs Number of observations
#' @param K Initial value of concentration matrix.
#' @param maxit Maximum number of iterations.
#' @param eps Convergence criterion.
#' @param convcrit Convergence criterions. See section `details`.
#' @param aux A list of form name=value.
#' @param method Either `"ips"` or `"fips"`.
#' @param print Should output from fitting be printed?
#' 
#' @details
#'
#' Convergence criterion:
#'
#' * 1: max absolute difference between S and Sigmahat on edges.
#'
#' * 2: difference in log likelihood divided by number of parameters in the
#' model (number of edges + number of nodes) between successive
#' iterations.
#' 
#' R-based / c++-based in combination with con / cov.
#'
#' @examples
#' options("digits"=3)
#' data(math)
#'
#' S <- cov(math)
#' nobs <- nrow(math)
#' gl <- list(1:3, 3:5)
#' em <- matrix(c(1,2, 2,3, 1,3, 3,4, 3,5, 4,5), nrow=2)
#'
#' EPS = 1e-2
#'
#' cg = fit_ggm(S, gl, nobs=nobs, eps=EPS, method="cov")
#' ci = fit_ggm(S, gl, nobs=nobs, eps=EPS, method="con")
#' rg = fit_ggm(S, gl, nobs=nobs, eps=EPS, method="cov", aux=list(engine="R"))
#' ri = fit_ggm(S, gl, nobs=nobs, eps=EPS, method="con", aux=list(engine="R"))
#'
#' K <- solve(S)
#' (ci$K - K)  %>% abs %>% max
#' (cg$K - K)  %>% abs %>% max
#' (ri$K - K)  %>% abs %>% max
#' (rg$K - K)  %>% abs %>% max
#'
NULL

#' @rdname fit_ggm
#' @export
fit_ggm <- function(S, edges=NULL, nobs, K=NULL, maxit=10000L, eps=1e-6, convcrit=1, aux=list(),
                    method="covips", print=0){

    t0 <- .get.time()
    method <- match.arg(tolower(method),
                        c("covips", "conips", "ncd", "sncd"))

    method_str <- method
    
    edges <- parse_edges(edges, nrow(S))
    elist <- .form2glist(edges)
    emat  <- as_elist2emat(elist)
    
    ig <- as_emat2igraph(emat, nrow(S))
    max_coreness <- max(coreness(ig))
    if (max_coreness > (nobs - 1)){  ## FIXME: evt nobs-2 (for det er antal frihedsgrader -1)
        stop(glue("Max coreness ({max_coreness}) is larger than nobs ({nobs}); mle may not exist.\n"))
    }


    switch(method,
           "sncd" = {ver=0; method="ncd"},
           "ncd"  = {ver=1},
           "covips"={ver=0},
           "conips"={ver=0}
           )

    
    aux0 <- list(method  = method,
                 version = ver,                 
                 engine  = "cpp")
    
    engine <- match.arg(tolower(aux0$engine), c("cpp", "r"))
    aux0[names(aux)] <- aux    
    
    if (is.null(K)) {
        if (identical(method, "ncd")){
            K <- diag(1, nrow(S))
        } else {
            K <- diag(1/diag(S))
        }
    }

    Ks <- .c_clone(K)
    t0 <- .get.time()
    comb <- paste0(engine, "_", method)
    switch(comb,
           "cpp_covips"     = {fitfun <- .c_covips_ggm_ },
           "cpp_conips"     = {fitfun <- .c_conips_ggm_ },
           "cpp_ncd"        = {fitfun <- .c_ncd_ggm_ },
           "r_covips"       = {fitfun <- .r_covips_ggm_ },           
           "r_conips"       = {fitfun <- .r_conips_ggm_ },
           "r_ncd"          = {fitfun <- .r_ncd_ggm_ },
           ## "r_cal"       = {fitfun <- .r_cal_ggm_  },           
           )    

    out <- fitfun(S=S, elist=elist, emat=emat,
                  nobs=nobs, K=Ks, maxit=maxit, eps=eps, convcrit=convcrit, print=print, aux=aux0)
    
    out <- c(out, list(edges=emat, nobs=nobs, eps=eps, max_coreness=max_coreness))
    out <- .finalize_fit(out, S=S, t0=t0, method=method_str, engine=engine)
    class(out) <- "gips_fit_class"
    out
}


## 1. edges is right hand sided formula -> returns list
## 2. edges is list -> returns list
## 3. edges is matrix -> returns matrix

parse_edges <- function(edges, nvar){

    if (is.null(edges))
        return(matrix(NA, nrow=2, ncol=0))
    
    if (inherits(edges, "matrix"))
        return(edges)
    
    if (inherits(edges, "list")) {
        if (length(edges) == 0)
            return(matrix(NA, nrow=2, ncol=0))
        else 
            return(edges)  
    }
            
    if (inherits(edges, "formula")){
        st <- gRbase::rhsf2vec(edges)

        if ((length(st) == 1) && (st %in% c(".^1", ".^."))){
            edges <-
                switch(st,
                       ".^1"={matrix(NA, nrow=2, ncol=0)},
                       ".^."={emat_saturated_model(1:nvar)}
                       )
        } else {
            st <- gRbase::rhsf2list(edges)
            edges <- lapply(st, as.numeric)
        }        
        return(edges)
    }
}

.form2glist <- function(glist){

    if (is.list(glist)){
        return(glist)
    }

    if (is.matrix(glist)){
        if (nrow(glist) == 2){
            return(colmat2list(glist))
        } 
        else stop("Need 2 x p or p x 2 matrix")      
    }
    else 
        stop("Need list or matrix")
}

.finalize_fit <- function(out, S=S, t0=NULL, method, engine, ...){
    ## cat(".finalize_fit\n")

    
    dots <- list(...)
    dimnames(out$K) <- dimnames(S)
    if (inherits(out$Sigma, "matrix"))
        dimnames(out$Sigma) <- dimnames(S)

    nparm <- ncol(out$edges) + nrow(out$K)
    
    if (!is.null(t0))
        out$time <- .get.diff.time(t0, units="millisecs")

    ## A HACK; converged is defined for NCD only
    if (is.null(out$converged))
        out$converged = TRUE
    
    trKS <- if (out$converged){
                sum(out$K * S)                
            } else {
                -1
            }

    
    
    out$details <- list(
        method = method,
        ## ver    = out$ver,
        eng    = engine,
        time   = unname(out$time),
        nobs   = out$nobs,
        iter   = out$iter,
        eps    = out$eps,
        dim    = nparm,
        idim   = nrow(out$K),
        trKS   = trKS,
        logL   = out$logL,
        ## For cov / con the lines below give the same, but for glasso there is no conv_check variable.
        ## made   = mean_abs_diff_on_emat_(out$Sigma, S, out$edges, 1)
        conv   = out$conv_check,
        dgap   = out$gap       
    )
    out$time <- out$iter <- out$eps <- NULL
    out$dim  <- out$diff <- NULL
    out$logL <- out$gap <- out$conv_check <- NULL
    
    out <- c(out, dots)
    class(out) <- "gips_fit_class"
    out    
}


#' @export
summary.gips_fit_class <- function(object, ...){
    unlist(object[c("dim", "time", "iter", "edgevisit", "diff")])
}


get_init_parm <- function(S, K){

    .parmInit <- function(S){
        list(K=diag(1/diag(S)), Sigma=diag(diag(S)))
    }
    
    if (is.null(K))
        parm <- .parmInit(S)
    else
        parm <- list(K=K, Sigma=solve_fun(K))
    parm
}

## .initK <- function(S){
##     diag(1/diag(S))
## }





