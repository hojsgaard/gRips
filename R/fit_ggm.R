#' @title Fit Gaussian graphical models
#' @description Fit Gaussian graphical models using various algorithms.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name fit_ggm
#'
#' @param S Sample covariance matrix.
#' @param formula Generators of model; a list of integer vectors or a 2
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
#' fit_cov = fit_ggm(S, gl, nobs=nobs, eps=EPS, method="cov")
#' fit_con = fit_ggm(S, gl, nobs=nobs, eps=EPS, method="con")
#' fit_ncd = fit_ggm(S, gl, nobs=nobs, eps=EPS, method="ncd")
#'
#' K <- solve(S)
#' (fit_con$K - K)  |> abs() |> max()
#' (fit_cov$K - K)  |> abs() |> max()
#' (fit_ncd$K - K)  |> abs() |> max()
#'
NULL

check_coreness <- function(emat, d, nobs){   
    ig <- as_emat2igraph(emat, d)
    max_coreness <- max(coreness(ig))
    if (max_coreness >= (nobs - 1)) {  ## FIXME: evt nobs-2 (for det er antal frihedsgrader -1)
        stop(glue("Max coreness ({max_coreness}) is larger than or equal nobs ({nobs}) minus one; mle may not exist.\n"))
    }
    max_coreness
}

#' @rdname fit_ggm
#' @export
fit_ggm <- function(S, formula=NULL, nobs, K=NULL, maxit=10000L, eps=1e-2, convcrit=1, aux=list(),
                    method="covips", print=0) {
    t0 <- .get.time()
    method <- match.arg(tolower(method),
                        c("covips", "conips", "ncd", "sncd"))

    method_str <- method
    
    formula <- parse_formula(formula, nrow(S))
    elst  <- formula2glist(formula)
    emat  <- as_elist2emat(elst)
    amat  <- as_emat2amat(emat, d=nrow(S))
    ## str(list(formula=formula, elst=elst, emat=emat, amat=amat))

    if (any(abs(diag(S)-1) > 1e-8)){
        S <- cov2cor(S)
    }

    max_coreness <- check_coreness(emat, nrow(S), nobs)
    deg <- rowSums(amat)
    maxdeg <- max(deg)
    
    if ((identical(method,"ncd"))&& (maxdeg>= nobs-1)){
       S <- sweep_start(S, amat,eps=eps,f=nobs-1) 
    }

    switch(method,
           "sncd" = {ver=0; method="ncd"},
           "ncd"  = {ver=1},
           "covips"={ver=0},
           "conips"={ver=0}
           )

    aux0 <- list(method  = method,
                 amat    = amat,
                 version = ver,                 
                 engine  = "cpp")
    
    engine <- match.arg(tolower(aux0$engine), c("cpp", "r"))
    aux0[names(aux)] <- aux    
    
    if (is.null(K)) {
        if (identical(method, "ncd")) {
            K <- diag(1, nrow(S))
        } else {
            K <- diag(1/diag(S))
        }
    }
    Ks <- .c_clone(K)

    t0 <- .get.time()
    comb <- paste0(engine, "_", method)
    switch(comb,
      "cpp_covips"     = {
        fitfun <- .c_covips_ggm_
      },
      "cpp_conips"     = {
        fitfun <- .c_conips_ggm_
      },
      "cpp_ncd"        = {
        fitfun <- .c_ncd_ggm_
      },
      "r_covips"       = {
        fitfun <- .r_covips_ggm_
      },
      "r_conips"       = {
        fitfun <- .r_conips_ggm_
      },
      "r_ncd"          = {
        fitfun <- .r_ncd_ggm_
      }
    )    

    out <- fitfun(S=S, elst=elst, emat=emat,
                  nobs=nobs, K=Ks, maxit=maxit, eps=eps, convcrit=convcrit, print=print, aux=aux0)
    
    out <- c(out, list(edges=emat, nobs=nobs, eps=eps, max_coreness=max_coreness))
    out <- .finalize_fit(out, S=S, t0=t0, method=method_str, engine=engine)
    class(out) <- "gips_fit_class"
    out
}


## 1. formula is right hand sided formula -> returns list
## 2. formula is list -> returns list
## 3. formula is matrix -> returns matrix

parse_formula <- function(formula, nvar) {

    if (is.null(formula))
        return(matrix(NA, nrow=2, ncol=0))
  
    if (inherits(formula, "matrix"))
        return(formula)
    
    if (inherits(formula, "list")) {
        if (length(formula) == 0)
            return(matrix(NA, nrow=2, ncol=0))
        else 
            return(formula)  
    }
            
  if (inherits(formula, "formula")) {
    st <- gRbase::rhsf2vec(formula)
    
    if ((length(st) == 1) && (st %in% c(".^1", ".^."))) {
      formula <-
        switch (st,
                ".^1" = {
                  matrix(NA, nrow = 2, ncol = 0)
                },
                ".^." = {
                  emat_saturated_model(1:nvar)
                })
    } else {
      st <- gRbase::rhsf2list(formula)
      formula <- lapply(st, as.numeric)
    }
    return(formula)
  }
}

formula2glist <- function(glist) {

    if (is.list(glist)) {
        return(glist)
    }

    if (is.matrix(glist)) {
        if (nrow(glist) == 2) {
            return(colmat2list(glist))
        } 
        else stop("Need 2 x p or p x 2 matrix")      
    }
    else 
        stop("Need list or matrix")
}

.finalize_fit <- function(out, S=S, t0=NULL, method, engine, ...) {
    ## cat(".finalize_fit\n")

    
    dots <- list(...)
    dimnames(out$K) <- dimnames(S)
    if (inherits(out$Sigma, "matrix"))
        dimnames(out$Sigma) <- dimnames(S)

    nparm <- ncol(out$edges) + nrow(out$K)
    
    if (!is.null(t0))
        out$time <- round(.get.diff.time(t0, units="secs"), 2)

    ## Variable 'converged' is defined for NCD only
    if (is.null(out$converged))
        out$converged = TRUE

    if (is.null(out$mev))
        out$mev = NA

    
    if (out$converged) {
        trKS = sum(out$K * S)
        logL = out$logL
        conv = out$conv_check
        dgap = out$gap
    } else {
        trKS <- logL <- conv <- dgap <- NA
    }
    
    if (!identical(method, "ncd"))
        dgap <- NA
    
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
        logL   = logL,
        conv   = conv,
        dgap   = dgap,
        mev    = out$mev
    )
    out$time <- out$iter <- out$eps <- NULL
    out$dim  <- out$diff <- NULL
    out$logL <- out$gap <- out$conv_check <- NULL
    
    out <- c(out, dots)
    class(out) <- "gips_fit_class"
    out    
}




get_init_parm <- function(S, K) {

    .parmInit <- function(S) {
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





