#' @title Fit Gaussian  graphical models
#' @description Fit Gaussian graphical models using various algorithms.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name fit-ggm
#'
#' @param S Sample covariance matrix.
#' @param edges Generators of model; a list of integer vectors or a 2
#'     x p matrix of integers.
#' @param nobs Number of observations
#' @param K Initial value of concentration matrix.
#' @param iter Maximum number of iterations.
#' @param eps Convergence criterion.
#' @param convcrit Convergence criterions. See section `details`.
#' @param aux A list of form name=value.
#' @param engine Either \code{"R"} or \code{"cpp"}.
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
#' R-based / c++-based in combination with ips / fips.
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
#' cg = fit_ggm(S, gl, nobs=nobs, eps=EPS, engine="cpp", method="fips")
#' ci = fit_ggm(S, gl, nobs=nobs, eps=EPS, engine="cpp", method="ips")
#' rg = fit_ggm(S, gl, nobs=nobs, eps=EPS, engine="R",   method="fips")
#' ri = fit_ggm(S, gl, nobs=nobs, eps=EPS, engine="R",   method="ips")
#'
#' K <- solve(S)
#' (ci$K - K)  %>% abs %>% max
#' (cg$K - K)  %>% abs %>% max
#' (ri$K - K)  %>% abs %>% max
#' (rg$K - K)  %>% abs %>% max
#'
#' cg = fit_ggm(S, em, nobs=nobs, eps=EPS, engine="cpp", method="fips")
#' ci = fit_ggm(S, em, nobs=nobs, eps=EPS, engine="cpp", method="ips")
#' rg = fit_ggm(S, em, nobs=nobs, eps=EPS, engine="R",   method="fips")
#' ri = fit_ggm(S, em, nobs=nobs, eps=EPS, engine="R",   method="ips")
#'
#' K <- solve(S)
#' (ci$K - K)  %>% abs %>% max
#' (cg$K - K)  %>% abs %>% max
#' (ri$K - K)  %>% abs %>% max
#' (rg$K - K)  %>% abs %>% max
#'
NULL


#' @rdname fit-ggm
#' @export
fit_ggm <- function(S, edges=NULL, nobs, K=NULL, iter=10000L, eps=1e-6, convcrit=1, aux=list(),
                    engine="cpp", method="fips", print=0){

    if (inherits(S, "data.frame")){
        nobs = nrow(S)
        S <- cov2cor(cov(S))      
    }
    
    ## t0 <- proc.time()
    t0 <- .get.time()

    if (is.null(K)) {
        ## cat("Generating initial K\n")
        K <- .initK(S)
    }

    
    ## cat("edges (in): \n"); print(edges)
    edges <- parse_edges(edges, nrow(S))
    ## cat("edges (tmp): \n"); print(edges)    
    Elist <- .form2glist(edges)
    ## cat("edges (to use): \n"); print(edges) ## A list
    ## Emat <- list2Emat_(edges, shift=0)

    Emat <- as_elist2emat(Elist)
    amat <- as_emat2amat(Emat, nrow(S))
    aux$amat <- amat ## FIXME: Hack?
    
    engine <- match.arg(tolower(engine), c("cpp", "r"))
    method <- match.arg(tolower(method), c("fips", "ips", "ncd", "cal", "glasso"))

    ## solve_fun <- solve_qr
    
    ## HACK - FIXME
    if (identical(method, "ncd")){
        ## cat("ncd; K is set to solve(S)\n")
        ## K <- solve_fun(S)
        K <- diag(1, nrow(S))
    }

    ## cat(sprintf("TIME: %f\n", .get.diff.time(t0, "millisecs")))

    
    
    if (print)
        cat("engine: ", engine, " method: ", method, "\n")
    
    comb <- paste0(engine, "_", method)
    switch(comb,
           "cpp_fips"    = {fitfun <- .c_fips_ggm_ },
           "r_fips"      = {fitfun <- .r_fips_ggm_ },           
           "cpp_ips"     = {fitfun <- .c_ips_ggm_  },
           "r_ips"       = {fitfun <- .r_ips_ggm_  },
           "cpp_ncd"     = {fitfun <- .c_ncd_ggm_  },
           "r_ncd"       = {fitfun <- .r_ncd_ggm_  },
           ## "r_cal"       = {fitfun <- .r_cal_ggm_  },
           )    
    
    out <- fitfun(S=S, Elist=Elist, Emat=Emat,
                    nobs=nobs, K=K, iter=iter, eps=eps, convcrit=convcrit, print=print, aux=aux)
    out <- c(out, list(edges=Emat, nobs=nobs, eps=eps))

    out   <- .finalize_fit(out, S=S, t0=t0, method=method, engine=engine)
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
    dimnames(out$K)     <- dimnames(S)
    if (inherits(out$Sigma, "matrix"))
        dimnames(out$Sigma) <- dimnames(S)

    nparm <- ncol(out$edges) + nrow(out$K)
    
    if (!is.null(t0))
        out$time <- .get.diff.time(t0, units="millisecs")
    
    out$details <- list(
        method =method,
        engine = engine,
        time   = unname(out$time),
        nobs   = out$nobs,
        iter   = out$iter,
        eps    = out$eps,
        dim    = nparm,
        idim   = nrow(out$K),
        trKS   = sum(out$K * S),
        logL   = out$logL,
        ## For fips / ips the lines below give the same, but for glasso there is no conv_check variable.
        ## made   = mean_abs_diff_on_Emat_(out$Sigma, S, out$edges, 1)
        conv_check  = out$conv_check
    )
    out$time <- out$iter <- out$eps <- NULL
    out$dim  <- out$diff <- NULL
    out <- c(out, dots)
    class(out) <- "gips_fit_class"
    out    
}


#' @export
summary.gips_fit_class <- function(object, ...){
    unlist(object[c("dim", "time", "iter", "edgevisit", "diff")])
}


get_init_parm <- function(S, K){
    solve_fun <- solve_qr
    .parmInit <- function(S){
        list(K=diag(1/diag(S)), Sigma=diag(diag(S)))
    }
    
    if (is.null(K))
        parm <- .parmInit(S)
    else
        parm <- list(K=K, Sigma=solve_fun(K))
    parm
}

## .initSigma  <- function(S){
    ## diag(diag(S))
## }

.initK <- function(S){
    diag(1/diag(S))
}





