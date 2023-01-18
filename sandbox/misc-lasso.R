#' @title Cheap graphical lasso
#'
#' @description Cheap graphical lasso by univariate regression models.
#'
#' @name glasso-misc
#'
#' @param data. A dataframe or matrix.
#' @param lambda Penalty parameter for lasso regressions.
#' @param fit Should model be fitted afterwards.
#' @param penalize.diagonal Should diagonal be penalized.
#' @param zero A 2xp matrix of edges that are not allowed to be in the model.
#' @param edges A 2xp matrix of edges that are allowed to be in the model.
#' @details
#'     Regressions are made in parallel on appropriate operating systems.
#'
#' @examples
#' options("digits"=3)
#' data(math)
#'
#' LAMBDA  <- .8
#' x1 <- glasso_ns(math, lambda=LAMBDA)
#'

## #' @export
## glasso_ns <- function(data., lambda=0, fit=TRUE, edges=NULL, zero=NULL){
##     if (!inherits(data., c("data.frame", "matrix")))
##         stop("data. must be data.frame or matrix\n")
##     if (inherits(data., "data.frame"))
##         data.  <- as.matrix(data.)

##     ## return <- match(return, c("list", "amat", "emat", "graphNEL"))

##     lambda <- lambda / 1    
##     nc <- ncol(data.)

##     out <- lapply(1:nc, .regress_one_variable, data.=data., lambda=lambda)
    
##     amat2 <- amat  <- .make_amat(out, nc)
##     amat2[upper.tri(amat2)] <- 0
##     ee  <- t.default(gRbase::which.arr.index(amat2))
##     em  <- emat_sort(ee)

##     print(em)
    
##     out <- list(amat=amat, emat=em, parents=out, graph=as(amat, "graphNEL"))
##     class(out) <- "ns_glasso_class"
    
##     if (fit) gips18(cov2cor(cov(data.)), edges=em)
##     else out
## }

## Lauritzen / Højsgaard

#' @export
#' @rdname glasso-misc
glasso_lh <- function(data., lambda=0, penalize.diagonal = TRUE, edges=NULL, zero=NULL){
    if (inherits(data., c("data.frame", "matrix"))){
        nobs <- nrow(data.)
        data. <- cov2cor(cov(data.))
    }

    if (!is.null(zero)){        
        esat <- emat_saturated_model(1:nrow(data.))
        emat <- emat_complement(edges, esat)
    } else {
        if (is.null(edges))
            emat <- emat_saturated_model(1:nrow(data.))
        else
            emat  <- edges
    }
    
    gips18(data., edges=emat, aux=list(method="lasso", lambda=lambda,
                                penalize.diagonal=penalize.diagonal))
}

## Friedman / Hastie / Tibshirani

#' @export
#' @rdname glasso-misc
glasso_fht <- function(data., lambda=0, penalize.diagonal = TRUE, edges=NULL, zero=NULL){
    if (inherits(data., c("data.frame", "matrix"))){
        nobs <- nrow(data.)
        data. <- cov2cor(cov(data.))
    }

    
    out <- glasso(data., rho=lambda, nobs=nobs, zero=zero,
                  penalize.diagonal=penalize.diagonal) %>% .iznogood()

    vn <- colnames(data.)
    dimnames(out$K) <- list(vn, vn)
    dimnames(out$S) <- list(vn, vn)
    
    class(out) <- "glasso_fit"
    out
}

.iznogood <- function(x){
  #x$emat <- emat_sort(as_amat2emat(x$wi))
  #x$graph <- as_K2graph(x$wi)
  idx <- match(c("w", "wi", "loglik"), names(x))
  names(x)[idx] <- c("S","K", "logL")
  x
}



.regress_one_variable <- function(j, data., lambda){    
    m <- glmnet::glmnet(x=data.[, -j], y=data.[, j], lambda=lambda,
                        intercept=FALSE,
                        standardize.response=FALSE,
                        standardize=FALSE)        
    eps <- 0.0001 #.Machine$double.eps
    b <- coef(m)[, 1]
    bn <- names(b[abs(b) > eps ])#[-1]
        idx <- match(bn, colnames(data.))
    idx
}

.make_amat <- function(out, nc){
    amat <- matrix(0, nrow=nc, ncol=nc)
    for (j in 1:nc){
        amat[j, out[[j]]] <- 1
        amat[out[[j]], j] <- 1
    }
    amat        
}





#' @export
print.ns_glasso_class <- function(x, ...){
    print(names(x))
    invisible(x)
}

#' @method plot gips_fit_class
#' @export
plot.gips_fit_class <- function(x, ...){
    x$K  %>% as_K2graph  %>% plot
}

#' @method plot glasso_fit
#' @export
plot.glasso_fit <- function(x, ...){
    x$K  %>% as_K2graph  %>% plot
}

